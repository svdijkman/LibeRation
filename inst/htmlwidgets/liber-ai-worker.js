/* LibeRation browser-local AI worker.
 *
 * The WebLLM runtime and selected weights are downloaded only when the user
 * first invokes AI. Once the engine reports ready, network-capable APIs in
 * this worker are disabled before any model prompt is accepted. The model is
 * deliberately given no tools or browser/DOM access.
 */
let runtime = null;
let engine = null;
let loadedModel = null;
let loadedContextWindow = 0;
let networkLocked = false;
let busy = false;
const nativeFetch = self.fetch.bind(self);
const downloadHosts = new Set([
  "esm.run", "cdn.jsdelivr.net", "huggingface.co", "cdn-lfs.huggingface.co",
  "cas-bridge.xethub.hf.co", "raw.githubusercontent.com", "github.com"
]);

// Install the gate before importing WebLLM. Even if the runtime retains this
// function, it observes `networkLocked` on every call. Only model/runtime
// artifact hosts are reachable during the pre-prompt initialization phase.
self.fetch = function (input, init) {
  if (networkLocked) {
    return Promise.reject(new Error("Network access is disabled for LibeRation local AI inference."));
  }
  let url;
  try { url = new URL(typeof input === "string" ? input : input.url); }
  catch (error) { return Promise.reject(new Error("Blocked non-URL AI download request.")); }
  if (url.protocol !== "https:" || !downloadHosts.has(url.hostname)) {
    return Promise.reject(new Error("Blocked AI download host: " + url.hostname));
  }
  return nativeFetch(input, init);
};

function send(type, detail) {
  self.postMessage(Object.assign({ type: type }, detail || {}));
}

function messageText(chunk) {
  try {
    return chunk.choices[0].delta.content || "";
  } catch (error) {
    return "";
  }
}

function lockNetwork() {
  if (networkLocked) return;
  const denied = function () {
    throw new Error("Network access is disabled for LibeRation local AI inference.");
  };
  try { self.XMLHttpRequest = denied; } catch (error) {}
  try { self.WebSocket = denied; } catch (error) {}
  try { self.EventSource = denied; } catch (error) {}
  try { self.navigator.sendBeacon = function () { return false; }; } catch (error) {}
  networkLocked = true;
  send("network_locked", { locked: true });
}

function normalizeContextWindow(value) {
  const numeric = Math.round(Number(value) || 4096);
  return Math.max(1024, Math.min(16384, numeric));
}

function contextFallbacks(requested) {
  return [requested, 12288, 8192, 6144, 4096, 1024].filter(function (size, index, values) {
    return size <= requested && values.indexOf(size) === index;
  });
}

function isMemoryAllocationFailure(error) {
  const message = error && error.message ? error.message : String(error || "");
  return /out of memory|memory allocation|insufficient memory|allocation failed|buffer.*(?:allocation|size)|exceeds.*(?:buffer|storage|memory)|GPUValidationError/i.test(message);
}

async function ensureEngine(model, requestedContextWindow) {
  const requested = normalizeContextWindow(requestedContextWindow);
  if (engine && loadedModel === model && loadedContextWindow === requested) return requested;
  if (networkLocked) {
    throw new Error("Changing the local model or context requires a new worker session.");
  }
  if (!self.isSecureContext || !self.navigator.gpu) {
    throw new Error("WebGPU is unavailable. Use a current WebGPU-enabled browser on localhost or HTTPS.");
  }
  send("status", { stage: "adapter", text: "Checking the browser WebGPU adapter" });
  const adapter = await self.navigator.gpu.requestAdapter({ powerPreference: "high-performance" });
  if (!adapter) {
    throw new Error("No usable WebGPU adapter is available in this browser session.");
  }
  send("status", { stage: "runtime", text: "Loading the local WebGPU runtime" });
  runtime = runtime || await import("https://esm.run/@mlc-ai/web-llm@0.2.84");
  const appConfig = Object.assign({}, runtime.prebuiltAppConfig, {
    cacheBackend: "indexeddb"
  });
  const engineConfig = {
    appConfig: appConfig,
    initProgressCallback: function (progress) {
      send("progress", {
        progress: Number(progress.progress || 0),
        text: String(progress.text || "Preparing browser-local model")
      });
    }
  };
  const candidates = contextFallbacks(requested);
  let failure = null;
  for (const contextWindow of candidates) {
    try {
      send("status", { stage: "context", text: "Allocating " + (contextWindow / 1024).toFixed(contextWindow % 1024 ? 1 : 0) + "K model context" });
      engine = await runtime.CreateMLCEngine(model, engineConfig, {
        context_window_size: contextWindow
      });
      loadedModel = model;
      loadedContextWindow = contextWindow;
      if (contextWindow !== requested) {
        send("context_fallback", {
          requested_context_window_size: requested,
          context_window_size: contextWindow
        });
      }
      break;
    } catch (error) {
      failure = error;
      const staleEngine = engine;
      engine = null;
      loadedModel = null;
      loadedContextWindow = 0;
      if (staleEngine && staleEngine.unload) {
        try { await staleEngine.unload(); } catch (unloadError) {}
      }
      if (!isMemoryAllocationFailure(error) || contextWindow === candidates[candidates.length - 1]) throw error;
      send("status", { stage: "fallback", text: "GPU memory was insufficient; trying a smaller context" });
    }
  }
  if (!engine) throw failure || new Error("The browser-local model could not be loaded.");
  lockNetwork();
  send("ready", { model: model, context_window_size: loadedContextWindow, network_locked: true });
  return loadedContextWindow;
}

function isDisposedFailure(error) {
  const message = error && error.message ? error.message : String(error || "");
  return /already been disposed|disposed object|device (?:was )?lost|DXGI_ERROR_DEVICE_(?:REMOVED|HUNG|RESET)|requestDevice|create command queue failed|GPUAdapter/i.test(message);
}

function isContextFailure(error) {
  const message = error && error.message ? error.message : String(error || "");
  return /prompt tokens exceed context window|context window size|maximum context|too many tokens/i.test(message);
}

function clipContext(text, limit) {
  text = String(text || "");
  if (text.length <= limit) return text;
  const marker = "\n[... context shortened after tokenizer preflight ...]\n";
  const available = Math.max(20, limit - marker.length);
  const head = Math.ceil(available * 0.72);
  return text.slice(0, head) + marker + text.slice(-(available - head));
}

function compactContextRequest(request, targetWindow) {
  const windowSize = normalizeContextWindow(targetWindow || request.context_window_size);
  const maxTokens = Math.min(Number(request.max_tokens) || 800, windowSize <= 1024 ? 160 : Math.max(384, Math.floor(windowSize * 0.2)));
  const charBudget = Math.max(600, Math.floor((windowSize - maxTokens - (windowSize <= 1024 ? 128 : 512)) * 2.7));
  const source = Array.isArray(request.messages) ? request.messages.map(function (item) {
    return { role: String(item.role || "user"), content: String(item.content || "") };
  }) : [];
  const first = source.length && source[0].role === "system" ? source.shift() : null;
  const last = source.length ? source.pop() : null;
  const lastBudget = last ? Math.min(last.content.length, Math.max(320, Math.floor(charBudget * 0.16))) : 0;
  const output = [];
  let used = 0;
  if (first) {
    const firstBudget = Math.min(first.content.length, Math.max(720, Math.floor(charBudget * 0.46)), charBudget - lastBudget);
    const content = clipContext(first.content, firstBudget);
    output.push({ role: first.role, content: content });
    used += content.length;
  }
  let middleBudget = Math.max(0, charBudget - used - lastBudget);
  const kept = [];
  for (let index = source.length - 1; index >= 0 && middleBudget > 80; index -= 1) {
    const content = clipContext(source[index].content, Math.min(3600, middleBudget));
    kept.unshift({ role: source[index].role, content: content });
    middleBudget -= content.length;
  }
  if (kept.length && kept[0].role === "assistant") kept.shift();
  output.push.apply(output, kept);
  if (last) output.push({ role: last.role, content: clipContext(last.content, lastBudget || 160) });
  return Object.assign({}, request, {
    messages: output,
    max_tokens: maxTokens,
    context_window_size: windowSize,
    context_compacted: true
  });
}

async function recoverEngine(model) {
  const staleEngine = engine;
  engine = null;
  loadedModel = null;
  loadedContextWindow = 0;
  // Re-open the tightly allow-listed download gate while CreateMLCEngine reads
  // the already cached runtime/model artifacts. It is locked again before any
  // project context is supplied to the replacement engine.
  networkLocked = false;
  if (staleEngine && staleEngine.unload) {
    try { await staleEngine.unload(); } catch (error) {}
  }
  await ensureEngine(model, 4096);
}

async function generateOnce(request, progress) {
    const actualContextWindow = await ensureEngine(request.model, request.context_window_size);
    if (actualContextWindow < normalizeContextWindow(request.context_window_size)) {
      request = compactContextRequest(request, actualContextWindow);
    }
    const messages = Array.isArray(request.messages) ? request.messages.map(function (item) {
      return { role: String(item.role || "user"), content: String(item.content || "") };
    }) : [];
    const stream = await engine.chat.completions.create({
      messages: messages,
      temperature: Number.isFinite(Number(request.temperature)) ? Number(request.temperature) : 0.1,
      top_p: Number.isFinite(Number(request.top_p)) ? Number(request.top_p) : 0.8,
      max_tokens: Math.max(64, Math.min(4096, Number(request.max_tokens) || 1200)),
      stream: true
    });
    let answer = "";
    for await (const chunk of stream) {
      const token = messageText(chunk);
      if (token) {
        answer += token;
        progress.emitted = true;
        send("token", { id: request.id, token: token });
      }
    }
    send("complete", { id: request.id, text: answer });
}

async function generate(request) {
  if (busy) throw new Error("The local model is already generating a response.");
  busy = true;
  const progress = { emitted: false };
  try {
    try {
      await generateOnce(request, progress);
    } catch (error) {
      if (!progress.emitted && isContextFailure(error) && !request.context_compacted) {
        send("status", { stage: "compacting", text: "Reducing project context to fit the selected local model" });
        try {
          await generateOnce(compactContextRequest(request), progress);
        } catch (compactError) {
          if (isContextFailure(compactError)) {
            throw new Error("The selected model could not fit even the reduced project context. Ask a narrower question or choose a model with a larger context window. Technical detail: " + (compactError.message || String(compactError)));
          }
          throw compactError;
        }
      } else if (!progress.emitted && isDisposedFailure(error)) {
        send("status", { stage: "recovering", text: "Refreshing the local GPU session after a graphics-device reset" });
        try {
          await recoverEngine(request.model);
          request = compactContextRequest(request, 4096);
          await generateOnce(request, progress);
        } catch (recoveryError) {
          // Never leave the temporary artifact-download gate open when an
          // automatic recovery itself fails.
          lockNetwork();
          throw recoveryError;
        }
      } else {
        throw error;
      }
    }
  } finally {
    busy = false;
  }
}

self.onmessage = function (event) {
  const request = event.data || {};
  if (request.type === "generate") {
    generate(request).catch(function (error) {
      busy = false;
      send("error", { id: request.id, message: error && error.message ? error.message : String(error) });
    });
  } else if (request.type === "interrupt") {
    try { if (engine && engine.interruptGenerate) engine.interruptGenerate(); } catch (error) {}
  } else if (request.type === "dispose") {
    Promise.resolve(engine && engine.unload ? engine.unload() : null).catch(function () {}).finally(function () { self.close(); });
  } else if (request.type === "capabilities") {
    send("capabilities", {
      secure_context: !!self.isSecureContext,
      webgpu: !!self.navigator.gpu,
      network_locked: networkLocked,
      model: loadedModel,
      context_window_size: loadedContextWindow
    });
  }
};
