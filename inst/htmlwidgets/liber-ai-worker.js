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

async function ensureEngine(model) {
  if (engine && loadedModel === model) return;
  if (networkLocked) {
    throw new Error("Changing the local model requires a new worker session.");
  }
  if (!self.isSecureContext || !self.navigator.gpu) {
    throw new Error("WebGPU is unavailable. Use a current WebGPU-enabled browser on localhost or HTTPS.");
  }
  send("status", { stage: "runtime", text: "Loading the local WebGPU runtime" });
  runtime = runtime || await import("https://esm.run/@mlc-ai/web-llm@0.2.84");
  const appConfig = Object.assign({}, runtime.prebuiltAppConfig, {
    cacheBackend: "indexeddb"
  });
  engine = await runtime.CreateMLCEngine(model, {
    appConfig: appConfig,
    initProgressCallback: function (progress) {
      send("progress", {
        progress: Number(progress.progress || 0),
        text: String(progress.text || "Preparing browser-local model")
      });
    }
  });
  loadedModel = model;
  lockNetwork();
  send("ready", { model: model, network_locked: true });
}

function isDisposedFailure(error) {
  const message = error && error.message ? error.message : String(error || "");
  return /already been disposed|disposed object|device (?:was )?lost/i.test(message);
}

async function recoverEngine(model) {
  const staleEngine = engine;
  engine = null;
  loadedModel = null;
  // Re-open the tightly allow-listed download gate while CreateMLCEngine reads
  // the already cached runtime/model artifacts. It is locked again before any
  // project context is supplied to the replacement engine.
  networkLocked = false;
  if (staleEngine && staleEngine.unload) {
    try { await staleEngine.unload(); } catch (error) {}
  }
  await ensureEngine(model);
}

async function generateOnce(request, progress) {
    await ensureEngine(request.model);
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
      if (!progress.emitted && isDisposedFailure(error)) {
        send("status", { stage: "recovering", text: "Refreshing the local GPU session after a WebGPU cache reset" });
        try {
          await recoverEngine(request.model);
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
      model: loadedModel
    });
  }
};
