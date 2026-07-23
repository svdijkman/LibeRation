(function () {
  "use strict";

  var e = React.createElement;
  var palette = ["#466d91", "#a36d3d", "#447963", "#745f86", "#a64e55", "#38777c", "#82713e"];

  function list(x) { return Array.isArray(x) ? x : []; }
  function value(x, fallback) { return x === null || x === undefined || x === "" ? fallback : x; }
  function number(x) { var n = Number(x); return isFinite(n) ? n : null; }
  function initialDarkTheme(legacyKey) {
    try {
      var shared = window.localStorage.getItem("liber.theme");
      if (shared === "dark" || shared === "light") return shared === "dark";
      var legacy = window.localStorage.getItem(legacyKey);
      if (legacy === "1" || legacy === "dark") return true;
      if (legacy === "0" || legacy === "light") return false;
    } catch (error) {}
    return !!(window.matchMedia && window.matchMedia("(prefers-color-scheme: dark)").matches);
  }
  function storeTheme(dark, legacyKey, numericLegacy) {
    try {
      window.localStorage.setItem("liber.theme", dark ? "dark" : "light");
      window.localStorage.setItem(legacyKey, numericLegacy ? (dark ? "1" : "0") : (dark ? "dark" : "light"));
      document.documentElement.setAttribute("data-liber-theme", dark ? "dark" : "light");
    } catch (error) {}
  }
  function useDialogFocus(open, onClose) {
    var dialog = React.useRef(null), close = React.useRef(onClose);
    close.current = onClose;
    React.useEffect(function () {
      if (!open) return;
      var prior = document.activeElement, node = dialog.current;
      function focusable() {
        return node ? Array.prototype.slice.call(node.querySelectorAll(
          'button:not([disabled]),input:not([disabled]),select:not([disabled]),textarea:not([disabled]),a[href],[tabindex]:not([tabindex="-1"])'
        )) : [];
      }
      function keydown(event) {
        if (event.key === "Escape") { event.preventDefault(); close.current(); return; }
        if (event.key !== "Tab" || !node) return;
        var items = focusable();
        if (!items.length) { event.preventDefault(); node.focus(); return; }
        var first = items[0], last = items[items.length - 1];
        if (event.shiftKey && document.activeElement === first) { event.preventDefault(); last.focus(); }
        else if (!event.shiftKey && document.activeElement === last) { event.preventDefault(); first.focus(); }
      }
      document.addEventListener("keydown", keydown);
      window.setTimeout(function () { var items = focusable(); (items[0] || node).focus(); }, 0);
      return function () {
        document.removeEventListener("keydown", keydown);
        if (prior && prior.focus) prior.focus();
      };
    }, [open]);
    return dialog;
  }
  function formatNumber(x) {
    var n = Number(x);
    if (!isFinite(n)) return "-";
    if (n !== 0 && (Math.abs(n) < 0.001 || Math.abs(n) >= 10000)) return n.toExponential(4);
    return n.toPrecision(6).replace(/\.?0+$/, "");
  }
  function emit(props, action, detail) {
    if (!window.Shiny || !window.Shiny.setInputValue) return false;
    window.Shiny.setInputValue(
      (props.inputId || "liber_workbench") + "_event",
      Object.assign({ action: action, nonce: Date.now() }, detail || {}),
      { priority: "event" }
    );
    return true;
  }
  function cloneRows(rows) { return list(rows).map(function (row) { return Object.assign({}, row); }); }
  function useSynced(initial, dependency) {
    var state = React.useState(initial);
    React.useEffect(function () { state[1](initial); }, dependency || []);
    return state;
  }

  var localAIContextPresets=[1024,4096,6144,8192,12288,16384];
  function localAIModelInfo(ai,model) {
    return list(ai&&ai.models).filter(function(item){return item.id===model;})[0]||{};
  }
  function localAIAutoContextWindow(ai,model,purpose) {
    if(/(?:^|-)1k(?:-|$)/i.test(String(model||"")))return 1024;
    var memory=Number(localAIModelInfo(ai,model).vram_mb)||0;
    if(memory>5500)return 4096;
    if(memory>3000)return purpose==="report"?6144:6144;
    return 8192;
  }
  function localAIContextWindow(ai,model,purpose,override) {
    var configured=override;
    if(configured===undefined||configured===null||configured==="")configured=purpose==="report"?ai&&ai.report_context:ai&&ai.help_context;
    if(!configured||String(configured).toLowerCase()==="auto")return localAIAutoContextWindow(ai,model,purpose);
    var numeric=Math.round(Number(configured));
    return isFinite(numeric)?Math.max(1024,Math.min(16384,numeric)):localAIAutoContextWindow(ai,model,purpose);
  }
  function localAISettingsDetail(ai) {
    return {activated:!!(ai&&ai.activated),consented:!!(ai&&ai.consented),help_model:value(ai&&ai.help_model,ai&&ai.model),report_model:value(ai&&ai.report_model,"same_as_help"),help_context:value(ai&&ai.help_context,"auto"),report_context:value(ai&&ai.report_context,"auto")};
  }
  function localAIClip(text, limit) {
    text=String(text||"");limit=Math.max(80,Number(limit)||80);
    if(text.length<=limit)return text;
    var marker="\n[... context shortened by LibeRation ...]\n",available=Math.max(20,limit-marker.length);
    var head=Math.ceil(available*.72),tail=Math.max(0,available-head);
    return text.slice(0,head)+marker+(tail?text.slice(-tail):"");
  }
  function localAIBudgetMessages(messages, model, requestedMaxTokens, contextWindow) {
    contextWindow=Math.max(1024,Number(contextWindow)||4096);
    var outputCeiling=contextWindow<=1024?256:Math.max(512,Math.floor(contextWindow*.25));
    var maxTokens=Math.max(64,Math.min(Number(requestedMaxTokens)||800,outputCeiling,2200));
    var safety=contextWindow<=1024?128:512;
    var inputTokenBudget=Math.max(256,contextWindow-maxTokens-safety);
    var inputBudget=Math.max(720,Math.floor(inputTokenBudget*3));
    var source=list(messages).map(function(item){return {role:String(item&&item.role||"user"),content:String(item&&item.content||"")};});
    var originalChars=source.reduce(function(total,item){return total+item.content.length;},0),originalCount=source.length;
    if(!source.length)return {messages:source,max_tokens:maxTokens,context_window_size:contextWindow,input_char_budget:inputBudget,prompt_tokens_estimated:0,original_message_count:0,retained_message_count:0,compacted:false};
    if(originalChars<=inputBudget)return {messages:source,max_tokens:maxTokens,context_window_size:contextWindow,input_char_budget:inputBudget,prompt_tokens_estimated:Math.ceil(originalChars/3),original_message_count:originalCount,retained_message_count:originalCount,compacted:false};
    var first=source[0].role==="system"?source.shift():null,last=source.length?source.pop():null;
    var reservedLast=last?Math.min(last.content.length,Math.max(480,Math.floor(inputBudget*.16))):0;
    var reservedFirst=first?Math.min(first.content.length,Math.max(900,Math.floor(inputBudget*.46))):0;
    if(first&&reservedFirst+reservedLast>inputBudget){reservedFirst=Math.max(300,inputBudget-reservedLast);}
    var output=[],used=0;
    if(first){first.content=localAIClip(first.content,reservedFirst);output.push(first);used+=first.content.length;}
    var middleBudget=Math.max(0,inputBudget-used-reservedLast),kept=[];
    for(var index=source.length-1;index>=0&&middleBudget>80;index--){
      var item=source[index],content=localAIClip(item.content,Math.min(4200,middleBudget));
      if(content.length>middleBudget)content=localAIClip(content,middleBudget);
      kept.unshift({role:item.role,content:content});middleBudget-=content.length;
    }
    if(kept.length&&kept[0].role==="assistant")kept.shift();
    output=output.concat(kept);
    if(last){last.content=localAIClip(last.content,Math.max(120,inputBudget-output.reduce(function(total,item){return total+item.content.length;},0)));output.push(last);}
    var retainedChars=output.reduce(function(total,item){return total+item.content.length;},0);
    return {messages:output,max_tokens:maxTokens,context_window_size:contextWindow,input_char_budget:inputBudget,prompt_tokens_estimated:Math.ceil(retainedChars/3),original_message_count:originalCount,retained_message_count:output.length,compacted:retainedChars<originalChars||output.length<originalCount};
  }

  var localAI = { worker:null, workerUrl:"", model:"", contextWindow:0, purpose:"", status:{stage:"idle",text:"Model not loaded",progress:0,locked:false,model:"",purpose:"",budget:null}, listeners:[], pending:{}, cooldownUntil:0 };
  function localAINotify() { localAI.listeners.slice().forEach(function (listener) { listener(Object.assign({}, localAI.status)); }); }
  function localAISetStatus(next) { localAI.status=Object.assign({},localAI.status,next);localAINotify(); }
  function localAIRejectPending(reason) {
    Object.keys(localAI.pending).forEach(function(id){localAI.pending[id].reject(reason);});
    localAI.pending={};
  }
  function localAIHasPendingPurpose(purpose) {
    return Object.keys(localAI.pending).some(function(id){return localAI.pending[id].purpose===purpose;});
  }
  function localAIRecoverableGPUFailure(reason) {
    var message=reason&&reason.message?reason.message:String(reason||"");
    return /already been disposed|disposed object|device (?:was )?lost|DXGI_ERROR_DEVICE_(?:REMOVED|HUNG|RESET)|requestDevice|create command queue failed|GPUAdapter/i.test(message);
  }
  function localAIUnavailableGPUFailure(reason) {
    var message=reason&&reason.message?reason.message:String(reason||"");
    return /unable to find a compatible GPU|no (?:usable|compatible) WebGPU adapter|WebGPU is unavailable|browser does not expose WebGPU/i.test(message);
  }
  function localAIGPUFailure(reason) {
    return localAIRecoverableGPUFailure(reason)||localAIUnavailableGPUFailure(reason);
  }
  function localAIFriendlyGPUFailure(reason) {
    var detail=reason&&reason.message?reason.message:String(reason||"WebGPU device failure");
    if(localAIUnavailableGPUFailure(reason))return new Error("No usable WebGPU adapter is available in this browser session. LibeRation remains fully usable with local AI switched off. If this followed a graphics-driver reset, fully exit every browser window and reopen the browser; if it persists, enable browser hardware acceleration and update the graphics driver. Technical detail: "+detail);
    return new Error("The browser's WebGPU device was reset and LibeRation could not recover it automatically. Close other GPU-heavy applications or tabs and try again; if it repeats, reload LibeRation or restart the browser. Technical detail: "+detail);
  }
  function localAICancelPurpose(purpose,reason) {
    var ids=Object.keys(localAI.pending).filter(function(id){return localAI.pending[id].purpose===purpose;});
    if(!ids.length)return;
    try{if(localAI.worker)localAI.worker.postMessage({type:"interrupt"});}catch(error){}
    localAI.cooldownUntil=Date.now()+600;
    ids.forEach(function(id){var pending=localAI.pending[id];delete localAI.pending[id];pending.reject(reason instanceof Error?reason:new Error("Local AI generation was stopped."));});
    localAISetStatus({stage:"ready",text:(purpose==="report"?"Report":"Help")+" model ready - local inference only",purpose:purpose});
  }
  function localAIShutdown(reason) {
    var worker=localAI.worker;
    if(worker){
      try{worker.postMessage({type:"dispose"});}catch(error){}
      setTimeout(function(){try{worker.terminate();}catch(error){}},750);
    }
    localAIRejectPending(reason instanceof Error?reason:new Error("Local AI was stopped."));
    localAI.worker=null;localAI.pending={};localAI.model="";localAI.contextWindow=0;localAI.purpose="";localAI.cooldownUntil=Date.now()+800;localAISetStatus({stage:"idle",text:"Model not loaded",progress:0,locked:false,model:"",purpose:"",budget:null});
  }
  function localAIFailPending(id,reason) {
    var pending=localAI.pending[id];if(!pending)return;
    delete localAI.pending[id];
    var gpuFailure=localAIGPUFailure(reason),failure=gpuFailure?localAIFriendlyGPUFailure(reason):(reason instanceof Error?reason:new Error(String(reason||"Local AI failed")));
    if(gpuFailure&&localAI.worker){try{localAI.worker.terminate();}catch(error){}localAI.worker=null;localAI.model="";localAI.contextWindow=0;localAI.purpose="";}
    localAISetStatus({stage:"error",text:failure.message,progress:0});pending.reject(failure);
  }
  function localAIPostPending(id,worker) {
    var pending=localAI.pending[id];if(!pending)return;
    var send=function(){
      if(!localAI.pending[id]||localAI.worker!==worker)return;
      try{localAISetStatus({stage:"generating",text:"Generating locally",purpose:pending.purpose});worker.postMessage(pending.request);}
      catch(error){localAIFailPending(id,error);}
    };
    var delay=Math.max(0,localAI.cooldownUntil-Date.now());
    if(delay)setTimeout(send,delay);else send();
  }
  function localAIRetryPending(id,reason,worker) {
    var pending=localAI.pending[id];
    if(!pending||pending.emitted||pending.retries>=2||!localAIRecoverableGPUFailure(reason))return false;
    var current=Number(pending.request.context_window_size)||4096,next=current;
    if(current>12288)next=12288;else if(current>8192)next=8192;else if(current>6144)next=6144;else if(current>4096)next=4096;
    pending.retries+=1;
    if(next<current){
      var rebudgeted=localAIBudgetMessages(pending.originalMessages,pending.request.model,pending.requestedMaxTokens,next);
      pending.request.messages=rebudgeted.messages;pending.request.max_tokens=rebudgeted.max_tokens;pending.request.context_window_size=next;pending.request.input_char_budget=rebudgeted.input_char_budget;
      pending.budget=rebudgeted;
    }
    localAISetStatus({stage:"recovering",text:"Restarting the local WebGPU session"+(next<current?" with a smaller "+Math.round(next/1024)+"K context":" after a graphics-device reset"),progress:0,purpose:pending.purpose,budget:pending.budget});
    try{if(worker)worker.terminate();}catch(error){}
    if(localAI.worker===worker){localAI.worker=null;localAI.model="";localAI.contextWindow=0;localAI.purpose="";}
    setTimeout(function(){
      if(!localAI.pending[id])return;
      try{var replacement=localAIWorker(pending.ai,pending.request.model,pending.purpose,pending.request.context_window_size);localAIPostPending(id,replacement);}
      catch(error){localAIFailPending(id,error);}
    },800);
    return true;
  }
  function localAIWorker(ai,model,purpose,contextWindow) {
    contextWindow=Math.max(1024,Number(contextWindow)||4096);
    if (localAI.worker && localAI.workerUrl===ai.worker_url && localAI.model===model && localAI.contextWindow===contextWindow) {
      localAI.purpose=purpose;
      if(localAI.status.stage==="ready")localAISetStatus({text:(purpose==="report"?"Report":"Help")+" model ready - local inference only",purpose:purpose});
      return localAI.worker;
    }
    if (localAI.worker && Object.keys(localAI.pending).length) throw new Error("Wait for the current local AI response before switching models.");
    if (localAI.worker) localAIShutdown();
    if (!ai.worker_url) throw new Error("The LibeRation AI worker is unavailable in this installation.");
    var worker=new Worker(ai.worker_url,{type:"module"});localAI.worker=worker;localAI.workerUrl=ai.worker_url;localAI.model=model;localAI.contextWindow=contextWindow;localAI.purpose=purpose;
    localAISetStatus({stage:"starting",text:"Starting "+(purpose==="report"?"Report":"Help")+" model with "+(contextWindow/1024).toFixed(contextWindow%1024?1:0)+"K context",progress:0,locked:false,model:model,purpose:purpose});
    worker.onmessage=function(event){if(localAI.worker!==worker)return;var message=event.data||{},pending=message.id&&localAI.pending[message.id];
      if(message.type==="progress")localAISetStatus({stage:"loading",text:value(message.text,"Downloading model"),progress:Number(message.progress)||0});
      else if(message.type==="status")localAISetStatus({stage:value(message.stage,"loading"),text:value(message.text,"Preparing model")});
      else if(message.type==="network_locked")localAISetStatus({locked:true});
      else if(message.type==="context_fallback"){
        localAI.contextWindow=Number(message.context_window_size)||localAI.contextWindow;
        if(pending){pending.request.context_window_size=localAI.contextWindow;pending.budget.context_window_size=localAI.contextWindow;pending.budget.compacted=true;}
        localAISetStatus({stage:"loading",text:"Using "+(localAI.contextWindow/1024).toFixed(localAI.contextWindow%1024?1:0)+"K context after a memory fallback",budget:pending&&pending.budget});
      }
      else if(message.type==="ready")localAISetStatus({stage:"ready",text:(localAI.purpose==="report"?"Report":"Help")+" model ready - local inference only",progress:1,locked:!!message.network_locked,model:localAI.model,purpose:localAI.purpose});
      else if(message.type==="token"&&pending){pending.emitted=true;pending.onToken(message.token||"");}
      else if(message.type==="complete"&&pending){delete localAI.pending[message.id];localAISetStatus({stage:"ready",text:(pending.purpose==="report"?"Report":"Help")+" model ready - local inference only",progress:1,purpose:pending.purpose});pending.resolve(message.text||"");}
      else if(message.type==="error"&&pending){var reason=new Error(message.message||"Local AI failed");if(!localAIRetryPending(message.id,reason,worker))localAIFailPending(message.id,reason);}
    };
    function workerFailure(event){
      if(localAI.worker!==worker)return;
      var reason=new Error(event&&event.message?event.message:"The local AI worker stopped unexpectedly.");
      var pendingId=Object.keys(localAI.pending)[0];
      if(pendingId&&localAIRetryPending(pendingId,reason,worker))return;
      try{worker.terminate();}catch(error){}
      var failure=localAIGPUFailure(reason)?localAIFriendlyGPUFailure(reason):reason;
      localAI.worker=null;localAI.model="";localAI.contextWindow=0;localAI.purpose="";localAIRejectPending(failure);
      localAISetStatus({stage:"error",text:failure.message,progress:0,locked:false,model:"",purpose:""});
    }
    worker.onerror=workerFailure;
    worker.onmessageerror=workerFailure;
    return worker;
  }
  function localAIGenerate(ai,messages,onToken,options){
    return new Promise(function(resolve,reject){
      if(!ai||!ai.activated){reject(new Error("Activate AI before requesting local generation."));return;}
      var purpose=options&&options.purpose==="report"?"report":"help",model=options&&options.model;
      if(!model)model=purpose==="report"?ai.report_model:ai.help_model;
      if(model==="same_as_help"||!model)model=ai.help_model||ai.model;
      if(!model){reject(new Error("No local AI model is configured for this function."));return;}
      var contextWindow=localAIContextWindow(ai,model,purpose,options&&options.context_window_size);
      if(localAI.worker&&localAI.model===model&&localAI.contextWindow&&localAI.contextWindow<contextWindow)contextWindow=localAI.contextWindow;
      var budgeted=localAIBudgetMessages(messages,model,options&&options.max_tokens,contextWindow);
      var worker;try{worker=localAIWorker(ai,model,purpose,contextWindow);}catch(error){reject(error);return;}
      var id="ai-"+Date.now()+"-"+Math.random().toString(16).slice(2),request={type:"generate",id:id,model:model,messages:budgeted.messages,temperature:options&&options.temperature,top_p:options&&options.top_p,max_tokens:budgeted.max_tokens,context_window_size:budgeted.context_window_size,input_char_budget:budgeted.input_char_budget};
      localAI.pending[id]={resolve:resolve,reject:reject,onToken:onToken||function(){},purpose:purpose,request:request,ai:ai,retries:0,emitted:false,originalMessages:messages,requestedMaxTokens:options&&options.max_tokens,budget:budgeted};
      localAISetStatus({budget:budgeted});
      localAIPostPending(id,worker);
    });
  }
  function useLocalAIStatus(){
    var state=React.useState(Object.assign({},localAI.status));
    React.useEffect(function(){var listener=function(next){state[1](next);};localAI.listeners.push(listener);return function(){localAI.listeners=localAI.listeners.filter(function(item){return item!==listener;});};},[]);
    return state[0];
  }

  function escapeCode(value) {
    return String(value || "").replace(/&/g, "&amp;").replace(/</g, "&lt;").replace(/>/g, "&gt;");
  }
  function highlightCode(source) {
    var text = String(source || ""), output = "", cursor = 0;
    var pattern = /(#.*$|\/\/.*$|"(?:\\.|[^"\\])*"|'(?:\\.|[^'\\])*'|\b(?:THETA|ETA|OMEGA|SIGMA|ERR|EPS)(?:_\d+|\s*\(\s*\d+(?:\s*,\s*\d+)?\s*\))|\$[A-Z][A-Z0-9_]*|\b(?:if|else|for|while|return|TRUE|FALSE|NA|NULL)\b|\b\d+(?:\.\d+)?(?:[eE][+-]?\d+)?\b|\b(?:exp|log|sqrt|sin|cos|tan|tanh|abs|expm1|log1p|ifelse|min|max|pow)\b(?=\s*\()|\b[A-Za-z_][A-Za-z0-9_.]*\b(?=\s*=)|\b(?:A|DADT)\s*\(\s*\d+\s*\)|\b(?:S\d+|F|Y|LIK|LOGLIK|IPRED|PRED|TIME|T|AMT|RATE|CMT|EVID|MDV|II|SS|DV|PREV_DV|PREV_TIME|DT|FIRST|DVID|LLOQ|BLQ|CENS|MIXNUM)\b)/gm;
    text.replace(pattern, function (token, match, offset) {
      output += escapeCode(text.slice(cursor, offset));
      var kind = /^#|^\/\//.test(token) ? "comment" :
        /^['\"]/.test(token) ? "string" :
        /^(?:THETA|ETA|OMEGA|SIGMA|ERR|EPS)(?:_\d+|\s*\()/.test(token) ? "parameter" :
        /^\$/.test(token) ? "block" :
        /^(?:if|else|for|while|return|TRUE|FALSE|NA|NULL)$/.test(token) ? "keyword" :
        /^\d/.test(token) ? "number" :
        /^(?:exp|log|sqrt|sin|cos|tan|tanh|abs|expm1|log1p|ifelse|min|max|pow)$/.test(token) ? "function" :
        /^\s*=/.test(text.slice(offset + token.length)) ? "definition" : "symbol";
      output += '<span class="lw-syntax-' + kind + '">' + escapeCode(token) + "</span>";
      cursor = offset + token.length;
      return token;
    });
    output += escapeCode(text.slice(cursor));
    return output + (text.slice(-1) === "\n" ? " " : "");
  }
  function CodeEditor(props) {
    var overlay = React.useRef(null);
    function scroll(event) {
      if (!overlay.current) return;
      overlay.current.scrollTop = event.target.scrollTop;
      overlay.current.scrollLeft = event.target.scrollLeft;
    }
    function keyDown(event) {
      if (event.key !== "Tab") return;
      event.preventDefault();
      var input = event.target, start = input.selectionStart, end = input.selectionEnd;
      var next = props.value.slice(0, start) + "  " + props.value.slice(end);
      props.onValue(next);
      window.setTimeout(function () { input.selectionStart = input.selectionEnd = start + 2; }, 0);
    }
    return e("div", { className: "lw-code-editor" },
      e("pre", { ref: overlay, "aria-hidden": "true", dangerouslySetInnerHTML: { __html: highlightCode(props.value) } }),
      e("textarea", { value: props.value, spellCheck: false, "aria-label": props.label,
        onScroll: scroll, onKeyDown: keyDown, onChange: function (event) { props.onValue(event.target.value); } }));
  }

  function StatusDot(props) { return e("span", { className: "lw-status-dot lw-status-" + value(props.status, "ready"), "aria-hidden": "true" }); }
  function Button(props) {
    return e("button", {
      type: "button", className: "lw-button " + value(props.className, ""), disabled: !!props.disabled,
      title: props.title, "aria-label": props.ariaLabel || props.title, onClick: props.onClick
    }, props.icon ? e("span", { className: "lw-button-icon", "aria-hidden": "true" }, props.icon) : null, props.children);
  }
  function Empty(props) {
    return e("div", { className: "lw-empty" }, e("strong", null, value(props.title, "Nothing to display")), e("span", null, value(props.detail, "")));
  }
  function Panel(props) {
    return e("section", { className: "lw-panel " + value(props.className, "") },
      props.title ? e("header", { className: "lw-panel-header" },
        e("div", null, e("strong", null, props.title), props.subtitle ? e("span", null, props.subtitle) : null), props.actions || null) : null,
      e("div", { className: "lw-panel-body " + value(props.bodyClass, "") }, props.children));
  }
  function Tabs(props) {
    return e("div", { className: "lw-tabs " + value(props.className, "") }, list(props.items).map(function (item) {
      return e("button", { type: "button", key: item.id, className: props.value === item.id ? "active" : "", onClick: function () { props.onChange(item.id); } }, item.label);
    }));
  }
  function Field(props) {
    return e("label", { className: "lw-field " + value(props.className, "") }, props.label ? e("span", null, props.label) : null, props.children);
  }
  function SimpleTable(props) {
    var rows = list(props.rows), columns = list(props.columns);
    if (!columns.length && rows.length) columns = Object.keys(rows[0]);
    if (!rows.length) return e(Empty, { title: value(props.empty, "No records"), detail: value(props.detail, "Nothing to display yet.") });
    return e("div", { className: "lw-table-wrap " + value(props.className, "") }, e("table", { className: "lw-table" },
      e("thead", null, e("tr", null, columns.map(function (column) { return e("th", { key: column }, column); }))),
      e("tbody", null, rows.map(function (row, index) {
        return e("tr", { key: value(row.id, index), className: props.selected === row.id ? "selected" : "", onClick: props.onRow ? function () { props.onRow(row); } : null },
          columns.map(function (column) {
            var cell = row[column];
            return e("td", { key: column, title: String(value(cell, "")) }, typeof cell === "number" ? formatNumber(cell) : value(cell, ""));
          }));
      }))));
  }
  function Modal(props) {
    var dialog = useDialogFocus(!!props.open, props.onClose);
    if (!props.open) return null;
    return e("div", { className: "lw-modal-backdrop", onMouseDown: function (event) { if (event.target === event.currentTarget) props.onClose(); } },
      e("div", { ref: dialog, tabIndex: -1, className: "lw-modal " + value(props.className, ""), role: "dialog", "aria-modal": "true", "aria-label": props.title },
        e("header", { className:"lw-modal-header" }, e("strong", null, props.title), e("button", { className:"lw-modal-close", type: "button", "aria-label":"Close", onClick: props.onClose }, "×")),
        e("div", { className: "lw-modal-body" }, props.children),
        props.footer ? e("footer", { className:"lw-modal-footer" }, props.footer) : null));
  }

  function extent(rows, key) {
    var values = rows.map(function (row) { return number(row[key]); }).filter(function (x) { return x !== null; });
    if (!values.length) return [0, 1];
    var lo = Math.min.apply(null, values), hi = Math.max.apply(null, values);
    if (lo === hi) { lo -= 0.5; hi += 0.5; }
    return [lo, hi];
  }
  function ScatterPlot(props) {
    var rows = list(props.rows).filter(function (row) { return number(row[props.x]) !== null && number(row[props.y]) !== null; });
    if (!rows.length) return e(Empty, { title: "No plottable records", detail: "Choose numeric axes or run the required diagnostic." });
    var width = value(props.width, 560), height = value(props.height, 250), margin = { l: 48, r: 16, t: 28, b: 40 };
    var overlayRows = list(props.overlayRows).filter(function (row) { return number(row[value(props.overlayX, props.x)]) !== null && number(row[value(props.overlayY, props.y)]) !== null; });
    var rangeRows = rows.concat(overlayRows.map(function (row) {
      var copy = {}; copy[props.x] = row[value(props.overlayX, props.x)]; copy[props.y] = row[value(props.overlayY, props.y)]; return copy;
    }));
    var xr = extent(rangeRows, props.x), yr = props.yRange ? props.yRange.slice() : extent(rangeRows, props.y);
    if (props.intervals) {
      var intervalValues = [];
      rows.forEach(function (row) { var lo=number(row.LOWER), hi=number(row.UPPER); if(lo!==null)intervalValues.push(lo);if(hi!==null)intervalValues.push(hi); });
      if (intervalValues.length) { yr=[Math.min(yr[0],Math.min.apply(null,intervalValues)),Math.max(yr[1],Math.max.apply(null,intervalValues))]; if(yr[0]===yr[1]){yr[0]-=0.5;yr[1]+=0.5;} }
    }
    var sx = function (x) { return margin.l + (number(x) - xr[0]) / (xr[1] - xr[0]) * (width - margin.l - margin.r); };
    var sy = function (y) { return height - margin.b - (number(y) - yr[0]) / (yr[1] - yr[0]) * (height - margin.t - margin.b); };
    var groupValues = [], groupIndex = {};
    rows.forEach(function (row) { var key = props.group ? String(value(row[props.group], "")) : "all"; if (groupIndex[key] === undefined) { groupIndex[key] = groupValues.length; groupValues.push(key); } });
    var points = rows.map(function (row, index) {
      var key = props.group ? String(value(row[props.group], "")) : "all";
      return { x: sx(row[props.x]), y: sy(row[props.y]), row: row, index: index, color: value(props.pointColor, palette[groupIndex[key] % palette.length]), group: key };
    });
    var overlayPoints = overlayRows.map(function (row, index) {
      var overlayX = value(props.overlayX, props.x), overlayY = value(props.overlayY, props.y), group = props.group ? String(value(row[props.group], "")) : "all";
      var phase = ((index * 37) % 101) / 100 - 0.5, scatter = Math.max(0, number(props.overlayScatter) || 0);
      return { x: sx(row[overlayX]) + phase * scatter * 18, y: sy(row[overlayY]), row: row, index: index, group: group, color: value(props.overlayColor, palette[value(groupIndex[group], 0) % palette.length]) };
    });
    var lines = [];
    if (number(props.intervalShade) !== null && number(props.intervalShade) > 0) {
      groupValues.forEach(function (group) {
        var sequence = points.filter(function (point) { return point.group === group && number(point.row.LOWER) !== null && number(point.row.UPPER) !== null; }).sort(function (a, b) { return number(a.row[props.x]) - number(b.row[props.x]); });
        if (sequence.length > 1) {
          var upper = sequence.map(function (point) { return point.x + "," + sy(point.row.UPPER); });
          var lower = sequence.slice().reverse().map(function (point) { return point.x + "," + sy(point.row.LOWER); });
          lines.push(e("polygon", { key: "band-" + group, points: upper.concat(lower).join(" "), fill: sequence[0].color, opacity: Math.max(0.05, Math.min(0.7, number(props.intervalShade))) }));
        }
      });
    }
    if (props.intervals) points.forEach(function (point) { var lo=number(point.row.LOWER),hi=number(point.row.UPPER); if(lo!==null&&hi!==null){lines.push(e("line",{key:"interval-"+point.index,x1:point.x,y1:sy(lo),x2:point.x,y2:sy(hi),stroke:point.color,strokeWidth:1.4,opacity:0.7}));lines.push(e("line",{key:"caplo-"+point.index,x1:point.x-3,y1:sy(lo),x2:point.x+3,y2:sy(lo),stroke:point.color}));lines.push(e("line",{key:"caphi-"+point.index,x1:point.x-3,y1:sy(hi),x2:point.x+3,y2:sy(hi),stroke:point.color}));} });
    if (props.lines) {
      var lineGroups = {}, lineKey = value(props.lineGroup, props.group);
      points.forEach(function (point) { var key = lineKey ? String(value(point.row[lineKey], "(missing)")) : point.group; if (!lineGroups[key]) lineGroups[key] = []; lineGroups[key].push(point); });
      Object.keys(lineGroups).forEach(function (group) {
        var sequence = lineGroups[group].sort(function (a, b) { return number(a.row[props.x]) - number(b.row[props.x]); });
        if (sequence.length > 1) lines.push(e("polyline", { key: "line-" + group, points: sequence.map(function (point) { return point.x + "," + point.y; }).join(" "), fill: "none", stroke: sequence[0].color, strokeWidth: 1.3, opacity: 0.7 }));
      });
    }
    if (props.overlayLines && overlayPoints.length > 1) {
      var overlayGroups = {};
      overlayPoints.forEach(function (point) { if (!overlayGroups[point.group]) overlayGroups[point.group] = []; overlayGroups[point.group].push(point); });
      Object.keys(overlayGroups).forEach(function (group) {
        var sequence = overlayGroups[group].sort(function (a,b) { return a.x-b.x; });
        lines.push(e("polyline", { key: "overlay-line-" + group, points: sequence.map(function (point) { return point.x + "," + point.y; }).join(" "), fill: "none", stroke: sequence[0].color, strokeWidth: 1.7, opacity: 0.92 }));
      });
    }
    if (props.identity) {
      var low = Math.max(xr[0], yr[0]), high = Math.min(xr[1], yr[1]);
      if (low < high) lines.push(e("line", { key: "identity", x1: sx(low), y1: sy(low), x2: sx(high), y2: sy(high), className: "lw-reference" }));
    }
    if (props.zero && yr[0] <= 0 && yr[1] >= 0) lines.push(e("line", { key: "zero", x1: margin.l, y1: sy(0), x2: width - margin.r, y2: sy(0), className: "lw-reference" }));
    list(props.referenceY).forEach(function(reference,index){var y=number(reference);if(y!==null&&y>=yr[0]&&y<=yr[1])lines.push(e("line",{key:"reference-y-"+index,x1:margin.l,y1:sy(y),x2:width-margin.r,y2:sy(y),stroke:value(props.referenceColor,"#b5484d"),strokeWidth:1.2,strokeDasharray:"3 4"}));});
    return e("svg", { className: "lw-chart", viewBox: "0 0 " + width + " " + height, role: "img", "aria-label": value(props.title, props.y + " versus " + props.x) },
      e("text", { x: margin.l, y: 17, className: "lw-chart-title" }, value(props.title, props.y + " vs " + props.x)),
      e("line", { x1: margin.l, y1: height - margin.b, x2: width - margin.r, y2: height - margin.b, className: "lw-axis" }),
      e("line", { x1: margin.l, y1: margin.t, x2: margin.l, y2: height - margin.b, className: "lw-axis" }),
      lines,
      overlayPoints.map(function (point) {
        return e("circle", { key: "overlay-" + point.index, cx: point.x, cy: point.y, r: Math.max(1.2, value(props.pointSize, 2.8) * 0.72), fill: point.color, opacity: props.overlayLines ? 0.9 : 0.32 });
      }),
      props.hidePoints ? null : points.map(function (point) {
        var radius = value(props.pointSize, 2.8), shape = String(value(props.pointShape, "16")), common = { key: "p" + point.index, opacity: 0.72 };
        if (shape === "1") return e("circle", Object.assign(common, { cx: point.x, cy: point.y, r: radius, fill: "none", stroke: point.color, strokeWidth: 1.2 }));
        if (shape === "15") return e("rect", Object.assign(common, { x: point.x - radius, y: point.y - radius, width: radius * 2, height: radius * 2, fill: point.color }));
        if (shape === "17") return e("polygon", Object.assign(common, { points: point.x + "," + (point.y-radius) + " " + (point.x-radius) + "," + (point.y+radius) + " " + (point.x+radius) + "," + (point.y+radius), fill: point.color }));
        if (shape === "18") return e("polygon", Object.assign(common, { points: point.x + "," + (point.y-radius) + " " + (point.x-radius) + "," + point.y + " " + point.x + "," + (point.y+radius) + " " + (point.x+radius) + "," + point.y, fill: point.color }));
        return e("circle", Object.assign(common, { cx: point.x, cy: point.y, r: radius, fill: point.color }));
      }),
      e("text", { x: (margin.l + width - margin.r) / 2, y: height - 8, textAnchor: "middle", className: "lw-axis-label" }, value(props.xLabel, props.x)),
      e("text", { x: 13, y: (margin.t + height - margin.b) / 2, textAnchor: "middle", transform: "rotate(-90 13 " + ((margin.t + height - margin.b) / 2) + ")", className: "lw-axis-label" }, value(props.yLabel, props.y)),
      e("text", { x: margin.l, y: height - margin.b + 14, className: "lw-tick" }, value(props.xTickStart, formatNumber(xr[0]))),
      e("text", { x: width - margin.r, y: height - margin.b + 14, textAnchor: "end", className: "lw-tick" }, value(props.xTickEnd, formatNumber(xr[1]))),
      e("text", { x: margin.l - 5, y: height - margin.b, textAnchor: "end", className: "lw-tick" }, formatNumber(yr[0])),
      e("text", { x: margin.l - 5, y: margin.t + 3, textAnchor: "end", className: "lw-tick" }, formatNumber(yr[1])));
  }
  function Histogram(props) {
    var values = list(props.values).map(number).filter(function (x) { return x !== null; });
    if (!values.length) return e(Empty, { title: "No values", detail: "Run this diagnostic first." });
    var bins = value(props.bins, 10), lo = Math.min.apply(null, values), hi = Math.max.apply(null, values);
    if (props.range) { lo = props.range[0]; hi = props.range[1]; }
    if (lo === hi) hi = lo + 1;
    var counts = Array(bins).fill(0);
    values.forEach(function (x) { var i = Math.min(bins - 1, Math.max(0, Math.floor((x - lo) / (hi - lo) * bins))); counts[i] += 1; });
    var rows = counts.map(function (count, i) { return { BIN: lo + (i + 0.5) * (hi - lo) / bins, COUNT: count }; });
    return e(ScatterPlot, { rows: rows, x: "BIN", y: "COUNT", title: props.title, lines: true, pointSize: 3.5, width: props.width, height: props.height });
  }
  function DistributionPlot(props) {
    var rows=list(props.rows).filter(function(row){return number(row.X)!==null&&number(row.Y)!==null;});
    if(!rows.length)return e(Empty,{title:"No plottable records",detail:"Choose numeric axes and a dataset."});
    var width=value(props.width,560),height=value(props.height,250),margin={l:48,r:16,t:28,b:40},xr=extent(rows,"X"), yvalues=[];
    rows.forEach(function(row){[row.LOWER,row.Q1,row.Y,row.Q3,row.UPPER].forEach(function(item){var n=number(item);if(n!==null)yvalues.push(n);});});
    var yr=yvalues.length?[Math.min.apply(null,yvalues),Math.max.apply(null,yvalues)]:[0,1];if(yr[0]===yr[1]){yr[0]-=0.5;yr[1]+=0.5;}
    var sx=function(x){return margin.l+(number(x)-xr[0])/(xr[1]-xr[0])*(width-margin.l-margin.r);},sy=function(y){return height-margin.b-(number(y)-yr[0])/(yr[1]-yr[0])*(height-margin.t-margin.b);};
    var groupIndex={},groupCount=0;
    rows.forEach(function(row){var key=props.group?String(value(row[props.group],"")):"all";if(groupIndex[key]===undefined)groupIndex[key]=groupCount++;});
    return e("svg",{className:"lw-chart",viewBox:"0 0 "+width+" "+height,role:"img","aria-label":value(props.title,"Distribution plot")},
      e("text",{x:margin.l,y:17,className:"lw-chart-title"},props.title),e("line",{x1:margin.l,y1:height-margin.b,x2:width-margin.r,y2:height-margin.b,className:"lw-axis"}),e("line",{x1:margin.l,y1:margin.t,x2:margin.l,y2:height-margin.b,className:"lw-axis"}),
      rows.map(function(row,index){var x=sx(row.X),color=palette[groupIndex[props.group?String(value(row[props.group],"")):"all"]%palette.length],w=Math.max(5,Math.min(14,(width-margin.l-margin.r)/(rows.length+2)*0.35));
        if(props.kind==="violin")return e("g",{key:index},e("polygon",{points:x+","+sy(row.LOWER)+" "+(x-w*0.55)+","+sy(row.Q1)+" "+(x-w)+","+sy(row.Y)+" "+(x-w*0.55)+","+sy(row.Q3)+" "+x+","+sy(row.UPPER)+" "+(x+w*0.55)+","+sy(row.Q3)+" "+(x+w)+","+sy(row.Y)+" "+(x+w*0.55)+","+sy(row.Q1),fill:color,opacity:0.42,stroke:color}),e("line",{x1:x-w,y1:sy(row.Y),x2:x+w,y2:sy(row.Y),stroke:color,strokeWidth:1.5}));
        return e("g",{key:index},e("line",{x1:x,y1:sy(row.LOWER),x2:x,y2:sy(row.UPPER),stroke:color}),e("rect",{x:x-w,y:sy(row.Q3),width:w*2,height:Math.max(1,sy(row.Q1)-sy(row.Q3)),fill:color,fillOpacity:0.28,stroke:color}),e("line",{x1:x-w,y1:sy(row.Y),x2:x+w,y2:sy(row.Y),stroke:color,strokeWidth:1.6}));}),
      e("text",{x:(margin.l+width-margin.r)/2,y:height-8,textAnchor:"middle",className:"lw-axis-label"},value(props.xLabel,"X")),e("text",{x:13,y:(margin.t+height-margin.b)/2,textAnchor:"middle",transform:"rotate(-90 13 "+((margin.t+height-margin.b)/2)+")",className:"lw-axis-label"},value(props.yLabel,"Y")),e("text",{x:margin.l,y:height-margin.b+14,className:"lw-tick"},formatNumber(xr[0])),e("text",{x:width-margin.r,y:height-margin.b+14,textAnchor:"end",className:"lw-tick"},formatNumber(xr[1])),e("text",{x:margin.l-5,y:height-margin.b,textAnchor:"end",className:"lw-tick"},formatNumber(yr[0])),e("text",{x:margin.l-5,y:margin.t+3,textAnchor:"end",className:"lw-tick"},formatNumber(yr[1])));
  }
  function QQPlot(props) {
    var values = list(props.values).map(number).filter(function (x) { return x !== null; }).sort(function (a, b) { return a - b; });
    var rows = values.map(function (x, i) { return { EXPECTED: normalQuantile((i + 0.5) / values.length), OBSERVED: x }; });
    return e(ScatterPlot, { rows: rows, x: "EXPECTED", y: "OBSERVED", title: props.title, identity: true, width: props.width, height: props.height });
  }
  function normalQuantile(p) {
    var a = [-39.6968302866538, 220.946098424521, -275.928510446969, 138.357751867269, -30.6647980661472, 2.50662827745924];
    var b = [-54.4760987982241, 161.585836858041, -155.698979859887, 66.8013118877197, -13.2806815528857];
    var c = [-0.00778489400243029, -0.322396458041136, -2.40075827716184, -2.54973253934373, 4.37466414146497, 2.93816398269878];
    var d = [0.00778469570904146, 0.32246712907004, 2.445134137143, 3.75440866190742];
    var q, r;
    if (p < 0.02425) { q = Math.sqrt(-2 * Math.log(p)); return (((((c[0]*q+c[1])*q+c[2])*q+c[3])*q+c[4])*q+c[5]) / ((((d[0]*q+d[1])*q+d[2])*q+d[3])*q+1); }
    if (p > 0.97575) { q = Math.sqrt(-2 * Math.log(1-p)); return -(((((c[0]*q+c[1])*q+c[2])*q+c[3])*q+c[4])*q+c[5]) / ((((d[0]*q+d[1])*q+d[2])*q+d[3])*q+1); }
    q = p - 0.5; r = q*q;
    return (((((a[0]*r+a[1])*r+a[2])*r+a[3])*r+a[4])*r+a[5])*q / (((((b[0]*r+b[1])*r+b[2])*r+b[3])*r+b[4])*r+1);
  }

  function VpcCharts(props) {
    var result = props.result || {}, observed = list(result.observed), simulated = list(result.simulated), rawPoints = list(result.points);
    if (!observed.length || !simulated.length) return e(Empty, { title: "No VPC summaries", detail: "The predictive simulation did not produce plottable bins." });
    var quantiles = Object.keys(observed[0]).filter(function (key) { return /^Q[0-9.]+$/.test(key) && Object.prototype.hasOwnProperty.call(simulated[0], key + "_median"); }).sort(function(a,b){return Number(a.slice(1))-Number(b.slice(1));});
    if (!quantiles.length) return e(Empty,{title:"No VPC quantiles",detail:"The saved VPC does not contain compatible quantile summaries."});
    var width=820,height=450,margin={l:68,r:22,t:42,b:56};
    function time(row,index){var found=number(row.TIME);return found===null?index+1:found;}
    var xValues=[];observed.forEach(function(row,index){xValues.push(time(row,index));});rawPoints.forEach(function(row){var x=number(row.TIME);if(x!==null)xValues.push(x);});
    var yValues=[];quantiles.forEach(function(key){observed.forEach(function(row){var y=number(row[key]);if(y!==null)yValues.push(y);});simulated.forEach(function(row){["_lo","_median","_hi"].forEach(function(suffix){var y=number(row[key+suffix]);if(y!==null)yValues.push(y);});});});rawPoints.forEach(function(row){var y=number(row.DV);if(y!==null)yValues.push(y);});
    var xr=[Math.min.apply(null,xValues),Math.max.apply(null,xValues)], positive=yValues.length&&yValues.every(function(y){return y>0;}), logScale=positive;
    var yr=[Math.min.apply(null,yValues),Math.max.apply(null,yValues)];if(xr[0]===xr[1])xr=[xr[0]-.5,xr[1]+.5];if(yr[0]===yr[1])yr=[yr[0]*.8,yr[1]*1.2];
    if(!logScale){var padding=(yr[1]-yr[0])*.06;yr=[yr[0]-padding,yr[1]+padding];}
    var transformY=function(y){return logScale?Math.log(Number(y))/Math.LN10:Number(y);},ty=[transformY(yr[0]),transformY(yr[1])];
    var sx=function(x){return margin.l+(Number(x)-xr[0])/(xr[1]-xr[0])*(width-margin.l-margin.r);},sy=function(y){return height-margin.b-(transformY(y)-ty[0])/(ty[1]-ty[0])*(height-margin.t-margin.b);};
    function band(key,color){var sequence=simulated.map(function(row,index){return{x:time(row,index),lo:number(row[key+"_lo"]),hi:number(row[key+"_hi"])};}).filter(function(item){return item.lo!==null&&item.hi!==null&&(!logScale||(item.lo>0&&item.hi>0));}).sort(function(a,b){return a.x-b.x;});if(sequence.length<2)return null;var upper=sequence.map(function(item){return sx(item.x)+","+sy(item.hi);}),lower=sequence.slice().reverse().map(function(item){return sx(item.x)+","+sy(item.lo);});return e("polygon",{key:"band-"+key,points:upper.concat(lower).join(" "),fill:color,opacity:.58});}
    var bandColors=quantiles.map(function(_,index){return index===Math.floor(quantiles.length/2)?"#ef8588":"#80b9ee";});
    var xTicks=Array.from({length:6},function(_,i){return xr[0]+i*(xr[1]-xr[0])/5;});
    var yTicks=[];if(logScale){for(var power=Math.floor(ty[0]);power<=Math.ceil(ty[1]);power+=1){var tick=Math.pow(10,power);if(tick>=yr[0]&&tick<=yr[1])yTicks.push(tick);}}else{yTicks=Array.from({length:5},function(_,i){return yr[0]+i*(yr[1]-yr[0])/4;});}
    return e("div",{className:"lw-vpc-figure"},e("svg",{className:"lw-chart lw-vpc-chart",viewBox:"0 0 "+width+" "+height,role:"img","aria-label":"Visual predictive check"},
      e("text",{x:margin.l,y:22,className:"lw-chart-title"},result.pc_correct?"Prediction-corrected visual predictive check":"Visual predictive check"),
      e("line",{x1:margin.l,y1:height-margin.b,x2:width-margin.r,y2:height-margin.b,className:"lw-axis"}),e("line",{x1:margin.l,y1:margin.t,x2:margin.l,y2:height-margin.b,className:"lw-axis"}),
      xTicks.map(function(tick,index){return e("g",{key:"xt"+index},e("line",{x1:sx(tick),y1:height-margin.b,x2:sx(tick),y2:height-margin.b+5,className:"lw-axis"}),e("text",{x:sx(tick),y:height-margin.b+18,textAnchor:"middle",className:"lw-tick"},formatNumber(tick)));}),
      yTicks.map(function(tick,index){return e("g",{key:"yt"+index},e("line",{x1:margin.l-5,y1:sy(tick),x2:margin.l,y2:sy(tick),className:"lw-axis"}),e("line",{x1:margin.l,y1:sy(tick),x2:width-margin.r,y2:sy(tick),className:"lw-vpc-gridline"}),e("text",{x:margin.l-9,y:sy(tick)+3,textAnchor:"end",className:"lw-tick"},formatNumber(tick)));}),
      quantiles.map(function(key,index){return band(key,bandColors[index]);}),
      quantiles.map(function(key,index){var sequence=observed.map(function(row,i){return{x:time(row,i),y:number(row[key])};}).filter(function(item){return item.y!==null&&(!logScale||item.y>0);}).sort(function(a,b){return a.x-b.x;});return e("polyline",{key:"observed-"+key,points:sequence.map(function(item){return sx(item.x)+","+sy(item.y);}).join(" "),fill:"none",stroke:"#111827",strokeWidth:index===Math.floor(quantiles.length/2)?2:1.6,strokeDasharray:"8 6",strokeLinecap:"round"});}),
      rawPoints.filter(function(row){return number(row.TIME)!==null&&number(row.DV)!==null&&(!logScale||Number(row.DV)>0);}).map(function(row,index){return e("circle",{key:"raw"+index,cx:sx(row.TIME),cy:sy(row.DV),r:2,fill:"#111827",opacity:.76});}),
      e("text",{x:(margin.l+width-margin.r)/2,y:height-14,textAnchor:"middle",className:"lw-axis-label"},"Time"),e("text",{x:17,y:(margin.t+height-margin.b)/2,textAnchor:"middle",transform:"rotate(-90 17 "+((margin.t+height-margin.b)/2)+")",className:"lw-axis-label"},"Observed concentration / response"),
      e("g",{className:"lw-vpc-legend",transform:"translate("+(width-310)+",14)"},e("rect",{x:0,y:0,width:14,height:8,fill:"#80b9ee",opacity:.7}),e("text",{x:19,y:8},"outer quantile intervals"),e("rect",{x:139,y:0,width:14,height:8,fill:"#ef8588",opacity:.7}),e("text",{x:158,y:8},"median interval"))));
  }

  function ParameterGrid(props) {
    var rows = list(props.rows);
    return e("div", { className: "lw-parameter-block" },
      e("div",{className:"lw-parameter-heading"},e("h5", null, props.title),e("div",null,
        props.onRemove?e(Button,{className:"lw-button-quiet lw-parameter-resize",disabled:!rows.length,onClick:props.onRemove,title:"Remove the last "+value(props.unitLabel,props.prefix)+" definition"},"−"):null,
        props.onAdd?e(Button,{className:"lw-button-quiet lw-parameter-resize",onClick:props.onAdd,title:"Add "+value(props.unitLabel,props.prefix)+" definition"},"+ Add "+value(props.unitLabel,props.prefix)):null)),
      rows.length ? e("table", { className: "lw-param-table" },
      e("thead", null, e("tr", null, e("th", null, props.matrix ? "Element" : "Name"), e("th", null, "Initial"), props.bounds?e("th",null,"Lower"):null, props.bounds?e("th",null,"Upper"):null, e("th", null, "Fixed"))),
      e("tbody", null, rows.map(function (row, index) { var name = props.matrix ? "OMEGA(" + value(row.ROW, index + 1) + "," + value(row.COL, index + 1) + ")" : props.prefix + value(row[props.indexName], index + 1); return e("tr", { key: name + "-" + index },
        e("td", null, name),
        e("td", null, e("input", { type: "number", step: "any", value: value(row.Value, 0), onChange: function (event) { props.onChange(index, "Value", Number(event.target.value)); } })),
        props.bounds?e("td",null,e("input",{type:"number",step:"any",value:row.LOWER===null||row.LOWER===undefined?"":row.LOWER,onChange:function(event){props.onChange(index,"LOWER",event.target.value===""?null:Number(event.target.value));}})):null,
        props.bounds?e("td",null,e("input",{type:"number",step:"any",value:row.UPPER===null||row.UPPER===undefined?"":row.UPPER,onChange:function(event){props.onChange(index,"UPPER",event.target.value===""?null:Number(event.target.value));}})):null,
        e("td", null, e("input", { type: "checkbox", checked: !!row.FIX, onChange: function (event) { props.onChange(index, "FIX", event.target.checked); } }))); }))) : e("span", { className: "lw-muted" }, "None"));
  }

  function PriorGrid(props) {
    var rows = list(props.rows), names = list(props.parameterNames);
    function update(index, field, nextValue) { var next=cloneRows(rows);next[index][field]=nextValue;props.onChange(next); }
    function add() { props.onChange(rows.concat([{parameter:value(names[0],"THETA1"),distribution:"normal",mean:0,sd:1,shape:null,rate:null}])); }
    function remove(index) { props.onChange(rows.filter(function (_,i) { return i!==index; })); }
    return e("div", { className:"lw-prior-block" },
      e("div", { className:"lw-prior-heading" }, e("div",null,e("h5",null,"Estimation priors"),e("small",null,"Saved with the model version and applied to every estimation method.")), e(Button,{className:"lw-button-quiet",disabled:!names.length,onClick:add},"+ Add prior")),
      rows.length ? e("div",{className:"lw-prior-table-wrap"},e("table",{className:"lw-param-table lw-prior-table"},
        e("thead",null,e("tr",null,e("th",null,"Parameter"),e("th",null,"Distribution"),e("th",null,"Mean / shape"),e("th",null,"SD / rate"),e("th",null,""))),
        e("tbody",null,rows.map(function(row,index){var inverse=row.distribution==="inverse_gamma";return e("tr",{key:index},
          e("td",null,e("select",{value:row.parameter,onChange:function(event){update(index,"parameter",event.target.value);}},names.map(function(name){return e("option",{key:name,value:name},name);}))),
          e("td",null,e("select",{value:row.distribution,onChange:function(event){update(index,"distribution",event.target.value);}},["normal","lognormal","half_normal","inverse_gamma"].map(function(name){return e("option",{key:name,value:name},name.replace("_"," "));}))),
          e("td",null,e("input",{type:"number",step:"any",value:inverse?value(row.shape,2):value(row.mean,0),onChange:function(event){update(index,inverse?"shape":"mean",Number(event.target.value));}})),
          e("td",null,e("input",{type:"number",step:"any",min:0,value:inverse?value(row.rate,1):value(row.sd,1),onChange:function(event){update(index,inverse?"rate":"sd",Number(event.target.value));}})),
          e("td",null,e(Button,{className:"lw-button-link lw-prior-remove",title:"Remove prior",onClick:function(){remove(index);}},"Remove")));})))) : e("p",{className:"lw-muted lw-prior-empty"},"No priors: estimation uses the likelihood only."));
  }

  function ModelEditor(props) {
    var model = props.model || {};
    var sourceState = useSynced({ pred: value(model.pred, ""), des: value(model.des, ""), alg: value(model.alg, ""), error: value(model.error, "Y=F") }, [model.pred, model.des, model.alg, model.error]);
    var source = sourceState[0], setSource = sourceState[1];
    var parameterState = useSynced({ theta: cloneRows(model.theta), omega: cloneRows(model.omega), sigma: cloneRows(model.sigma) }, [model.theta, model.omega, model.sigma]);
    var parameters = parameterState[0], setParameters = parameterState[1];
    var priorState = useSynced(cloneRows(model.priors), [model.priors]), priors = priorState[0], setPriors = priorState[1];
    var omegaStructureState = useSynced(value(model.omega_structure,"diagonal"), [model.omega_structure]), omegaStructure = omegaStructureState[0];
    var advanState = useSynced(String(value(model.advan, 4)), [model.advan]);
    var transState = useSynced(String(value(model.trans, 2)), [model.trans]);
    var problemState = useSynced(value(model.name, "Untitled model"), [model.name]);
    var nState = useSynced(value(model.n_state, 2), [model.n_state]);
    var inputState = useSynced(list(model.input), [model.input]), selectedInput = inputState[0], setSelectedInput = inputState[1];
    var outputState = useSynced(list(model.output), [model.output]), selectedOutput = outputState[0], setSelectedOutput = outputState[1];
    var columnModal = React.useState(false), transHelp = React.useState(false);
    var dirty = source.pred !== value(model.pred, "") || source.des !== value(model.des, "") || source.alg !== value(model.alg, "") || source.error !== value(model.error, "Y=F") || advanState[0] !== String(value(model.advan, 4)) || transState[0] !== String(value(model.trans, 2)) || problemState[0] !== value(model.name, "Untitled model") || JSON.stringify(parameters) !== JSON.stringify({ theta: cloneRows(model.theta), omega: cloneRows(model.omega), sigma: cloneRows(model.sigma) }) || JSON.stringify(priors) !== JSON.stringify(cloneRows(model.priors)) || omegaStructure !== value(model.omega_structure,"diagonal") || JSON.stringify(selectedInput) !== JSON.stringify(list(model.input)) || JSON.stringify(selectedOutput) !== JSON.stringify(list(model.output));
    var validationNonce=props.result&&props.result.kind==="model_validation"?props.result.nonce:null;
    React.useEffect(function(){
      if(!validationNonce||!props.result.parameters)return;
      setParameters({theta:cloneRows(props.result.parameters.theta),omega:cloneRows(props.result.parameters.omega),sigma:cloneRows(props.result.parameters.sigma)});
      if(props.result.parameters.priors)setPriors(cloneRows(props.result.parameters.priors));
    },[validationNonce]);
    function updateParameter(kind, index, field, nextValue) {
      var next = Object.assign({}, parameters); next[kind] = cloneRows(parameters[kind]); next[kind][index][field] = nextValue; setParameters(next);
    }
    function omegaDimension(rows){var dimension=0;list(rows).forEach(function(row,index){dimension=Math.max(dimension,Number(value(row.ROW,row.OMEGA||index+1))||0,Number(value(row.COL,row.OMEGA||index+1))||0);});return dimension;}
    function omegaRowsForDimension(count,current){
      count=Math.max(0,Number(count)||0);current=cloneRows(current);
      function find(row,column){return current.filter(function(item,index){return Number(value(item.ROW,item.OMEGA||index+1))===row&&Number(value(item.COL,item.OMEGA||index+1))===column;})[0];}
      var rows=[];
      if(omegaStructure==="full"){
        for(var row=1;row<=count;row+=1)for(var column=1;column<=row;column+=1){var existing=find(row,column);rows.push({OMEGA:rows.length+1,ROW:row,COL:column,Value:existing?Number(existing.Value):(row===column?0.1:0),FIX:existing?!!existing.FIX:false});}
      }else{
        for(var diagonal=1;diagonal<=count;diagonal+=1){var item=find(diagonal,diagonal);rows.push({OMEGA:diagonal,ROW:diagonal,COL:diagonal,Value:item?Number(item.Value):0.1,FIX:item?!!item.FIX:false});}
      }
      return rows;
    }
    function simpleRowsForCount(kind,count,current){
      count=Math.max(0,Number(count)||0);var rows=cloneRows(current).slice(0,count),indexName=kind==="theta"?"THETA":"SIGMA";
      while(rows.length<count)rows.push(kind==="theta"?{THETA:rows.length+1,Value:1,LOWER:null,UPPER:null,FIX:false}:{SIGMA:rows.length+1,Value:0.1,FIX:false});
      rows.forEach(function(row,index){row[indexName]=index+1;});return rows;
    }
    function parameterNames(next){return next.theta.map(function(_,i){return "THETA"+(i+1);}).concat(next.omega.map(function(_,i){return "OMEGA"+(i+1);})).concat(next.sigma.map(function(_,i){return "SIGMA"+(i+1);}));}
    function commitParameterRows(next){var valid=parameterNames(next);setParameters(next);setPriors(cloneRows(priors).filter(function(prior){return valid.indexOf(String(prior.parameter))>=0;}));return next;}
    function resizeParameterKind(kind,delta){
      var next={theta:cloneRows(parameters.theta),omega:cloneRows(parameters.omega),sigma:cloneRows(parameters.sigma)};
      if(kind==="omega")next.omega=omegaRowsForDimension(Math.max(0,omegaDimension(next.omega)+delta),next.omega);
      else next[kind]=simpleRowsForCount(kind,Math.max(0,next[kind].length+delta),next[kind]);
      commitParameterRows(next);
    }
    function maximumCodeReference(names){
      var code=[source.pred,source.des,source.alg,source.error].join("\n").replace(/\/\*[\s\S]*?\*\//g,"").replace(/\/\/.*$/gm,"").replace(/#.*$/gm,"");
      var pattern=new RegExp("\\b(?:"+names.join("|")+")\\s*\\(\\s*([1-9][0-9]*)\\s*\\)","gi"),match,maximum=0;
      while((match=pattern.exec(code))!==null)maximum=Math.max(maximum,Number(match[1])||0);
      return maximum;
    }
    function synchronizedParameters(){
      var thetaCount=Math.max(parameters.theta.length,maximumCodeReference(["THETA"])),etaCount=Math.max(omegaDimension(parameters.omega),maximumCodeReference(["ETA"])),sigmaCount=Math.max(parameters.sigma.length,maximumCodeReference(["ERR","EPS","SIGMA"]));
      return {theta:simpleRowsForCount("theta",thetaCount,parameters.theta),omega:omegaRowsForDimension(etaCount,parameters.omega),sigma:simpleRowsForCount("sigma",sigmaCount,parameters.sigma)};
    }
    function changeOmegaStructure(nextStructure) {
      var nEta=omegaDimension(parameters.omega), current=cloneRows(parameters.omega), nextRows=[];
      function find(row,column){return current.filter(function(item,index){var r=Number(value(item.ROW,index+1)),c=Number(value(item.COL,index+1));return r===row&&c===column;})[0];}
      if(nextStructure==="full"){
        for(var row=1;row<=nEta;row+=1)for(var column=1;column<=row;column+=1){var existing=find(row,column);nextRows.push({OMEGA:nextRows.length+1,ROW:row,COL:column,Value:existing?Number(existing.Value):(row===column?0.1:0),FIX:existing?!!existing.FIX:false});}
      }else{
        for(var diagonal=1;diagonal<=nEta;diagonal+=1){var item=find(diagonal,diagonal);nextRows.push({OMEGA:diagonal,ROW:diagonal,COL:diagonal,Value:item?Number(item.Value):0.1,FIX:item?!!item.FIX:false});}
      }
      var remapped=cloneRows(priors).map(function(prior){if(!/^OMEGA\d+$/.test(prior.parameter))return prior;var old=current[Number(prior.parameter.replace("OMEGA",""))-1];if(!old)return null;var oldRow=Number(value(old.ROW,old.OMEGA)),oldCol=Number(value(old.COL,old.OMEGA));var nextIndex=nextRows.findIndex(function(item){return Number(item.ROW)===oldRow&&Number(item.COL)===oldCol;});if(nextIndex<0)return null;prior.parameter="OMEGA"+(nextIndex+1);return prior;}).filter(function(prior){return !!prior;});
      setPriors(remapped);setParameters(Object.assign({},parameters,{omega:nextRows}));omegaStructureState[1](nextStructure);
    }
    var priorParameterNames = parameterNames(parameters);
    function draftPayload(mode,parameterRows) {
      parameterRows=parameterRows||parameters;var validParameters=parameterNames(parameterRows),draftPriors=cloneRows(priors).filter(function(prior){return validParameters.indexOf(String(prior.parameter))>=0;});
      return {
        pred: source.pred, des: source.des, alg: source.alg, error: source.error, advan: Number(advanState[0]), trans: Number(transState[0]),
        n_state: Number(nState[0]), problem: problemState[0], theta: parameterRows.theta, omega: parameterRows.omega, sigma: parameterRows.sigma,
        omega_structure: omegaStructure, priors: draftPriors, input: selectedInput, output: selectedOutput, save_mode: mode
      };
    }
    function save(mode) {
      var next=synchronizedParameters();commitParameterRows(next);emit(props, "update_model", draftPayload(mode,next));
    }
    function validateDraft(){
      var next=synchronizedParameters();commitParameterRows(next);emit(props,"validate",draftPayload("validate",next));
    }
    return e("div", { className: "lw-code-workspace" },
      dirty ? e("div", { className: "lw-dirty-banner" }, "Unsaved editor changes — apply changes before saving or running the model.") : null,
      model.experimental ? e("div",{className:"lw-info-banner lw-outcome-banner"},
        e("strong",null,"Experimental engine"+(model.experimental.strict?" · strict":"")),
        e("span",null,list(model.experimental.features).join(" · ")),
        e("small",null,value(model.experimental.label,"This model records experimental solver provenance with every run."))) : null,
      model.dde ? e("div",{className:"lw-info-banner lw-outcome-banner"},e("strong",null,"Delay differential equations"),e("span",null,value(model.dde.lag_count,0)+" lag(s) · step "+value(model.dde.step,"-")+" · "+list(model.dde.delays).join(", ")),e("small",null,"Fixed-step method-of-steps with differentiable linear history interpolation.")) : null,
      model.dae ? e("div",{className:"lw-info-banner lw-outcome-banner"},e("strong",null,"Index-1 DAE"+(model.dae.sparse?" · block sparse":"")),e("span",null,list(model.dae.variables).join(" · ")),e("small",null,"Algebraic residuals are edited in $ALG and solved inside C++/CppAD.")) : null,
      model.qsp ? e("div",{className:"lw-info-banner lw-outcome-banner"},e("strong",null,"QSP reaction network"),e("span",null,list(model.qsp.species).join(" · ")),e("small",null,value(model.qsp.reactions,0)+" stoichiometric reaction(s) compiled into $DES.")) : null,
      list(model.components).length ? e("div",{className:"lw-info-banner lw-outcome-banner"},e("strong",null,"Offline hybrid components"),e("span",null,list(model.components).map(function(item){return item.name+" ["+item.type+" / "+item.scope+"]";}).join(" · ")),e("small",null,"Immutable hashed payloads run within the differentiable C++ objective.")) : null,
      list(model.outcomes).length ? e("div",{className:"lw-info-banner lw-outcome-banner"},
        e("strong",null,"First-class outcomes"),
        e("span",null,list(model.outcomes).map(function(item){return item.name+" ["+item.family.replace(/_/g," ")+"]"+(item.dvid===null||item.dvid===undefined?"":" · DVID "+item.dvid);}).join("; ")),
        e("small",null,"The editable $ERROR likelihood, stochastic outcome generator and family-specific diagnostics share this declaration.")) : null,
      model.hmm ? e("div",{className:"lw-info-banner lw-outcome-banner"},
        e("strong",null,value(model.hmm.model_type,"Hidden Markov")+" model"),
        e("span",null,list(model.hmm.states).join(" · ")),
        e("small",null,"Filtering, retrospective smoothing and Viterbi decoding use the compiled sequence likelihood.")) : null,
      model.kalman ? e("div",{className:"lw-info-banner lw-outcome-banner"},
        e("strong",null,String(value(model.kalman.filter,"linear")).toUpperCase()+" "+value(model.kalman.model_type,"state-space")+" model"),
        e("span",null,list(model.kalman.states).join(" · ")+" · "+value(model.kalman.dynamics,"discrete")+" dynamics"+(list(model.kalman.regimes).length?" · regimes: "+list(model.kalman.regimes).join(", "):"")),
        e("small",null,"State likelihood, filtering, smoothing and stochastic simulation run in the C++ engine.")) : null,
      model.random_effects ? e("div",{className:"lw-info-banner lw-outcome-banner"},
        e("strong",null,"Generalized random-effect design"),
        e("span",null,list(model.random_effects.blocks).map(function(block){return block.name+" ["+block.column+" → ETA("+list(block.etas).join(",")+")]";}).join(" · ")),
        e("small",null,"Independent objective clusters: "+value(model.random_effects.cluster,"automatic connected components"))) : null,
      e("div", { className: "lw-control-row" },
        e(Field, { label: "$PROBLEM", className: "lw-grow" }, e("input", { value: problemState[0], onChange: function (event) { problemState[1](event.target.value); } })),
        e(Field, { label: "ADVAN" }, e("select", { value: advanState[0], onChange: function (event) { advanState[1](event.target.value); } }, [1,2,3,4,6,11,12,13].map(function (x) { return e("option", { key: x, value: x }, x); }))),
        Number(advanState[0]) === 6 || Number(advanState[0]) === 13 ? e(Field, { label: "Compartments" }, e("input", { type: "number", min: 1, max: 20, value: nState[0], onChange: function (event) { nState[1](Number(event.target.value)); } })) :
          e(Field, { label: "TRANS" }, e("select", { value: transState[0], onChange: function (event) { transState[1](event.target.value); } }, [1,2,3,4,5,6].map(function (x) { return e("option", { key: x, value: x }, x); }))),
        Number(advanState[0]) === 6 || Number(advanState[0]) === 13 ? null : e(Button, { className: "lw-button-link lw-help-button", onClick: function () { transHelp[1](true); } }, "?"),
        e(Field, { label: "Dataset" }, e("select", { value: props.dataset.loaded ? "current" : "" }, e("option", { value: props.dataset.loaded ? "current" : "" }, props.dataset.loaded ? value(props.dataset.name, "Current dataset") : "No dataset"))),
        e(Button, { className: "lw-button-quiet lw-columns-button", onClick: function () { columnModal[1](true); } }, "Columns...")),
      e("div", { className: "lw-editor-grid " + (model.dae ? "lw-editor-grid-four" : ((Number(advanState[0]) === 6 || Number(advanState[0]) === 13) ? "lw-editor-grid-three" : "")) },
        e("div", { className: "lw-editor-box" }, e("h5", null, "$PK / $PRED"), e(CodeEditor,{label:"PK or PRED model code",value:source.pred,onValue:function(next){setSource(Object.assign({},source,{pred:next}));}})),
        Number(advanState[0]) === 6 || Number(advanState[0]) === 13 ? e("div", { className: "lw-editor-box" }, e("h5", null, "$DES"), e(CodeEditor,{label:"DES differential equation code",value:source.des,onValue:function(next){setSource(Object.assign({},source,{des:next}));}})) : null,
        model.dae ? e("div", { className: "lw-editor-box" }, e("h5", null, "$ALG"), e(CodeEditor,{label:"Algebraic residual equation code",value:source.alg,onValue:function(next){setSource(Object.assign({},source,{alg:next}));}})) : null,
        e("div", { className: "lw-editor-box" }, e("h5", null, "$ERROR"), e(CodeEditor,{label:"ERROR model code",value:source.error,onValue:function(next){setSource(Object.assign({},source,{error:next}));}}))),
      e("div", { className: "lw-parameter-grid" },
        e(ParameterGrid, { title: "THETA", bounds:true, prefix: "THETA", indexName: "THETA", rows: parameters.theta, unitLabel:"THETA", onAdd:function(){resizeParameterKind("theta",1);}, onRemove:function(){resizeParameterKind("theta",-1);}, onChange: function (i,f,v) { updateParameter("theta",i,f,v); } }),
        e("div",{className:"lw-omega-block"},e(ParameterGrid, { title: omegaStructure==="full"?"OMEGA lower triangle":"OMEGA", matrix:omegaStructure==="full", prefix: "OMEGA", indexName: "OMEGA", rows: parameters.omega, unitLabel:"ETA", onAdd:function(){resizeParameterKind("omega",1);}, onRemove:function(){resizeParameterKind("omega",-1);}, onChange: function (i,f,v) { updateParameter("omega",i,f,v); } }),e("label",{className:"lw-check lw-omega-matrix-toggle"},e("input",{type:"checkbox",checked:omegaStructure==="full",onChange:function(event){changeOmegaStructure(event.target.checked?"full":"diagonal");}})," OMEGA matrix")),
        e(ParameterGrid, { title: "SIGMA", prefix: "SIGMA", indexName: "SIGMA", rows: parameters.sigma, unitLabel:"SIGMA", onAdd:function(){resizeParameterKind("sigma",1);}, onRemove:function(){resizeParameterKind("sigma",-1);}, onChange: function (i,f,v) { updateParameter("sigma",i,f,v); } })),
      e(PriorGrid,{rows:priors,parameterNames:priorParameterNames,onChange:setPriors}),
      e("div", { className: "lw-inline-actions lw-editor-actions" },
        e(Button, { className: "lw-button-quiet", onClick: validateDraft }, "Validate"),
        e(Button, { className: "lw-button-primary", onClick: function () { save("current"); } }, "Apply changes")),
      e(Modal, { open: columnModal[0], onClose: function () { columnModal[1](false); }, title: "$INPUT / generated OUTPUT columns", footer: e(Button, { className: "lw-button-primary", onClick: function () { columnModal[1](false); } }, "Done") },
        e("div", { className: "lw-column-sections" },
          e("section", null,
            e("h4", null, "$INPUT dataset columns"),
            e("div", { className: "lw-column-list" }, list(props.dataset.columns).map(function (column) {
              return e("label", { key: column }, e("input", { type: "checkbox", checked: selectedInput.indexOf(column) >= 0, onChange: function (event) { var next=selectedInput.filter(function(item){return item!==column;});if(event.target.checked)next.push(column);setSelectedInput(next); } }), column);
            }))),
          e("section", null,
            e("h4", null, "Generated run columns"),
            e("p", { className: "lw-help-text" }, "Validate the current draft to refresh variables assigned in $PK/$PRED."),
            e("div", { className: "lw-column-list lw-output-column-list" }, list(model.outputs).map(function (item) {
              var name=String(item.name),disabled=item.selectable===false||item.selectable===0;
              return e("label", { key:name,className:disabled?"disabled":"",title:value(item.description,"") }, e("input", { type: "checkbox",disabled:disabled, checked:selectedOutput.indexOf(name)>=0,onChange:function(event){var next=selectedOutput.filter(function(itemName){return itemName!==name;});if(event.target.checked)next.push(name);setSelectedOutput(next);} }),e("span",null,name),e("small",null,value(item.source,"generated")+(item.availability==="estimation"?" (estimation only)":"")));
            }))))),
      e(Modal, { open: transHelp[0], onClose: function () { transHelp[1](false); }, title: "TRANS parameterization", footer: e(Button, { className: "lw-button-primary", onClick: function () { transHelp[1](false); } }, "Close") },
        e("div", { className: "lw-form-stack" }, e("p", null, "TRANS selects the NONMEM-compatible micro- or macro-parameterization used by the chosen ADVAN model."), e("p", { className: "lw-muted" }, "Common choices: TRANS2 for one-compartment CL/V models; TRANS4 for multi-compartment CL/V/Q parameterizations. The engine validates the required symbols when the model is compiled."))));
  }

  function DiagnosticsPane(props) {
    var fit = props.fit || {}, tab = props.tab, diagnostics = props.diagnostics || {};
    var result = diagnostics[tab] || {};
    if (tab === "gof" && !fit.available) return e(Empty, { title: "No fitted model", detail: "Run or open an estimation to populate diagnostics." });
    if (tab === "gof" && !fit.gof_loaded) return e(Empty, { title: "Loading GOF plots", detail: "The selected run data is loaded only for this tab and then cached." });
    if (tab !== "gof" && diagnostics.available && diagnostics.available[tab] && !diagnostics[tab]) return e(Empty, { title: "Loading saved diagnostic", detail: "The plot payload is being loaded once for this tab." });
    var gof = list(fit.gof), observed = gof.filter(function (row) { return number(row.DV) !== null; });
    if (tab === "gof") return e("div", { className: "lw-diagnostic-grid" },
      e(ScatterPlot, { rows: observed, x: "PRED", y: "DV", identity: true, title: "DV vs PRED" }),
      e(ScatterPlot, { rows: observed, x: "IPRED", y: "DV", identity: true, title: "DV vs IPRED" }),
      e(ScatterPlot, { rows: observed, x: "TIME", y: "CWRES", zero: true, yRange:[-5,5], referenceY:[-2,2], pointColor:"#b5484d", referenceColor:"#c45b61", title: "CWRES vs time" }),
      e(ScatterPlot, { rows: observed, x: "PRED", y: "CWRES", zero: true, yRange:[-5,5], referenceY:[-2,2], pointColor:"#b5484d", referenceColor:"#c45b61", title: "CWRES vs PRED" }),
      e(QQPlot, { values: observed.map(function (row) { return row.CWRES; }), title: "Normal Q-Q of CWRES" }));
    if (tab === "npc") {
      return e("div", { className: "lw-diagnostic-grid" }, e(Histogram, { values: list(result.table).map(function (row) { return row.PERCENTILE; }), range: [0,1], title: "Predictive percentiles" }), e(ScatterPlot, { rows: result.table, x: "TIME", y: "PERCENTILE", group: "ID", title: "Predictive percentile vs time" }));
    }
    if (tab === "npde") {
      return e("div", { className: "lw-diagnostic-grid" }, e(QQPlot, { values: list(result.table).map(function (row) { return row.NPDE; }), title: "Normal Q-Q of NPDE" }), e(ScatterPlot, { rows: result.table, x: "TIME", y: "NPDE", group: "ID", zero: true, title: "NPDE vs time" }));
    }
    if (tab === "vpc") {
      var strata = list(result.stratified);
      return e("div", null,
        e("div", { className: "lw-vpc-options" }, e("span", null, result.pc_correct ? "Prediction-corrected VPC" : "Standard VPC"), e("span", null, result.nsim + " simulations (saved)")),
        strata.length ? e("h4", { className: "lw-vpc-section-title" }, "Overall population") : null,
        e(VpcCharts, { result:result }),
        strata.length ? e("section", { className: "lw-vpc-stratified" },
          e("h4", { className: "lw-vpc-section-title" }, "Stratified by " + value(result.stratify, "variable")),
          e("div", { className: "lw-vpc-strata-grid" }, strata.map(function (item, index) {
            return e("article", { className: "lw-vpc-stratum", key: value(item.level, index) },
              e("h5", null, value(result.stratify, "Stratum") + " = " + value(item.level, "(missing)")),
              e(VpcCharts, { result:item }));
          }))) : null,
        e("details", { className:"lw-vpc-table" }, e("summary", null, "Overall simulation interval table"), e(SimpleTable, { rows: result.simulated })));
    }
    if (tab === "vpc_categorical") {
      var categoricalRows=list(result.simulated).map(function(row){return {TIME:row.TIME,Y:row.median,LOWER:row.lower,UPPER:row.upper,CATEGORY:row.CATEGORY};});
      return e("div",{className:"lw-diagnostic-grid"},e(ScatterPlot,{rows:categoricalRows,x:"TIME",y:"Y",group:"CATEGORY",lineGroup:"CATEGORY",lines:true,intervals:true,intervalShade:0.18,overlayRows:result.observed,overlayX:"TIME",overlayY:"PROPORTION",overlayColor:"#17202a",title:"Categorical VPC: observed proportions",xLabel:"Time",yLabel:"Proportion"}),e(SimpleTable,{rows:result.simulated}));
    }
    if (tab === "vpc_count") {
      var countRows=list(result.simulated).map(function(row){return {TIME:row.TIME,Y:row.MEAN_median,LOWER:row.MEAN_lower,UPPER:row.MEAN_upper};});
      return e("div",{className:"lw-diagnostic-grid"},e(ScatterPlot,{rows:countRows,x:"TIME",y:"Y",lines:true,intervals:true,intervalShade:0.25,overlayRows:result.observed,overlayX:"TIME",overlayY:"MEAN",overlayLines:true,overlayColor:"#17202a",title:"Count VPC: mean response",xLabel:"Time",yLabel:"Mean count"}),e(SimpleTable,{rows:result.simulated}));
    }
    if (tab === "vpc_tte") {
      var tteRows=list(result.simulated).map(function(row){return {TIME:row.TIME,Y:row.median,LOWER:row.lower,UPPER:row.upper};});
      return e("div",{className:"lw-diagnostic-grid"},e(ScatterPlot,{rows:tteRows,x:"TIME",y:"Y",lines:true,intervalShade:0.28,overlayRows:result.observed,overlayX:"TIME",overlayY:"SURVIVAL",overlayLines:true,overlayColor:"#17202a",hidePoints:true,title:"Time-to-event VPC",xLabel:"Time",yLabel:"Event-free survival"}),e(SimpleTable,{rows:result.simulated}));
    }
    if (tab === "vpc_competing") {
      var competingRows=list(result.simulated).map(function(row){return {TIME:row.TIME,Y:row.median,LOWER:row.lower,UPPER:row.upper,CAUSE:row.CAUSE};});
      return e("div",{className:"lw-diagnostic-grid"},e(ScatterPlot,{rows:competingRows,x:"TIME",y:"Y",group:"CAUSE",lineGroup:"CAUSE",lines:true,intervalShade:0.18,overlayRows:result.observed,overlayX:"TIME",overlayY:"CIF",overlayLines:true,overlayColor:"#17202a",hidePoints:true,title:"Competing-risk VPC",xLabel:"Time",yLabel:"Cumulative incidence"}),e(SimpleTable,{rows:result.simulated}));
    }
    if (tab === "vpc_recurrent") {
      var recurrentRows=list(result.simulated).map(function(row){return {TIME:row.TIME,Y:row.median,LOWER:row.lower,UPPER:row.upper};});
      return e("div",{className:"lw-diagnostic-grid"},e(ScatterPlot,{rows:recurrentRows,x:"TIME",y:"Y",lines:true,intervalShade:0.28,overlayRows:result.observed,overlayX:"TIME",overlayY:"MEAN_CUMULATIVE",overlayLines:true,overlayColor:"#17202a",hidePoints:true,title:"Recurrent-event VPC",xLabel:"Time",yLabel:"Mean cumulative events"}),e(SimpleTable,{rows:result.simulated}));
    }
    if (tab === "bootstrap") return e("div",{className:"lw-diagnostic-run"},e("p",null,value(result.successful,0)+" of "+value(result.n,0)+" bootstrap fits completed."),e(SimpleTable,{rows:result.summary}));
    if (tab === "profile") return e("div",{className:"lw-diagnostic-run"},e(SimpleTable,{rows:result.intervals}),e(ScatterPlot,{rows:result.grid,x:"value",y:"delta",group:"parameter",lines:true,title:"Profile likelihood",xLabel:"Fixed parameter value",yLabel:"Objective difference"}));
    if (tab === "scm") return e("div",{className:"lw-diagnostic-run"},e("p",null,"Objective: "+formatNumber(result.base_objective)+" to "+formatNumber(result.final_objective)),e(SimpleTable,{rows:result.selected}),e(SimpleTable,{rows:result.steps}));
    return null;
  }
  var CachedDiagnosticsPane = React.memo(DiagnosticsPane, function(previous,next) {
    function signature(props) {
      var tab=props.tab,fit=props.fit||{},diagnostics=props.diagnostics||{},result=diagnostics[tab]||{};
      return [tab,value(props.workspace&&props.workspace.current_run,""),!!fit.available,!!fit.gof_loaded,list(fit.gof).length,!!(diagnostics.available&&diagnostics.available[tab]),value(result.nsim,""),list(result.table).length,list(result.observed).length,list(result.simulated).length,list(result.summary).length,list(result.grid).length,list(result.intervals).length,list(result.steps).length].join("|");
    }
    return signature(previous)===signature(next);
  });

  function HmmPane(props) {
    var hmm=props.hmm||{},rows=list(hmm.rows),states=list(hmm.states),summaries=list(hmm.sequence_summary);
    var decoderState=React.useState("smoothed"),subjectState=React.useState(""),sequenceState=React.useState(""),hiddenState=React.useState("");
    if (!hmm.loaded) return e(Empty,{title:"Loading HMM results",detail:"Filtering, retrospective smoothing and Viterbi decoding are calculated once for the selected fitted run."});
    if (!rows.length) return e(Empty,{title:"No HMM observations",detail:"The fitted run did not contain decodable observation records."});
    function uniqueFor(source,key){var seen={},output=[];source.forEach(function(row){var item=String(value(row[key],""));if(item&&!seen[item]){seen[item]=true;output.push(item);}});return output;}
    var subjects=uniqueFor(rows,"SUBJECT"),selectedSubject=subjects.indexOf(subjectState[0])>=0?subjectState[0]:subjects[0];
    var subjectRows=rows.filter(function(row){return String(value(row.SUBJECT,""))===selectedSubject;});
    var sequences=uniqueFor(subjectRows,"SEQUENCE"),selectedSequence=sequences.indexOf(sequenceState[0])>=0?sequenceState[0]:sequences[0];
    var selectedRows=subjectRows.filter(function(row){return String(value(row.SEQUENCE,""))===selectedSequence;});
    var stateKeys=states.map(function(item){return String(item.key);}),selectedStateKey=stateKeys.indexOf(hiddenState[0])>=0?hiddenState[0]:stateKeys[0];
    var stateInfo=states.filter(function(item){return String(item.key)===selectedStateKey;})[0]||states[0]||{label:"state",key:"",index:1};
    var decoder=decoderState[0],probabilityRows=[],pathRows=[];
    selectedRows.forEach(function(row){
      if(decoder==="combined"){
        probabilityRows.push({TIME:row.TIME,VALUE:row["HMM_FILTER_PROB_"+stateInfo.key],METHOD:"Filtered"});
        probabilityRows.push({TIME:row.TIME,VALUE:row["HMM_SMOOTH_PROB_"+stateInfo.key],METHOD:"Smoothed"});
        pathRows.push({TIME:row.TIME,VALUE:row.HMM_FILTER_STATE_INDEX,METHOD:"Filtered"});
        pathRows.push({TIME:row.TIME,VALUE:row.HMM_SMOOTH_STATE_INDEX,METHOD:"Smoothed"});
        pathRows.push({TIME:row.TIME,VALUE:row.HMM_VITERBI_STATE_INDEX,METHOD:"Viterbi"});
      }else{
        var prefix=decoder==="filtered"?"FILTER":decoder==="smoothed"?"SMOOTH":"VITERBI";
        if(decoder!=="viterbi")probabilityRows.push({TIME:row.TIME,VALUE:row["HMM_"+prefix+"_PROB_"+stateInfo.key],METHOD:decoder.charAt(0).toUpperCase()+decoder.slice(1)});
        pathRows.push({TIME:row.TIME,VALUE:row["HMM_"+prefix+"_STATE_INDEX"],METHOD:decoder.charAt(0).toUpperCase()+decoder.slice(1)});
      }
    });
    var summaryRows=summaries.filter(function(row){return String(value(row.SUBJECT,""))===selectedSubject&&String(value(row.SEQUENCE,""))===selectedSequence;});
    var probabilityColumns=["SUBJECT","SEQUENCE","TIME","HMM_FILTER_STATE","HMM_SMOOTH_STATE","HMM_VITERBI_STATE","HMM_ROW_NLL"];
    if(decoder==="filtered"||decoder==="combined")probabilityColumns.push("HMM_FILTER_PROB_"+stateInfo.key);
    if(decoder==="smoothed"||decoder==="combined")probabilityColumns.push("HMM_SMOOTH_PROB_"+stateInfo.key);
    return e("div",{className:"lw-hmm-view"},
      e("div",{className:"lw-hmm-controls"},
        e(Field,{label:"Decoder"},e("select",{value:decoder,onChange:function(event){decoderState[1](event.target.value);}},
          e("option",{value:"filtered"},"Filtered"),e("option",{value:"smoothed"},"Retrospective smoothed"),e("option",{value:"viterbi"},"Viterbi path"),e("option",{value:"combined"},"Combined"))),
        e(Field,{label:"Subject"},e("select",{value:selectedSubject,onChange:function(event){subjectState[1](event.target.value);sequenceState[1]("");}},subjects.map(function(item){return e("option",{key:item,value:item},item);}))),
        sequences.length>1?e(Field,{label:"Sequence / DVID"},e("select",{value:selectedSequence,onChange:function(event){sequenceState[1](event.target.value);}},sequences.map(function(item){return e("option",{key:item,value:item},item);}))) : null,
        e(Field,{label:"State probability"},e("select",{value:stateInfo.key,onChange:function(event){hiddenState[1](event.target.value);}},states.map(function(item){return e("option",{key:item.key,value:item.key},item.label);})))),
      e("div",{className:"lw-hmm-meta"},
        e("span",null,"Log likelihood ",e("strong",null,formatNumber(hmm.log_likelihood))),
        e("span",null,value(hmm.observations,0)+" observations"),e("span",null,value(hmm.sequences,0)+" sequences"),
        e("span",null,"ETA: "+value(hmm.eta_type,"individual")),hmm.truncated?e("strong",{className:"lw-warning-text"},"Display limited to the first 50,000 observations") : null),
      e("div",{className:"lw-hmm-state-key"},states.map(function(item){return e("span",{key:item.key},item.index+" = "+item.label);})),
      e("div",{className:"lw-diagnostic-grid lw-hmm-grid"},
        decoder==="viterbi"?e(Empty,{title:"Viterbi has no marginal state probability",detail:"Select Filtered, Retrospective smoothed, or Combined to inspect uncertainty for a state."}):
          e(ScatterPlot,{rows:probabilityRows,x:"TIME",y:"VALUE",group:"METHOD",lines:true,yRange:[0,1],title:(decoder==="combined"?"Filtered and smoothed probability of ":"State probability: ")+stateInfo.label,xLabel:"Time",yLabel:"Probability"}),
        e(ScatterPlot,{rows:pathRows,x:"TIME",y:"VALUE",group:"METHOD",lines:true,yRange:[0.5,Math.max(1.5,states.length+0.5)],title:decoder==="viterbi"?"Most probable Viterbi path":"Decoded state over time",xLabel:"Time",yLabel:"State index"})),
      e("section",{className:"lw-hmm-summary"},e("h4",null,"Sequence likelihood and Viterbi path evidence"),e(SimpleTable,{rows:summaryRows})),
      e("details",{className:"lw-hmm-table"},e("summary",null,"Observation-level decoded values"),e(SimpleTable,{rows:selectedRows.slice(0,500),columns:probabilityColumns}),selectedRows.length>500?e("p",{className:"lw-muted"},"Showing the first 500 rows for this sequence."):null));
  }

  function KalmanPane(props) {
    var kalman=props.kalman||{},rows=list(kalman.rows),states=list(kalman.states);
    var estimateState=React.useState("combined"),subjectState=React.useState(""),sequenceState=React.useState(""),latentState=React.useState("");
    if (!kalman.loaded) return e(Empty,{title:"Loading state estimates",detail:"Kalman filtering and retrospective smoothing are calculated once for the selected fitted run."});
    if (!rows.length) return e(Empty,{title:"No state-space observations",detail:"The fitted run did not contain filterable observation records."});
    function uniqueFor(source,key){var seen={},output=[];source.forEach(function(row){var item=String(value(row[key],""));if(item&&!seen[item]){seen[item]=true;output.push(item);}});return output;}
    var subjects=uniqueFor(rows,"SUBJECT"),selectedSubject=subjects.indexOf(subjectState[0])>=0?subjectState[0]:subjects[0];
    var subjectRows=rows.filter(function(row){return String(value(row.SUBJECT,""))===selectedSubject;});
    var sequences=uniqueFor(subjectRows,"SEQUENCE"),selectedSequence=sequences.indexOf(sequenceState[0])>=0?sequenceState[0]:sequences[0];
    var selectedRows=subjectRows.filter(function(row){return String(value(row.SEQUENCE,""))===selectedSequence;});
    var stateKeys=states.map(function(item){return String(item.key);}),selectedStateKey=stateKeys.indexOf(latentState[0])>=0?latentState[0]:stateKeys[0];
    var stateInfo=states.filter(function(item){return String(item.key)===selectedStateKey;})[0]||states[0]||{label:"state",key:"",index:1};
    var estimate=estimateState[0],stateRows=[];
    selectedRows.forEach(function(row){
      function add(method,prefix,sdPrefix){
        var mean=Number(row[prefix+stateInfo.key]),sd=Number(row[sdPrefix+stateInfo.key]);
        if(Number.isFinite(mean))stateRows.push({TIME:row.TIME,Y:mean,LOWER:Number.isFinite(sd)?mean-1.96*sd:null,UPPER:Number.isFinite(sd)?mean+1.96*sd:null,METHOD:method});
      }
      if(estimate==="filtered"||estimate==="combined")add("Filtered","KF_FILTER_","KF_FILTER_SD_");
      if(estimate==="smoothed"||estimate==="combined")add("Smoothed","KF_SMOOTH_","KF_SMOOTH_SD_");
    });
    var innovationRows=selectedRows.map(function(row){return {TIME:row.TIME,VALUE:row.KF_STANDARDIZED_INNOVATION,METHOD:"Standardized innovation"};});
    var tableColumns=["SUBJECT","SEQUENCE","TIME","DV","KF_PRED_"+stateInfo.key,"KF_FILTER_"+stateInfo.key,"KF_FILTER_SD_"+stateInfo.key,"KF_SMOOTH_"+stateInfo.key,"KF_SMOOTH_SD_"+stateInfo.key,"KF_STANDARDIZED_INNOVATION","KF_ROW_NLL"];
    return e("div",{className:"lw-hmm-view"},
      e("div",{className:"lw-hmm-controls"},
        e(Field,{label:"Estimate"},e("select",{value:estimate,onChange:function(event){estimateState[1](event.target.value);}},e("option",{value:"filtered"},"Filtered"),e("option",{value:"smoothed"},"Retrospective smoothed"),e("option",{value:"combined"},"Combined"))),
        e(Field,{label:"Subject"},e("select",{value:selectedSubject,onChange:function(event){subjectState[1](event.target.value);sequenceState[1]("");}},subjects.map(function(item){return e("option",{key:item,value:item},item);}))),
        sequences.length>1?e(Field,{label:"Sequence / DVID"},e("select",{value:selectedSequence,onChange:function(event){sequenceState[1](event.target.value);}},sequences.map(function(item){return e("option",{key:item,value:item},item);}))) : null,
        e(Field,{label:"Latent state"},e("select",{value:stateInfo.key,onChange:function(event){latentState[1](event.target.value);}},states.map(function(item){return e("option",{key:item.key,value:item.key},item.label);})))),
      e("div",{className:"lw-hmm-meta"},e("span",null,"Log likelihood ",e("strong",null,formatNumber(kalman.log_likelihood))),e("span",null,value(kalman.observations,0)+" observations"),e("span",null,value(kalman.sequences,0)+" sequences"),e("span",null,"Filter: "+String(value(kalman.filter,"linear")).toUpperCase()),e("span",null,"Smoother: "+value(kalman.smoother,"RTS")),e("span",null,"ETA: "+value(kalman.eta_type,"individual")),kalman.truncated?e("strong",{className:"lw-warning-text"},"Display limited to the first 50,000 observations") : null),
      e("div",{className:"lw-hmm-state-key"},states.map(function(item){return e("span",{key:item.key},item.index+" = "+item.label);})),
      e("div",{className:"lw-diagnostic-grid lw-hmm-grid"},
        e(ScatterPlot,{rows:stateRows,x:"TIME",y:"Y",group:"METHOD",lineGroup:"METHOD",lines:true,intervals:true,intervalShade:0.16,title:(estimate==="combined"?"Filtered and smoothed ":estimate.charAt(0).toUpperCase()+estimate.slice(1)+" ")+stateInfo.label,xLabel:"Time",yLabel:"Latent state"}),
        e(ScatterPlot,{rows:innovationRows,x:"TIME",y:"VALUE",group:"METHOD",lines:false,yRange:[-5,5],title:"Standardized innovations",xLabel:"Time",yLabel:"Innovation / SD"})),
      e("details",{className:"lw-hmm-table"},e("summary",null,"Observation-level state estimates"),e(SimpleTable,{rows:selectedRows.slice(0,500),columns:tableColumns}),selectedRows.length>500?e("p",{className:"lw-muted"},"Showing the first 500 rows for this sequence."):null));
  }

  function defaultTemplateTrans(advan) { return [3,4,11,12].indexOf(Number(advan)) >= 0 ? "4" : Number(advan) === 6 || Number(advan) === 13 ? "1" : "2"; }
  function TemplateFields(props) {
    var advan = Number(props.advan[0]), ode = advan === 6 || advan === 13;
    var labels = {1:"One compartment, IV",2:"One compartment, absorption",3:"Two compartment, IV",4:"Two compartment, absorption",6:"General ODE",11:"Three compartment, IV",12:"Three compartment, absorption",13:"General stiff ODE"};
    return e("div", { className:"lw-template-fields" },
      e("div", { className:"lw-form-grid lw-form-grid-two" },
        e(Field, { label:"Version label (optional)" }, e("input", { value:props.label[0], onChange:function(event){props.label[1](event.target.value);} })),
        e(Field, { label:"Problem statement" }, e("input", { value:props.problem[0], onChange:function(event){props.problem[1](event.target.value);} }))),
      e("div", { className:"lw-form-grid" },
        e(Field, { label:"ADVAN" }, e("select", { value:props.advan[0], onChange:function(event){props.advan[1](event.target.value);props.trans[1](defaultTemplateTrans(event.target.value));} }, [1,2,3,4,6,11,12,13].map(function(item){return e("option",{key:item,value:item},"ADVAN"+item+" - "+labels[item]);}))),
        ode ? e(Field, { label:"Compartments" }, e("input", { type:"number", min:1, max:20, value:props.nState[0], onChange:function(event){props.nState[1](Number(event.target.value));} })) :
          e(Field, { label:"TRANS" }, e("select", { value:props.trans[0], onChange:function(event){props.trans[1](event.target.value);} }, [1,2,3,4,5,6].map(function(item){return e("option",{key:item,value:item},"TRANS"+item);}))),
        e("div", { className:"lw-template-summary" }, e("strong",null,labels[advan]), e("span",null,ode?"C++ ODE solver with an editable $DES block.":"Analytical C++ ADVAN solution with NONMEM-style parameters."))),
      e("p", { className:"lw-help-text" }, "Creates standard THETA, OMEGA, SIGMA, PK/PRED and ERROR blocks. The generated version remains fully editable."));
  }

  function ProjectRunRow(options) {
    var run=options.run,workspace=options.workspace,queued=!!run.queued_job;
    var estimation=run.result_type==="estimation",runSelected=!queued&&workspace.current_run===run.id;
    var status=value(run.job_status,"queued");
    return e("div",{className:"lw-run-row "+(runSelected?"selected ":"")+(queued?"lw-run-pending lw-run-pending-"+status:""),key:run.id},
      estimation&&!queued?e("label",{className:"lw-compare-check",title:"Select estimation run for comparison",onClick:function(event){event.stopPropagation();}},e("input",{type:"checkbox",checked:options.comparison.indexOf(run.id)>=0,onChange:function(event){options.toggleComparison(run.id,event.target.checked);}})):null,
      e("button",{type:"button",className:"lw-run-main",onClick:function(){if(queued)emit(options.workbench,"job_select",{id:run.job_id,queueId:run.queue_id});else emit(options.workbench,"run_open",{id:workspace.current,run:run.id});}},
        e("span",{className:"lw-run-number"},queued?"Job":(run.result_type==="simulation"?"Sim":"Run")+String(value(run.run_number,"?")).padStart(3,"0")),
        e("span",{className:"lw-run-copy"},e("strong",null,run.label),e("small",null,queued?status:value(run.method,run.result_type))),
        queued?e("span",{className:"lw-run-job-status"},status):null,
        run.has_vpc?e("span",{className:"lw-run-flag lw-run-flag-vpc"},"VPC"):null,
        run.has_npde?e("span",{className:"lw-run-flag lw-run-flag-npde"},"NPDE"):null,
        run.has_npc?e("span",{className:"lw-run-flag lw-run-flag-npc"},"NPC"):null,
        run.has_covariance?e("span",{className:"lw-run-flag lw-run-flag-cov"},"COV"):null));
  }

  function ProjectTree(props) {
    var workspace = props.workspace || {}, projects = list(workspace.projects), versions = list(workspace.versions);
    var newProject = React.useState(false), projectName = React.useState(""), projectDescription = React.useState(""), projectMode = React.useState("empty"), projectDataSource = React.useState("synthetic"), projectExample = React.useState("theophylline"), projectSubjects = React.useState(10), projectFile = React.useState(null);
    var projectAdvan = React.useState("2"), projectTrans = React.useState("2"), projectNState = React.useState(2), projectLabel = React.useState(""), projectProblem = React.useState("Synthetic demo");
    var comparison = React.useState([]), templateModal = React.useState(false), copyModal = React.useState(false), copyUpdateInits = React.useState(true), deleteModal = React.useState(null), deleteConfirmation = React.useState("");
    var expandedVersions = React.useState({});
    var templateAdvan = React.useState("4"), templateTrans = React.useState("4"), templateNState = React.useState(2), templateStructural = React.useState("standard"), templateLabel = React.useState(""), templateProblem = React.useState("Template model");
    var estimationModal = React.useState(false), estimationMethod = React.useState("FOCEI"), estimationLabel = React.useState(""), estimationMaxit = React.useState(200), etaMaxit = React.useState(100), tolerance = React.useState(0.000001), estimationCores = React.useState(1), printEvery = React.useState(0), methodSeed = React.useState(20260713), nImp = React.useState(200), gqOrder = React.useState(5), gqGrid = React.useState("auto"), gqLevel = React.useState(3), gqAdaptive = React.useState(true), gqMaxPoints = React.useState(100000), nIter = React.useState(200), burn = React.useState(60), mcmcSteps = React.useState(2), nBurn = React.useState(500), nSample = React.useState(1000), nThin = React.useState(1), nChains = React.useState(4), targetAcceptance = React.useState(0.8), maxTreeDepth = React.useState(10), nLeapfrog = React.useState(10), npPoints = React.useState(25), npCycles = React.useState(3), npMaxSupport = React.useState(100), npGridStep = React.useState(1);
    var estimationPreStages = React.useState([]);
    var covarianceStep = React.useState(false), covarianceType = React.useState("hessian"), covarianceTolerance = React.useState(0.00000001), covarianceSamples = React.useState(200);
    var simulationModal = React.useState(false), simulationLabel = React.useState("Simulation"), simulationSeed = React.useState(Math.floor(Math.random() * 99999) + 1), simulationCores = React.useState(1);
    var simulationSubjects = React.useState(value(props.dataset.subjects, 10)), simulationReplicates = React.useState(1), simulationDays = React.useState(1), simulationUseDesign = React.useState(false);
    var diagnosticModal = React.useState(false), diagnosticVpc = React.useState(true), diagnosticNpc = React.useState(false), diagnosticNpde = React.useState(false), diagnosticCategorical = React.useState(false), diagnosticCount = React.useState(false), diagnosticTte = React.useState(false), diagnosticCompeting = React.useState(false), diagnosticRecurrent = React.useState(false), diagnosticOutcome = React.useState("DV"), diagnosticEvent = React.useState("DV"), diagnosticDvid = React.useState(""), diagnosticNsim = React.useState(200), diagnosticSeed = React.useState(20260713), diagnosticPc = React.useState(false), diagnosticStratify = React.useState("");
    var uncertaintyModal = React.useState(false), uncertaintyBootstrap = React.useState(true), uncertaintyProfile = React.useState(false), uncertaintyReplicates = React.useState(100), uncertaintyPoints = React.useState(9), uncertaintySpan = React.useState(3), uncertaintyLevel = React.useState(0.95), uncertaintyParameters = React.useState(""), uncertaintyMaxit = React.useState(100);
    var scmModal = React.useState(false), scmCandidates = React.useState("CL,WT,power\nV,WT,power"), scmDirection = React.useState("both"), scmForward = React.useState(0.05), scmBackward = React.useState(0.01), scmMaxSteps = React.useState(20), scmMaxit = React.useState(100), scmLabel = React.useState("SCM model");
    var controlModal = React.useState(false), controlFile = React.useState(null), controlData = React.useState(null), controlNewProject = React.useState(!workspace.current), controlProjectName = React.useState("NONMEM import"), controlLabel = React.useState("NONMEM import"), exportModal = React.useState(false), exportName = React.useState("model.ctl"), exportDataPath = React.useState("data.csv");
    var libraryModal = React.useState(false), libraryQuery = React.useState(""), librarySelected = React.useState(""), libraryNewProject = React.useState(!workspace.current), libraryProjectName = React.useState(""), libraryLabel = React.useState("");
    var libraryInfo=props.library||{},libraryEntries=list(libraryInfo.entries),libraryNeedle=libraryQuery[0].trim().toLowerCase();
    var libraryFiltered=libraryEntries.filter(function(item){return !libraryNeedle||[item.library_id,item.title,item.compound,item.population,item.status,"ADVAN"+value(item.advan,"")].join(" ").toLowerCase().indexOf(libraryNeedle)>=0;});
    var doseMode = React.useState("single"), doseAmount = React.useState(320), doseCmt = React.useState(1), doseN = React.useState(3), doseII = React.useState(12), doseTable = React.useState("0 320"), obsPerDay = React.useState(8), simulationUseFit = React.useState(true);
    var userLikelihood = value(props.model.likelihood_type, "none") === "likelihood";
    var outcomeFamilies=list(props.model.outcomes).map(function(item){return item.family;}),outcomeDvids=list(props.model.outcomes).filter(function(item){return item.dvid!==null&&item.dvid!==undefined;});
    var estimationMethods = ["FO","FOCE","FOCEI","LAPLACE","ITS","GQ","IMP","SAEM","BAYES","HMC","NUTS","NPML","NPAG"].filter(function(method){return !userLikelihood||["FO","FOCE","FOCEI"].indexOf(method)<0;});
    var sequenceMethods = ["FO","FOCE","FOCEI","LAPLACE","ITS","GQ","IMP","SAEM"].filter(function(method){return !userLikelihood||["FO","FOCE","FOCEI"].indexOf(method)<0;});
    React.useEffect(function(){if(workspace.current_version){var next=Object.assign({},expandedVersions[0]);next[workspace.current_version]=true;expandedVersions[1](next);}},[workspace.current_version]);
    React.useEffect(function(){if(userLikelihood){if(["FO","FOCE","FOCEI"].indexOf(estimationMethod[0])>=0)estimationMethod[1]("LAPLACE");estimationPreStages[1](function(stages){return stages.filter(function(stage){return ["FO","FOCE","FOCEI"].indexOf(stage.method)<0;});});}},[userLikelihood,estimationMethod[0]]);
    function toggleExpanded(id) { var next=Object.assign({},expandedVersions[0]);next[id]=!next[id];expandedVersions[1](next); }
    function toggleComparison(id, checked) { var next=comparison[0].filter(function(item){return item!==id;});if(checked)next=next.concat([id]).slice(-2);comparison[1](next); }
    function readProjectDataset(event) { var file=event.target.files&&event.target.files[0];if(!file){projectFile[1](null);return;}var reader=new FileReader();reader.onload=function(){projectFile[1]({name:file.name,text:String(reader.result||"")});};reader.readAsText(file); }
    function readControlFile(event) { var file=event.target.files&&event.target.files[0];if(!file){controlFile[1](null);return;}var reader=new FileReader();reader.onload=function(){controlFile[1]({name:file.name,text:String(reader.result||"")});};reader.readAsText(file); }
    function readControlData(event) { var file=event.target.files&&event.target.files[0];if(!file){controlData[1](null);return;}var reader=new FileReader();reader.onload=function(){controlData[1]({name:file.name,text:String(reader.result||"")});};reader.readAsText(file); }
    function submitNewProject() { emit(props,"project_create",{name:projectName[0],description:projectDescription[0],mode:projectMode[0],dataSource:projectDataSource[0],example:projectExample[0],nSubjects:Number(projectSubjects[0]),fileName:projectFile[0]&&projectFile[0].name,text:projectFile[0]&&projectFile[0].text,advan:Number(projectAdvan[0]),trans:Number(projectTrans[0]),nState:Number(projectNState[0]),label:projectLabel[0],problem:projectProblem[0]});newProject[1](false); }
    function submitTemplate() { emit(props,"model_template",{advan:Number(templateAdvan[0]),trans:Number(templateTrans[0]),nState:Number(templateNState[0]),structuralTemplate:templateStructural[0],label:templateLabel[0],problem:templateProblem[0]});templateModal[1](false); }
    var covarianceSupported = ["FO","FOCE","FOCEI","LAPLACE","ITS","GQ","IMP","SAEM"].indexOf(estimationMethod[0]) >= 0;
    function updatePreStage(index,field,nextValue){var next=estimationPreStages[0].map(function(stage){return Object.assign({},stage);});next[index][field]=nextValue;estimationPreStages[1](next);}
    function movePreStage(index,direction){var next=estimationPreStages[0].slice(),target=index+direction;if(target<0||target>=next.length)return;var item=next[index];next[index]=next[target];next[target]=item;estimationPreStages[1](next);}
    function submitEstimate() { var covType=estimationMethod[0]==="FO"?"hessian":covarianceType[0],covSamples=estimationMethod[0]==="IMP"?Number(nImp[0]):Number(covarianceSamples[0]);var payload={label:estimationLabel[0],method:estimationMethod[0],maxit:Number(estimationMaxit[0]),etaMaxit:Number(etaMaxit[0]),tolerance:Number(tolerance[0]),nCores:Number(estimationCores[0]),printEvery:Number(printEvery[0]),methodSeed:Number(methodSeed[0]),nImp:Number(nImp[0]),gqOrder:Number(gqOrder[0]),gqGrid:gqGrid[0],gqLevel:Number(gqLevel[0]),gqAdaptive:!!gqAdaptive[0],gqMaxPoints:Number(gqMaxPoints[0]),nIter:Number(nIter[0]),burn:Number(burn[0]),mcmcSteps:Number(mcmcSteps[0]),nBurn:Number(nBurn[0]),nSample:Number(nSample[0]),nThin:Number(nThin[0]),nChains:Number(nChains[0]),targetAcceptance:Number(targetAcceptance[0]),maxTreeDepth:Number(maxTreeDepth[0]),nLeapfrog:Number(nLeapfrog[0]),npPoints:Number(npPoints[0]),npCycles:Number(npCycles[0]),npMaxSupport:Number(npMaxSupport[0]),npGridStep:Number(npGridStep[0]),covariance:covarianceStep[0]&&covarianceSupported,covarianceType:covType,covarianceTolerance:Number(covarianceTolerance[0]),covarianceSamples:covSamples,covarianceSeed:Number(methodSeed[0])};payload.stages=estimationPreStages[0].map(function(stage){return {method:stage.method,maxit:Number(stage.maxit),etaMaxit:Number(stage.etaMaxit),tolerance:Number(tolerance[0]),nCores:Number(estimationCores[0]),printEvery:Number(printEvery[0])};}).concat([Object.assign({},payload)]);emit(props,"estimate",payload);estimationModal[1](false); }
    function submitSimulation() { emit(props,"simulate",{label:simulationLabel[0],seed:Number(simulationSeed[0]),nCores:Number(simulationCores[0]),nSubjects:Number(simulationSubjects[0]),replicates:Number(simulationReplicates[0]),days:Number(simulationDays[0]),useDesign:simulationUseDesign[0],doseMode:doseMode[0],doseAmt:Number(doseAmount[0]),doseCmt:Number(doseCmt[0]),doseN:Number(doseN[0]),doseII:Number(doseII[0]),doseTable:doseTable[0],obsPerDay:Number(obsPerDay[0]),useFit:simulationUseFit[0]});simulationModal[1](false); }
    function submitDiagnostic() { var types=[];if(diagnosticVpc[0])types.push("vpc");if(diagnosticNpc[0])types.push("npc");if(diagnosticNpde[0])types.push("npde");if(diagnosticCategorical[0])types.push("vpc_categorical");if(diagnosticCount[0])types.push("vpc_count");if(diagnosticTte[0])types.push("vpc_tte");if(diagnosticCompeting[0])types.push("vpc_competing");if(diagnosticRecurrent[0])types.push("vpc_recurrent");emit(props,"run_diagnostic",{types:types,nsim:Number(diagnosticNsim[0]),seed:Number(diagnosticSeed[0]),pcCorrect:diagnosticPc[0],stratify:diagnosticVpc[0]&&diagnosticStratify[0]?diagnosticStratify[0]:null,categoricalOutcome:diagnosticOutcome[0],countOutcome:diagnosticOutcome[0],countDvid:diagnosticDvid[0],tteEvent:diagnosticEvent[0],competingDvid:diagnosticDvid[0],recurrentDvid:diagnosticDvid[0]});diagnosticModal[1](false); }
    function submitUncertainty() { var types=[];if(uncertaintyBootstrap[0])types.push("bootstrap");if(uncertaintyProfile[0])types.push("profile");emit(props,"run_uncertainty",{types:types,replicates:Number(uncertaintyReplicates[0]),points:Number(uncertaintyPoints[0]),span:Number(uncertaintySpan[0]),level:Number(uncertaintyLevel[0]),parameters:uncertaintyParameters[0],maxit:Number(uncertaintyMaxit[0]),seed:Number(diagnosticSeed[0])});uncertaintyModal[1](false); }
    var projectUploadMissing=projectMode[0]==="template"&&projectDataSource[0]==="upload"&&!projectFile[0];
    return e("div", { className:"lw-project-sidebar" },
      e("div",{className:"lw-tree-title"},"Projects"),
      e("div",{className:"lw-tree-list lw-project-list"},projects.length?projects.map(function(project){return e("button",{type:"button",key:project.id,className:workspace.current===project.id?"selected":"",title:value(project.description,""),onClick:function(){emit(props,"project_open",{id:project.id});}},e("strong",null,project.name),e("span",null,value(project.versions,project.snapshots)+" versions"));}):e(Empty,{title:"No projects",detail:"Create a project below."})),
      e("div",{className:"lw-sidebar-actions lw-action-grid"},e(Button,{className:"lw-button-primary",icon:"+",title:"Create a new project",disabled:!workspace.enabled,onClick:function(){newProject[1](true);}},"New project"),e(Button,{className:"lw-button-quiet lw-action-button",icon:"NM",title:"Load a NONMEM control stream",disabled:!workspace.enabled,onClick:function(){controlNewProject[1](!workspace.current);controlModal[1](true);}},"Load .ctl"),e(Button,{className:"lw-button-quiet lw-action-button lw-action-span",icon:"L",title:"Browse and import from the LibeRary model catalogue",disabled:!workspace.enabled,onClick:function(){libraryNewProject[1](!workspace.current);libraryModal[1](true);}},"Model library")),
      e("div",{className:"lw-tree-title lw-version-title"},"Model versions"),
      e("div",{className:"lw-tree-list lw-version-list"},versions.length?versions.map(function(version){var runs=list(version.runs),expanded=!!expandedVersions[0][version.id],selected=workspace.current_version===version.id;return e("div",{className:"lw-version-group",key:version.id},e("div",{className:"lw-version-row "+(selected?"selected":"")},e("button",{type:"button",className:"lw-version-toggle",title:expanded?"Collapse runs":"Expand runs",onClick:function(){toggleExpanded(version.id);}},expanded?"▼":"▶"),e("button",{type:"button",className:"lw-version-main",onClick:function(){var next=Object.assign({},expandedVersions[0]);next[version.id]=true;expandedVersions[1](next);emit(props,"project_open",{id:workspace.current,snapshot:version.id});}},e("span",{className:"lw-version-number"},"v"+value(version.version,"?")),e("span",{className:"lw-version-copy"},e("strong",null,version.label),e("small",null,runs.length+" run"+(runs.length===1?"":"s"))))),expanded?e("div",{className:"lw-run-list"},runs.length?runs.map(function(run){return e(ProjectRunRow,{key:run.id,run:run,workspace:workspace,comparison:comparison[0],toggleComparison:toggleComparison,workbench:props});}):e("div",{className:"lw-run-empty"},"No runs yet")):null);}):e(Empty,{title:"No versions",detail:"Create from a template or save the current model."})),
      e("div",{className:"lw-compare-runs"},e(Button,{className:"lw-button-quiet lw-action-button",icon:"=",title:"Compare two selected estimation runs",disabled:comparison[0].length!==2,onClick:function(){emit(props,"project_compare",{runs:comparison[0]});}},comparison[0].length?"Compare selected ("+comparison[0].length+"/2)":"Compare runs")),
      e("div",{className:"lw-action-group"},e("span",{className:"lw-action-group-label"},"Model workflow"),e("div",{className:"lw-sidebar-actions lw-action-grid"},
        e(Button,{className:"lw-button-quiet lw-action-button",icon:"C",title:"Copy the selected model version",disabled:!workspace.current_version,onClick:function(){copyUpdateInits[1](!!props.fit.available);copyModal[1](true);}},"Copy version"),
        e(Button,{className:"lw-button-quiet lw-action-button",icon:"{}",title:"Create a model version from a template",disabled:!workspace.current,onClick:function(){templateModal[1](true);}},"From template"),
        e(Button,{className:"lw-button-quiet lw-action-button",icon:">",title:"Run a model estimation",disabled:!props.model.loaded||!props.dataset.loaded,onClick:function(){estimationModal[1](true);}},"Estimate"),
        e(Button,{className:"lw-button-quiet lw-action-button",icon:"D",title:"Run VPC, NPDE or NPC diagnostics",disabled:workspace.current_result_type!=="estimation"||!workspace.current_run,onClick:function(){diagnosticModal[1](true);}},"Diagnostics"),
        e(Button,{className:"lw-button-quiet lw-action-button",icon:"U",title:"Run bootstrap or profile-likelihood uncertainty",disabled:workspace.current_result_type!=="estimation"||!workspace.current_run,onClick:function(){uncertaintyModal[1](true);}},"Uncertainty"),
        e(Button,{className:"lw-button-quiet lw-action-button",icon:"SCM",title:"Run stepwise covariate modelling",disabled:workspace.current_result_type!=="estimation"||!workspace.current_run,onClick:function(){scmModal[1](true);}},"Covariates"),
        e(Button,{className:"lw-button-quiet lw-action-button",icon:".ctl",title:"Export the current model as a NONMEM control stream",disabled:!props.model.loaded,onClick:function(){exportModal[1](true);}},"Export .ctl"),
        e(Button,{className:"lw-button-quiet lw-action-button lw-action-span",icon:"~",title:"Create a simulation from the selected model",disabled:!props.model.loaded||!props.dataset.loaded,onClick:function(){simulationModal[1](true);}},"Simulate"))),
      e("div",{className:"lw-action-group lw-danger-zone"},e("span",{className:"lw-action-group-label"},"Delete"),e("div",{className:"lw-sidebar-actions lw-action-grid"},
        e(Button,{className:"lw-button-danger-ghost",title:"Delete the current project",disabled:!workspace.current,onClick:function(){deleteConfirmation[1]("");deleteModal[1]("project");}},"Project"),
        e(Button,{className:"lw-button-danger-ghost",title:"Delete the selected model version",disabled:!workspace.current_version,onClick:function(){deleteModal[1]("version");}},"Version"),
        e(Button,{className:"lw-button-danger-ghost lw-action-span",title:"Delete the selected estimation or simulation run",disabled:!workspace.current_run,onClick:function(){deleteModal[1]("run");}},"Selected run"))),

      e(Modal,{open:newProject[0],className:"lw-modal-wide",onClose:function(){newProject[1](false);},title:"New project",footer:e(React.Fragment,null,e(Button,{className:"lw-button-quiet",onClick:function(){newProject[1](false);}},"Cancel"),e(Button,{className:"lw-button-primary",disabled:!projectName[0].trim()||projectUploadMissing,onClick:submitNewProject},"Create project"))},
        e("div",{className:"lw-modal-section"},e("div",{className:"lw-form-grid lw-form-grid-two"},e(Field,{label:"Project name"},e("input",{autoFocus:true,value:projectName[0],onChange:function(event){projectName[1](event.target.value);}})),e(Field,{label:"Description (optional)"},e("textarea",{rows:2,value:projectDescription[0],onChange:function(event){projectDescription[1](event.target.value);}})))),
        e("div",{className:"lw-choice-cards"},[["empty","Empty project"],["template","Create from template"]].map(function(item){return e("label",{key:item[0],className:"lw-choice-card "+(projectMode[0]===item[0]?"selected":"")},e("input",{type:"radio",name:"project-mode",checked:projectMode[0]===item[0],onChange:function(){projectMode[1](item[0]);}}),e("span",null,e("strong",null,item[1])));})),
        projectMode[0]==="template"?e("div",{className:"lw-modal-section lw-modal-section-tinted"},e("h4",null,"Initial model version"),e("p",{className:"lw-help-text"},"Choose a built-in example or import an existing NONMEM-style dataset."),e("div",{className:"lw-choice-row"},e("label",{className:"lw-check"},e("input",{type:"radio",name:"project-data",checked:projectDataSource[0]==="synthetic",onChange:function(){projectDataSource[1]("synthetic");}})," Built-in synthetic example"),e("label",{className:"lw-check"},e("input",{type:"radio",name:"project-data",checked:projectDataSource[0]==="upload",onChange:function(){projectDataSource[1]("upload");}})," Upload dataset")),projectDataSource[0]==="synthetic"?e(React.Fragment,null,e("div",{className:"lw-form-grid lw-form-grid-two"},e(Field,{label:"Example"},e("select",{value:projectExample[0],onChange:function(event){projectExample[1](event.target.value);}},e("option",{value:"theophylline"},"Theophylline-style oral PK"),e("option",{value:"sparse"},"Sparse oral PK"),e("option",{value:"rich"},"Rich sampling oral PK"))),e(Field,{label:"Number of subjects"},e("input",{type:"number",min:1,max:500,value:projectSubjects[0],onChange:function(event){projectSubjects[1](Number(event.target.value));}}))),e("div",{className:"lw-form-grid lw-form-grid-two"},e(Field,{label:"Version label (optional)"},e("input",{value:projectLabel[0],onChange:function(event){projectLabel[1](event.target.value);}})),e(Field,{label:"Problem statement"},e("input",{value:projectProblem[0],onChange:function(event){projectProblem[1](event.target.value);}})))):e(React.Fragment,null,e(Field,{label:"Dataset file (.csv, .txt, .dat, .tsv)"},e("input",{type:"file",accept:".csv,.txt,.dat,.tsv,text/csv,text/plain",onChange:readProjectDataset})),e("p",{className:"lw-help-text"},projectFile[0]?"Loaded "+projectFile[0].name:"Expected NONMEM-style ID, TIME, DV, AMT, EVID, CMT and MDV columns."),e(TemplateFields,{advan:projectAdvan,trans:projectTrans,nState:projectNState,label:projectLabel,problem:projectProblem}))):null),

      e(Modal,{open:copyModal[0],onClose:function(){copyModal[1](false);},title:"Copy to new model version",footer:e(React.Fragment,null,e(Button,{className:"lw-button-quiet",onClick:function(){copyModal[1](false);}},"Cancel"),e(Button,{className:"lw-button-primary",onClick:function(){emit(props,"project_copy",{id:workspace.current,snapshot:workspace.current_snapshot,updateInits:copyUpdateInits[0]});copyModal[1](false);}},"Copy"))},e("p",{className:"lw-help-text"},props.fit.available?"A fitted run is loaded; its final estimates can become the new version's initial values.":"No fitted run is loaded, so initials will match the source version."),e("label",{className:"lw-check"},e("input",{type:"checkbox",disabled:!props.fit.available,checked:copyUpdateInits[0]&&props.fit.available,onChange:function(event){copyUpdateInits[1](event.target.checked);}})," Update THETA / OMEGA / SIGMA initials from current fit")),

      e(Modal,{open:templateModal[0],onClose:function(){templateModal[1](false);},title:"New version from template",footer:e(React.Fragment,null,e(Button,{className:"lw-button-quiet",onClick:function(){templateModal[1](false);}},"Cancel"),e(Button,{className:"lw-button-primary",disabled:!props.dataset.loaded,onClick:submitTemplate},"Create version"))},e(Field,{label:"Dataset"},e("select",{disabled:!props.dataset.loaded,value:props.dataset.loaded?"current":""},e("option",{value:props.dataset.loaded?"current":""},props.dataset.loaded?value(props.dataset.name,"Current dataset"):"No dataset loaded"))),e(Field,{label:"Model family"},e("select",{value:templateStructural[0],onChange:function(event){templateStructural[1](event.target.value);}},e("option",{value:"standard"},"Standard ADVAN template"),[["nonlinear_elimination","Nonlinear elimination"],["transit_absorption","Transit absorption"],["dual_absorption","Dual absorption"],["parent_metabolite","Parent–metabolite"],["effect_compartment","Effect compartment"],["indirect_response","Indirect response"],["tumour_growth","Tumour growth"],["tmdd","Target-mediated disposition"]].map(function(item){return e("option",{key:item[0],value:item[0]},item[1]);}))),templateStructural[0]==="standard"?e(TemplateFields,{advan:templateAdvan,trans:templateTrans,nState:templateNState,label:templateLabel,problem:templateProblem}):e("div",{className:"lw-modal-section lw-modal-section-tinted"},e("p",{className:"lw-help-text"},"Creates a complete editable ADVAN13 $PK/$PRED and $DES model. Initial-state requirements are documented with the template."),e("div",{className:"lw-form-grid lw-form-grid-two"},e(Field,{label:"Version label (optional)"},e("input",{value:templateLabel[0],onChange:function(event){templateLabel[1](event.target.value);}})),e(Field,{label:"Problem statement"},e("input",{value:templateProblem[0],onChange:function(event){templateProblem[1](event.target.value);}}))))),

      e(Modal,{open:estimationModal[0],className:"lw-modal-wide",onClose:function(){estimationModal[1](false);},title:"Run estimation",footer:e(React.Fragment,null,e(Button,{className:"lw-button-quiet",onClick:function(){estimationModal[1](false);}},"Cancel"),e(Button,{className:"lw-button-primary",onClick:submitEstimate},"Submit estimation"))},
        e("div",{className:"lw-modal-section"},userLikelihood?e("div",{className:"lw-info-banner"},"User-defined likelihood detected. LAPLACE is the NONMEM-like default; Gaussian FO/FOCE/FOCEI linearizations are not applicable."):null,e("div",{className:"lw-form-grid"},e(Field,{label:"Run on"},e("select",{value:value(props.server.queue_id,"local"),onChange:function(event){emit(props,"queue_select",{id:event.target.value});}},list(props.server.queues).map(function(queue){return e("option",{key:queue.id,value:queue.id},queue.name);}))),e(Field,{label:"Method"},e("select",{value:estimationMethod[0],onChange:function(event){estimationMethod[1](event.target.value);}},estimationMethods.map(function(method){var label=method==="GQ"?"GQ (adaptive Gauss-Hermite)":method==="NPML"?"NPML (fixed support)":method==="NPAG"?"NPAG (adaptive grid)":method;return e("option",{key:method,value:method},label);}))),e(Field,{label:"Job label (optional)"},e("input",{value:estimationLabel[0],onChange:function(event){estimationLabel[1](event.target.value);}}))),e("div",{className:"lw-form-grid"},["HMC","NUTS"].indexOf(estimationMethod[0])<0?e(Field,{label:"Outer iterations"},e("input",{type:"number",min:1,value:estimationMaxit[0],onChange:function(event){estimationMaxit[1](Number(event.target.value));}})):null,["BAYES","HMC","NUTS"].indexOf(estimationMethod[0])<0?e(Field,{label:"ETA iterations"},e("input",{type:"number",min:1,value:etaMaxit[0],onChange:function(event){etaMaxit[1](Number(event.target.value));}})):null,e(Field,{label:"Tolerance"},e("input",{type:"number",min:1e-12,step:"any",value:tolerance[0],onChange:function(event){tolerance[1](Number(event.target.value));}})),["HMC","NUTS","NPML","NPAG"].indexOf(estimationMethod[0])<0?e(Field,{label:"Parallel cores"},e("input",{type:"number",min:1,max:64,value:estimationCores[0],onChange:function(event){estimationCores[1](Number(event.target.value));}})):null,e(Field,{label:"Print gradients every N (0 = off)"},e("input",{type:"number",min:0,value:printEvery[0],onChange:function(event){printEvery[1](Number(event.target.value));}})))),
        e("div",{className:"lw-modal-section lw-modal-section-tinted lw-estimation-sequence"},e("div",{className:"lw-prior-heading"},e("div",null,e("h4",null,"Sequential estimation"),e("small",null,"Completed steps pass THETA, OMEGA, SIGMA and compatible ETA starts to the next step.")),e(Button,{className:"lw-button-quiet",onClick:function(){estimationPreStages[1](estimationPreStages[0].concat([{method:userLikelihood?"LAPLACE":"FOCE",maxit:100,etaMaxit:100}]));}},"+ Add preceding step")),estimationPreStages[0].map(function(stage,index){return e("div",{className:"lw-sequence-row",key:index},e("strong",null,"Step "+(index+1)),e("select",{value:stage.method,onChange:function(event){updatePreStage(index,"method",event.target.value);}},sequenceMethods.map(function(method){return e("option",{key:method,value:method},method);})),e("label",null,"Iterations",e("input",{type:"number",min:1,value:stage.maxit,onChange:function(event){updatePreStage(index,"maxit",Number(event.target.value));}})),e("label",null,"ETA iterations",e("input",{type:"number",min:1,value:stage.etaMaxit,onChange:function(event){updatePreStage(index,"etaMaxit",Number(event.target.value));}})),e(Button,{className:"lw-button-link",disabled:index===0,onClick:function(){movePreStage(index,-1);}},"Up"),e(Button,{className:"lw-button-link",disabled:index===estimationPreStages[0].length-1,onClick:function(){movePreStage(index,1);}},"Down"),e(Button,{className:"lw-button-link",onClick:function(){estimationPreStages[1](estimationPreStages[0].filter(function(_,item){return item!==index;}));}},"Remove"));}),e("div",{className:"lw-sequence-final"},e("strong",null,"Final step "+(estimationPreStages[0].length+1)),e("span",null,estimationMethod[0]),e("small",null,"The method-specific controls below apply to this final step."))),
        estimationMethod[0]==="GQ"?e("div",{className:"lw-modal-section lw-modal-section-tinted"},e("h4",null,"Gaussian quadrature"),e("div",{className:"lw-form-grid"},e(Field,{label:"Grid strategy"},e("select",{value:gqGrid[0],onChange:function(event){gqGrid[1](event.target.value);}},e("option",{value:"auto"},"Automatic"),e("option",{value:"tensor"},"Tensor product"),e("option",{value:"smolyak"},"Smolyak sparse grid"))),gqGrid[0]!=="smolyak"?e(Field,{label:"Tensor nodes / ETA"},e("input",{type:"number",min:1,max:50,value:gqOrder[0],onChange:function(event){gqOrder[1](Number(event.target.value));}})):null,gqGrid[0]!=="tensor"?e(Field,{label:"Smolyak level"},e("input",{type:"number",min:1,max:25,value:gqLevel[0],onChange:function(event){gqLevel[1](Number(event.target.value));}})):null,e(Field,{label:"Maximum total grid points"},e("input",{type:"number",min:1,value:gqMaxPoints[0],onChange:function(event){gqMaxPoints[1](Number(event.target.value));}}))),e("label",{className:"lw-check"},e("input",{type:"checkbox",checked:gqAdaptive[0],onChange:function(event){gqAdaptive[1](event.target.checked);}})," Adapt nodes to each subject's conditional ETA mode and curvature"),e("p",{className:"lw-help-text"},gqGrid[0]==="auto"?"Automatic uses the tensor rule for up to three ETAs and a Smolyak sparse grid above that. Tensor order 5 and sparse level 3 are practical starting points.":gqGrid[0]==="smolyak"?"Smolyak grids retain low-order interactions while avoiding order^ETAs growth. Increase the level to assess quadrature convergence.":"The tensor grid uses order^number-of-ETAs points per subject and is most effective for low-dimensional ETA models.")):null,
        estimationMethod[0]==="IMP"?e("div",{className:"lw-modal-section lw-modal-section-tinted"},e("h4",null,"Importance sampling"),e("div",{className:"lw-form-grid lw-form-grid-two"},e(Field,{label:"Importance samples"},e("input",{type:"number",min:5,value:nImp[0],onChange:function(event){nImp[1](Number(event.target.value));}})),e(Field,{label:"Random seed"},e("input",{type:"number",min:1,value:methodSeed[0],onChange:function(event){methodSeed[1](Number(event.target.value));}})))):null,
        estimationMethod[0]==="SAEM"?e("div",{className:"lw-modal-section lw-modal-section-tinted"},e("h4",null,"SAEM controls"),e("div",{className:"lw-form-grid"},e(Field,{label:"SAEM iterations"},e("input",{type:"number",min:2,value:nIter[0],onChange:function(event){nIter[1](Number(event.target.value));}})),e(Field,{label:"Burn-in"},e("input",{type:"number",min:0,value:burn[0],onChange:function(event){burn[1](Number(event.target.value));}})),e(Field,{label:"MCMC steps / subject"},e("input",{type:"number",min:1,value:mcmcSteps[0],onChange:function(event){mcmcSteps[1](Number(event.target.value));}})),e(Field,{label:"Random seed"},e("input",{type:"number",min:1,value:methodSeed[0],onChange:function(event){methodSeed[1](Number(event.target.value));}})))):null,
        estimationMethod[0]==="BAYES"?e("div",{className:"lw-modal-section lw-modal-section-tinted"},e("h4",null,"Bayesian random-walk sampling"),e("div",{className:"lw-form-grid"},e(Field,{label:"Burn-in"},e("input",{type:"number",min:0,value:nBurn[0],onChange:function(event){nBurn[1](Number(event.target.value));}})),e(Field,{label:"Posterior samples"},e("input",{type:"number",min:1,value:nSample[0],onChange:function(event){nSample[1](Number(event.target.value));}})),e(Field,{label:"Thinning interval"},e("input",{type:"number",min:1,value:nThin[0],onChange:function(event){nThin[1](Number(event.target.value));}})),e(Field,{label:"Random seed"},e("input",{type:"number",min:1,value:methodSeed[0],onChange:function(event){methodSeed[1](Number(event.target.value));}})))):null,
        ["HMC","NUTS"].indexOf(estimationMethod[0])>=0?e("div",{className:"lw-modal-section lw-modal-section-tinted"},e("h4",null,estimationMethod[0]+" sampling"),e("div",{className:"lw-form-grid"},e(Field,{label:"Warmup iterations"},e("input",{type:"number",min:0,value:nBurn[0],onChange:function(event){nBurn[1](Number(event.target.value));}})),e(Field,{label:"Samples / chain"},e("input",{type:"number",min:1,value:nSample[0],onChange:function(event){nSample[1](Number(event.target.value));}})),e(Field,{label:"Chains"},e("input",{type:"number",min:1,max:16,value:nChains[0],onChange:function(event){nChains[1](Number(event.target.value));}})),e(Field,{label:"Thinning interval"},e("input",{type:"number",min:1,value:nThin[0],onChange:function(event){nThin[1](Number(event.target.value));}})),e(Field,{label:"Target acceptance"},e("input",{type:"number",min:0.5,max:0.99,step:0.01,value:targetAcceptance[0],onChange:function(event){targetAcceptance[1](Number(event.target.value));}})),estimationMethod[0]==="NUTS"?e(Field,{label:"Maximum tree depth"},e("input",{type:"number",min:1,max:15,value:maxTreeDepth[0],onChange:function(event){maxTreeDepth[1](Number(event.target.value));}})):e(Field,{label:"Leapfrog steps"},e("input",{type:"number",min:1,value:nLeapfrog[0],onChange:function(event){nLeapfrog[1](Number(event.target.value));}})),e(Field,{label:"Random seed"},e("input",{type:"number",min:1,value:methodSeed[0],onChange:function(event){methodSeed[1](Number(event.target.value));}}))),e("p",{className:"lw-help-text"},"Exact joint CppAD gradients, dual-averaged step size, diagonal mass adaptation, divergence diagnostics, R-hat and ESS are retained with the run.")):null,
        ["NPML","NPAG"].indexOf(estimationMethod[0])>=0?e("div",{className:"lw-modal-section lw-modal-section-tinted"},e("h4",null,"Nonparametric population distribution"),e("div",{className:"lw-form-grid"},e(Field,{label:"Initial support points"},e("input",{type:"number",min:1,value:npPoints[0],onChange:function(event){npPoints[1](Number(event.target.value));}})),e(Field,{label:"Alternating cycles"},e("input",{type:"number",min:1,value:npCycles[0],onChange:function(event){npCycles[1](Number(event.target.value));}})),e(Field,{label:"Maximum retained support"},e("input",{type:"number",min:1,value:npMaxSupport[0],onChange:function(event){npMaxSupport[1](Number(event.target.value));}})),estimationMethod[0]==="NPAG"?e(Field,{label:"Initial grid step"},e("input",{type:"number",min:0.001,step:0.1,value:npGridStep[0],onChange:function(event){npGridStep[1](Number(event.target.value));}})):null),e("p",{className:"lw-help-text"},estimationMethod[0]==="NPAG"?"NPAG expands and prunes the ETA support grid while optimizing mixture weights.":"NPML optimizes weights on a fixed ETA support initialized from conditional modes. Custom support matrices are available through nm_est().")):null,
        e("div",{className:"lw-modal-section"},e("h4",null,"Estimation priors"),list(props.model.priors).length?e(React.Fragment,null,e(SimpleTable,{rows:list(props.model.priors),columns:["parameter","distribution","mean","sd","shape","rate"],className:"lw-active-priors"}),e("p",{className:"lw-help-text"},list(props.model.priors).length+" prior"+(list(props.model.priors).length===1?" is":"s are")+" active for this run. Edit priors in the Code tab before submitting to change them.")):e("p",{className:"lw-help-text"},"No parameter priors are active. Add them under Estimation priors in the Code tab to make them part of the reproducible model version.")),
        e("div",{className:"lw-modal-section lw-modal-section-tinted"},e("h4",null,["BAYES","HMC","NUTS"].indexOf(estimationMethod[0])>=0?"Posterior uncertainty":"Covariance step"),["BAYES","HMC","NUTS"].indexOf(estimationMethod[0])>=0?e("p",{className:"lw-help-text"},"Posterior SDs, posterior CVs and 95% credible intervals are calculated from the saved samples automatically."):e(React.Fragment,null,e("label",{className:"lw-check"},e("input",{type:"checkbox",checked:covarianceStep[0]&&covarianceSupported,disabled:!covarianceSupported,onChange:function(event){covarianceStep[1](event.target.checked);}})," Run covariance step after estimation"),!covarianceSupported?e("p",{className:"lw-help-text"},["NPML","NPAG"].indexOf(estimationMethod[0])>=0?"Nonparametric support uncertainty should be assessed by bootstrap; the regular Hessian covariance step is not applicable.":"Covariance is available for FO, FOCE, FOCEI, LAPLACE, ITS, GQ, IMP and SAEM."):covarianceStep[0]?e(React.Fragment,null,e("div",{className:"lw-form-grid"},e(Field,{label:"Estimator"},e("select",{value:estimationMethod[0]==="FO"?"hessian":covarianceType[0],onChange:function(event){covarianceType[1](event.target.value);}},e("option",{value:"hessian"},"Hessian (R matrix)"),estimationMethod[0]!=="FO"?e("option",{value:"opg"},"Gradient outer product (S matrix)"):null)),e(Field,{label:"Regularization tolerance"},e("input",{type:"number",min:1e-14,step:"any",value:covarianceTolerance[0],onChange:function(event){covarianceTolerance[1](Number(event.target.value));}})),estimationMethod[0]==="SAEM"?e(Field,{label:"Marginal samples"},e("input",{type:"number",min:5,value:covarianceSamples[0],onChange:function(event){covarianceSamples[1](Number(event.target.value));}})):null),e("p",{className:"lw-help-text"},estimationMethod[0]==="GQ"?"The covariance step reuses the fitted quadrature grid and its adaptive or fixed proposal.":estimationMethod[0]==="IMP"?"The covariance calculation uses deterministic Gauss-Hermite integration when feasible and otherwise reuses the IMP sample budget and seed.":estimationMethod[0]==="SAEM"?"Observed marginal information uses deterministic Gauss-Hermite integration when feasible, with common-random-number importance sampling as the high-dimensional fallback.":"Standard errors, RSEs, covariance and correlation matrices will be saved with the estimation run.")):e("p",{className:"lw-help-text"},"Enable this to calculate parameter uncertainty after the fit."))),
        e("p",{className:"lw-help-text"},"All likelihood, automatic differentiation, ADVAN, matrix-exponential and ODE calculations run in the C++ engine. Queued runs execute in an isolated worker.")),

      e(Modal,{open:simulationModal[0],className:"lw-modal-wide",onClose:function(){simulationModal[1](false);},title:"Create simulation - "+value(props.model.name,"model"),footer:e(React.Fragment,null,e(Button,{className:"lw-button-quiet",onClick:function(){simulationModal[1](false);}},"Cancel"),e(Button,{className:"lw-button-primary",onClick:submitSimulation},"Run simulation"))},
        e("div",{className:"lw-modal-section"},e("div",{className:"lw-form-grid"},e(Field,{label:"Run on"},e("select",{value:value(props.server.queue_id,"local"),onChange:function(event){emit(props,"queue_select",{id:event.target.value});}},list(props.server.queues).map(function(queue){return e("option",{key:queue.id,value:queue.id},queue.name);}))),e(Field,{label:"Label"},e("input",{value:simulationLabel[0],onChange:function(event){simulationLabel[1](event.target.value);}})),e(Field,{label:"Random seed"},e("input",{type:"number",min:1,value:simulationSeed[0],onChange:function(event){simulationSeed[1](Number(event.target.value));}})),e(Field,{label:"Parallel cores"},e("input",{type:"number",min:1,value:simulationCores[0],onChange:function(event){simulationCores[1](Number(event.target.value));}}))),e("div",{className:"lw-form-grid"},e(Field,{label:"Individuals"},e("input",{type:"number",min:1,max:10000,value:simulationSubjects[0],onChange:function(event){simulationSubjects[1](Number(event.target.value));}})),e(Field,{label:"Replications"},e("input",{type:"number",min:1,max:1000,value:simulationReplicates[0],onChange:function(event){simulationReplicates[1](Number(event.target.value));}})),e(Field,{label:"Days (TIME horizon)"},e("input",{type:"number",min:1,max:365,value:simulationDays[0],onChange:function(event){simulationDays[1](Number(event.target.value));}})))),
        e("div",{className:"lw-modal-section lw-modal-section-tinted"},e("h4",null,"Parameters"),e("label",{className:"lw-check"},e("input",{type:"checkbox",checked:simulationUseFit[0]&&props.fit.available,disabled:!props.fit.available,onChange:function(event){simulationUseFit[1](event.target.checked);}})," Use fitted THETA / OMEGA / SIGMA"),!props.fit.available?e("p",{className:"lw-help-text"},"No estimation is loaded; simulation will use the model's initial parameter values."):e("p",{className:"lw-help-text"},"Diagnostics are run separately from the selected estimation run.")),
        e("label",{className:"lw-check lw-design-toggle"},e("input",{type:"checkbox",checked:simulationUseDesign[0],onChange:function(event){simulationUseDesign[1](event.target.checked);}})," Custom dosing / sampling design"),
        simulationUseDesign[0]?e("div",{className:"lw-modal-section lw-simulation-design"},e("h4",null,"Dosing and sampling design"),e("div",{className:"lw-form-grid"},e(Field,{label:"Dosing"},e("select",{value:doseMode[0],onChange:function(event){doseMode[1](event.target.value);}},e("option",{value:"single"},"Single dose"),e("option",{value:"repeat"},"Repeat doses"),e("option",{value:"steady_state"},"Steady state"))),e(Field,{label:"Default dose amount"},e("input",{type:"number",min:0,value:doseAmount[0],onChange:function(event){doseAmount[1](Number(event.target.value));}})),e(Field,{label:"Dose CMT"},e("input",{type:"number",min:1,value:doseCmt[0],onChange:function(event){doseCmt[1](Number(event.target.value));}})),doseMode[0]==="repeat"?e(Field,{label:"Number of doses"},e("input",{type:"number",min:1,max:100,value:doseN[0],onChange:function(event){doseN[1](Number(event.target.value));}})):null,doseMode[0]!=="single"?e(Field,{label:"Dosing interval (h)"},e("input",{type:"number",min:.1,step:.5,value:doseII[0],onChange:function(event){doseII[1](Number(event.target.value));}})):null,e(Field,{label:"Observations / day"},e("input",{type:"number",min:3,max:48,value:obsPerDay[0],onChange:function(event){obsPerDay[1](Number(event.target.value));}}))),e(Field,{label:"Dose amounts (TIME AMT per line, or AMT only)"},e("textarea",{rows:3,value:doseTable[0],placeholder:"0 320\n12 320",onChange:function(event){doseTable[1](event.target.value);}}))):e("p",{className:"lw-help-text"},"The linked dataset structure is retained and resampled to the requested number of individuals.")),

      e(Modal,{open:diagnosticModal[0],className:"lw-modal-wide",onClose:function(){diagnosticModal[1](false);},title:"Run diagnostic",footer:e(React.Fragment,null,e(Button,{className:"lw-button-quiet",onClick:function(){diagnosticModal[1](false);}},"Cancel"),e(Button,{className:"lw-button-primary",disabled:!diagnosticVpc[0]&&!diagnosticNpc[0]&&!diagnosticNpde[0]&&!diagnosticCategorical[0]&&!diagnosticCount[0]&&!diagnosticTte[0]&&!diagnosticCompeting[0]&&!diagnosticRecurrent[0],onClick:submitDiagnostic},"Run selected"))},
        e("p",{className:"lw-help-text"},"Diagnostics are calculated for the selected estimation run and saved with it."),
        e("div",{className:"lw-choice-row"},
          e("label",{className:"lw-check"},e("input",{type:"checkbox",checked:diagnosticVpc[0],onChange:function(event){diagnosticVpc[1](event.target.checked);}})," Continuous VPC"),
          e("label",{className:"lw-check"},e("input",{type:"checkbox",checked:diagnosticCategorical[0],onChange:function(event){diagnosticCategorical[1](event.target.checked);}})," Categorical VPC"),
          outcomeFamilies.some(function(item){return ["poisson","negative_binomial","binomial","zero_inflated_poisson","hurdle_poisson"].indexOf(item)>=0;})?e("label",{className:"lw-check"},e("input",{type:"checkbox",checked:diagnosticCount[0],onChange:function(event){diagnosticCount[1](event.target.checked);}})," Count VPC"):null,
          e("label",{className:"lw-check"},e("input",{type:"checkbox",checked:diagnosticTte[0],onChange:function(event){diagnosticTte[1](event.target.checked);}})," Time-to-event VPC"),
          outcomeFamilies.indexOf("competing_risks")>=0?e("label",{className:"lw-check"},e("input",{type:"checkbox",checked:diagnosticCompeting[0],onChange:function(event){diagnosticCompeting[1](event.target.checked);}})," Competing-risk VPC"):null,
          outcomeFamilies.indexOf("recurrent_event")>=0?e("label",{className:"lw-check"},e("input",{type:"checkbox",checked:diagnosticRecurrent[0],onChange:function(event){diagnosticRecurrent[1](event.target.checked);}})," Recurrent-event VPC"):null,
          e("label",{className:"lw-check"},e("input",{type:"checkbox",checked:diagnosticNpde[0],onChange:function(event){diagnosticNpde[1](event.target.checked);}})," NPDE"),
          e("label",{className:"lw-check"},e("input",{type:"checkbox",checked:diagnosticNpc[0],onChange:function(event){diagnosticNpc[1](event.target.checked);}})," NPC")),
        e("div",{className:"lw-form-grid"},e(Field,{label:"Simulations"},e("input",{type:"number",min:20,max:10000,value:diagnosticNsim[0],onChange:function(event){diagnosticNsim[1](Number(event.target.value));}})),e(Field,{label:"Random seed"},e("input",{type:"number",min:1,value:diagnosticSeed[0],onChange:function(event){diagnosticSeed[1](Number(event.target.value));}}))),
        diagnosticVpc[0]?e("div",{className:"lw-modal-section lw-modal-section-tinted"},e("h4",null,"Continuous VPC options"),e("label",{className:"lw-check"},e("input",{type:"checkbox",checked:diagnosticPc[0],onChange:function(event){diagnosticPc[1](event.target.checked);}})," Prediction-corrected VPC (DV x PRED / IPRED)"),e(Field,{label:"Stratify by"},e("select",{value:diagnosticStratify[0],onChange:function(event){diagnosticStratify[1](event.target.value);}},[e("option",{key:"none",value:""},"(no stratification)")].concat(list(props.dataset.columns).map(function(column){return e("option",{key:column,value:column},column);})))),e("p",{className:"lw-help-text"},"The VPC tab retains the overall population plot and adds one saved plot per stratum.")):null,
        diagnosticCategorical[0]?e("div",{className:"lw-modal-section lw-modal-section-tinted"},e("h4",null,"Categorical VPC options"),e(Field,{label:"Outcome column"},e("select",{value:diagnosticOutcome[0],onChange:function(event){diagnosticOutcome[1](event.target.value);}},list(props.dataset.columns).map(function(column){return e("option",{key:column,value:column},column);}))),e("p",{className:"lw-help-text"},outcomeFamilies.some(function(item){return ["categorical","ordinal","markov"].indexOf(item)>=0;})?"Declared probability vectors are used for every category.":"Legacy binary models interpret F/IPRED as the non-reference probability.")):null,
        (diagnosticCount[0]||diagnosticCompeting[0]||diagnosticRecurrent[0])&&outcomeDvids.length>1?e(Field,{label:"Endpoint (DVID)"},e("select",{value:diagnosticDvid[0],onChange:function(event){diagnosticDvid[1](event.target.value);}},e("option",{value:""},"Select endpoint"),outcomeDvids.map(function(item){return e("option",{key:item.dvid,value:item.dvid},item.name+" (DVID "+item.dvid+")");}))):null,
        diagnosticTte[0]?e("div",{className:"lw-modal-section lw-modal-section-tinted"},e("h4",null,"Time-to-event VPC options"),e(Field,{label:"Event indicator column"},e("select",{value:diagnosticEvent[0],onChange:function(event){diagnosticEvent[1](event.target.value);}},list(props.dataset.columns).map(function(column){return e("option",{key:column,value:column},column);}))),e("p",{className:"lw-help-text"},"F/IPRED is interpreted as a non-negative hazard on the observation-time grid.")):null),

      e(Modal,{open:uncertaintyModal[0],className:"lw-modal-wide",onClose:function(){uncertaintyModal[1](false);},title:"Parameter uncertainty",footer:e(React.Fragment,null,e(Button,{className:"lw-button-quiet",onClick:function(){uncertaintyModal[1](false);}},"Cancel"),e(Button,{className:"lw-button-primary",disabled:!uncertaintyBootstrap[0]&&!uncertaintyProfile[0],onClick:submitUncertainty},"Run selected"))},
        e("div",{className:"lw-choice-row"},e("label",{className:"lw-check"},e("input",{type:"checkbox",checked:uncertaintyBootstrap[0],onChange:function(event){uncertaintyBootstrap[1](event.target.checked);}})," Subject bootstrap"),e("label",{className:"lw-check"},e("input",{type:"checkbox",checked:uncertaintyProfile[0],onChange:function(event){uncertaintyProfile[1](event.target.checked);}})," Profile likelihood")),
        e("div",{className:"lw-form-grid"},uncertaintyBootstrap[0]?e(Field,{label:"Bootstrap replicates"},e("input",{type:"number",min:1,value:uncertaintyReplicates[0],onChange:function(event){uncertaintyReplicates[1](Number(event.target.value));}})):null,e(Field,{label:"Confidence level"},e("input",{type:"number",min:0.5,max:0.999,step:0.01,value:uncertaintyLevel[0],onChange:function(event){uncertaintyLevel[1](Number(event.target.value));}})),e(Field,{label:"Maximum fit iterations"},e("input",{type:"number",min:1,value:uncertaintyMaxit[0],onChange:function(event){uncertaintyMaxit[1](Number(event.target.value));}}))),
        uncertaintyProfile[0]?e("div",{className:"lw-modal-section lw-modal-section-tinted"},e("div",{className:"lw-form-grid"},e(Field,{label:"Grid points / parameter"},e("input",{type:"number",min:3,step:2,value:uncertaintyPoints[0],onChange:function(event){uncertaintyPoints[1](Number(event.target.value));}})),e(Field,{label:"Grid half-width (SE)"},e("input",{type:"number",min:0.1,step:0.5,value:uncertaintySpan[0],onChange:function(event){uncertaintySpan[1](Number(event.target.value));}}))),e(Field,{label:"Parameters (blank = all free)"},e("input",{placeholder:"THETA1, THETA2, SIGMA1",value:uncertaintyParameters[0],onChange:function(event){uncertaintyParameters[1](event.target.value);}})),e("p",{className:"lw-help-text"},"Each grid point fixes one parameter and re-estimates the remaining free parameters.")):null),

      e(Modal,{open:scmModal[0],className:"lw-modal-wide",onClose:function(){scmModal[1](false);},title:"Stepwise covariate modelling",footer:e(React.Fragment,null,e(Button,{className:"lw-button-quiet",onClick:function(){scmModal[1](false);}},"Cancel"),e(Button,{className:"lw-button-primary",onClick:function(){emit(props,"run_scm",{candidates:scmCandidates[0],direction:scmDirection[0],pForward:Number(scmForward[0]),pBackward:Number(scmBackward[0]),maxSteps:Number(scmMaxSteps[0]),maxit:Number(scmMaxit[0]),label:scmLabel[0]});scmModal[1](false);}},"Run SCM"))},e(Field,{label:"Candidate relationships (parameter,covariate,form,reference,category)"},e("textarea",{rows:6,value:scmCandidates[0],onChange:function(event){scmCandidates[1](event.target.value);}})),e("div",{className:"lw-form-grid"},e(Field,{label:"Direction"},e("select",{value:scmDirection[0],onChange:function(event){scmDirection[1](event.target.value);}},e("option",{value:"forward"},"Forward"),e("option",{value:"backward"},"Backward"),e("option",{value:"both"},"Forward + backward"))),e(Field,{label:"Forward p"},e("input",{type:"number",min:0.0001,max:0.5,step:0.01,value:scmForward[0],onChange:function(event){scmForward[1](Number(event.target.value));}})),e(Field,{label:"Backward p"},e("input",{type:"number",min:0.0001,max:0.5,step:0.01,value:scmBackward[0],onChange:function(event){scmBackward[1](Number(event.target.value));}})),e(Field,{label:"Maximum steps"},e("input",{type:"number",min:1,value:scmMaxSteps[0],onChange:function(event){scmMaxSteps[1](Number(event.target.value));}})),e(Field,{label:"Fit iterations"},e("input",{type:"number",min:1,value:scmMaxit[0],onChange:function(event){scmMaxit[1](Number(event.target.value));}})),e(Field,{label:"New version label"},e("input",{value:scmLabel[0],onChange:function(event){scmLabel[1](event.target.value);}}))),e("p",{className:"lw-help-text"},"Forms are continuous, power, or categorical. The accepted SCM model is saved as a new version and estimation run.")),

      e(Modal,{open:controlModal[0],className:"lw-modal-wide",onClose:function(){controlModal[1](false);},title:"Load NONMEM control stream",footer:e(React.Fragment,null,e(Button,{className:"lw-button-quiet",onClick:function(){controlModal[1](false);}},"Cancel"),e(Button,{className:"lw-button-primary",disabled:!controlFile[0]||(controlNewProject[0]&&!controlProjectName[0].trim()),onClick:function(){emit(props,"control_import",{text:controlFile[0]&&controlFile[0].text,fileName:controlFile[0]&&controlFile[0].name,dataText:controlData[0]&&controlData[0].text,dataName:controlData[0]&&controlData[0].name,newProject:controlNewProject[0],projectName:controlProjectName[0],label:controlLabel[0]});controlModal[1](false);}},"Import"))},e(Field,{label:"Control stream (.ctl, .mod)"},e("input",{type:"file",accept:".ctl,.mod,.txt,text/plain",onChange:readControlFile})),e("p",{className:"lw-help-text"},controlFile[0]?"Loaded "+controlFile[0].name:"Unsupported records are preserved and reported instead of silently discarded."),e(Field,{label:"Dataset (optional; otherwise keep current dataset)"},e("input",{type:"file",accept:".csv,.txt,.dat,.tsv,text/csv,text/plain",onChange:readControlData})),e("label",{className:"lw-check"},e("input",{type:"checkbox",checked:controlNewProject[0],onChange:function(event){controlNewProject[1](event.target.checked);}})," Create a new project"),controlNewProject[0]?e(Field,{label:"Project name"},e("input",{value:controlProjectName[0],onChange:function(event){controlProjectName[1](event.target.value);}})):null,e(Field,{label:"Model version label"},e("input",{value:controlLabel[0],onChange:function(event){controlLabel[1](event.target.value);}}))),

      e(Modal,{open:libraryModal[0],className:"lw-modal-wide",onClose:function(){libraryModal[1](false);},title:"LibeRary model catalogue",footer:e(React.Fragment,null,e(Button,{className:"lw-button-quiet",onClick:function(){libraryModal[1](false);}},"Cancel"),e(Button,{className:"lw-button-primary",disabled:!libraryInfo.available||!librarySelected[0]||(libraryNewProject[0]&&!libraryProjectName[0].trim()),onClick:function(){emit(props,"library_import",{libraryId:librarySelected[0],newProject:libraryNewProject[0],projectName:libraryProjectName[0],label:libraryLabel[0]});libraryModal[1](false);}},"Import model"))},
        libraryInfo.available?e(React.Fragment,null,
          e("div",{className:"lw-library-toolbar"},e(Field,{label:"Search catalogue"},e("input",{value:libraryQuery[0],placeholder:"Compound, population, ADVAN, title...",onChange:function(event){libraryQuery[1](event.target.value);}})),e("span",{className:"lw-help-text"},libraryFiltered.length+" of "+libraryEntries.length+" models")),
          e("div",{className:"lw-library-list"},libraryFiltered.length?libraryFiltered.map(function(item){var selected=librarySelected[0]===item.library_id;return e("label",{key:item.library_id,className:"lw-library-entry "+(selected?"selected":"")},e("input",{type:"radio",name:"library-entry",checked:selected,onChange:function(){librarySelected[1](item.library_id);if(!libraryProjectName[0])libraryProjectName[1](String(item.library_id).indexOf("lib_")===0?item.library_id:"lib_"+item.library_id);if(!libraryLabel[0])libraryLabel[1](value(item.title,item.library_id));}}),e("span",{className:"lw-library-copy"},e("strong",null,value(item.title,item.library_id)),e("small",null,[value(item.compound,"Unspecified compound"),value(item.population,"Unspecified population"),item.advan?"ADVAN"+item.advan:"",value(item.status,"")].filter(Boolean).join(" · "))),e("span",{className:"lw-library-confidence"},number(item.confidence_overall)!==null?Math.round(Number(item.confidence_overall)*100)+"%":"—"));}):e(Empty,{title:"No matching models",detail:"Change the search terms or add models through LibeRary::ingest_shiny()."})),
          e("div",{className:"lw-modal-section lw-modal-section-tinted"},e("label",{className:"lw-check"},e("input",{type:"checkbox",checked:libraryNewProject[0],onChange:function(event){libraryNewProject[1](event.target.checked);}})," Create a new project"),libraryNewProject[0]?e(Field,{label:"Project name"},e("input",{value:libraryProjectName[0],onChange:function(event){libraryProjectName[1](event.target.value);}})):e("p",{className:"lw-help-text"},"The model will be added as a new version under the current project."),e(Field,{label:"Model version label"},e("input",{value:libraryLabel[0],onChange:function(event){libraryLabel[1](event.target.value);}})),e("p",{className:"lw-help-text"},"Catalogue evidence, assessment and qualification provenance are retained with the model version.")))
        :e(Empty,{title:"LibeRary is unavailable",detail:value(libraryInfo.message,"Install LibeRary to browse the pharmacometric model catalogue.")})),

      e(Modal,{open:exportModal[0],onClose:function(){exportModal[1](false);},title:"Export NONMEM control stream",footer:e(React.Fragment,null,e(Button,{className:"lw-button-quiet",onClick:function(){exportModal[1](false);}},"Cancel"),e(Button,{className:"lw-button-primary",onClick:function(){emit(props,"control_export",{name:exportName[0],dataPath:exportDataPath[0]});exportModal[1](false);}},"Export"))},e(Field,{label:"File name"},e("input",{value:exportName[0],onChange:function(event){exportName[1](event.target.value);}})),e(Field,{label:"$DATA path"},e("input",{value:exportDataPath[0],onChange:function(event){exportDataPath[1](event.target.value);}})),e("p",{className:"lw-help-text"},"The file is written to the workspace exports directory. Preserved records from an imported stream are retained.")),

      e(Modal,{open:!!deleteModal[0],onClose:function(){deleteConfirmation[1]("");deleteModal[1](null);},title:deleteModal[0]==="project"?"Delete project":deleteModal[0]==="run"?"Delete model run":"Delete model version",footer:e(React.Fragment,null,e(Button,{className:"lw-button-quiet",onClick:function(){deleteConfirmation[1]("");deleteModal[1](null);}},"Cancel"),e(Button,{className:"lw-button-danger",disabled:deleteModal[0]==="project"&&deleteConfirmation[0]!=="YES",onClick:function(){if(deleteModal[0]==="project")emit(props,"project_delete",{id:workspace.current,confirmation:deleteConfirmation[0]});else if(deleteModal[0]==="run")emit(props,"project_delete_run",{id:workspace.current,run:workspace.current_run});else emit(props,"project_delete_snapshot",{id:workspace.current,snapshot:workspace.current_version});deleteConfirmation[1]("");deleteModal[1](null);}},"Delete"))},e("div",{className:"lw-destructive-note"},e("strong",null,"This action cannot be undone."),e("p",null,deleteModal[0]==="project"?"The project and all of its model versions and runs will be removed.":deleteModal[0]==="run"?"The selected estimation or simulation run and its saved diagnostics will be removed.":"The selected model version and all of its runs will be removed.")),deleteModal[0]==="project"?e("div",{className:"lw-modal-section"},e(Field,{label:'Type "YES" to confirm'},e("input",{autoFocus:true,value:deleteConfirmation[0],autoComplete:"off",onChange:function(event){deleteConfirmation[1](event.target.value);}})),e("p",{className:"lw-help-text"},"Confirmation is case-sensitive.")):null));
  }

  function AIStatusLine(props) {
    var status=useLocalAIStatus(),percent=Math.round((Number(status.progress)||0)*100),budget=status.budget;
    return e("div",{className:"lw-ai-status lw-ai-status-"+status.stage},
      e("div",null,e("span",null,status.text),status.locked?e("strong",null,"Network locked"):null),
      budget?e("small",{className:budget.compacted?"lw-ai-budget-compact":""},"Approx. "+Number(budget.prompt_tokens_estimated||0).toLocaleString()+" prompt tokens / "+Number(budget.context_window_size||0).toLocaleString()+" context; "+budget.retained_message_count+" of "+budget.original_message_count+" messages retained; "+Number(budget.max_tokens||0).toLocaleString()+" output tokens reserved"+(budget.compacted?" (context compacted)":"")):null,
      status.stage==="loading"?e("progress",{max:100,value:percent},percent+"%"):null,
      status.stage==="error"?e(Button,{className:"lw-button-quiet",onClick:function(){localAIShutdown();}},"Reset local AI"):null);
  }

  function AIModelSelect(props) {
    var ai=props.ai||{},models=list(ai.models),className=value(props.className,""),purpose=props.purpose==="report"?"report":"help";
    var configured=purpose==="report"?value(ai.report_model,"same_as_help"):value(ai.help_model,ai.model),resolved=configured==="same_as_help"?value(ai.help_model,ai.model):configured;
    var selected=models.filter(function(model){return model.id===resolved;})[0];
    function setModel(event){var detail=localAISettingsDetail(ai);detail[purpose+"_model"]=event.target.value;localAIShutdown();emit(props,"ai_settings",detail);}
    if(!models.length)return null;
    var options=models.map(function(model){return e("option",{key:model.id,value:model.id},model.label);});
    if(purpose==="report")options=[e("option",{key:"same_as_help",value:"same_as_help"},"Same as Help model")].concat(options);
    var description=configured==="same_as_help"?"Uses the Help model and avoids a model switch.":selected&&selected.description;
    return e("label",{className:"lw-ai-model-select "+value(className,""),title:description||"Choose the browser-local language model"},
      e("span",null,props.label||((purpose==="report"?"Report":"Help")+" LLM")),
      e("select",{value:configured,onChange:setModel,"aria-label":purpose==="report"?"Report builder language model":"Help language model"},options),
      description&&className==="lw-ai-model-panel"?e("small",null,description):null);
  }

  function AIContextSelect(props) {
    var ai=props.ai||{},purpose=props.purpose==="report"?"report":"help",model=purpose==="report"?value(ai.report_model,"same_as_help"):value(ai.help_model,ai.model);
    if(model==="same_as_help")model=value(ai.help_model,ai.model);
    var configured=String(value(ai[purpose+"_context"],"auto")),presetStrings=localAIContextPresets.map(String),isPreset=configured==="auto"||presetStrings.indexOf(configured)>=0;
    var mode=useSynced(isPreset?configured:"custom",[configured]),custom=useSynced(isPreset?8192:Number(configured)||8192,[configured]);
    var resolved=localAIContextWindow(ai,model,purpose),className=value(props.className,"");
    function commit(next){var detail=localAISettingsDetail(ai);detail[purpose+"_context"]=String(next);localAIShutdown();emit(props,"ai_settings",detail);}
    function change(event){var next=event.target.value;mode[1](next);if(next!=="custom")commit(next);}
    function commitCustom(){var next=Math.max(1024,Math.min(16384,Math.round((Number(custom[0])||8192)/512)*512));custom[1](next);commit(next);}
    var warning=resolved>8192?"Larger contexts use substantially more GPU memory and may be unstable on an 8 GB GPU.":"Changing context reloads the model lazily on its next use.";
    return e("label",{className:"lw-ai-context-select "+className,title:warning},e("span",null,props.label||((purpose==="report"?"Report":"Help")+" context")),e("div",null,e("select",{value:mode[0],onChange:change,"aria-label":purpose+" AI context window"},e("option",{value:"auto"},"Auto ("+(resolved/1024).toFixed(resolved%1024?1:0)+"K)"),localAIContextPresets.map(function(size){return e("option",{key:size,value:String(size)},(size/1024).toFixed(size%1024?1:0)+"K"+(size>8192?" - higher memory":""));}),e("option",{value:"custom"},"Custom")),mode[0]==="custom"?e("input",{type:"number",min:1024,max:16384,step:512,value:custom[0],"aria-label":purpose+" custom context tokens",onChange:function(event){custom[1](event.target.value);},onBlur:commitCustom,onKeyDown:function(event){if(event.key==="Enter"){event.preventDefault();commitCustom();}}}):null),className==="lw-ai-context-panel"?e("small",null,warning):null);
  }

  function helpQuestionScope(question) {
    var text=String(question||"").toLowerCase();
    var results=/(?:estimate|result|objective|\bofv\b|converg|covariance|\bse\b|\brse\b|\bgof\b|\bvpc\b|\bnpde\b|\bnpc\b|cwres|residual|diagnostic|compare|difference|best run|successful|failed|timing|iteration|posthoc|posterior|support point)/.test(text);
    var model=/(?:\$(?:pk|pred|des|error)|\bcode\b|equation|advan|trans|compartment|\btheta\b|\bomega\b|\bsigma\b|\beta\b|structural model|error model|parameter definition)/.test(text);
    return results&&model?"full":results?"results":model?"model":"overview";
  }
  function helpRunLine(run, detailed) {
    var parts=[value(run.model_version,"Model"),value(run.label,run.id||"Run")].join(" / ");
    var facts=[value(run.result_type,"run")];
    if(run.method)facts.push("method="+run.method);
    if(run.objective!==null&&run.objective!==undefined)facts.push("OFV="+formatNumber(run.objective));
    if(run.convergence!==null&&run.convergence!==undefined)facts.push("convergence="+run.convergence);
    if(run.iterations!==null&&run.iterations!==undefined&&!isNaN(Number(run.iterations)))facts.push("iterations="+run.iterations);
    var diagnostics=Object.keys(run.diagnostics||{}).filter(function(name){return !!run.diagnostics[name];});
    if(diagnostics.length)facts.push("diagnostics="+diagnostics.join(","));
    if(detailed&&list(run.parameters).length){
      facts.push("parameters="+list(run.parameters).map(function(parameter){return value(parameter.name,"parameter")+"="+formatNumber(parameter.estimate);}).join(","));
    }
    return "- "+parts+" ("+facts.join("; ")+")";
  }
  function helpProjectEvidence(projectEvidence, detailed) {
    if(!projectEvidence)return "Saved project evidence: unavailable in this view.";
    var heading="Saved project evidence: "+value(projectEvidence.project_name,projectEvidence.project||"project")+" has "+value(projectEvidence.run_count,0)+" completed run(s). "+value(projectEvidence.included_runs,0)+" are represented"+(Number(projectEvidence.omitted_runs)>0?"; "+projectEvidence.omitted_runs+" older run(s) are omitted.":".");
    var rows=list(projectEvidence.runs).map(function(run){return helpRunLine(run,detailed);});
    return heading+(rows.length?"\n"+rows.join("\n"):" "+value(projectEvidence.message,"No completed run summaries are available."));
  }

  function AIHelpPanel(props) {
    var ai=props.ai||{},welcome={role:"assistant",content:"Ask me about the current model, estimation workflow, or diagnostics. I run locally in this browser."};
    var messages=React.useState([welcome]),prompt=React.useState(""),busy=React.useState(false),error=React.useState("");
    var workspace=props.workspace||{},projectId=value(workspace.current,""),versionId=value(workspace.current_version,""),runId=value(workspace.current_run,"");
    var project=list(workspace.projects).filter(function(item){return String(item.id)===String(projectId);})[0];
    var version=list(workspace.versions).filter(function(item){return String(item.id)===String(versionId);})[0];
    var run=version?list(version.runs).filter(function(item){return String(item.id)===String(runId);})[0]:null;
    var aiContext=props.ai_context||{},contextKey=[projectId,versionId,runId].join("|");
    var contextRef=React.useRef(contextKey),pendingContext=React.useRef(null);
    React.useEffect(function(){
      contextRef.current=contextKey;
      pendingContext.current=null;
      if(localAIHasPendingPurpose("help"))localAICancelPurpose("help",new Error("The Help request was stopped because the selected project context changed."));
      messages[1]([welcome]);prompt[1]("");busy[1](false);error[1]("");
    },[contextKey]);
    React.useEffect(function(){
      var waiting=pendingContext.current,requestId=String(aiContext.request_id||"");
      if(!waiting||!requestId||requestId!==waiting.requestId)return;
      pendingContext.current=null;
      if(contextRef.current!==waiting.sentContext)return;
      runGeneration(waiting.question,waiting.history,waiting.sentContext,aiContext);
    },[String(aiContext.request_id||"")]);
    function runGeneration(question,history,sentContext,projectEvidence){
      var scope=helpQuestionScope(question),needsResults=scope==="results"||scope==="full",needsModel=scope==="model"||scope==="full";
      var selection=projectId?["Selected project: "+value(project&&project.name,project&&project.label||projectId)+" ["+projectId+"]","Selected model version: "+value(version&&version.label,versionId||"none")+(versionId?" ["+versionId+"]":""),runId?"Selected model run: "+value(run&&run.label,runId)+" ["+runId+"]":"Selected model run: none"]:["Selected project: none. The displayed model is not currently associated with an opened project in this Help context.","Selected model version: none","Selected model run: none"];
      var dataset=props.dataset||{};
      var compactProject=projectEvidence&&String(projectEvidence.project||"")===String(projectId)?{project_name:projectEvidence.project_name,run_count:projectEvidence.run_count,included_runs:projectEvidence.included_runs,omitted_runs:projectEvidence.omitted_runs,message:projectEvidence.message,runs:list(projectEvidence.runs)}:null;
      var selectedFit=props.fit&&props.fit.available?{method:props.fit.method,method_sequence:props.fit.method_sequence,objective:props.fit.objective,convergence:props.fit.convergence,parameters:props.fit.parameters,covariance:props.fit.covariance,run_info:props.fit.run_info}:null;
      var context=["You are LibeRation's browser-local pharmacometric modelling assistant.","Evidence rules: use the supplied project context as the only source for claims about this model, dataset, run, or result. If a requested fact is absent, say exactly that it is not available in the supplied context. Never invent or approximate parameter values, run results, diagnostics, validation status, dataset characteristics, or code behaviour. Clearly label general PK/PD knowledge as general guidance rather than a fact about the current project. State uncertainty and ask the user to verify consequential modelling decisions. You have no tools and no network access.","Evidence scope selected for this question: "+scope].concat(selection,["Displayed model: "+value(props.model&&props.model.name,"none")+"; ADVAN/TRANS "+value(props.model&&props.model.advan,"-")+"/"+value(props.model&&props.model.trans,"-"),dataset.loaded?"Dataset metadata: "+value(dataset.name,"Current dataset")+"; "+value(dataset.records,0)+" records, "+value(dataset.subjects,0)+" subjects, "+value(dataset.observations,0)+" observations; columns "+list(dataset.columns).join(", "):"Dataset metadata: none loaded",helpProjectEvidence(compactProject,needsResults)]);
      if(needsResults)context.push(selectedFit?"Selected fit detail: "+localAIClip(JSON.stringify(selectedFit),1600):"Selected fit detail: no estimation run is selected.","Selected-run diagnostic availability: "+JSON.stringify(props.diagnostics&&props.diagnostics.available||{}));
      if(needsModel)context=context.concat(["THETA definitions: "+localAIClip(JSON.stringify(list(props.model&&props.model.theta)),800),"OMEGA definitions: "+localAIClip(JSON.stringify(list(props.model&&props.model.omega)),800),"SIGMA definitions: "+localAIClip(JSON.stringify(list(props.model&&props.model.sigma)),600),"$PK/$PRED:\n"+localAIClip(value(props.model&&props.model.pred,""),1200),"$DES:\n"+localAIClip(value(props.model&&props.model.des,""),1200),"$ERROR:\n"+localAIClip(value(props.model&&props.model.error,""),800)]);
      context=context.join("\n\n");
      localAIGenerate(ai,[{role:"system",content:context}].concat(history,[{role:"user",content:question}]),function(token){if(contextRef.current!==sentContext)return;messages[1](function(current){var next=current.slice(),last=Object.assign({},next[next.length-1]);last.content+=token;next[next.length-1]=last;return next;});},{purpose:"help",max_tokens:1000,temperature:.1,top_p:.8}).then(function(answer){if(contextRef.current!==sentContext)return;messages[1](function(current){var next=current.slice(),last=Object.assign({},next[next.length-1]);if(!last.content)last.content=answer||"No response was generated.";next[next.length-1]=last;return next;});busy[1](false);}).catch(function(reason){if(contextRef.current!==sentContext)return;messages[1](function(current){var next=current.slice(),last=Object.assign({},next[next.length-1]);if(!last.content)last.content="Generation stopped before a response was produced.";last.failed=true;next[next.length-1]=last;return next;});busy[1](false);error[1](reason.message);});
    }
    function send(){var question=prompt[0].trim();if(!question||busy[0]||!ai.activated)return;
      var sentContext=contextKey,scope=helpQuestionScope(question),projectScope=scope==="results"||scope==="full"?"results":"index",history=messages[0].filter(function(item,index){return index>0&&!item.failed;}).map(function(item){return {role:item.role,content:item.content};});
      var conversation=messages[0].concat([{role:"user",content:question},{role:"assistant",content:""}]);messages[1](conversation);prompt[1]("");busy[1](true);error[1]("");
      var contextReady=projectId&&String(aiContext.project||"")===String(projectId)&&!!aiContext.request_id&&(String(aiContext.scope||"results")===projectScope||String(aiContext.scope||"")==="results");
      if(projectId&&!contextReady){
        var requestId="ai-context-"+Date.now()+"-"+Math.random().toString(16).slice(2);
        pendingContext.current={requestId:requestId,question:question,history:history,sentContext:sentContext};
        if(emit(props,"ai_context_request",{project:projectId,requestId:requestId,scope:projectScope}))return;
        pendingContext.current=null;
      }
      runGeneration(question,history,sentContext,contextReady?aiContext:null);
    }
    if(!ai.activated)return e(Empty,{title:"Local AI is off",detail:"Use Activate AI in the header to enable browser-local help."});
    return e("div",{className:"lw-ai-help"},
      e("div",{className:"lw-ai-toolbar"},e("div",{className:"lw-ai-config-row"},e(AIModelSelect,Object.assign({},props,{className:"lw-ai-model-panel",purpose:"help",label:"Help model"})),e(AIContextSelect,Object.assign({},props,{className:"lw-ai-context-panel",purpose:"help",label:"Context"}))),e(AIStatusLine,props)),
      e("div",{className:"lw-ai-messages"},messages[0].map(function(item,index){return e("div",{key:index,className:"lw-ai-message lw-ai-"+item.role+(item.failed?" lw-ai-failed":"")},e("strong",null,item.role==="user"?"You":"Local AI"),e("div",null,item.content||e("span",{className:"lw-ai-cursor"},busy[0]&&index===messages[0].length-1?"Generating...":"No response was generated.")));})),
      error[0]?e("div",{className:"lw-destructive-note"},error[0]):null,
      e("div",{className:"lw-ai-compose"},e("textarea",{rows:3,value:prompt[0],placeholder:"Ask about the model or workflow...",onChange:function(event){prompt[1](event.target.value);},onKeyDown:function(event){if(event.key==="Enter"&&!event.shiftKey){event.preventDefault();send();}}}),e(Button,{className:"lw-button-primary",disabled:busy[0]||!prompt[0].trim(),onClick:send},busy[0]?"Working...":"Send")),
      e("p",{className:"lw-help-text"},"Model details and compact saved-run summaries are loaded only when Help needs them. Row-level result data are excluded. The worker has no tools, DOM access, or network capability during inference."));
  }

  var reportBlockLabels={title:"Title",introduction:"Introduction",methods:"Methods",run:"Model run",comparison:"Model comparison",discussion:"Discussion",conclusion:"Conclusion",appendix:"Appendix",text:"Text",page_break:"Page break"};
  function reportBlock(type){return {id:"block-"+Date.now()+"-"+Math.random().toString(16).slice(2),type:type,title:reportBlockLabels[type],source:["run","comparison"].indexOf(type)>=0?"run":"user",text:"",run_ids:[],elements:["run","comparison"].indexOf(type)>=0?["summary","parameters","gof"]:[],options:{instruction:"",template:"",source_name:"",source_text:""}};}
  function reportRuns(workspace){var rows=[];list(workspace&&workspace.versions).forEach(function(version){list(version.runs).forEach(function(run){if(!run.queued_job)rows.push({id:run.id,label:value(run.label,run.id),version:value(version.label,version.id),type:run.result_type});});});return rows;}
  function reportSelections(blocks){
    var selected={};
    list(blocks).forEach(function(block){
      list(block.run_ids).forEach(function(id){
        id=String(id);if(!selected[id])selected[id]={id:id,elements:[]};
        list(block.elements).forEach(function(element){if(selected[id].elements.indexOf(element)<0)selected[id].elements.push(element);});
      });
    });
    return Object.keys(selected).map(function(id){if(!selected[id].elements.length)selected[id].elements=["summary","parameters"];return selected[id];});
  }
  function reportSelectionKey(selections){return list(selections).map(function(item){return item.id+":"+list(item.elements).slice().sort().join(",");}).sort().join("|");}
  function reportEvidenceText(context,selections){
    if(!context||!context.available)return "Selected report-run evidence is unavailable: "+value(context&&context.message,"No saved runs were loaded.");
    var selected={};list(selections).forEach(function(item){selected[item.id]=list(item.elements);});
    var output=["Selected report evidence from "+value(context.project_name,context.project||"the project")+": "+list(context.runs).length+" completed run(s)."];
    list(context.runs).forEach(function(run){
      var elements=selected[String(run.id)]||["summary","parameters"],section=[helpRunLine(run,elements.indexOf("parameters")>=0)];
      if(elements.indexOf("summary")>=0){
        var model=run.model||{},data=run.data||{},timing=run.timing_seconds||{};
        section.push("Model: "+value(model.name,"unnamed")+"; ADVAN/TRANS "+value(model.advan,"-")+"/"+value(model.trans,"-")+"; solver="+value(model.solver,"-")+"; language="+value(model.language,"-")+"; OMEGA structure="+value(model.omega_structure,"-")+".");
        section.push("Data: "+value(data.records,0)+" records; "+value(data.subjects,0)+" subjects; "+value(data.observations,0)+" observations; columns="+list(data.columns).join(",")+".");
        section.push("Timing (seconds): fit="+formatNumber(timing.model_fit)+", covariance="+formatNumber(timing.covariance)+", total="+formatNumber(timing.total)+".");
      }
      if(elements.indexOf("parameters")>=0||elements.indexOf("covariance")>=0){
        section.push("Parameter estimates: "+list(run.parameters).map(function(parameter){var text=value(parameter.name,"parameter")+"="+formatNumber(parameter.estimate);if(parameter.se!==null&&parameter.se!==undefined)text+=" (SE "+formatNumber(parameter.se)+", RSE "+formatNumber(parameter.rse)+")";return text;}).join("; ")+". Covariance: "+JSON.stringify(run.covariance||{})+".");
      }
      if(elements.indexOf("gof")>=0)section.push("GOF summary: "+JSON.stringify(run.gof_summary||"not available")+".");
      var diagnosticNames=["vpc","vpc_categorical","vpc_count","vpc_tte","vpc_competing","vpc_recurrent","npde","npc"];
      diagnosticNames.forEach(function(name){if(elements.indexOf(name)>=0){var detail=list(run.diagnostic_details).filter(function(item){return item.type===name;})[0];section.push(name.toUpperCase()+": "+JSON.stringify(detail||"not available")+".");}});
      if(elements.indexOf("run_info")>=0)section.push("Run information: iterations="+value(run.iterations,"not available")+", sequence="+list(run.method_sequence).join(" -> ")+".");
      if(elements.indexOf("code")>=0){var modelCode=run.model||{};section.push("$PK/$PRED:\n"+localAIClip(modelCode.pred,900)+"\n$DES:\n"+localAIClip(modelCode.des,900)+"\n$ERROR:\n"+localAIClip(modelCode.error,600));}
      output.push(section.join("\n"));
    });
    return output.join("\n\n");
  }

  function ReportDesigner(props) {
    var open=props.open,onClose=props.onClose,initial=[reportBlock("introduction"),reportBlock("methods"),reportBlock("run"),reportBlock("discussion"),reportBlock("conclusion")];
    var blocks=React.useState(initial),selected=React.useState(null),title=React.useState("LibeRation modelling report"),name=React.useState("liberation-report"),directory=React.useState(value(props.report_directory,"")),formats=React.useState({docx:true,pdf:true}),dragged=React.useRef(null),runs=reportRuns(props.workspace),ai=props.ai||{},reportContext=props.report_ai_context||{},pendingDraft=React.useRef(null),drafting=React.useState("");
    function update(id,changes){blocks[1](function(current){return current.map(function(block){return block.id===id?Object.assign({},block,changes):block;});});}
    function drop(event,index){event.preventDefault();var type=event.dataTransfer.getData("liber/block-type"),id=event.dataTransfer.getData("liber/block-id"),next=blocks[0].slice();if(type){next.splice(index,0,reportBlock(type));}else if(id){var old=next.findIndex(function(block){return block.id===id;});if(old>=0){var item=next.splice(old,1)[0];if(old<index)index-=1;next.splice(index,0,item);}}blocks[1](next);}
    function payload(){return {id:"report-main",title:title[0],name:name[0],directory:directory[0],formats:Object.keys(formats[0]).filter(function(key){return formats[0][key];}),blocks:blocks[0]};}
    function readFile(block,event){var file=event.target.files&&event.target.files[0];if(!file)return;var ext=file.name.split(".").pop().toLowerCase();if(["txt","md","markdown","rmd"].indexOf(ext)>=0){var reader=new FileReader();reader.onload=function(){update(block.id,{text:String(reader.result||""),options:Object.assign({},block.options,{source_name:file.name,source_text:String(reader.result||"")})});};reader.readAsText(file);}else{var binary=new FileReader();binary.onload=function(){emit(props,"report_document",{blockId:block.id,name:file.name,data:String(binary.result||""),nonce:Date.now()});};binary.readAsDataURL(file);}}
    React.useEffect(function(){var handler=function(event){var result=event.detail||{};if(result.input_id&&result.input_id!==props.inputId)return;blocks[1](function(current){return current.map(function(block){return block.id===result.block_id?Object.assign({},block,{text:block.source==="user"?result.text:block.text,options:Object.assign({},block.options,{source_name:result.name,source_text:result.text})}):block;});});};window.addEventListener("liber-report-document",handler);return function(){window.removeEventListener("liber-report-document",handler);};},[props.inputId]);
    React.useEffect(function(){var handler=function(event){directory[1](value(event&&event.detail&&event.detail.path,directory[0]));};window.addEventListener("liber-report-directory",handler);return function(){window.removeEventListener("liber-report-directory",handler);};},[]);
    React.useEffect(function(){var design=props.report_design;if(!design||!list(design.blocks).length)return;blocks[1](list(design.blocks).map(function(block){return Object.assign(reportBlock(block.type),block,{run_ids:list(block.run_ids),elements:list(block.elements),options:Object.assign({instruction:"",template:"",source_name:"",source_text:""},block.options||{})});}));title[1](value(design.title,"LibeRation modelling report"));name[1](value(design.name,"liberation-report"));directory[1](value(design.directory,props.report_directory||""));var selectedFormats={docx:false,pdf:false};list(design.formats).forEach(function(format){selectedFormats[format]=true;});formats[1](selectedFormats);},[props.report_design&&props.report_design.updated]);
    React.useEffect(function(){if(!(props.report_design&&props.report_design.directory))directory[1](value(props.report_directory,""));},[props.report_directory]);
    React.useEffect(function(){var waiting=pendingDraft.current,requestId=String(reportContext.request_id||"");if(!waiting||!requestId||requestId!==waiting.requestId)return;pendingDraft.current=null;runDraft(waiting.block,waiting.selections,reportContext);},[String(reportContext.request_id||"")]);
    function runDraft(block,selections,context){
      if(list(selections).length&&!context.available){update(block.id,{text:"[AI drafting failed: "+value(context.message,"Selected run evidence could not be loaded.")+"]"});drafting[1]("");return;}
      var fallback=props.fit&&props.fit.available?"Currently open fit (used only because the workflow has no selected runs): "+props.fit.method+", objective "+formatNumber(props.fit.objective)+"; parameters "+JSON.stringify(props.fit.parameters):"The workflow has no selected model runs.";
      var evidence=["Draft the report section titled: "+block.title,"Synthesize a coherent account across every selected run. Compare methods, estimates, uncertainty and diagnostics where relevant to this section. Do not produce a generic checklist of missing facts. Mention an unavailable fact only when it is material to the requested section.","Report workflow: "+blocks[0].map(function(item){return item.title+" ["+item.type+"]";}).join(" -> "),"Instruction: "+value(block.options.instruction,"None"),"Template: "+value(block.options.template,"None"),"Source material: "+localAIClip(value(block.options.source_text,"None"),1200),list(selections).length?reportEvidenceText(context,selections):fallback].join("\n\n");
      update(block.id,{text:""});
      localAIGenerate(ai,[{role:"system",content:"You draft concise, connected pharmacometric report prose from the supplied evidence. Treat every listed saved run as available evidence. Integrate the evidence into a narrative; do not invent results and do not repeat a missing-information checklist unless explicitly requested."},{role:"user",content:evidence}],function(token){blocks[1](function(current){return current.map(function(item){return item.id===block.id?Object.assign({},item,{text:item.text+token}):item;});});},{purpose:"report",max_tokens:1800,temperature:.1,top_p:.8}).then(function(){drafting[1]("");}).catch(function(error){drafting[1]("");update(block.id,{text:"[AI drafting failed: "+error.message+"]"});});
    }
    function draft(block){
      if(!ai.activated||drafting[0])return;
      var selections=reportSelections(blocks[0]),runIds=selections.map(function(item){return item.id;}),projectId=String(props.workspace&&props.workspace.current||""),contextIds=list(reportContext.run_ids).map(String).slice().sort(),requestedIds=runIds.slice().sort();
      drafting[1](block.id);
      if(runIds.length&&(!reportContext.available||String(reportContext.project||"")!==projectId||contextIds.join("|")!==requestedIds.join("|"))){
        var requestId="report-ai-context-"+Date.now()+"-"+Math.random().toString(16).slice(2);pendingDraft.current={requestId:requestId,block:block,selections:selections};
        if(emit(props,"report_ai_context_request",{project:projectId,runs:runIds,requestId:requestId}))return;
        pendingDraft.current=null;
      }
      runDraft(block,selections,runIds.length?reportContext:null);
    }
    var current=blocks[0].filter(function(block){return block.id===selected[0];})[0],optionBody=null;
    if(current){
      var common=e(Field,{label:"Heading"},e("input",{value:current.title,onChange:function(event){update(current.id,{title:event.target.value});}}));
      if(["run","comparison"].indexOf(current.type)>=0){
        optionBody=e("div",{className:"lw-form-stack"},common,e("strong",null,"Model runs"),runs.length?runs.map(function(run){return e("label",{className:"lw-check",key:run.id},e("input",{type:"checkbox",checked:current.run_ids.indexOf(run.id)>=0,onChange:function(event){var next=current.run_ids.filter(function(id){return id!==run.id;});if(event.target.checked)next.push(run.id);update(current.id,{run_ids:next});}})," "+run.version+" / "+run.label);}):e("p",{className:"lw-help-text"},"No completed runs are available."),e("strong",null,"Evidence"),["summary","parameters","code","gof","vpc","vpc_categorical","vpc_count","vpc_tte","vpc_competing","vpc_recurrent","npde","npc","covariance","run_info"].map(function(item){return e("label",{className:"lw-check",key:item},e("input",{type:"checkbox",checked:current.elements.indexOf(item)>=0,onChange:function(event){var next=current.elements.filter(function(value){return value!==item;});if(event.target.checked)next.push(item);update(current.id,{elements:next});}})," "+item.replace(/_/g," "));}));
      }else{
        optionBody=e("div",{className:"lw-form-stack"},common,current.source==="ai"?e(Field,{label:"AI instruction"},e("textarea",{rows:4,value:value(current.options.instruction,""),onChange:function(event){update(current.id,{options:Object.assign({},current.options,{instruction:event.target.value})});}})):null,e(Field,{label:current.source==="ai"?"Template text":"Section text"},e("textarea",{rows:8,value:current.source==="ai"?value(current.options.template,""):current.text,onChange:function(event){if(current.source==="ai")update(current.id,{options:Object.assign({},current.options,{template:event.target.value})});else update(current.id,{text:event.target.value});}})),e(Field,{label:"Source document (TXT, Markdown, DOCX, PDF)"},e("input",{type:"file",accept:".txt,.md,.docx,.pdf,text/plain,application/pdf,application/vnd.openxmlformats-officedocument.wordprocessingml.document",onChange:function(event){readFile(current,event);}})),current.options.source_name?e("p",{className:"lw-help-text"},"Loaded: "+current.options.source_name):null,current.source==="ai"?e(Button,{className:"lw-button-primary",disabled:!ai.activated||!!drafting[0],onClick:function(){draft(current);}},drafting[0]===current.id?"Loading evidence...":"Draft with local AI"):null);
      }
    }
    return e(Modal,{open:open,className:"lw-modal-report-designer",onClose:onClose,title:"Visual report designer",footer:e(React.Fragment,null,e(Button,{className:"lw-button-quiet",onClick:function(){emit(props,"report_design_save",payload());}},"Save workflow"),e(Button,{className:"lw-button-primary",disabled:!formats[0].docx&&!formats[0].pdf,onClick:function(){emit(props,"report_design_render",payload());}},"Generate DOCX / PDF"),e(Button,{className:"lw-button-quiet",onClick:onClose},"Close"))},
      e("div",{className:"lw-report-designer-toolbar"},e(Field,{label:"Report title"},e("input",{value:title[0],onChange:function(event){title[1](event.target.value);}})),e(Field,{label:"Filename"},e("input",{value:name[0],onChange:function(event){name[1](event.target.value);}})),e(Field,{label:"Save location",className:"lw-report-location"},e("div",null,e("input",{value:directory[0],placeholder:"Report output folder",onChange:function(event){directory[1](event.target.value);}}),e(Button,{className:"lw-button-quiet",title:"Choose report output folder",onClick:function(){emit(props,"report_directory_choose",{directory:directory[0]});}},"Browse..."))),e(AIModelSelect,Object.assign({},props,{className:"lw-ai-model-report",purpose:"report",label:"Report AI model"})),e(AIContextSelect,Object.assign({},props,{className:"lw-ai-context-report",purpose:"report",label:"Context"})),e("label",{className:"lw-check"},e("input",{type:"checkbox",checked:formats[0].docx,onChange:function(event){formats[1](Object.assign({},formats[0],{docx:event.target.checked}));}})," DOCX"),e("label",{className:"lw-check"},e("input",{type:"checkbox",checked:formats[0].pdf,onChange:function(event){formats[1](Object.assign({},formats[0],{pdf:event.target.checked}));}})," PDF")),
      e("div",{className:"lw-report-designer"},
        e("aside",{className:"lw-report-palette"},e("strong",null,"Blocks"),Object.keys(reportBlockLabels).map(function(type){return e("button",{type:"button",key:type,draggable:true,onDragStart:function(event){event.dataTransfer.setData("liber/block-type",type);}},reportBlockLabels[type]);}),e("p",{className:"lw-help-text"},"Drag blocks into the workflow. It remains a single top-to-bottom line.")),
        e("main",{className:"lw-report-lane",onDragOver:function(event){event.preventDefault();},onDrop:function(event){drop(event,blocks[0].length);}},blocks[0].map(function(block,index){return e(React.Fragment,{key:block.id},e("div",{className:"lw-report-drop",onDragOver:function(event){event.preventDefault();},onDrop:function(event){event.stopPropagation();drop(event,index);}}),e("article",{className:"lw-report-block",draggable:true,onDragStart:function(event){event.dataTransfer.setData("liber/block-id",block.id);}},e("span",{className:"lw-report-order"},index+1),e("div",null,e("strong",null,block.title),e("small",null,block.type.replace(/_/g," "))),!["run","comparison","page_break","title"].includes(block.type)?e("select",{value:block.source,onChange:function(event){update(block.id,{source:event.target.value});}},e("option",{value:"user"},"User generated"),e("option",{value:"ai"},"AI generated")):e("span",{className:"lw-report-source"},block.source),e(Button,{className:"lw-button-quiet",title:"Block details",onClick:function(){selected[1](block.id);}},"..."),block.source==="ai"?e(Button,{className:"lw-button-primary",disabled:!ai.activated||!!drafting[0],onClick:function(){draft(block);}},drafting[0]===block.id?"Loading...":"Draft"):null,e(Button,{className:"lw-button-link lw-prior-remove",title:"Remove block",onClick:function(){blocks[1](blocks[0].filter(function(item){return item.id!==block.id;}));}},"x"),block.text?e("p",null,block.text.slice(0,180)+(block.text.length>180?"...":"")):null));}))),
      e(Modal,{open:!!current,onClose:function(){selected[1](null);},title:current?current.title+" options":"Block options",footer:e(Button,{className:"lw-button-primary",onClick:function(){selected[1](null);}},"Done")},optionBody),
      e(AIStatusLine,props),props.report&&(props.report.docx||props.report.pdf)?e("div",{className:"lw-report-status"},e("strong",null,"Latest report"),props.report.docx?e("span",null,props.report.docx):null,props.report.pdf?e("span",null,props.report.pdf):null):null);
  }

  function ResultsPanel(props) {
    var tabState = React.useState("parameters"), tab = tabState[0];
    var reportDesigner = React.useState(false);
    var parameters = props.fit.available ? props.fit.parameters : [];
    var covariance = props.fit.covariance || {requested:false,status:"not_requested"};
    var posterior = props.fit.posterior || {available:false,parameters:[]};
    var nonparametric = props.fit.nonparametric || {available:false,supports:[]};
    var resultTabs = [{id:"parameters",label:"Parameters"}];
    if (covariance.requested) resultTabs.push({id:"covariance",label:"Covariance"});
    if (posterior.available) resultTabs.push({id:"posterior",label:"Posterior"});
    if (nonparametric.available) resultTabs.push({id:"support",label:"Support distribution"});
    resultTabs.push({id:"run",label:"Run info"},{id:"help",label:"Help"},{id:"report",label:"Report"});
    React.useEffect(function(){if(tab==="covariance"&&!covariance.requested)tabState[1]("parameters");},[covariance.requested]);
    React.useEffect(function(){if(tab==="posterior"&&!posterior.available)tabState[1]("parameters");},[posterior.available]);
    React.useEffect(function(){if(tab==="support"&&!nonparametric.available)tabState[1]("parameters");},[nonparametric.available]);
    return e(Panel, { title: "Results", className: "lw-results-panel", bodyClass: "lw-results-body" },
      e(Tabs, { value: tab, onChange: tabState[1], items: resultTabs }),
      tab === "parameters" ? e("div", { className: "lw-results-tab" },
        props.result && props.result.kind === "comparison" ? e(React.Fragment, null,
          e("div", { className: "lw-fit-summary" }, e("strong", null, "Run comparison"), e("span", null, "Side-by-side estimates")),
          e(SimpleTable, { rows: props.result.parameters })) : props.fit.available ? e(React.Fragment, null,
          e("div", { className: "lw-fit-summary" }, e("strong", null, props.fit.method + " fit"), e("span", null, "OFV " + formatNumber(props.fit.objective)), e("span", null, props.fit.convergence === 0 ? "Converged" : "Code " + props.fit.convergence)),
          e(SimpleTable, { rows: parameters, columns: posterior.available?["name","value","posterior_sd","posterior_cv","median","lower_95","upper_95"]:covariance.status==="completed"?["name","value","se","rse"]:["name","value"] })) : e(Empty, { title: "No estimates", detail: "Open or run an estimation." })) : null,
      tab === "covariance" ? e("div", { className: "lw-results-tab lw-covariance-tab" }, covariance.status==="failed"?e("div",{className:"lw-destructive-note"},e("strong",null,"Covariance step failed"),e("p",null,value(covariance.error,"No covariance result was produced."))):e(React.Fragment,null,e("div",{className:"lw-fit-summary"},e("strong",null,(covariance.type||"covariance").toUpperCase()+" covariance"),e("span",null,"Condition "+formatNumber(covariance.condition)),e("span",null,"Regularization "+formatNumber(covariance.regularization))),e("h4",null,"Covariance matrix"),e(SimpleTable,{rows:covariance.covariance,empty:"No free parameters"}),e("h4",null,"Correlation matrix"),e(SimpleTable,{rows:covariance.correlation,empty:"No correlations"}))) : null,
      tab === "posterior" ? e("div", { className: "lw-results-tab lw-covariance-tab" }, e("div",{className:"lw-fit-summary"},e("strong",null,"Bayesian posterior uncertainty"),e("span",null,formatNumber(posterior.samples)+" saved samples"),posterior.mean_acceptance!=null?e("span",null,"Acceptance "+formatNumber(posterior.mean_acceptance)):e("span",null,"Outer acceptance "+formatNumber(posterior.outer_acceptance)),posterior.divergences!=null?e("span",null,"Divergences "+formatNumber(posterior.divergences)):e("span",null,"ETA acceptance "+formatNumber(posterior.eta_acceptance))),e(SimpleTable,{rows:posterior.parameters,columns:["name","mean","posterior_sd","posterior_cv","median","lower_95","upper_95","rhat","ess"]}),e("h4",null,"Posterior covariance"),e(SimpleTable,{rows:posterior.covariance,empty:"No posterior covariance"}),e("h4",null,"Posterior correlation"),e(SimpleTable,{rows:posterior.correlation,empty:"No posterior correlations"})) : null,
      tab === "support" ? e("div", { className: "lw-results-tab lw-covariance-tab" },e("div",{className:"lw-fit-summary"},e("strong",null,props.fit.method+" discrete population distribution"),e("span",null,formatNumber(nonparametric.support_count)+" retained support points"),e("span",null,"Log-likelihood "+formatNumber(nonparametric.log_likelihood))),e("p",{className:"lw-help-text"},nonparametric.interpretation),e(SimpleTable,{rows:nonparametric.supports,empty:"No retained support points"})) : null,
      tab === "run" ? e("div", { className: "lw-results-tab lw-run-info" }, props.result && props.result.kind === "comparison" ? e(SimpleTable, { rows: props.result.runs }) : props.fit.available ?
        e(SimpleTable, { rows: Object.keys(props.fit.run_info || {}).map(function (key) { return { Item: key, Value: props.fit.run_info[key] }; }), columns: ["Item","Value"] }) : e(Empty, { title: "No run information", detail: "Run an estimation first." })) : null,
      tab === "help" ? e("div",{className:"lw-results-tab"},e(AIHelpPanel,props)) : null,
      tab === "report" ? e("div", { className: "lw-results-tab lw-report-controls" },
        e("p", null, "Build a top-to-bottom workflow from narrative, model-run, diagnostic, and comparison blocks. Narrative blocks can be user-authored or drafted by the browser-local AI."),
        e(Button,{className:"lw-button-primary",disabled:!props.workspace.current,onClick:function(){reportDesigner[1](true);}},"Open report designer"),
        !props.workspace.current?e("p",{className:"lw-help-text"},"Open a project to design a report."):null,
        props.report&&(props.report.docx||props.report.pdf)?e("div", { className: "lw-report-status" }, e("strong", null, "Report created"),props.report.docx?e("span", null, props.report.docx):null,props.report.pdf?e("span", null, props.report.pdf):null, props.report.json ? e("span", null, props.report.json) : null) : null) : null,
      e(ReportDesigner,Object.assign({},props,{open:reportDesigner[0],onClose:function(){reportDesigner[1](false);}})));
  }

  function ComparisonPlots(props) {
    var plots = props.plots || {};
    var labels = { gof:"Goodness-of-fit plots", vpc:"Visual predictive checks", vpc_categorical:"Categorical VPCs", vpc_count:"Count VPCs", vpc_tte:"Time-to-event VPCs", vpc_competing:"Competing-risk VPCs", vpc_recurrent:"Recurrent-event VPCs", npde:"NPDE plots", npc:"NPC plots" };
    var kinds = ["gof","vpc","vpc_categorical","vpc_count","vpc_tte","vpc_competing","vpc_recurrent","npde","npc"].filter(function(kind){return list(plots[kind]).length === 2;});
    if (!kinds.length) return null;
    return e("div", { className:"lw-comparison-plot-sections" }, kinds.map(function(kind){
      return e("section", { key:kind }, e("h4",null,labels[kind]), e("div",{className:"lw-comparison-plot-grid"},list(plots[kind]).map(function(item,index){
        var diagnosticProps = { tab:kind, fit:{available:false}, diagnostics:{} };
        if (kind === "gof") diagnosticProps.fit = item.fit || {available:false};
        else { diagnosticProps.diagnostics.available = {}; diagnosticProps.diagnostics.available[kind] = true; diagnosticProps.diagnostics[kind] = item.result; }
        return e("article",{className:"lw-comparison-plot-card",key:value(item.label,index)},e("h5",null,value(item.label,"Run "+(index+1))),e(DiagnosticsPane,diagnosticProps));
      })));
    }));
  }

  function diagramClone(value){return JSON.parse(JSON.stringify(value));}
  function diagramStarter(model){return {schema:"liber.model-diagram/1",version:1,title:value(model&&model.name,"Visual model"),advan:6,residual:"additive",covariates:[],compartments:[{id:1,name:"CENTRAL",kind:"amount",volume_parameter:"V",scale_parameter:"V",dose:true,observe:true,x:280,y:190}],flows:[{id:"flow-cl",from:1,to:0,type:"clearance",parameter:"CL",secondary_parameter:"",expression:"",label:"Elimination"}],parameters:[{name:"V",initial:20,lower:null,upper:null,fixed:false,iiv:true,eta_variance:.1},{name:"CL",initial:2,lower:null,upper:null,fixed:false,iiv:true,eta_variance:.1}]};}
  function diagramSymbols(graph){var names=[];function add(name){name=String(name||"").trim().toUpperCase();if(/^[A-Z][A-Z0-9_]*$/.test(name)&&names.indexOf(name)<0)names.push(name);}list(graph.compartments).forEach(function(item){add(item.volume_parameter);add(item.scale_parameter);});list(graph.flows).forEach(function(flow){add(flow.parameter);add(flow.secondary_parameter);if(flow.type==="custom")String(flow.expression||"").match(/\b[A-Za-z][A-Za-z0-9_]*\b/g)?.forEach(function(token){var upper=token.toUpperCase();if(["A","C","DADT","EXP","LOG","SQRT","TIME","T","IFELSE","MIN","MAX"].indexOf(upper)<0)add(upper);});});return names;}
  function diagramSyncParameters(graph){var required=diagramSymbols(graph),existing=list(graph.parameters),next=required.map(function(name){return existing.filter(function(item){return String(item.name).toUpperCase()===name;})[0]||{name:name,initial:1,lower:null,upper:null,fixed:false,iiv:true,eta_variance:.1};});existing.forEach(function(item){if(next.indexOf(item)<0&&required.indexOf(String(item.name).toUpperCase())<0)next.push(item);});graph.parameters=next;return graph;}
  function diagramRenameParameter(graph,oldName,newName){
    oldName=String(oldName||"").trim().toUpperCase();newName=String(newName||"").trim().toUpperCase();
    if(!oldName||oldName===newName||!/^[A-Z][A-Z0-9_]*$/.test(newName)||diagramSymbols(graph).indexOf(oldName)>=0)return graph;
    var oldIndex=list(graph.parameters).findIndex(function(item){return String(item.name||"").toUpperCase()===oldName;}),newIndex=list(graph.parameters).findIndex(function(item){return String(item.name||"").toUpperCase()===newName;});
    if(oldIndex<0)return graph;
    if(newIndex>=0)graph.parameters.splice(oldIndex,1);else graph.parameters[oldIndex].name=newName;
    return graph;
  }
  function VisualModelEditor(props){
    var supplied=props.model&&props.model.diagram&&props.model.diagram.available?props.model.diagram.graph:null;
    var graphState=React.useState(function(){return diagramClone(supplied||diagramStarter(props.model));}),graph=graphState[0],selected=React.useState(null),connect=React.useState(null),preview=React.useState(null),canvas=React.useRef(null),drag=React.useRef(null);
    React.useEffect(function(){if(supplied)graphState[1](diagramClone(supplied));},[props.model&&props.model.diagram&&props.model.diagram.graph&&props.model.diagram.graph.generated&&props.model.diagram.graph.generated.generated_at]);
    React.useEffect(function(){if(props.result&&props.result.kind==="diagram_preview")preview[1](props.result);},[props.result&&props.result.nonce]);
    function setGraph(next){graphState[1](diagramSyncParameters(diagramClone(next)));}
    function addCompartment(kind,x,y){var next=diagramClone(graph),id=Math.max.apply(null,[0].concat(next.compartments.map(function(item){return Number(item.id)||0;})))+1,name=kind==="response"?"RESPONSE"+id:"COMP"+id;next.compartments.push({id:id,name:name,kind:kind,volume_parameter:kind==="amount"?"V"+id:"",scale_parameter:kind==="amount"?"V"+id:"",dose:next.compartments.length===0,observe:kind==="response",x:x,y:y});setGraph(next);selected[1]("c:"+id);}
    function addFlow(from,to,type){var next=diagramClone(graph),id="flow-"+Date.now()+"-"+Math.random().toString(16).slice(2),parameter=type==="bidirectional_clearance"?"Q"+(next.flows.length+1):type==="michaelis_menten"?"VMAX"+(next.flows.length+1):type==="zero_order"?"K0"+(next.flows.length+1):type==="clearance"?"CL"+(next.flows.length+1):"K"+(next.flows.length+1);next.flows.push({id:id,from:from,to:to,type:type,parameter:parameter,secondary_parameter:type==="michaelis_menten"?"KM"+(next.flows.length+1):"",expression:type==="custom"?parameter+" * C("+from+")":"",label:""});setGraph(next);selected[1]("f:"+id);}
    function nodeClick(id){if(connect[0]&&connect[0].from){addFlow(connect[0].from,id,connect[0].type);connect[1](null);}else if(connect[0])connect[1](Object.assign({},connect[0],{from:id}));else selected[1]("c:"+id);}
    function pointerDown(event,item){if(connect[0])return;event.preventDefault();var rect=canvas.current.getBoundingClientRect();drag.current={id:item.id,dx:event.clientX-rect.left-Number(item.x),dy:event.clientY-rect.top-Number(item.y)};window.addEventListener("pointermove",pointerMove);window.addEventListener("pointerup",pointerUp,{once:true});}
    function pointerMove(event){if(!drag.current)return;var rect=canvas.current.getBoundingClientRect(),next=diagramClone(graph),item=next.compartments.filter(function(compartment){return compartment.id===drag.current.id;})[0];if(item){item.x=Math.max(60,Math.min(rect.width-60,event.clientX-rect.left-drag.current.dx));item.y=Math.max(45,Math.min(rect.height-45,event.clientY-rect.top-drag.current.dy));graphState[1](next);}}
    function pointerUp(){drag.current=null;window.removeEventListener("pointermove",pointerMove);}
    function dropPalette(event){event.preventDefault();var kind=event.dataTransfer.getData("liber/compartment"),rect=canvas.current.getBoundingClientRect();if(kind)addCompartment(kind,event.clientX-rect.left,event.clientY-rect.top);}
    function updateComp(id,field,value){var next=diagramClone(graph),item=next.compartments.filter(function(row){return row.id===id;})[0],oldValue=item[field];item[field]=value;if(field==="volume_parameter"&&(!item.scale_parameter||String(item.scale_parameter).startsWith("V")))item.scale_parameter=value;if(field==="volume_parameter"||field==="scale_parameter")diagramRenameParameter(next,oldValue,value);setGraph(next);}
    function updateFlow(id,field,value){var next=diagramClone(graph),item=next.flows.filter(function(row){return row.id===id;})[0],oldValue=item[field];item[field]=["from","to"].indexOf(field)>=0?Number(value):value;if(field==="parameter"||field==="secondary_parameter")diagramRenameParameter(next,oldValue,value);setGraph(next);}
    function renumberCompartment(oldId,newId){
      newId=Number(newId);var next=diagramClone(graph),count=next.compartments.length;if(!Number.isInteger(newId)||newId<1||newId>count||newId===oldId)return;
      var replacement=next.compartments.filter(function(item){return item.id===newId;})[0];next.compartments.forEach(function(item){if(item.id===oldId)item.id=newId;else if(replacement&&item.id===newId)item.id=oldId;});
      function swapId(id){return id===oldId?newId:(replacement&&id===newId?oldId:id);}
      next.flows.forEach(function(flow){flow.from=swapId(flow.from);flow.to=swapId(flow.to);if(flow.type==="custom")flow.expression=String(flow.expression||"").replace(/\b(C|A|DADT)\s*\(\s*([0-9]+)\s*\)/g,function(_,symbol,id){return symbol+"("+swapId(Number(id))+")";});});
      next.compartments.sort(function(a,b){return a.id-b.id;});selected[1]("c:"+newId);setGraph(next);
    }
    function updateParameter(index,field,value){var next=diagramClone(graph);next.parameters[index][field]=["initial","lower","upper","eta_variance"].indexOf(field)>=0?(value===""?null:Number(value)):value;graphState[1](next);}
    function removeParameter(index){var next=diagramClone(graph),item=next.parameters[index],required=diagramSymbols(next);if(!item||required.indexOf(String(item.name||"").toUpperCase())>=0)return;next.parameters.splice(index,1);graphState[1](next);}
    function removeSelected(){if(!selected[0])return;var next=diagramClone(graph);if(selected[0].startsWith("c:")){var id=Number(selected[0].slice(2));if(next.compartments.length<=1)return;next.compartments=next.compartments.filter(function(item){return item.id!==id;});next.flows=next.flows.filter(function(item){return item.from!==id&&item.to!==id;});}else{var flowId=selected[0].slice(2);next.flows=next.flows.filter(function(item){return item.id!==flowId;});}selected[1](null);setGraph(next);}
    var selectedComp=selected[0]&&selected[0].startsWith("c:")?graph.compartments.filter(function(item){return item.id===Number(selected[0].slice(2));})[0]:null,selectedFlow=selected[0]&&selected[0].startsWith("f:")?graph.flows.filter(function(item){return item.id===selected[0].slice(2);})[0]:null;
    function line(flow){var from=flow.from===0?{x:20,y:50}:graph.compartments.filter(function(item){return item.id===flow.from;})[0],to=flow.to===0?{x:760,y:50}:graph.compartments.filter(function(item){return item.id===flow.to;})[0];if(!from||!to)return null;var label=flow.type==="custom"?flow.expression:[flow.parameter,flow.secondary_parameter].filter(Boolean).join(" / ");return e("g",{key:flow.id,className:selected[0]==="f:"+flow.id?"selected":"",onClick:function(event){event.stopPropagation();selected[1]("f:"+flow.id);}},e("line",{x1:from.x,y1:from.y,x2:to.x,y2:to.y,markerEnd:"url(#lw-arrow)",markerStart:flow.type==="bidirectional_clearance"?"url(#lw-arrow-start)":null}),e("text",{x:(Number(from.x)+Number(to.x))/2,y:(Number(from.y)+Number(to.y))/2-7},label));}
    return e("div",{className:"lw-diagram-workspace"},
      e("div",{className:"lw-diagram-toolbar"},e("div",{className:"lw-diagram-palette"},[["amount","Amount compartment"],["response","Response compartment"]].map(function(item){return e("button",{key:item[0],type:"button",draggable:true,onDragStart:function(event){event.dataTransfer.setData("liber/compartment",item[0]);}},item[1]);})),e(Field,{label:"ADVAN"},e("select",{value:graph.advan,onChange:function(event){setGraph(Object.assign({},graph,{advan:Number(event.target.value)}));}},e("option",{value:6},"6 (general ODE)"),e("option",{value:13},"13 (stiff ODE)"))),e(Field,{label:"Residual error"},e("select",{value:graph.residual,onChange:function(event){setGraph(Object.assign({},graph,{residual:event.target.value}));}},e("option",{value:"additive"},"Additive"),e("option",{value:"proportional"},"Proportional"),e("option",{value:"combined"},"Combined"))),e(Button,{className:"lw-button-quiet",onClick:function(){connect[1]({type:"rate",from:null});}},"Connect"),e(Button,{className:"lw-button-quiet",onClick:function(){connect[1]({type:"bidirectional_clearance",from:null});}},"Bidirectional"),e(Button,{className:"lw-button-quiet",onClick:function(){connect[1]({type:"michaelis_menten",from:null});}},"Nonlinear"),e(Button,{className:"lw-button-quiet",disabled:!selectedComp,onClick:function(){if(selectedComp)addFlow(selectedComp.id,0,"clearance");}},"Elimination"),connect[0]?e("span",{className:"lw-diagram-mode"},connect[0].from?"Select target compartment":"Select source compartment"):null),
      e("div",{className:"lw-diagram-main"},e("div",{className:"lw-diagram-canvas",ref:canvas,onDragOver:function(event){event.preventDefault();},onDrop:dropPalette,onClick:function(){selected[1](null);}},e("svg",{className:"lw-diagram-links"},e("defs",null,e("marker",{id:"lw-arrow",viewBox:"0 0 10 10",refX:9,refY:5,markerWidth:7,markerHeight:7,orient:"auto-start-reverse"},e("path",{d:"M 0 0 L 10 5 L 0 10 z"})),e("marker",{id:"lw-arrow-start",viewBox:"0 0 10 10",refX:1,refY:5,markerWidth:7,markerHeight:7,orient:"auto-start-reverse"},e("path",{d:"M 10 0 L 0 5 L 10 10 z"}))),graph.flows.map(line)),graph.compartments.map(function(item){return e("button",{type:"button",key:item.id,className:"lw-diagram-node lw-diagram-"+item.kind+(selected[0]==="c:"+item.id?" selected":""),style:{left:item.x,top:item.y},onPointerDown:function(event){pointerDown(event,item);},onClick:function(event){event.stopPropagation();nodeClick(item.id);}},e("strong",null,item.name),e("small",null,item.kind==="amount"?(item.volume_parameter?"A / "+item.volume_parameter:"Amount"):"Response"),item.dose?e("i",null,"Dose"):null,item.observe?e("i",null,"Obs"):null);})),
        e("aside",{className:"lw-diagram-inspector"},selectedComp?e("div",{className:"lw-form-stack"},e("strong",null,"Compartment "+selectedComp.id),e(Field,{label:"Compartment number"},e("select",{value:selectedComp.id,onChange:function(event){renumberCompartment(selectedComp.id,event.target.value);}},graph.compartments.map(function(_,index){return e("option",{key:index+1,value:index+1},index+1);}))),e(Field,{label:"Name"},e("input",{value:selectedComp.name,onChange:function(event){updateComp(selectedComp.id,"name",event.target.value);}})),e(Field,{label:"Kind"},e("select",{value:selectedComp.kind,onChange:function(event){updateComp(selectedComp.id,"kind",event.target.value);}},e("option",{value:"amount"},"Amount"),e("option",{value:"response"},"Response"))),e(Field,{label:"Volume / scale parameter"},e("input",{value:selectedComp.volume_parameter,onChange:function(event){updateComp(selectedComp.id,"volume_parameter",event.target.value.toUpperCase());}})),e("label",{className:"lw-check"},e("input",{type:"checkbox",checked:!!selectedComp.dose,onChange:function(event){updateComp(selectedComp.id,"dose",event.target.checked);}})," Dose compartment"),e("label",{className:"lw-check"},e("input",{type:"checkbox",checked:!!selectedComp.observe,onChange:function(event){updateComp(selectedComp.id,"observe",event.target.checked);}})," Observation compartment"),e(Button,{className:"lw-button-danger",disabled:graph.compartments.length<=1,onClick:removeSelected},"Delete")):selectedFlow?e("div",{className:"lw-form-stack"},e("strong",null,"Flow"),e(Field,{label:"From"},e("select",{value:selectedFlow.from,onChange:function(event){updateFlow(selectedFlow.id,"from",event.target.value);}},[e("option",{key:0,value:0},"Outside")].concat(graph.compartments.map(function(item){return e("option",{key:item.id,value:item.id},item.name);})))),e(Field,{label:"To"},e("select",{value:selectedFlow.to,onChange:function(event){updateFlow(selectedFlow.id,"to",event.target.value);}},[e("option",{key:0,value:0},"Outside")].concat(graph.compartments.map(function(item){return e("option",{key:item.id,value:item.id},item.name);})))),e(Field,{label:"Flow type"},e("select",{value:selectedFlow.type,onChange:function(event){updateFlow(selectedFlow.id,"type",event.target.value);}},[["rate","First-order rate"],["clearance","Clearance"],["bidirectional_clearance","Bidirectional clearance"],["michaelis_menten","Michaelis-Menten"],["zero_order","Zero-order"],["custom","Custom nonlinear"]].map(function(item){return e("option",{key:item[0],value:item[0]},item[1]);}))),selectedFlow.type!=="custom"?e(Field,{label:selectedFlow.type==="michaelis_menten"?"VMAX":"Parameter"},e("input",{value:selectedFlow.parameter,onChange:function(event){updateFlow(selectedFlow.id,"parameter",event.target.value.toUpperCase());}})):null,selectedFlow.type==="michaelis_menten"?e(Field,{label:"KM"},e("input",{value:selectedFlow.secondary_parameter,onChange:function(event){updateFlow(selectedFlow.id,"secondary_parameter",event.target.value.toUpperCase());}})):null,selectedFlow.type==="custom"?e(Field,{label:"Flux expression (C(1) allowed)"},e("textarea",{rows:5,value:selectedFlow.expression,onChange:function(event){updateFlow(selectedFlow.id,"expression",event.target.value);}})):null,e(Button,{className:"lw-button-danger",onClick:removeSelected},"Delete")):e("p",{className:"lw-help-text"},"Drag compartments onto the canvas. Select a connection tool, then its source and target. Select any item to edit its model semantics."))),
      e("details",{className:"lw-diagram-parameters",open:true},
        e("summary",null,"Structural parameters - log-normal ETA by default"),
        e("div",{className:"lw-table-wrap"},
          e("table",{className:"lw-param-table"},
            e("thead",null,e("tr",null,["Parameter","Initial","Lower","Upper","IIV","OMEGA",""].map(function(label,index){
              return e("th",{key:label||"action-"+index},label);
            }))),
            e("tbody",null,graph.parameters.map(function(item,index){
              var required=diagramSymbols(graph).indexOf(String(item.name||"").toUpperCase())>=0;
              return e("tr",{key:item.name+"-"+index},
                e("td",null,item.name),
                ["initial","lower","upper"].map(function(field){
                  return e("td",{key:field},e("input",{type:"number",value:item[field]==null?"":item[field],onChange:function(event){updateParameter(index,field,event.target.value);}}));
                }),
                e("td",null,e("input",{type:"checkbox",checked:!!item.iiv,onChange:function(event){updateParameter(index,"iiv",event.target.checked);}})),
                e("td",null,e("input",{type:"number",min:.000001,step:.01,value:item.eta_variance,onChange:function(event){updateParameter(index,"eta_variance",event.target.value);}})),
                e("td",null,e(Button,{className:"lw-button-danger-ghost lw-diagram-param-remove",disabled:required,title:required?"This parameter is used by a compartment or flow. Remove or rename that reference first.":"Delete "+item.name,onClick:function(){removeParameter(index);}},"Delete"))
              );
            }))
          )
        )
      ),
      props.model&&props.model.diagram&&props.model.diagram.code_changed?e("div",{className:"lw-dirty-banner"},"The editable $DES differs from the last graph-generated version. Previewing is safe; applying the diagram will explicitly replace the current model code."):null,
      e("div",{className:"lw-inline-actions lw-editor-actions"},e(Button,{className:"lw-button-quiet",onClick:function(){graphState[1](diagramStarter(props.model));selected[1](null);}},"New diagram"),e(Button,{className:"lw-button-primary",onClick:function(){emit(props,"diagram_preview",{graph:diagramSyncParameters(diagramClone(graph)),nonce:Date.now()});}},"Preview generated code")),
      e(Modal,{open:!!preview[0],className:"lw-modal-wide",onClose:function(){preview[1](null);},title:"Apply diagram-generated model code?",footer:e(React.Fragment,null,e(Button,{className:"lw-button-quiet",onClick:function(){preview[1](null);}},"Keep current code"),e(Button,{className:"lw-button-primary",onClick:function(){emit(props,"diagram_apply",{graph:diagramSyncParameters(diagramClone(graph))});preview[1](null);}},"Apply generated code"))},preview[0]?e("div",{className:"lw-diagram-preview"},preview[0].code_changed?e("div",{className:"lw-destructive-note"},e("strong",null,"Manual code changes detected"),e("p",null,"Applying replaces $PK/$PRED, $DES, $ERROR, THETA, OMEGA and SIGMA. The visual graph and previous generated baseline remain stored with the model.")):null,e("div",{className:"lw-editor-grid lw-editor-grid-three"},e("div",{className:"lw-editor-box"},e("h5",null,"$PK / $PRED"),e("pre",null,preview[0].preview.pred)),e("div",{className:"lw-editor-box"},e("h5",null,"$DES"),e("pre",null,preview[0].preview.des)),e("div",{className:"lw-editor-box"},e("h5",null,"$ERROR"),e("pre",null,preview[0].preview.error))),e("p",{className:"lw-help-text"},"All generated structural parameters use PARAM = THETA(i) * exp(ETA(i)) unless IIV is unchecked. You can continue editing $DES in the Code tab after applying.")):null));
  }

  function HomePage(props) {
    var centerTab = React.useState("code"), tab = centerTab[0];
    var saveModal = React.useState(false), saveLabel = React.useState("Model revision");
    var comparisonModal = React.useState(false);
    var workspace = props.workspace || {};
    var diagnostics = props.diagnostics || {};
    var diagnosticAvailability = diagnostics.available || {};
    var centerTabs = [{id:"diagram",label:"Visual model"},{id:"code",label:"Code"},{id:"gof",label:"GOF"}];
    if (props.hmm&&props.hmm.available) centerTabs.push({id:"hmm",label:"HMM"});
    if (props.kalman&&props.kalman.available) centerTabs.push({id:"kalman",label:"States"});
    if (diagnosticAvailability.vpc) centerTabs.push({id:"vpc",label:"VPC"});
    if (diagnosticAvailability.vpc_categorical) centerTabs.push({id:"vpc_categorical",label:"Cat VPC"});
    if (diagnosticAvailability.vpc_count) centerTabs.push({id:"vpc_count",label:"Count VPC"});
    if (diagnosticAvailability.vpc_tte) centerTabs.push({id:"vpc_tte",label:"TTE VPC"});
    if (diagnosticAvailability.vpc_competing) centerTabs.push({id:"vpc_competing",label:"Risk VPC"});
    if (diagnosticAvailability.vpc_recurrent) centerTabs.push({id:"vpc_recurrent",label:"Recurrent VPC"});
    if (diagnosticAvailability.npde) centerTabs.push({id:"npde",label:"NPDE"});
    if (diagnosticAvailability.npc) centerTabs.push({id:"npc",label:"NPC"});
    if (diagnosticAvailability.bootstrap) centerTabs.push({id:"bootstrap",label:"Bootstrap"});
    if (diagnosticAvailability.profile) centerTabs.push({id:"profile",label:"Profile"});
    if (diagnosticAvailability.scm) centerTabs.push({id:"scm",label:"SCM"});
    React.useEffect(function(){if(!centerTabs.some(function(item){return item.id===centerTab[0];}))centerTab[1]("code");},[!!(props.hmm&&props.hmm.available),!!(props.kalman&&props.kalman.available),!!diagnosticAvailability.vpc,!!diagnosticAvailability.vpc_categorical,!!diagnosticAvailability.vpc_count,!!diagnosticAvailability.vpc_tte,!!diagnosticAvailability.vpc_competing,!!diagnosticAvailability.vpc_recurrent,!!diagnosticAvailability.npde,!!diagnosticAvailability.npc,!!diagnosticAvailability.bootstrap,!!diagnosticAvailability.profile,!!diagnosticAvailability.scm]);
    React.useEffect(function(){if(props.result&&props.result.kind==="comparison")comparisonModal[1](true);},[props.result&&props.result.comparison_id]);
    function closeComparison(){comparisonModal[1](false);emit(props,"comparison_close");}
    function selectCenterTab(next) {
      centerTab[1](next);
      if (next === "gof" && props.fit.available && !props.fit.gof_loaded) emit(props,"load_payload",{kind:"gof"});
      if (next === "hmm" && props.hmm&&props.hmm.available&&!props.hmm.loaded) emit(props,"load_payload",{kind:"hmm"});
      if (next === "kalman" && props.kalman&&props.kalman.available&&!props.kalman.loaded) emit(props,"load_payload",{kind:"kalman"});
      if (["vpc","npde","npc","vpc_categorical","vpc_count","vpc_tte","vpc_competing","vpc_recurrent","bootstrap","profile","scm"].indexOf(next) >= 0 && diagnosticAvailability[next] && !diagnostics[next]) emit(props,"load_payload",{kind:next});
    }
    var actions = e("div", { className: "lw-header-actions" },
      e(Button, { className: "lw-button-quiet", disabled: !workspace.current_version, onClick: function () { if(workspace.current_run)emit(props,"run_open",{id:workspace.current,run:workspace.current_run});else emit(props, "project_open", { id: workspace.current, snapshot: workspace.current_version }); } }, "Reload"),
      e(Button, { className: "lw-button-primary", disabled: !workspace.current, onClick: function () { emit(props, "project_save", { label: "" }); } }, "Save new version"),
      e(Button, { className: "lw-button-quiet", disabled: !workspace.current, onClick: function () { saveModal[1](true); } }, "Save as new"));
    return e("div", { className: "lw-home-grid" },
      e(ProjectTree, props),
      e("div", { className: "lw-center-column" },
        e(Panel, { title: props.model.loaded ? value(props.model.name, "Model version") : "Model version", subtitle: workspace.current_run ? "Selected run " + workspace.current_run : workspace.current_version ? "Selected model version" : "Unsaved", actions: actions, className: "lw-center-panel", bodyClass: "lw-center-body" },
          e(Tabs, { value: tab, onChange: selectCenterTab, items: centerTabs }),
          centerTabs.map(function(item){return e("div",{key:item.id,className:"lw-center-tab "+(tab===item.id?"":"lw-center-tab-hidden"),"aria-hidden":tab===item.id?"false":"true"},item.id==="diagram"?e(VisualModelEditor,props):item.id==="code"?e(ModelEditor,props):item.id==="hmm"?e(HmmPane,props):item.id==="kalman"?e(KalmanPane,props):e(CachedDiagnosticsPane,Object.assign({},props,{tab:item.id})));}))),
      e(ResultsPanel, props),
      e(Modal,{open:comparisonModal[0]&&props.result&&props.result.kind==="comparison",className:"lw-modal-comparison",onClose:closeComparison,title:"Compare estimation runs",footer:e(Button,{className:"lw-button-primary",onClick:closeComparison},"Close")},
        e("div",{className:"lw-comparison-content"},
          e("section",null,e("h4",null,"Run summary"),e(SimpleTable,{rows:list(props.result&&props.result.runs)})),
          e("section",null,e("h4",null,"Goodness-of-fit and information criteria"),e(SimpleTable,{rows:list(props.result&&props.result.gof)})),
          e("section",null,e("h4",null,"Parameter estimates"),e(SimpleTable,{rows:list(props.result&&props.result.parameters)})),
          e(ComparisonPlots,{plots:props.result&&props.result.plots}))),
      e(Modal, { open: saveModal[0], onClose: function () { saveModal[1](false); }, title: "Save model as new version", footer: e(React.Fragment, null,
        e(Button, { className: "lw-button-quiet", onClick: function () { saveModal[1](false); } }, "Cancel"),
        e(Button, { className: "lw-button-primary", onClick: function () { emit(props, "project_save", { label: saveLabel[0] }); saveModal[1](false); } }, "Save")) },
        e(Field, { label: "Version label" }, e("input", { value: saveLabel[0], onChange: function (event) { saveLabel[1](event.target.value); } }))));
  }

  function JobsPage(props) {
    var jobs = list(props.jobs), selected = React.useState(null), poll = React.useState(5), remoteModal = React.useState(false);
    var remoteName = React.useState("Remote server"), remoteUrl = React.useState("https://"), remoteToken = React.useState(""), remoteEditId = React.useState(null);
    var queues = list(props.server.queues), selectedQueue = queues.filter(function (queue) { return queue.id === props.server.queue_id; })[0] || {};
    React.useEffect(function () { var timer = window.setInterval(function () { emit(props, "jobs_refresh"); }, Math.max(1, poll[0]) * 1000); return function () { window.clearInterval(timer); }; }, [poll[0], props.inputId]);
    return e("div", { className: "lw-ribbon-page lw-jobs-page" },
      e("div", { className: "lw-jobs-toolbar" },
        e(Field, { label: "Queue" }, e("select", { value: value(props.server.queue_id, "local"), onChange: function (event) { emit(props, "queue_select", { id: event.target.value }); } }, queues.map(function (queue) { return e("option", { key: queue.id, value: queue.id }, queue.name); }))),
        e(Field, { label: "Hub poll interval (seconds)" }, e("input", { type: "number", min: 1, max: 120, value: poll[0], onChange: function (event) { poll[1](Number(event.target.value)); } })),
        e("div", { className: "lw-inline-actions" },
          e(Button, { className: "lw-button-quiet", onClick: function () { emit(props, "jobs_refresh"); } }, "Refresh"),
          e(Button, { className: "lw-button-warning", disabled: !selected[0], onClick: function () { emit(props, "job_cancel", { id: selected[0] }); } }, "Cancel selected"),
          e(Button, { className: "lw-button-quiet", onClick: function () { emit(props, "jobs_clear"); } }, "Clear finished"),
          e(Button, { className: "lw-button-quiet", onClick: function () { remoteEditId[1](null); remoteName[1]("Remote server"); remoteUrl[1]("https://"); remoteToken[1](""); remoteModal[1](true); } }, "Add server"),
          e(Button, { className: "lw-button-quiet", disabled: value(props.server.queue_id, "local") === "local", onClick: function () { remoteEditId[1](props.server.queue_id); remoteName[1](value(selectedQueue.name, "Remote server")); remoteUrl[1](value(selectedQueue.url, "https://")); remoteToken[1](""); remoteModal[1](true); } }, "Edit server"),
          e(Button, { className: "lw-button-danger", disabled: value(props.server.queue_id, "local") === "local", onClick: function () { emit(props, "queue_remove", { id: props.server.queue_id }); } }, "Remove server"))),
      e("div", { className: "lw-jobs-clock" }, value(props.server.refreshed, "Queue ready"), " | LibeRation ", value(props.server.package_version,"unknown"), props.server.queue_root ? " | "+props.server.queue_root : ""),
      e(Panel, { title: "Jobs", subtitle: jobs.length + " jobs", className: "lw-jobs-table-panel" },
        jobs.length ? e("div", { className: "lw-job-tree" }, jobs.map(function (job, index) { var id = value(job.id, String(index)); return e("div", { key: id, className: "lw-job-row " + (selected[0] === id ? "selected" : ""), onClick: function () { selected[1](id); emit(props, "job_select", { id: id }); } },
          e(StatusDot, { status: job.status === "failed" ? "error" : job.status === "completed" ? "ready" : "running" }),
          e("div", { className: "lw-job-main" }, e("strong", null, value(job.label, id)), e("span", null, id)),
          e("span", null, value(job.type, "job")), e("span", { className: "lw-job-state" }, value(job.status, "unknown")),
          job.status === "completed" ? e(Button, { className: "lw-button-quiet", title:"Open the saved model run and its results", onClick: function (event) { event.stopPropagation(); emit(props, "job_result", { id: id, queueId:props.server.queue_id }); } }, "Open run") : null); })) : e(Empty, { title: "No jobs", detail: "Run an estimation or simulation to populate the queue." })),
      e(Panel, { title: "Worker log", className: "lw-worker-log-panel" }, e("pre", { className: "lw-worker-log" }, list(props.job_log).join("\n") || "Select a job to view its worker log.")),
      e("p", { className: "lw-muted lw-worker-location" }, value(props.server.worker, "in-process"), " | ", value(props.server.isolation, "current R session")),
      e(Modal, {
        open: remoteModal[0], onClose: function () { remoteModal[1](false); },
        title: remoteEditId[0] ? "Edit remote LibeRties server" : "Add remote LibeRties server",
        footer: e(Button, { className: "lw-button-primary", onClick: function () {
          emit(props, "queue_add", { id: remoteEditId[0], name: remoteName[0], url: remoteUrl[0], token: remoteToken[0] });
          remoteToken[1](""); remoteModal[1](false);
        } }, remoteEditId[0] ? "Save and reconnect" : "Connect")
      }, e("div", { className: "lw-form-stack" },
        e(Field, { label: "Name" }, e("input", { value: remoteName[0], onChange: function (event) { remoteName[1](event.target.value); } })),
        e(Field, { label: "Server URL" }, e("input", { value: remoteUrl[0], onChange: function (event) { remoteUrl[1](event.target.value); } })),
        e(Field, { label: "Bearer token" }, e("input", { type: "password", value: remoteToken[0], placeholder:remoteEditId[0]?"Leave blank to keep saved token":"Required", onChange: function (event) { remoteToken[1](event.target.value); } })),
        e("p",{className:"lw-help-text"},"Server definitions and tokens are stored in the workspace settings directory so they survive package upgrades. The settings file is restricted to the current OS user where supported.")
      ))
    );
  }

  function parseBreaks(text) {
    return String(text || "").split(/[ ,;\n\t]+/).map(Number).filter(function (item) { return isFinite(item); }).sort(function (a,b) { return a-b; }).filter(function (item,index,array) { return !index || item !== array[index-1]; });
  }
  function summarizeData(rows, x, y, bins, mode, manualBreaks, groupKeys, quantileInterval, binPosition) {
    var filtered = rows.filter(function (row) { return number(row[x]) !== null && number(row[y]) !== null; });
    if (!filtered.length) return [];
    var range = extent(filtered, x), breaks = list(manualBreaks);
    if (breaks.length < 2) breaks = Array.from({length: bins + 1}, function (_, i) { return range[0] + i * (range[1] - range[0]) / bins; });
    var groups = {}, keys = list(groupKeys).filter(function (key) { return !!key; });
    filtered.forEach(function (row) {
      var xv = number(row[x]), index = -1;
      for (var i=0; i<breaks.length-1; i += 1) if (xv >= breaks[i] && (xv < breaks[i+1] || (i === breaks.length-2 && xv <= breaks[i+1]))) { index = i; break; }
      if (index < 0) return;
      var groupValues = keys.map(function (key) { return String(value(row[key], "(missing)")); });
      var id = [index].concat(groupValues).join("\u001f");
      if (!groups[id]) groups[id] = { index: index, values: [], groupValues: groupValues };
      groups[id].values.push(number(row[y]));
    });
    return Object.keys(groups).map(function (id) {
      var group = groups[id], values = group.values.sort(function (a,b) { return a-b; }), mean = values.reduce(function (a,b) { return a+b; },0)/values.length;
      var q = function (p) { var position = (values.length-1)*p, lo=Math.floor(position), hi=Math.ceil(position); return values[lo] + (values[hi]-values[lo])*(position-lo); };
      var variance = values.length > 1 ? values.reduce(function (total,item) { return total + Math.pow(item-mean,2); },0)/(values.length-1) : 0;
      var interval = Math.max(50, Math.min(99.9, number(quantileInterval) || 95)) / 100, tail = (1 - interval) / 2;
      var row = { X: binPosition === "midpoint" ? (breaks[group.index] + breaks[group.index+1])/2 : group.index + 1, X_LO: breaks[group.index], X_HI: breaks[group.index+1], X_LABEL: formatNumber(breaks[group.index]) + "-" + formatNumber(breaks[group.index+1]), Y: mode === "mean_se" ? mean : q(0.5), N: values.length, LOWER: mode === "mean_se" ? mean-Math.sqrt(variance/values.length) : q(tail), UPPER: mode === "mean_se" ? mean+Math.sqrt(variance/values.length) : q(1-tail), Q1: q(0.25), Q3: q(0.75) };
      keys.forEach(function (key,index) { row[key] = group.groupValues[index]; });
      return row;
    });
  }
  function binData(rows, key, bins, manualBreaks, position, outputKey) {
    var filtered = rows.filter(function (row) { return number(row[key]) !== null; }), range = extent(filtered, key), breaks = list(manualBreaks);
    if (breaks.length < 2) breaks = Array.from({length: bins + 1}, function (_, i) { return range[0] + i * (range[1] - range[0]) / bins; });
    var maximumBreak = breaks[breaks.length - 1];
    return rows.map(function (row) {
      var xv = number(row[key]), index = -1;
      for (var i=0; i<breaks.length-1; i += 1) if (xv !== null && xv >= breaks[i] && (xv < breaks[i+1] || (xv === maximumBreak && xv === breaks[i+1]))) { index=i; break; }
      if (index < 0) return null;
      var copy = Object.assign({}, row);
      copy[outputKey] = position === "midpoint" ? (breaks[index] + breaks[index+1])/2 : index+1;
      copy[outputKey + "_LABEL"] = formatNumber(breaks[index]) + "-" + formatNumber(breaks[index+1]);
      return copy;
    }).filter(function (row) { return !!row; });
  }
  function regressionData(rows, x, y, group, split) {
    var keys = {}, output = [];
    rows.forEach(function (row) { var groupValue=group?String(value(row[group],"(missing)")):"all",splitValue=split?String(value(row[split],"(missing)")):"all",key=groupValue+"\u001f"+splitValue; if (!keys[key]) keys[key]={rows:[],group:groupValue,split:splitValue}; if (number(row[x]) !== null && number(row[y]) !== null) keys[key].rows.push(row); });
    Object.keys(keys).forEach(function (key) { var entry=keys[key],values=entry.rows; if (values.length<2) return; var mx=values.reduce(function(a,r){return a+number(r[x]);},0)/values.length, my=values.reduce(function(a,r){return a+number(r[y]);},0)/values.length; var den=values.reduce(function(a,r){return a+Math.pow(number(r[x])-mx,2);},0); var slope=den ? values.reduce(function(a,r){return a+(number(r[x])-mx)*(number(r[y])-my);},0)/den : 0; var xr=extent(values,x); var a={X:xr[0],Y:my+slope*(xr[0]-mx)}, b={X:xr[1],Y:my+slope*(xr[1]-mx)}; if(group){a[group]=entry.group;b[group]=entry.group;}if(split){a[split]=entry.split;b[split]=entry.split;} output.push(a,b); });
    return output;
  }
  function DataPage(props) {
    var dataset = props.dataset || {}, rows = list(dataset.plot_rows);
    var numeric = list(dataset.numeric_columns), columns = list(dataset.columns);
    var x = React.useState(numeric.indexOf("TIME") >= 0 ? "TIME" : value(numeric[0], ""));
    var y = React.useState(numeric.indexOf("DV") >= 0 ? "DV" : value(numeric[1], value(numeric[0], "")));
    var strat = React.useState(""), split = React.useState(""), plotType = React.useState("points"), bins = React.useState(10);
    var binX = React.useState(false), binY = React.useState(false), showTable = React.useState(false), allRows = React.useState(false), adjust = React.useState(false);
    var manualX = React.useState(false), manualY = React.useState(false), xBreaks = React.useState(""), yBreaks = React.useState(""), xBinPosition = React.useState("equal"), yBinPosition = React.useState("equal");
    var showPoints = React.useState(false), lineMode = React.useState("individual"), pointShape = React.useState("16");
    var title = React.useState(""), xlab = React.useState(""), ylab = React.useState(""), pointSize = React.useState(0.85), quantile = React.useState(95), shade = React.useState(25), scatter = React.useState(25);
    function load(event) { var file = event.target.files && event.target.files[0]; if (!file) return; var reader = new FileReader(); reader.onload = function () { emit(props, "load_csv", { name: file.name, text: String(reader.result || "") }); }; reader.readAsText(file); }
    if (dataset.loaded && !dataset.payload_loaded) {
      return e("div", { className:"lw-ribbon-page lw-data-page lw-data-lazy" },
        e("aside",{className:"lw-data-controls"},
          e("div",{className:"lw-data-summary"},e("strong",null,dataset.records+" records / "+dataset.subjects+" subjects"),e("span",null,dataset.observations+" observations")),
          e(Field,{label:"Dataset"},e("select",{value:"current",disabled:true},e("option",{value:"current"},value(dataset.name,"Selected model dataset")))),
          e(Button,{className:"lw-button-primary",title:"Load this dataset into the Data explorer",onClick:function(){emit(props,"load_payload",{kind:"data"});}},"Load selected dataset")),
        e("main",{className:"lw-data-canvas"},e(Empty,{title:"Dataset is not loaded into the browser",detail:"Choose Load selected dataset when you want to explore it. Model and run selection stay lightweight."})));
    }
    var observations = rows.filter(function (row) { return allRows[0] || (Number(value(row.EVID,0)) === 0 && Number(value(row.MDV,0)) === 0); });
    var chartRows = observations, chartX = x[0], chartY = y[0], lines = ["lines","both"].indexOf(plotType[0]) >= 0, hidePoints = plotType[0] === "lines";
    var summaryType = ["mean_se","median_q","boxplot","violin"].indexOf(plotType[0]) >= 0;
    var aggregateBins = binX[0] && ["lines","smooth","regression"].indexOf(plotType[0]) >= 0;
    if (summaryType || plotType[0] === "smooth" || aggregateBins) { chartRows = summarizeData(observations, x[0], y[0], plotType[0] === "smooth" ? Math.max(6,bins[0]) : bins[0], plotType[0] === "median_q" ? "median_q" : "mean_se", manualX[0] ? parseBreaks(xBreaks[0]) : [], [strat[0],split[0]], quantile[0], xBinPosition[0]); chartX = "X"; chartY = "Y"; lines = plotType[0] !== "boxplot" && plotType[0] !== "violin"; hidePoints = ["boxplot","violin"].indexOf(plotType[0]) >= 0 ? false : summaryType ? false : plotType[0] === "smooth"; }
    if (binX[0] && ["points","jitter","both"].indexOf(plotType[0]) >= 0) { chartRows=binData(observations,x[0],bins[0],manualX[0]?parseBreaks(xBreaks[0]):[],xBinPosition[0],"__BINX");chartX="__BINX"; }
    if (plotType[0] === "regression") { chartRows = aggregateBins ? regressionData(chartRows, "X", "Y", strat[0], split[0]) : regressionData(observations, x[0], y[0], strat[0], split[0]); chartX="X"; chartY="Y"; lines=true; hidePoints=true; }
    if (plotType[0] === "jitter") { var xr=extent(chartRows,chartX), yr=extent(chartRows,chartY); chartRows=chartRows.map(function(row,index){var next=Object.assign({},row), phase=((index*37)%101)/100-0.5; next.__X=number(row[chartX])+phase*(xr[1]-xr[0])*0.012; next.__Y=number(row[chartY])-phase*(yr[1]-yr[0])*0.012; return next;}); chartX="__X";chartY="__Y"; }
    if (binY[0] && chartY === y[0]) { var yrange=extent(chartRows,y[0]), ybreak=manualY[0]?parseBreaks(yBreaks[0]):[]; if(ybreak.length<2)ybreak=Array.from({length:bins[0]+1},function(_,i){return yrange[0]+i*(yrange[1]-yrange[0])/bins[0];}); chartRows=chartRows.map(function(row){var next=Object.assign({},row),yv=number(row[y[0]]);next.__BINY=yv;for(var i=0;i<ybreak.length-1;i+=1)if(yv>=ybreak[i]&&(yv<ybreak[i+1]||(i===ybreak.length-2&&yv<=ybreak[i+1]))){next.__BINY=yBinPosition[0] === "midpoint" ? (ybreak[i]+ybreak[i+1])/2 : i+1;break;}return next;});chartY="__BINY"; }
    if (["lines","both"].indexOf(plotType[0]) >= 0 && lineMode[0] !== "individual") {
      if (lineMode[0] === "none") lines = false;
      else { chartRows=summarizeData(observations,x[0],y[0],Math.max(3,bins[0]),lineMode[0] === "mean" ? "mean_se" : "median_q",manualX[0]?parseBreaks(xBreaks[0]):[],[strat[0],split[0]],quantile[0],xBinPosition[0]);chartX="X";chartY="Y";lines=true;hidePoints=plotType[0] === "lines"; }
    }
    function renderDataChart(plotRows, plotTitle, width, height) {
      if (["boxplot","violin"].indexOf(plotType[0]) >= 0) return e(DistributionPlot, { rows: plotRows, kind: plotType[0] === "boxplot" ? "box" : "violin", group: strat[0], title: plotTitle, xLabel: xlab[0] || x[0], yLabel: ylab[0] || y[0], width: width, height: height });
      var isSummary = ["mean_se","median_q"].indexOf(plotType[0]) >= 0, first = plotRows[0], last = plotRows[plotRows.length-1];
      var overlayData = [];
      if (isSummary && showPoints[0]) {
        var candidates = observations;
        if (split[0] && plotRows.length) candidates = candidates.filter(function (row) { return String(value(row[split[0]], "(missing)")) === String(value(plotRows[0][split[0]], "(missing)")); });
        var maximumBreak = plotRows.length ? Math.max.apply(null, plotRows.map(function (item) { return item.X_HI; })) : null;
        overlayData = candidates.map(function (row) {
          var xv = number(row[x[0]]), bin = plotRows.filter(function (item) { return xv !== null && xv >= item.X_LO && (xv < item.X_HI || (xv === maximumBreak && xv === item.X_HI)); })[0];
          if (!bin) return null;
          var copy = Object.assign({}, row); copy.__SUMMARY_X = bin.X; return copy;
        }).filter(function (row) { return !!row; });
      }
      return e(ScatterPlot, { rows: plotRows, x: chartX, y: chartY, group: strat[0] || ((["lines","both"].indexOf(plotType[0]) >= 0 && lineMode[0] === "individual" && columns.indexOf("ID") >= 0) ? "ID" : ""), lineGroup: (["lines","both"].indexOf(plotType[0]) >= 0 && lineMode[0] === "individual" && columns.indexOf("ID") >= 0) ? "ID" : strat[0], lines: lines, intervals: isSummary, intervalShade: plotType[0] === "median_q" ? shade[0] / 100 : 0, overlayRows: overlayData, overlayX: "__SUMMARY_X", overlayY: y[0], overlayScatter: scatter[0] / 100, hidePoints: hidePoints, pointShape: pointShape[0], title: plotTitle, pointSize: pointSize[0] * 3, xLabel: xlab[0] || x[0], yLabel: ylab[0] || y[0], xTickStart: chartX === "X" && xBinPosition[0] === "equal" && first ? first.X_LABEL : null, xTickEnd: chartX === "X" && xBinPosition[0] === "equal" && last ? last.X_LABEL : null, width: width, height: height });
    }
    return e("div", { className: "lw-ribbon-page lw-data-page" },
      e("aside", { className: "lw-data-controls" },
        e("div", { className: "lw-data-summary" }, e("strong", null, dataset.loaded ? dataset.records + " records / " + dataset.subjects + " subjects" : "No dataset"), e("span", null, dataset.loaded ? dataset.observations + " observations" : "Load a CSV to begin")),
        e("label", { className: "lw-button lw-button-primary lw-file-button" }, "Import dataset", e("input", { type: "file", accept: ".csv,text/csv", onChange: load })),
        e("label", { className: "lw-check" }, e("input", { type: "checkbox", checked: showTable[0], onChange: function (event) { showTable[1](event.target.checked); } }), " Show dataset table"),
        showTable[0] ? e("label", { className: "lw-check" }, e("input", { type: "checkbox", checked: allRows[0], onChange: function (event) { allRows[1](event.target.checked); } }), " Show all rows (including doses)") : null,
        e("h4", null, "Plot"),
        e(Field, { label: "X axis" }, e("select", { value: x[0], onChange: function (event) { x[1](event.target.value); } }, numeric.map(function (column) { return e("option", { key: column, value: column }, column); }))),
        e(Field, { label: "Y axis" }, e("select", { value: y[0], onChange: function (event) { y[1](event.target.value); } }, numeric.map(function (column) { return e("option", { key: column, value: column }, column); }))),
        e(Field, { label: "Stratify / colour by" }, e("select", { value: strat[0], onChange: function (event) { strat[1](event.target.value); } }, [e("option", { key: "none", value: "" }, "(none)")].concat(columns.map(function (column) { return e("option", { key: column, value: column }, column); })))),
        e(Field, { label: "Split by (facets)" }, e("select", { value: split[0], onChange: function (event) { split[1](event.target.value); } }, [e("option", { key: "none", value: "" }, "(none)")].concat(columns.map(function (column) { return e("option", { key: column, value: column }, column); })))),
        e("label", { className: "lw-check" }, e("input", { type: "checkbox", checked: binX[0], onChange: function (event) { binX[1](event.target.checked); } }), " Bin X (continuous)"),
        binX[0] ? e("div", { className: "lw-bin-controls" },
          e(Field, { label: "X bin position" }, e("select", { value: xBinPosition[0], onChange: function (event) { xBinPosition[1](event.target.value); } }, e("option", { value: "equal" }, "Equidistant bins"), e("option", { value: "midpoint" }, "At bin midpoints"))),
          e("label", { className: "lw-check" }, e("input", { type: "checkbox", checked: manualX[0], onChange: function (event) { manualX[1](event.target.checked); } }), " Manual X breaks"),
          manualX[0] ? e(Field, { label: "X break values" }, e("input", { value: xBreaks[0], placeholder: "0, 2, 5, 10", onChange: function (event) { xBreaks[1](event.target.value); } })) : null) : null,
        e("label", { className: "lw-check" }, e("input", { type: "checkbox", checked: binY[0], onChange: function (event) { binY[1](event.target.checked); } }), " Bin Y (continuous)"),
        binY[0] ? e("div", { className: "lw-bin-controls" },
          e(Field, { label: "Y bin position" }, e("select", { value: yBinPosition[0], onChange: function (event) { yBinPosition[1](event.target.value); } }, e("option", { value: "equal" }, "Equidistant bins"), e("option", { value: "midpoint" }, "At bin midpoints"))),
          e("label", { className: "lw-check" }, e("input", { type: "checkbox", checked: manualY[0], onChange: function (event) { manualY[1](event.target.checked); } }), " Manual Y breaks"),
          manualY[0] ? e(Field, { label: "Y break values" }, e("input", { value: yBreaks[0], placeholder: "0, 1, 5, 20", onChange: function (event) { yBreaks[1](event.target.value); } })) : null) : null,
        (binX[0] || binY[0] || ["mean_se","median_q","boxplot","violin"].indexOf(plotType[0]) >= 0) ? e(Field, { label: "Number of bins" }, e("input", { type: "number", min: 3, max: 50, value: bins[0], onChange: function (event) { bins[1](Number(event.target.value)); } })) : null,
        e(Field, { label: "Plot type" }, e("select", { value: plotType[0], onChange: function (event) { plotType[1](event.target.value); } }, [
          ["points","Points"],["jitter","Jittered points"],["lines","Lines"],["both","Points + lines"],["smooth","Smooth (moving mean)"],["regression","Linear regression"],["boxplot","Box plot"],["violin","Violin"],["mean_se","Mean +/- SE"],["median_q","Median + quantiles"]
        ].map(function (item) { return e("option", { key: item[0], value: item[0] }, item[1]); }))),
        ["mean_se","median_q"].indexOf(plotType[0]) >= 0 ? e("label", { className: "lw-check" }, e("input", { type: "checkbox", checked: showPoints[0], onChange: function (event) { showPoints[1](event.target.checked); } }), " Show individual points") : null,
        ["both","lines"].indexOf(plotType[0]) >= 0 ? e(Field, { label: "Line mode" }, e("select", { value: lineMode[0], onChange: function (event) { lineMode[1](event.target.value); } }, e("option", { value: "individual" }, "Each individual"), e("option", { value: "mean" }, "Mean across subjects"), e("option", { value: "median" }, "Median across subjects"), e("option", { value: "none" }, "No line"))) : null,
        e("label", { className: "lw-check" }, e("input", { type: "checkbox", checked: adjust[0], onChange: function (event) { adjust[1](event.target.checked); } }), " Adjust plot options"),
        adjust[0] ? e("div", { className: "lw-adjust-controls" },
          e(Field, { label: "Plot title" }, e("input", { value: title[0], onChange: function (event) { title[1](event.target.value); } })),
          e(Field, { label: "X axis label" }, e("input", { value: xlab[0], onChange: function (event) { xlab[1](event.target.value); } })),
          e(Field, { label: "Y axis label" }, e("input", { value: ylab[0], onChange: function (event) { ylab[1](event.target.value); } })),
          e(Field, { label: "Point size" }, e("input", { type: "range", min: 0.4, max: 1.8, step: 0.05, value: pointSize[0], onChange: function (event) { pointSize[1](Number(event.target.value)); } })),
          e(Field, { label: "Point shape" }, e("select", { value: pointShape[0], onChange: function (event) { pointShape[1](event.target.value); } }, e("option", { value: "16" }, "Circle"), e("option", { value: "1" }, "Open circle"), e("option", { value: "15" }, "Square"), e("option", { value: "17" }, "Triangle"), e("option", { value: "18" }, "Diamond"))),
          plotType[0] === "median_q" ? e(Field, { label: "Quantile interval (%)" }, e("input", { type: "number", min: 50, max: 99.9, step: 0.5, value: quantile[0], onChange: function (event) { quantile[1](Number(event.target.value)); } })) : null,
          plotType[0] === "median_q" ? e(Field, { label: "Quantile shade intensity" }, e("input", { type: "range", min: 5, max: 70, step: 5, value: shade[0], onChange: function (event) { shade[1](Number(event.target.value)); } })) : null,
          e(Field, { label: "Point scatter" }, e("input", { type: "range", min: 0, max: 100, step: 5, value: scatter[0], onChange: function (event) { scatter[1](Number(event.target.value)); } }))) : null),
      e("main", { className: "lw-data-canvas" },
        showTable[0] ? e(Panel, { title: "Dataset", className: "lw-data-table-panel" }, e(SimpleTable, { rows: allRows[0] ? list(dataset.preview_all) : list(dataset.preview), columns: columns })) : null,
        split[0] ? e("div", { className: "lw-facet-grid" }, Array.from(new Set(chartRows.map(function (row) { return String(value(row[split[0]], "(missing)")); }))).slice(0,9).map(function (facet) { return e(React.Fragment, { key: facet }, renderDataChart(chartRows.filter(function (row) { return String(value(row[split[0]], "(missing)")) === facet; }), facet, 500, 250)); })) :
          e("div",{className:"lw-data-chart-single"},renderDataChart(chartRows, title[0] || y[0] + " vs " + x[0], 700, 400))));
  }

  function LogBanner(props) {
    var open = React.useState(false), log = props.log || {};
    return e("div", { className: "lw-log-wrap" },
      e("div", { className: "lw-log-banner lw-log-" + value(log.level, "info") }, e(StatusDot, { status: log.level === "error" ? "error" : "ready" }), e("span", null, value(log.current, "Ready")),
        e("div", null, e(Button, { className: "lw-button-link", onClick: function () { open[1](!open[0]); } }, open[0] ? "Hide history" : "Show history"), e(Button, { className: "lw-button-link", onClick: function () { emit(props, "clear_log"); } }, "Clear log"))),
      open[0] ? e("div", { className: "lw-log-history" }, list(log.history).map(function (line, index) { return e("div", { key: index }, line); })) : null);
  }

  function AIActivation(props) {
    var ai=props.ai||{},active=useSynced(!!ai.activated,[!!ai.activated]),consent=useSynced(!!ai.consented,[!!ai.consented]),consentModal=React.useState(false),settingsModal=React.useState(false);
    function save(next,agreed){active[1](next);consent[1](agreed);if(!next)localAIShutdown();var detail=localAISettingsDetail(ai);detail.activated=next;detail.consented=agreed;emit(props,"ai_settings",detail);}
    function toggle(event){var next=event.target.checked;if(next&&!consent[0]){consentModal[1](true);return;}save(next,consent[0]);}
    return e(React.Fragment,null,
      e("label",{className:"lw-ai-toggle",title:"Enable optional browser-local WebGPU assistance"},e("span",null,"Activate AI"),e("input",{type:"checkbox",checked:active[0],onChange:toggle}),e("i",null)),
      e(Button,{className:"lw-ai-settings-button",title:"Local AI settings",onClick:function(){settingsModal[1](true);}},"..."),
      e(Modal,{open:settingsModal[0],className:"lw-modal-ai-settings",onClose:function(){settingsModal[1](false);},title:"Local AI settings",footer:e(Button,{className:"lw-button-primary",onClick:function(){settingsModal[1](false);}},"Done")},
        e("div",{className:"lw-ai-settings-grid"},
          e("section",{className:"lw-modal-section lw-modal-section-tinted"},e("h4",null,"Help assistant"),e("p",{className:"lw-help-text"},"Optimised for model code, syntax and workflow questions."),e(AIModelSelect,Object.assign({},props,{className:"lw-ai-model-settings",purpose:"help",label:"Model"})),e(AIContextSelect,Object.assign({},props,{className:"lw-ai-context-settings",purpose:"help",label:"Context window"}))),
          e("section",{className:"lw-modal-section lw-modal-section-tinted"},e("h4",null,"Report builder"),e("p",{className:"lw-help-text"},"A separate, larger model can synthesize selected model runs and diagnostics."),e(AIModelSelect,Object.assign({},props,{className:"lw-ai-model-settings",purpose:"report",label:"Model"})),e(AIContextSelect,Object.assign({},props,{className:"lw-ai-context-settings",purpose:"report",label:"Context window"})))),
        e("p",{className:"lw-help-text lw-ai-settings-note"},"Selections are saved immediately. Models remain lazy-loaded and only one model occupies GPU memory at a time.")),
      e(Modal,{open:consentModal[0],onClose:function(){consentModal[1](false);},title:"Activate browser-local AI",footer:e(React.Fragment,null,e(Button,{className:"lw-button-quiet",onClick:function(){consentModal[1](false);}},"Cancel"),e(Button,{className:"lw-button-primary",onClick:function(){consentModal[1](false);save(true,true);}},"Activate AI"))},
        e("div",{className:"lw-ai-consent"},e("p",null,"LibeRation uses WebGPU language models that run entirely in a dedicated worker inside this browser session."),e("ul",null,e("li",null,"Each selected open model is downloaded on first actual use and cached for later sessions."),e("li",null,"Only one model is held in GPU memory. Switching model or context unloads the resident model and lazily loads the new configuration."),e("li",null,"Activation, model choices, and Help/Report context settings are remembered, but no model loads until you ask for help or draft report text."),e("li",null,"Auto uses a model-aware context size and falls back to a smaller window if browser GPU allocation fails."),e("li",null,"After loading, network APIs in the AI worker are disabled before model context is supplied."),e("li",null,"The model receives no tools, DOM access, or ability to change a model or report.")),e("p",{className:"lw-help-text"},"A compromised browser, extension, or operating system is outside this isolation boundary. AI output can be wrong and must be reviewed as modelling assistance, not clinical advice."))));
  }

  function LibeRWorkbench(props) {
    var ribbon = React.useState("home"), theme = React.useState(function () { return initialDarkTheme("liberationDarkTheme"); });
    React.useEffect(function(){
      if(localAI.status.stage==="error"&&!Object.keys(localAI.pending).length)localAIShutdown();
      if(window.Shiny&&window.Shiny.addCustomMessageHandler){window.Shiny.addCustomMessageHandler("liber-report-document",function(message){window.dispatchEvent(new CustomEvent("liber-report-document",{detail:message||{}}));});window.Shiny.addCustomMessageHandler("liber-report-directory",function(message){window.dispatchEvent(new CustomEvent("liber-report-directory",{detail:message||{}}));});}
    },[]);
    React.useEffect(function () { storeTheme(theme[0], "liberationDarkTheme", true); }, [theme[0]]);
    function toggleTheme() { theme[1](!theme[0]); }
    function changeRibbon(next){ribbon[1](next);emit(props,"page_change",{page:next});}
    var ribbonItems = [{id:"home",label:"Home"},{id:"jobs",label:"Jobs"},{id:"data",label:"Data"}];
    return e("div", { className: "lw-legacy-shell " + (theme[0] ? "theme-dark" : "theme-light") },
      e("header", { className: "lw-app-header" }, e("div", {className:"lw-app-brand"}, props.server&&props.server.icon?e("img",{className:"lw-app-icon",src:props.server.icon,alt:""}):null, e("div",{className:"lw-app-title"},e("strong", null, "LibeRation"), e("span", null, "Population PK/PD modelling"))),
        e("div", { className: "lw-app-header-right" },e(AIActivation,props), e("span", {className:"lw-app-version"}, "v"+value(props.server&&props.server.package_version,"unknown")), e("label", { className: "lw-theme-toggle" }, e("span", null, theme[0] ? "Dark" : "Light"), e("input", { type: "checkbox", checked: theme[0], onChange: toggleTheme }), e("i", null)), e("span", { className: "lw-workspace-path" }, value(props.workspace && props.workspace.path, "No workspace selected")))),
      e("div", { className: "lw-ribbon" }, e(Tabs, { value: ribbon[0], onChange: changeRibbon, items: ribbonItems, className: "lw-ribbon-tabs" })),
      e(LogBanner, props),
      e("div", { className: "lw-page-host" }, ribbon[0] === "home" ? e(HomePage, props) : ribbon[0] === "jobs" ? e(JobsPage, props) : e(DataPage, props)),
      e("footer", { className: "lw-statusbar" }, e("span", null, e(StatusDot, { status: props.result && props.result.status === "error" ? "error" : "ready" }), value(props.result && props.result.message, "Engine ready")), e("span", null, value(props.server && props.server.mode, "local") + " | C++17 | Exact AD")));
  }

  reactR.reactWidget("liberWorkbench", "output", { LibeRWorkbench: LibeRWorkbench }, {});
}());
