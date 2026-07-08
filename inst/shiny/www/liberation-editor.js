/* global CodeMirror, Shiny, $ */
(function () {
  "use strict";

  var BUILTIN = { THETA: 1, OMEGA: 1, SIGMA: 1, ETA: 1, ERR: 1 };
  var FUNCS = { EXP: 1, LOG: 1, SQRT: 1, ABS: 1 };

  function upperSet(arr) {
    var o = {};
    if (!arr) return o;
    for (var i = 0; i < arr.length; i++) {
      if (arr[i]) o[String(arr[i]).toUpperCase()] = true;
    }
    return o;
  }

  function editorMode(config, stream) {
    var hl = config.liberationHighlight || {};
    if (hl.variant) {
      return hl.variant;
    }
    if (stream && stream.cm && stream.cm.getTextArea) {
      var ta = stream.cm.getTextArea();
      if (ta && ta.getAttribute) {
        return ta.getAttribute("data-editor-mode") || "pk";
      }
    }
    return "pk";
  }

  CodeMirror.defineMode("liberation-pk", function (config) {
    function hlSets() {
      var hl = config.liberationHighlight || {};
      return {
        pk: upperSet(hl.pk),
        flow: upperSet(hl.flows),
        variant: hl.variant || "pk"
      };
    }

    return {
      startState: function () {
        return { lhsDone: false };
      },
      token: function (stream, state) {
        var sets = hlSets();

        if (stream.sol()) {
          state.lhsDone = false;
        }
        if (stream.eatSpace()) {
          return null;
        }
        if (stream.match(/^#.*/)) {
          return "pk-comment";
        }
        if (sets.variant === "des") {
          if (stream.match(/^DADT\s*\(\s*\d+\s*\)/i)) {
            return "pk-flow";
          }
          if (stream.match(/^A\s*\(\s*\d+\s*\)/i)) {
            return "pk-flow";
          }
        }
        if (stream.match(/^(<-|<=|>=|==|!=)/) || stream.match(/^(=|\^|\*\*|\+|-|\*|\/|\(|\)|,|<|>)/)) {
          if (stream.current() === "=" || stream.current() === "<-") {
            state.lhsDone = true;
          }
          return "pk-op";
        }
        if (stream.match(/^\d+(\.\d+)?([eE][+-]?\d+)?/)) {
          return "pk-num";
        }
        if (stream.match(/^[A-Za-z_][A-Za-z0-9_]*/)) {
          var up = stream.current().toUpperCase();
          if (!state.lhsDone && stream.match(/^\s*(<-|=)/, false)) {
            state.lhsDone = true;
            return "pk-def";
          }
          if (BUILTIN[up]) {
            return "pk-builtin";
          }
          if (FUNCS[up]) {
            return "pk-func";
          }
          if (sets.flow[up]) {
            return "pk-flow";
          }
          if (sets.pk[up]) {
            return "pk-param";
          }
          return "pk-var";
        }
        stream.next();
        return null;
      }
    };
  });

  CodeMirror.defineMIME("text/x-liberation-pk", "liberation-pk");

  function rowsToHeight(rows) {
    var r = parseInt(rows, 10);
    if (isNaN(r) || r < 3) r = 6;
    return Math.max(80, r * 18 + 16);
  }

  function modeOptions(textarea, highlight) {
    var variant = (highlight && highlight.variant) ||
      textarea.getAttribute("data-editor-mode") || "pk";
    var hl = highlight || {};
    hl.variant = variant;
    return { name: "liberation-pk", liberationHighlight: hl };
  }

  function initEditor(textarea) {
    if (textarea._cm) {
      return textarea._cm;
    }
    var rows = textarea.getAttribute("data-rows") || "6";
    var cm = CodeMirror.fromTextArea(textarea, {
      mode: modeOptions(textarea, { variant: textarea.getAttribute("data-editor-mode") || "pk" }),
      lineNumbers: true,
      lineWrapping: false,
      indentUnit: 2,
      tabSize: 2
    });
    cm.setSize(null, rowsToHeight(rows));
    textarea._cm = cm;

    cm.on("change", function () {
      cm.save();
      $(textarea).trigger("change");
    });

    return cm;
  }

  if (typeof Shiny !== "undefined") {
    Shiny.inputBindings.register({
      find: function (scope) {
        return $(scope).find("textarea.liberation-code-editor");
      },
      getId: function (el) {
        return el.id;
      },
      getType: function () {
        return "liberation.codeEditor";
      },
      getValue: function (el) {
        if (el._cm) {
          return el._cm.getValue();
        }
        return el.value;
      },
      setValue: function (el, value) {
        var cm = initEditor(el);
        cm.setValue(value == null ? "" : value);
      },
      receiveMessage: function (el, data) {
        var cm = initEditor(el);
        if (Object.prototype.hasOwnProperty.call(data, "value")) {
          cm.setValue(data.value == null ? "" : data.value);
        }
        if (data.highlight) {
          cm.setOption("mode", modeOptions(el, data.highlight));
        }
      },
      subscribe: function (el, callback) {
        $(el).on("change.liberationCodeEditor", function () {
          callback(true);
        });
      },
      unsubscribe: function (el) {
        $(el).off(".liberationCodeEditor");
      },
      initialize: function (el) {
        initEditor(el);
      },
      getRatePolicy: function () {
        return { policy: "debounce", delay: 250 };
      }
    });
  }
})();
