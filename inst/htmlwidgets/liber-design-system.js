(function (global) {
  "use strict";
  var api = global.LibeRDesign || {};
  api.theme = api.theme || {
    initialDark: function (legacyKey) {
      try {
        var shared = global.localStorage.getItem("liber.theme");
        if (shared === "dark" || shared === "light") return shared === "dark";
        var legacy = legacyKey ? global.localStorage.getItem(legacyKey) : null;
        if (legacy === "dark" || legacy === "1") return true;
        if (legacy === "light" || legacy === "0") return false;
      } catch (_) {}
      return !!(global.matchMedia &&
        global.matchMedia("(prefers-color-scheme: dark)").matches);
    },
    store: function (dark, legacyKey, numericLegacy) {
      var value = dark ? "dark" : "light";
      try {
        global.localStorage.setItem("liber.theme", value);
        if (legacyKey) {
          global.localStorage.setItem(
            legacyKey,
            numericLegacy ? (dark ? "1" : "0") : value
          );
        }
      } catch (_) {}
      global.document.documentElement.setAttribute("data-liber-theme", value);
      return dark;
    },
    bootstrap: function (legacyKey, numericLegacy) {
      var dark = this.initialDark(legacyKey);
      global.document.documentElement.setAttribute(
        "data-liber-theme", dark ? "dark" : "light"
      );
      return dark;
    }
  };
  api.messageClass = function (level, prefix) {
    var allowed = ["info", "working", "success", "warning", "error"];
    level = allowed.indexOf(level) >= 0 ? level : "info";
    return (prefix || "liber") + "-message " + (prefix || "liber") + "-" + level;
  };
  api.taskState = api.taskState || {
    eventName: "liber:task-state",
    register: function () {
      if (this.registered || !global.Shiny ||
          !global.Shiny.addCustomMessageHandler) return;
      this.registered = true;
      global.Shiny.addCustomMessageHandler("liber-task-state", function (message) {
        global.dispatchEvent(new CustomEvent(api.taskState.eventName, {
          detail: message || {}
        }));
      });
    },
    use: function (React, inputId, initial) {
      var state = React.useState(initial || {
        running: false, id: "", label: "", cancellable: false
      });
      React.useEffect(function () {
        api.taskState.register();
        function update(event) {
          var detail = event.detail || {};
          if (detail.inputId && detail.inputId !== inputId) return;
          state[1](detail);
        }
        global.addEventListener(api.taskState.eventName, update);
        return function () {
          global.removeEventListener(api.taskState.eventName, update);
        };
      }, [inputId]);
      return state[0];
    }
  };
  if (global.document) {
    global.document.addEventListener("shiny:connected", function () {
      api.taskState.register();
    });
  }
  global.LibeRDesign = api;
})(window);
