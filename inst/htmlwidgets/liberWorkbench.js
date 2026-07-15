(function () {
  "use strict";

  var e = React.createElement;
  var palette = ["#466d91", "#a36d3d", "#447963", "#745f86", "#a64e55", "#38777c", "#82713e"];

  function list(x) { return Array.isArray(x) ? x : []; }
  function value(x, fallback) { return x === null || x === undefined || x === "" ? fallback : x; }
  function number(x) { var n = Number(x); return isFinite(n) ? n : null; }
  function formatNumber(x) {
    var n = Number(x);
    if (!isFinite(n)) return "-";
    if (n !== 0 && (Math.abs(n) < 0.001 || Math.abs(n) >= 10000)) return n.toExponential(4);
    return n.toPrecision(6).replace(/\.?0+$/, "");
  }
  function emit(props, action, detail) {
    if (!window.Shiny || !window.Shiny.setInputValue) return;
    window.Shiny.setInputValue(
      (props.inputId || "liber_workbench") + "_event",
      Object.assign({ action: action, nonce: Date.now() }, detail || {}),
      { priority: "event" }
    );
  }
  function cloneRows(rows) { return list(rows).map(function (row) { return Object.assign({}, row); }); }
  function useSynced(initial, dependency) {
    var state = React.useState(initial);
    React.useEffect(function () { state[1](initial); }, dependency || []);
    return state;
  }

  function escapeCode(value) {
    return String(value || "").replace(/&/g, "&amp;").replace(/</g, "&lt;").replace(/>/g, "&gt;");
  }
  function highlightCode(source) {
    var text = String(source || ""), output = "", cursor = 0;
    var pattern = /(#.*$|\/\/.*$|"(?:\\.|[^"\\])*"|'(?:\\.|[^'\\])*'|\b(?:THETA|ETA|OMEGA|SIGMA|ERR|EPS)(?:_\d+|\s*\(\s*\d+(?:\s*,\s*\d+)?\s*\))|\$[A-Z][A-Z0-9_]*|\b(?:if|else|for|while|return|TRUE|FALSE|NA|NULL)\b|\b\d+(?:\.\d+)?(?:[eE][+-]?\d+)?\b|\b(?:exp|log|sqrt|sin|cos|tan|tanh|abs|expm1|log1p|ifelse|min|max|pow)\b(?=\s*\()|\b[A-Za-z_][A-Za-z0-9_.]*\b(?=\s*=)|\b(?:A|DADT)\s*\(\s*\d+\s*\)|\b(?:S\d+|F|Y|IPRED|PRED|TIME|T|AMT|RATE|CMT|EVID|MDV|II|SS|DV|DVID|LLOQ|BLQ|CENS|MIXNUM)\b)/gm;
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

  function StatusDot(props) { return e("span", { className: "lw-status-dot lw-status-" + value(props.status, "ready") }); }
  function Button(props) {
    return e("button", {
      type: "button", className: "lw-button " + value(props.className, ""), disabled: !!props.disabled,
      title: props.title, "aria-label": props.ariaLabel || props.title, onClick: props.onClick
    }, props.icon ? e("span", { className: "lw-button-icon" }, props.icon) : null, props.children);
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
    if (!props.open) return null;
    return e("div", { className: "lw-modal-backdrop", onMouseDown: function (event) { if (event.target === event.currentTarget) props.onClose(); } },
      e("div", { className: "lw-modal " + value(props.className, ""), role: "dialog", "aria-modal": "true" },
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
    return e("div", { className: "lw-parameter-block" }, e("h5", null, props.title), rows.length ? e("table", { className: "lw-param-table" },
      e("thead", null, e("tr", null, e("th", null, props.matrix ? "Element" : "Name"), e("th", null, "Initial"), props.bounds?e("th",null,"Lower"):null, props.bounds?e("th",null,"Upper"):null, e("th", null, "Fixed"))),
      e("tbody", null, rows.map(function (row, index) { var name = props.matrix ? "ETA" + value(row.ROW, index + 1) + " / ETA" + value(row.COL, index + 1) : props.prefix + value(row[props.indexName], index + 1); return e("tr", { key: name + "-" + index },
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
    var sourceState = useSynced({ pred: value(model.pred, ""), des: value(model.des, ""), error: value(model.error, "Y=F") }, [model.pred, model.des, model.error]);
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
    var columnModal = React.useState(false), transHelp = React.useState(false);
    var dirty = source.pred !== value(model.pred, "") || source.des !== value(model.des, "") || source.error !== value(model.error, "Y=F") || advanState[0] !== String(value(model.advan, 4)) || transState[0] !== String(value(model.trans, 2)) || problemState[0] !== value(model.name, "Untitled model") || JSON.stringify(parameters) !== JSON.stringify({ theta: cloneRows(model.theta), omega: cloneRows(model.omega), sigma: cloneRows(model.sigma) }) || JSON.stringify(priors) !== JSON.stringify(cloneRows(model.priors)) || omegaStructure !== value(model.omega_structure,"diagonal") || JSON.stringify(selectedInput) !== JSON.stringify(list(model.input));
    function updateParameter(kind, index, field, nextValue) {
      var next = Object.assign({}, parameters); next[kind] = cloneRows(parameters[kind]); next[kind][index][field] = nextValue; setParameters(next);
    }
    function changeOmegaStructure(nextStructure) {
      var nEta=Number(value(model.n_eta,parameters.omega.length)), current=cloneRows(parameters.omega), nextRows=[];
      function find(row,column){return current.filter(function(item,index){var r=Number(value(item.ROW,index+1)),c=Number(value(item.COL,index+1));return r===row&&c===column;})[0];}
      if(nextStructure==="full"){
        for(var row=1;row<=nEta;row+=1)for(var column=1;column<=row;column+=1){var existing=find(row,column);nextRows.push({OMEGA:nextRows.length+1,ROW:row,COL:column,Value:existing?Number(existing.Value):(row===column?0.1:0),FIX:existing?!!existing.FIX:false});}
      }else{
        for(var diagonal=1;diagonal<=nEta;diagonal+=1){var item=find(diagonal,diagonal);nextRows.push({OMEGA:diagonal,ROW:diagonal,COL:diagonal,Value:item?Number(item.Value):0.1,FIX:item?!!item.FIX:false});}
      }
      var remapped=cloneRows(priors).map(function(prior){if(!/^OMEGA\d+$/.test(prior.parameter))return prior;var old=current[Number(prior.parameter.replace("OMEGA",""))-1];if(!old)return null;var oldRow=Number(value(old.ROW,old.OMEGA)),oldCol=Number(value(old.COL,old.OMEGA));var nextIndex=nextRows.findIndex(function(item){return Number(item.ROW)===oldRow&&Number(item.COL)===oldCol;});if(nextIndex<0)return null;prior.parameter="OMEGA"+(nextIndex+1);return prior;}).filter(function(prior){return !!prior;});
      setPriors(remapped);setParameters(Object.assign({},parameters,{omega:nextRows}));omegaStructureState[1](nextStructure);
    }
    var priorParameterNames = parameters.theta.map(function(_,i){return "THETA"+(i+1);}).concat(parameters.sigma.map(function(_,i){return "SIGMA"+(i+1);})).concat(parameters.omega.map(function(_,i){return "OMEGA"+(i+1);}));
    function save(mode) {
      emit(props, "update_model", {
        pred: source.pred, des: source.des, error: source.error, advan: Number(advanState[0]), trans: Number(transState[0]),
        n_state: Number(nState[0]), problem: problemState[0], theta: parameters.theta, omega: parameters.omega, sigma: parameters.sigma,
        omega_structure: omegaStructure, priors: priors, input: selectedInput, save_mode: mode
      });
    }
    return e("div", { className: "lw-code-workspace" },
      dirty ? e("div", { className: "lw-dirty-banner" }, "Unsaved editor changes — apply changes before saving or running the model.") : null,
      e("div", { className: "lw-control-row" },
        e(Field, { label: "$PROBLEM", className: "lw-grow" }, e("input", { value: problemState[0], onChange: function (event) { problemState[1](event.target.value); } })),
        e(Field, { label: "ADVAN" }, e("select", { value: advanState[0], onChange: function (event) { advanState[1](event.target.value); } }, [1,2,3,4,6,11,12,13].map(function (x) { return e("option", { key: x, value: x }, x); }))),
        Number(advanState[0]) === 6 || Number(advanState[0]) === 13 ? e(Field, { label: "Compartments" }, e("input", { type: "number", min: 1, max: 20, value: nState[0], onChange: function (event) { nState[1](Number(event.target.value)); } })) :
          e(Field, { label: "TRANS" }, e("select", { value: transState[0], onChange: function (event) { transState[1](event.target.value); } }, [1,2,3,4,5,6].map(function (x) { return e("option", { key: x, value: x }, x); }))),
        Number(advanState[0]) === 6 || Number(advanState[0]) === 13 ? null : e(Button, { className: "lw-button-link lw-help-button", onClick: function () { transHelp[1](true); } }, "?"),
        e(Field, { label: "Dataset" }, e("select", { value: props.dataset.loaded ? "current" : "" }, e("option", { value: props.dataset.loaded ? "current" : "" }, props.dataset.loaded ? value(props.dataset.name, "Current dataset") : "No dataset"))),
        e(Button, { className: "lw-button-quiet lw-columns-button", onClick: function () { columnModal[1](true); } }, "Columns...")),
      e("div", { className: "lw-editor-grid " + ((Number(advanState[0]) === 6 || Number(advanState[0]) === 13) ? "lw-editor-grid-three" : "") },
        e("div", { className: "lw-editor-box" }, e("h5", null, "$PK / $PRED"), e(CodeEditor,{label:"PK or PRED model code",value:source.pred,onValue:function(next){setSource(Object.assign({},source,{pred:next}));}})),
        Number(advanState[0]) === 6 || Number(advanState[0]) === 13 ? e("div", { className: "lw-editor-box" }, e("h5", null, "$DES"), e(CodeEditor,{label:"DES differential equation code",value:source.des,onValue:function(next){setSource(Object.assign({},source,{des:next}));}})) : null,
        e("div", { className: "lw-editor-box" }, e("h5", null, "$ERROR"), e(CodeEditor,{label:"ERROR model code",value:source.error,onValue:function(next){setSource(Object.assign({},source,{error:next}));}}))),
      e("div", { className: "lw-parameter-grid" },
        e(ParameterGrid, { title: "THETA", bounds:true, prefix: "THETA", indexName: "THETA", rows: parameters.theta, onChange: function (i,f,v) { updateParameter("theta",i,f,v); } }),
        e("div",{className:"lw-omega-block"},e(ParameterGrid, { title: omegaStructure==="full"?"OMEGA lower triangle":"OMEGA", matrix:omegaStructure==="full", prefix: "OMEGA", indexName: "OMEGA", rows: parameters.omega, onChange: function (i,f,v) { updateParameter("omega",i,f,v); } }),e("label",{className:"lw-check lw-omega-matrix-toggle"},e("input",{type:"checkbox",checked:omegaStructure==="full",onChange:function(event){changeOmegaStructure(event.target.checked?"full":"diagonal");}})," OMEGA matrix")),
        e(ParameterGrid, { title: "SIGMA", prefix: "SIGMA", indexName: "SIGMA", rows: parameters.sigma, onChange: function (i,f,v) { updateParameter("sigma",i,f,v); } })),
      e(PriorGrid,{rows:priors,parameterNames:priorParameterNames,onChange:setPriors}),
      e("div", { className: "lw-inline-actions lw-editor-actions" },
        e(Button, { className: "lw-button-quiet", onClick: function () { emit(props, "validate"); } }, "Validate"),
        e(Button, { className: "lw-button-primary", onClick: function () { save("current"); } }, "Apply changes")),
      e(Modal, { open: columnModal[0], onClose: function () { columnModal[1](false); }, title: "$INPUT / $OUTPUT columns", footer: e(Button, { className: "lw-button-primary", onClick: function () { columnModal[1](false); } }, "Done") },
        e("div", { className: "lw-column-list" }, list(props.dataset.columns).map(function (column) { return e("label", { key: column }, e("input", { type: "checkbox", checked: selectedInput.indexOf(column) >= 0, onChange: function (event) { var next=selectedInput.filter(function(item){return item!==column;});if(event.target.checked)next.push(column);setSelectedInput(next); } }), column); }))),
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
      var categoricalRows=list(result.simulated).map(function(row){return {TIME:row.TIME,Y:row.median,LOWER:row.lower,UPPER:row.upper};});
      return e("div",{className:"lw-diagnostic-grid"},e(ScatterPlot,{rows:categoricalRows,x:"TIME",y:"Y",lines:true,intervals:true,intervalShade:0.25,overlayRows:result.observed,overlayX:"TIME",overlayY:"PROPORTION",overlayColor:"#17202a",title:"Categorical VPC: observed proportion",xLabel:"Time",yLabel:"Proportion"}),e(SimpleTable,{rows:result.simulated}));
    }
    if (tab === "vpc_tte") {
      var tteRows=list(result.simulated).map(function(row){return {TIME:row.TIME,Y:row.median,LOWER:row.lower,UPPER:row.upper};});
      return e("div",{className:"lw-diagnostic-grid"},e(ScatterPlot,{rows:tteRows,x:"TIME",y:"Y",lines:true,intervalShade:0.28,overlayRows:result.observed,overlayX:"TIME",overlayY:"SURVIVAL",overlayLines:true,overlayColor:"#17202a",hidePoints:true,title:"Time-to-event VPC",xLabel:"Time",yLabel:"Event-free survival"}),e(SimpleTable,{rows:result.simulated}));
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
    var comparison = React.useState([]), templateModal = React.useState(false), copyModal = React.useState(false), copyUpdateInits = React.useState(true), deleteModal = React.useState(null);
    var expandedVersions = React.useState({});
    var templateAdvan = React.useState("4"), templateTrans = React.useState("4"), templateNState = React.useState(2), templateLabel = React.useState(""), templateProblem = React.useState("Template model");
    var estimationModal = React.useState(false), estimationMethod = React.useState("FOCEI"), estimationLabel = React.useState(""), estimationMaxit = React.useState(200), etaMaxit = React.useState(100), tolerance = React.useState(0.000001), estimationCores = React.useState(1), printEvery = React.useState(0), methodSeed = React.useState(20260713), nImp = React.useState(200), nIter = React.useState(200), burn = React.useState(60), mcmcSteps = React.useState(2), nBurn = React.useState(500), nSample = React.useState(1000), nThin = React.useState(1);
    var covarianceStep = React.useState(false), covarianceType = React.useState("hessian"), covarianceTolerance = React.useState(0.00000001), covarianceSamples = React.useState(200);
    var simulationModal = React.useState(false), simulationLabel = React.useState("Simulation"), simulationSeed = React.useState(Math.floor(Math.random() * 99999) + 1), simulationCores = React.useState(1);
    var simulationSubjects = React.useState(value(props.dataset.subjects, 10)), simulationReplicates = React.useState(1), simulationDays = React.useState(1), simulationUseDesign = React.useState(false);
    var diagnosticModal = React.useState(false), diagnosticVpc = React.useState(true), diagnosticNpc = React.useState(false), diagnosticNpde = React.useState(false), diagnosticCategorical = React.useState(false), diagnosticTte = React.useState(false), diagnosticOutcome = React.useState("DV"), diagnosticEvent = React.useState("DV"), diagnosticNsim = React.useState(200), diagnosticSeed = React.useState(20260713), diagnosticPc = React.useState(false), diagnosticStratify = React.useState("");
    var uncertaintyModal = React.useState(false), uncertaintyBootstrap = React.useState(true), uncertaintyProfile = React.useState(false), uncertaintyReplicates = React.useState(100), uncertaintyPoints = React.useState(9), uncertaintySpan = React.useState(3), uncertaintyLevel = React.useState(0.95), uncertaintyParameters = React.useState(""), uncertaintyMaxit = React.useState(100);
    var scmModal = React.useState(false), scmCandidates = React.useState("CL,WT,power\nV,WT,power"), scmDirection = React.useState("both"), scmForward = React.useState(0.05), scmBackward = React.useState(0.01), scmMaxSteps = React.useState(20), scmMaxit = React.useState(100), scmLabel = React.useState("SCM model");
    var controlModal = React.useState(false), controlFile = React.useState(null), controlData = React.useState(null), controlNewProject = React.useState(!workspace.current), controlProjectName = React.useState("NONMEM import"), controlLabel = React.useState("NONMEM import"), exportModal = React.useState(false), exportName = React.useState("model.ctl"), exportDataPath = React.useState("data.csv");
    var doseMode = React.useState("single"), doseAmount = React.useState(320), doseCmt = React.useState(1), doseN = React.useState(3), doseII = React.useState(12), doseTable = React.useState("0 320"), obsPerDay = React.useState(8), simulationUseFit = React.useState(true);
    React.useEffect(function(){if(workspace.current_version){var next=Object.assign({},expandedVersions[0]);next[workspace.current_version]=true;expandedVersions[1](next);}},[workspace.current_version]);
    function toggleExpanded(id) { var next=Object.assign({},expandedVersions[0]);next[id]=!next[id];expandedVersions[1](next); }
    function toggleComparison(id, checked) { var next=comparison[0].filter(function(item){return item!==id;});if(checked)next=next.concat([id]).slice(-2);comparison[1](next); }
    function readProjectDataset(event) { var file=event.target.files&&event.target.files[0];if(!file){projectFile[1](null);return;}var reader=new FileReader();reader.onload=function(){projectFile[1]({name:file.name,text:String(reader.result||"")});};reader.readAsText(file); }
    function readControlFile(event) { var file=event.target.files&&event.target.files[0];if(!file){controlFile[1](null);return;}var reader=new FileReader();reader.onload=function(){controlFile[1]({name:file.name,text:String(reader.result||"")});};reader.readAsText(file); }
    function readControlData(event) { var file=event.target.files&&event.target.files[0];if(!file){controlData[1](null);return;}var reader=new FileReader();reader.onload=function(){controlData[1]({name:file.name,text:String(reader.result||"")});};reader.readAsText(file); }
    function submitNewProject() { emit(props,"project_create",{name:projectName[0],description:projectDescription[0],mode:projectMode[0],dataSource:projectDataSource[0],example:projectExample[0],nSubjects:Number(projectSubjects[0]),fileName:projectFile[0]&&projectFile[0].name,text:projectFile[0]&&projectFile[0].text,advan:Number(projectAdvan[0]),trans:Number(projectTrans[0]),nState:Number(projectNState[0]),label:projectLabel[0],problem:projectProblem[0]});newProject[1](false); }
    function submitTemplate() { emit(props,"model_template",{advan:Number(templateAdvan[0]),trans:Number(templateTrans[0]),nState:Number(templateNState[0]),label:templateLabel[0],problem:templateProblem[0]});templateModal[1](false); }
    var covarianceSupported = ["FO","FOCE","FOCEI","LAPLACE","ITS","IMP","SAEM"].indexOf(estimationMethod[0]) >= 0;
    function submitEstimate() { var covType=estimationMethod[0]==="FO"?"hessian":covarianceType[0],covSamples=estimationMethod[0]==="IMP"?Number(nImp[0]):Number(covarianceSamples[0]);emit(props,"estimate",{label:estimationLabel[0],method:estimationMethod[0],maxit:Number(estimationMaxit[0]),etaMaxit:Number(etaMaxit[0]),tolerance:Number(tolerance[0]),nCores:Number(estimationCores[0]),printEvery:Number(printEvery[0]),methodSeed:Number(methodSeed[0]),nImp:Number(nImp[0]),nIter:Number(nIter[0]),burn:Number(burn[0]),mcmcSteps:Number(mcmcSteps[0]),nBurn:Number(nBurn[0]),nSample:Number(nSample[0]),nThin:Number(nThin[0]),covariance:covarianceStep[0]&&covarianceSupported,covarianceType:covType,covarianceTolerance:Number(covarianceTolerance[0]),covarianceSamples:covSamples,covarianceSeed:Number(methodSeed[0])});estimationModal[1](false); }
    function submitSimulation() { emit(props,"simulate",{label:simulationLabel[0],seed:Number(simulationSeed[0]),nCores:Number(simulationCores[0]),nSubjects:Number(simulationSubjects[0]),replicates:Number(simulationReplicates[0]),days:Number(simulationDays[0]),useDesign:simulationUseDesign[0],doseMode:doseMode[0],doseAmt:Number(doseAmount[0]),doseCmt:Number(doseCmt[0]),doseN:Number(doseN[0]),doseII:Number(doseII[0]),doseTable:doseTable[0],obsPerDay:Number(obsPerDay[0]),useFit:simulationUseFit[0]});simulationModal[1](false); }
    function submitDiagnostic() { var types=[];if(diagnosticVpc[0])types.push("vpc");if(diagnosticNpc[0])types.push("npc");if(diagnosticNpde[0])types.push("npde");if(diagnosticCategorical[0])types.push("vpc_categorical");if(diagnosticTte[0])types.push("vpc_tte");emit(props,"run_diagnostic",{types:types,nsim:Number(diagnosticNsim[0]),seed:Number(diagnosticSeed[0]),pcCorrect:diagnosticPc[0],stratify:diagnosticVpc[0]&&diagnosticStratify[0]?diagnosticStratify[0]:null,categoricalOutcome:diagnosticOutcome[0],tteEvent:diagnosticEvent[0]});diagnosticModal[1](false); }
    function submitUncertainty() { var types=[];if(uncertaintyBootstrap[0])types.push("bootstrap");if(uncertaintyProfile[0])types.push("profile");emit(props,"run_uncertainty",{types:types,replicates:Number(uncertaintyReplicates[0]),points:Number(uncertaintyPoints[0]),span:Number(uncertaintySpan[0]),level:Number(uncertaintyLevel[0]),parameters:uncertaintyParameters[0],maxit:Number(uncertaintyMaxit[0]),seed:Number(diagnosticSeed[0])});uncertaintyModal[1](false); }
    var projectUploadMissing=projectMode[0]==="template"&&projectDataSource[0]==="upload"&&!projectFile[0];
    return e("div", { className:"lw-project-sidebar" },
      e("div",{className:"lw-tree-title"},"Projects"),
      e("div",{className:"lw-tree-list lw-project-list"},projects.length?projects.map(function(project){return e("button",{type:"button",key:project.id,className:workspace.current===project.id?"selected":"",title:value(project.description,""),onClick:function(){emit(props,"project_open",{id:project.id});}},e("strong",null,project.name),e("span",null,value(project.versions,project.snapshots)+" versions"));}):e(Empty,{title:"No projects",detail:"Create a project below."})),
      e("div",{className:"lw-sidebar-actions lw-action-grid"},e(Button,{className:"lw-button-primary",icon:"+",title:"Create a new project",disabled:!workspace.enabled,onClick:function(){newProject[1](true);}},"New project"),e(Button,{className:"lw-button-quiet lw-action-button",icon:"NM",title:"Load a NONMEM control stream",disabled:!workspace.enabled,onClick:function(){controlNewProject[1](!workspace.current);controlModal[1](true);}},"Load .ctl")),
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
        e(Button,{className:"lw-button-danger-ghost",title:"Delete the current project",disabled:!workspace.current,onClick:function(){deleteModal[1]("project");}},"Project"),
        e(Button,{className:"lw-button-danger-ghost",title:"Delete the selected model version",disabled:!workspace.current_version,onClick:function(){deleteModal[1]("version");}},"Version"),
        e(Button,{className:"lw-button-danger-ghost lw-action-span",title:"Delete the selected estimation or simulation run",disabled:!workspace.current_run,onClick:function(){deleteModal[1]("run");}},"Selected run"))),

      e(Modal,{open:newProject[0],className:"lw-modal-wide",onClose:function(){newProject[1](false);},title:"New project",footer:e(React.Fragment,null,e(Button,{className:"lw-button-quiet",onClick:function(){newProject[1](false);}},"Cancel"),e(Button,{className:"lw-button-primary",disabled:!projectName[0].trim()||projectUploadMissing,onClick:submitNewProject},"Create project"))},
        e("div",{className:"lw-modal-section"},e("div",{className:"lw-form-grid lw-form-grid-two"},e(Field,{label:"Project name"},e("input",{autoFocus:true,value:projectName[0],onChange:function(event){projectName[1](event.target.value);}})),e(Field,{label:"Description (optional)"},e("textarea",{rows:2,value:projectDescription[0],onChange:function(event){projectDescription[1](event.target.value);}})))),
        e("div",{className:"lw-choice-cards"},[["empty","Empty project"],["template","Create from template"]].map(function(item){return e("label",{key:item[0],className:"lw-choice-card "+(projectMode[0]===item[0]?"selected":"")},e("input",{type:"radio",name:"project-mode",checked:projectMode[0]===item[0],onChange:function(){projectMode[1](item[0]);}}),e("span",null,e("strong",null,item[1])));})),
        projectMode[0]==="template"?e("div",{className:"lw-modal-section lw-modal-section-tinted"},e("h4",null,"Initial model version"),e("p",{className:"lw-help-text"},"Choose a built-in example or import an existing NONMEM-style dataset."),e("div",{className:"lw-choice-row"},e("label",{className:"lw-check"},e("input",{type:"radio",name:"project-data",checked:projectDataSource[0]==="synthetic",onChange:function(){projectDataSource[1]("synthetic");}})," Built-in synthetic example"),e("label",{className:"lw-check"},e("input",{type:"radio",name:"project-data",checked:projectDataSource[0]==="upload",onChange:function(){projectDataSource[1]("upload");}})," Upload dataset")),projectDataSource[0]==="synthetic"?e(React.Fragment,null,e("div",{className:"lw-form-grid lw-form-grid-two"},e(Field,{label:"Example"},e("select",{value:projectExample[0],onChange:function(event){projectExample[1](event.target.value);}},e("option",{value:"theophylline"},"Theophylline-style oral PK"),e("option",{value:"sparse"},"Sparse oral PK"),e("option",{value:"rich"},"Rich sampling oral PK"))),e(Field,{label:"Number of subjects"},e("input",{type:"number",min:1,max:500,value:projectSubjects[0],onChange:function(event){projectSubjects[1](Number(event.target.value));}}))),e("div",{className:"lw-form-grid lw-form-grid-two"},e(Field,{label:"Version label (optional)"},e("input",{value:projectLabel[0],onChange:function(event){projectLabel[1](event.target.value);}})),e(Field,{label:"Problem statement"},e("input",{value:projectProblem[0],onChange:function(event){projectProblem[1](event.target.value);}})))):e(React.Fragment,null,e(Field,{label:"Dataset file (.csv, .txt, .dat, .tsv)"},e("input",{type:"file",accept:".csv,.txt,.dat,.tsv,text/csv,text/plain",onChange:readProjectDataset})),e("p",{className:"lw-help-text"},projectFile[0]?"Loaded "+projectFile[0].name:"Expected NONMEM-style ID, TIME, DV, AMT, EVID, CMT and MDV columns."),e(TemplateFields,{advan:projectAdvan,trans:projectTrans,nState:projectNState,label:projectLabel,problem:projectProblem}))):null),

      e(Modal,{open:copyModal[0],onClose:function(){copyModal[1](false);},title:"Copy to new model version",footer:e(React.Fragment,null,e(Button,{className:"lw-button-quiet",onClick:function(){copyModal[1](false);}},"Cancel"),e(Button,{className:"lw-button-primary",onClick:function(){emit(props,"project_copy",{id:workspace.current,snapshot:workspace.current_snapshot,updateInits:copyUpdateInits[0]});copyModal[1](false);}},"Copy"))},e("p",{className:"lw-help-text"},props.fit.available?"A fitted run is loaded; its final estimates can become the new version's initial values.":"No fitted run is loaded, so initials will match the source version."),e("label",{className:"lw-check"},e("input",{type:"checkbox",disabled:!props.fit.available,checked:copyUpdateInits[0]&&props.fit.available,onChange:function(event){copyUpdateInits[1](event.target.checked);}})," Update THETA / OMEGA / SIGMA initials from current fit")),

      e(Modal,{open:templateModal[0],onClose:function(){templateModal[1](false);},title:"New version from template",footer:e(React.Fragment,null,e(Button,{className:"lw-button-quiet",onClick:function(){templateModal[1](false);}},"Cancel"),e(Button,{className:"lw-button-primary",disabled:!props.dataset.loaded,onClick:submitTemplate},"Create version"))},e(Field,{label:"Dataset"},e("select",{disabled:!props.dataset.loaded,value:props.dataset.loaded?"current":""},e("option",{value:props.dataset.loaded?"current":""},props.dataset.loaded?value(props.dataset.name,"Current dataset"):"No dataset loaded"))),e(TemplateFields,{advan:templateAdvan,trans:templateTrans,nState:templateNState,label:templateLabel,problem:templateProblem})),

      e(Modal,{open:estimationModal[0],className:"lw-modal-wide",onClose:function(){estimationModal[1](false);},title:"Run estimation",footer:e(React.Fragment,null,e(Button,{className:"lw-button-quiet",onClick:function(){estimationModal[1](false);}},"Cancel"),e(Button,{className:"lw-button-primary",onClick:submitEstimate},"Submit estimation"))},
        e("div",{className:"lw-modal-section"},e("div",{className:"lw-form-grid"},e(Field,{label:"Run on"},e("select",{value:value(props.server.queue_id,"local"),onChange:function(event){emit(props,"queue_select",{id:event.target.value});}},list(props.server.queues).map(function(queue){return e("option",{key:queue.id,value:queue.id},queue.name);}))),e(Field,{label:"Method"},e("select",{value:estimationMethod[0],onChange:function(event){estimationMethod[1](event.target.value);}},["FO","FOCE","FOCEI","LAPLACE","ITS","IMP","SAEM","BAYES"].map(function(method){return e("option",{key:method,value:method},method);}))),e(Field,{label:"Job label (optional)"},e("input",{value:estimationLabel[0],onChange:function(event){estimationLabel[1](event.target.value);}}))),e("div",{className:"lw-form-grid"},e(Field,{label:"Outer iterations"},e("input",{type:"number",min:1,value:estimationMaxit[0],onChange:function(event){estimationMaxit[1](Number(event.target.value));}})),estimationMethod[0]!=="BAYES"?e(Field,{label:"ETA iterations"},e("input",{type:"number",min:1,value:etaMaxit[0],onChange:function(event){etaMaxit[1](Number(event.target.value));}})):null,e(Field,{label:"Tolerance"},e("input",{type:"number",min:1e-12,step:"any",value:tolerance[0],onChange:function(event){tolerance[1](Number(event.target.value));}})),e(Field,{label:"Parallel cores"},e("input",{type:"number",min:1,max:64,value:estimationCores[0],onChange:function(event){estimationCores[1](Number(event.target.value));}})),e(Field,{label:"Print gradients every N (0 = off)"},e("input",{type:"number",min:0,value:printEvery[0],onChange:function(event){printEvery[1](Number(event.target.value));}})))),
        estimationMethod[0]==="IMP"?e("div",{className:"lw-modal-section lw-modal-section-tinted"},e("h4",null,"Importance sampling"),e("div",{className:"lw-form-grid lw-form-grid-two"},e(Field,{label:"Importance samples"},e("input",{type:"number",min:5,value:nImp[0],onChange:function(event){nImp[1](Number(event.target.value));}})),e(Field,{label:"Random seed"},e("input",{type:"number",min:1,value:methodSeed[0],onChange:function(event){methodSeed[1](Number(event.target.value));}})))):null,
        estimationMethod[0]==="SAEM"?e("div",{className:"lw-modal-section lw-modal-section-tinted"},e("h4",null,"SAEM controls"),e("div",{className:"lw-form-grid"},e(Field,{label:"SAEM iterations"},e("input",{type:"number",min:2,value:nIter[0],onChange:function(event){nIter[1](Number(event.target.value));}})),e(Field,{label:"Burn-in"},e("input",{type:"number",min:0,value:burn[0],onChange:function(event){burn[1](Number(event.target.value));}})),e(Field,{label:"MCMC steps / subject"},e("input",{type:"number",min:1,value:mcmcSteps[0],onChange:function(event){mcmcSteps[1](Number(event.target.value));}})),e(Field,{label:"Random seed"},e("input",{type:"number",min:1,value:methodSeed[0],onChange:function(event){methodSeed[1](Number(event.target.value));}})))):null,
        estimationMethod[0]==="BAYES"?e("div",{className:"lw-modal-section lw-modal-section-tinted"},e("h4",null,"Bayesian sampling"),e("div",{className:"lw-form-grid"},e(Field,{label:"Burn-in"},e("input",{type:"number",min:0,value:nBurn[0],onChange:function(event){nBurn[1](Number(event.target.value));}})),e(Field,{label:"Posterior samples"},e("input",{type:"number",min:1,value:nSample[0],onChange:function(event){nSample[1](Number(event.target.value));}})),e(Field,{label:"Thinning interval"},e("input",{type:"number",min:1,value:nThin[0],onChange:function(event){nThin[1](Number(event.target.value));}})),e(Field,{label:"Random seed"},e("input",{type:"number",min:1,value:methodSeed[0],onChange:function(event){methodSeed[1](Number(event.target.value));}})))):null,
        e("div",{className:"lw-modal-section"},e("h4",null,"Estimation priors"),list(props.model.priors).length?e(React.Fragment,null,e(SimpleTable,{rows:list(props.model.priors),columns:["parameter","distribution","mean","sd","shape","rate"],className:"lw-active-priors"}),e("p",{className:"lw-help-text"},list(props.model.priors).length+" prior"+(list(props.model.priors).length===1?" is":"s are")+" active for this run. Edit priors in the Code tab before submitting to change them.")):e("p",{className:"lw-help-text"},"No parameter priors are active. Add them under Estimation priors in the Code tab to make them part of the reproducible model version.")),
        e("div",{className:"lw-modal-section lw-modal-section-tinted"},e("h4",null,estimationMethod[0]==="BAYES"?"Posterior uncertainty":"Covariance step"),estimationMethod[0]==="BAYES"?e("p",{className:"lw-help-text"},"Posterior SDs, posterior CVs and 95% credible intervals are calculated from the saved samples automatically."):e(React.Fragment,null,e("label",{className:"lw-check"},e("input",{type:"checkbox",checked:covarianceStep[0]&&covarianceSupported,disabled:!covarianceSupported,onChange:function(event){covarianceStep[1](event.target.checked);}})," Run covariance step after estimation"),!covarianceSupported?e("p",{className:"lw-help-text"},"Covariance is available for FO, FOCE, FOCEI, LAPLACE, ITS, IMP and SAEM."):covarianceStep[0]?e(React.Fragment,null,e("div",{className:"lw-form-grid"},e(Field,{label:"Estimator"},e("select",{value:estimationMethod[0]==="FO"?"hessian":covarianceType[0],onChange:function(event){covarianceType[1](event.target.value);}},e("option",{value:"hessian"},"Hessian (R matrix)"),estimationMethod[0]!=="FO"?e("option",{value:"opg"},"Gradient outer product (S matrix)"):null)),e(Field,{label:"Regularization tolerance"},e("input",{type:"number",min:1e-14,step:"any",value:covarianceTolerance[0],onChange:function(event){covarianceTolerance[1](Number(event.target.value));}})),estimationMethod[0]==="SAEM"?e(Field,{label:"Marginal samples"},e("input",{type:"number",min:5,value:covarianceSamples[0],onChange:function(event){covarianceSamples[1](Number(event.target.value));}})):null),e("p",{className:"lw-help-text"},estimationMethod[0]==="IMP"?"The covariance calculation uses deterministic Gauss-Hermite integration when feasible and otherwise reuses the IMP sample budget and seed.":estimationMethod[0]==="SAEM"?"Observed marginal information uses deterministic Gauss-Hermite integration when feasible, with common-random-number importance sampling as the high-dimensional fallback.":"Standard errors, RSEs, covariance and correlation matrices will be saved with the estimation run.")):e("p",{className:"lw-help-text"},"Enable this to calculate parameter uncertainty after the fit."))),
        e("p",{className:"lw-help-text"},"All likelihood, automatic differentiation, ADVAN, matrix-exponential and ODE calculations run in the C++ engine. Queued runs execute in an isolated worker.")),

      e(Modal,{open:simulationModal[0],className:"lw-modal-wide",onClose:function(){simulationModal[1](false);},title:"Create simulation - "+value(props.model.name,"model"),footer:e(React.Fragment,null,e(Button,{className:"lw-button-quiet",onClick:function(){simulationModal[1](false);}},"Cancel"),e(Button,{className:"lw-button-primary",onClick:submitSimulation},"Run simulation"))},
        e("div",{className:"lw-modal-section"},e("div",{className:"lw-form-grid"},e(Field,{label:"Run on"},e("select",{value:value(props.server.queue_id,"local"),onChange:function(event){emit(props,"queue_select",{id:event.target.value});}},list(props.server.queues).map(function(queue){return e("option",{key:queue.id,value:queue.id},queue.name);}))),e(Field,{label:"Label"},e("input",{value:simulationLabel[0],onChange:function(event){simulationLabel[1](event.target.value);}})),e(Field,{label:"Random seed"},e("input",{type:"number",min:1,value:simulationSeed[0],onChange:function(event){simulationSeed[1](Number(event.target.value));}})),e(Field,{label:"Parallel cores"},e("input",{type:"number",min:1,value:simulationCores[0],onChange:function(event){simulationCores[1](Number(event.target.value));}}))),e("div",{className:"lw-form-grid"},e(Field,{label:"Individuals"},e("input",{type:"number",min:1,max:10000,value:simulationSubjects[0],onChange:function(event){simulationSubjects[1](Number(event.target.value));}})),e(Field,{label:"Replications"},e("input",{type:"number",min:1,max:1000,value:simulationReplicates[0],onChange:function(event){simulationReplicates[1](Number(event.target.value));}})),e(Field,{label:"Days (TIME horizon)"},e("input",{type:"number",min:1,max:365,value:simulationDays[0],onChange:function(event){simulationDays[1](Number(event.target.value));}})))),
        e("div",{className:"lw-modal-section lw-modal-section-tinted"},e("h4",null,"Parameters"),e("label",{className:"lw-check"},e("input",{type:"checkbox",checked:simulationUseFit[0]&&props.fit.available,disabled:!props.fit.available,onChange:function(event){simulationUseFit[1](event.target.checked);}})," Use fitted THETA / OMEGA / SIGMA"),!props.fit.available?e("p",{className:"lw-help-text"},"No estimation is loaded; simulation will use the model's initial parameter values."):e("p",{className:"lw-help-text"},"Diagnostics are run separately from the selected estimation run.")),
        e("label",{className:"lw-check lw-design-toggle"},e("input",{type:"checkbox",checked:simulationUseDesign[0],onChange:function(event){simulationUseDesign[1](event.target.checked);}})," Custom dosing / sampling design"),
        simulationUseDesign[0]?e("div",{className:"lw-modal-section lw-simulation-design"},e("h4",null,"Dosing and sampling design"),e("div",{className:"lw-form-grid"},e(Field,{label:"Dosing"},e("select",{value:doseMode[0],onChange:function(event){doseMode[1](event.target.value);}},e("option",{value:"single"},"Single dose"),e("option",{value:"repeat"},"Repeat doses"),e("option",{value:"steady_state"},"Steady state"))),e(Field,{label:"Default dose amount"},e("input",{type:"number",min:0,value:doseAmount[0],onChange:function(event){doseAmount[1](Number(event.target.value));}})),e(Field,{label:"Dose CMT"},e("input",{type:"number",min:1,value:doseCmt[0],onChange:function(event){doseCmt[1](Number(event.target.value));}})),doseMode[0]==="repeat"?e(Field,{label:"Number of doses"},e("input",{type:"number",min:1,max:100,value:doseN[0],onChange:function(event){doseN[1](Number(event.target.value));}})):null,doseMode[0]!=="single"?e(Field,{label:"Dosing interval (h)"},e("input",{type:"number",min:.1,step:.5,value:doseII[0],onChange:function(event){doseII[1](Number(event.target.value));}})):null,e(Field,{label:"Observations / day"},e("input",{type:"number",min:3,max:48,value:obsPerDay[0],onChange:function(event){obsPerDay[1](Number(event.target.value));}}))),e(Field,{label:"Dose amounts (TIME AMT per line, or AMT only)"},e("textarea",{rows:3,value:doseTable[0],placeholder:"0 320\n12 320",onChange:function(event){doseTable[1](event.target.value);}}))):e("p",{className:"lw-help-text"},"The linked dataset structure is retained and resampled to the requested number of individuals.")),

      e(Modal,{open:diagnosticModal[0],className:"lw-modal-wide",onClose:function(){diagnosticModal[1](false);},title:"Run diagnostic",footer:e(React.Fragment,null,e(Button,{className:"lw-button-quiet",onClick:function(){diagnosticModal[1](false);}},"Cancel"),e(Button,{className:"lw-button-primary",disabled:!diagnosticVpc[0]&&!diagnosticNpc[0]&&!diagnosticNpde[0]&&!diagnosticCategorical[0]&&!diagnosticTte[0],onClick:submitDiagnostic},"Run selected"))},
        e("p",{className:"lw-help-text"},"Diagnostics are calculated for the selected estimation run and saved with it."),
        e("div",{className:"lw-choice-row"},
          e("label",{className:"lw-check"},e("input",{type:"checkbox",checked:diagnosticVpc[0],onChange:function(event){diagnosticVpc[1](event.target.checked);}})," Continuous VPC"),
          e("label",{className:"lw-check"},e("input",{type:"checkbox",checked:diagnosticCategorical[0],onChange:function(event){diagnosticCategorical[1](event.target.checked);}})," Categorical VPC"),
          e("label",{className:"lw-check"},e("input",{type:"checkbox",checked:diagnosticTte[0],onChange:function(event){diagnosticTte[1](event.target.checked);}})," Time-to-event VPC"),
          e("label",{className:"lw-check"},e("input",{type:"checkbox",checked:diagnosticNpde[0],onChange:function(event){diagnosticNpde[1](event.target.checked);}})," NPDE"),
          e("label",{className:"lw-check"},e("input",{type:"checkbox",checked:diagnosticNpc[0],onChange:function(event){diagnosticNpc[1](event.target.checked);}})," NPC")),
        e("div",{className:"lw-form-grid"},e(Field,{label:"Simulations"},e("input",{type:"number",min:20,max:10000,value:diagnosticNsim[0],onChange:function(event){diagnosticNsim[1](Number(event.target.value));}})),e(Field,{label:"Random seed"},e("input",{type:"number",min:1,value:diagnosticSeed[0],onChange:function(event){diagnosticSeed[1](Number(event.target.value));}}))),
        diagnosticVpc[0]?e("div",{className:"lw-modal-section lw-modal-section-tinted"},e("h4",null,"Continuous VPC options"),e("label",{className:"lw-check"},e("input",{type:"checkbox",checked:diagnosticPc[0],onChange:function(event){diagnosticPc[1](event.target.checked);}})," Prediction-corrected VPC (DV x PRED / IPRED)"),e(Field,{label:"Stratify by"},e("select",{value:diagnosticStratify[0],onChange:function(event){diagnosticStratify[1](event.target.value);}},[e("option",{key:"none",value:""},"(no stratification)")].concat(list(props.dataset.columns).map(function(column){return e("option",{key:column,value:column},column);})))),e("p",{className:"lw-help-text"},"The VPC tab retains the overall population plot and adds one saved plot per stratum.")):null,
        diagnosticCategorical[0]?e("div",{className:"lw-modal-section lw-modal-section-tinted"},e("h4",null,"Categorical VPC options"),e(Field,{label:"Binary outcome column"},e("select",{value:diagnosticOutcome[0],onChange:function(event){diagnosticOutcome[1](event.target.value);}},list(props.dataset.columns).map(function(column){return e("option",{key:column,value:column},column);}))),e("p",{className:"lw-help-text"},"F/IPRED must be the conditional probability of the non-reference category.")):null,
        diagnosticTte[0]?e("div",{className:"lw-modal-section lw-modal-section-tinted"},e("h4",null,"Time-to-event VPC options"),e(Field,{label:"Event indicator column"},e("select",{value:diagnosticEvent[0],onChange:function(event){diagnosticEvent[1](event.target.value);}},list(props.dataset.columns).map(function(column){return e("option",{key:column,value:column},column);}))),e("p",{className:"lw-help-text"},"F/IPRED is interpreted as a non-negative hazard on the observation-time grid.")):null),

      e(Modal,{open:uncertaintyModal[0],className:"lw-modal-wide",onClose:function(){uncertaintyModal[1](false);},title:"Parameter uncertainty",footer:e(React.Fragment,null,e(Button,{className:"lw-button-quiet",onClick:function(){uncertaintyModal[1](false);}},"Cancel"),e(Button,{className:"lw-button-primary",disabled:!uncertaintyBootstrap[0]&&!uncertaintyProfile[0],onClick:submitUncertainty},"Run selected"))},
        e("div",{className:"lw-choice-row"},e("label",{className:"lw-check"},e("input",{type:"checkbox",checked:uncertaintyBootstrap[0],onChange:function(event){uncertaintyBootstrap[1](event.target.checked);}})," Subject bootstrap"),e("label",{className:"lw-check"},e("input",{type:"checkbox",checked:uncertaintyProfile[0],onChange:function(event){uncertaintyProfile[1](event.target.checked);}})," Profile likelihood")),
        e("div",{className:"lw-form-grid"},uncertaintyBootstrap[0]?e(Field,{label:"Bootstrap replicates"},e("input",{type:"number",min:1,value:uncertaintyReplicates[0],onChange:function(event){uncertaintyReplicates[1](Number(event.target.value));}})):null,e(Field,{label:"Confidence level"},e("input",{type:"number",min:0.5,max:0.999,step:0.01,value:uncertaintyLevel[0],onChange:function(event){uncertaintyLevel[1](Number(event.target.value));}})),e(Field,{label:"Maximum fit iterations"},e("input",{type:"number",min:1,value:uncertaintyMaxit[0],onChange:function(event){uncertaintyMaxit[1](Number(event.target.value));}}))),
        uncertaintyProfile[0]?e("div",{className:"lw-modal-section lw-modal-section-tinted"},e("div",{className:"lw-form-grid"},e(Field,{label:"Grid points / parameter"},e("input",{type:"number",min:3,step:2,value:uncertaintyPoints[0],onChange:function(event){uncertaintyPoints[1](Number(event.target.value));}})),e(Field,{label:"Grid half-width (SE)"},e("input",{type:"number",min:0.1,step:0.5,value:uncertaintySpan[0],onChange:function(event){uncertaintySpan[1](Number(event.target.value));}}))),e(Field,{label:"Parameters (blank = all free)"},e("input",{placeholder:"THETA1, THETA2, SIGMA1",value:uncertaintyParameters[0],onChange:function(event){uncertaintyParameters[1](event.target.value);}})),e("p",{className:"lw-help-text"},"Each grid point fixes one parameter and re-estimates the remaining free parameters.")):null),

      e(Modal,{open:scmModal[0],className:"lw-modal-wide",onClose:function(){scmModal[1](false);},title:"Stepwise covariate modelling",footer:e(React.Fragment,null,e(Button,{className:"lw-button-quiet",onClick:function(){scmModal[1](false);}},"Cancel"),e(Button,{className:"lw-button-primary",onClick:function(){emit(props,"run_scm",{candidates:scmCandidates[0],direction:scmDirection[0],pForward:Number(scmForward[0]),pBackward:Number(scmBackward[0]),maxSteps:Number(scmMaxSteps[0]),maxit:Number(scmMaxit[0]),label:scmLabel[0]});scmModal[1](false);}},"Run SCM"))},e(Field,{label:"Candidate relationships (parameter,covariate,form,reference,category)"},e("textarea",{rows:6,value:scmCandidates[0],onChange:function(event){scmCandidates[1](event.target.value);}})),e("div",{className:"lw-form-grid"},e(Field,{label:"Direction"},e("select",{value:scmDirection[0],onChange:function(event){scmDirection[1](event.target.value);}},e("option",{value:"forward"},"Forward"),e("option",{value:"backward"},"Backward"),e("option",{value:"both"},"Forward + backward"))),e(Field,{label:"Forward p"},e("input",{type:"number",min:0.0001,max:0.5,step:0.01,value:scmForward[0],onChange:function(event){scmForward[1](Number(event.target.value));}})),e(Field,{label:"Backward p"},e("input",{type:"number",min:0.0001,max:0.5,step:0.01,value:scmBackward[0],onChange:function(event){scmBackward[1](Number(event.target.value));}})),e(Field,{label:"Maximum steps"},e("input",{type:"number",min:1,value:scmMaxSteps[0],onChange:function(event){scmMaxSteps[1](Number(event.target.value));}})),e(Field,{label:"Fit iterations"},e("input",{type:"number",min:1,value:scmMaxit[0],onChange:function(event){scmMaxit[1](Number(event.target.value));}})),e(Field,{label:"New version label"},e("input",{value:scmLabel[0],onChange:function(event){scmLabel[1](event.target.value);}}))),e("p",{className:"lw-help-text"},"Forms are continuous, power, or categorical. The accepted SCM model is saved as a new version and estimation run.")),

      e(Modal,{open:controlModal[0],className:"lw-modal-wide",onClose:function(){controlModal[1](false);},title:"Load NONMEM control stream",footer:e(React.Fragment,null,e(Button,{className:"lw-button-quiet",onClick:function(){controlModal[1](false);}},"Cancel"),e(Button,{className:"lw-button-primary",disabled:!controlFile[0]||(controlNewProject[0]&&!controlProjectName[0].trim()),onClick:function(){emit(props,"control_import",{text:controlFile[0]&&controlFile[0].text,fileName:controlFile[0]&&controlFile[0].name,dataText:controlData[0]&&controlData[0].text,dataName:controlData[0]&&controlData[0].name,newProject:controlNewProject[0],projectName:controlProjectName[0],label:controlLabel[0]});controlModal[1](false);}},"Import"))},e(Field,{label:"Control stream (.ctl, .mod)"},e("input",{type:"file",accept:".ctl,.mod,.txt,text/plain",onChange:readControlFile})),e("p",{className:"lw-help-text"},controlFile[0]?"Loaded "+controlFile[0].name:"Unsupported records are preserved and reported instead of silently discarded."),e(Field,{label:"Dataset (optional; otherwise keep current dataset)"},e("input",{type:"file",accept:".csv,.txt,.dat,.tsv,text/csv,text/plain",onChange:readControlData})),e("label",{className:"lw-check"},e("input",{type:"checkbox",checked:controlNewProject[0],onChange:function(event){controlNewProject[1](event.target.checked);}})," Create a new project"),controlNewProject[0]?e(Field,{label:"Project name"},e("input",{value:controlProjectName[0],onChange:function(event){controlProjectName[1](event.target.value);}})):null,e(Field,{label:"Model version label"},e("input",{value:controlLabel[0],onChange:function(event){controlLabel[1](event.target.value);}}))),

      e(Modal,{open:exportModal[0],onClose:function(){exportModal[1](false);},title:"Export NONMEM control stream",footer:e(React.Fragment,null,e(Button,{className:"lw-button-quiet",onClick:function(){exportModal[1](false);}},"Cancel"),e(Button,{className:"lw-button-primary",onClick:function(){emit(props,"control_export",{name:exportName[0],dataPath:exportDataPath[0]});exportModal[1](false);}},"Export"))},e(Field,{label:"File name"},e("input",{value:exportName[0],onChange:function(event){exportName[1](event.target.value);}})),e(Field,{label:"$DATA path"},e("input",{value:exportDataPath[0],onChange:function(event){exportDataPath[1](event.target.value);}})),e("p",{className:"lw-help-text"},"The file is written to the workspace exports directory. Preserved records from an imported stream are retained.")),

      e(Modal,{open:!!deleteModal[0],onClose:function(){deleteModal[1](null);},title:deleteModal[0]==="project"?"Delete project":deleteModal[0]==="run"?"Delete model run":"Delete model version",footer:e(React.Fragment,null,e(Button,{className:"lw-button-quiet",onClick:function(){deleteModal[1](null);}},"Cancel"),e(Button,{className:"lw-button-danger",onClick:function(){if(deleteModal[0]==="project")emit(props,"project_delete",{id:workspace.current});else if(deleteModal[0]==="run")emit(props,"project_delete_run",{id:workspace.current,run:workspace.current_run});else emit(props,"project_delete_snapshot",{id:workspace.current,snapshot:workspace.current_version});deleteModal[1](null);}},"Delete"))},e("div",{className:"lw-destructive-note"},e("strong",null,"This action cannot be undone."),e("p",null,deleteModal[0]==="project"?"The project and all of its model versions and runs will be removed.":deleteModal[0]==="run"?"The selected estimation or simulation run and its saved diagnostics will be removed.":"The selected model version and all of its runs will be removed."))));
  }

  function ResultsPanel(props) {
    var tabState = React.useState("parameters"), tab = tabState[0];
    var sections = React.useState({ summary: true, parameters: true, gof: true, eta: true, vpc: true, narrative_stub: true });
    var reportName = React.useState("report_" + new Date().toISOString().slice(0,10).replace(/-/g,""));
    var parameters = props.fit.available ? props.fit.parameters : [];
    var covariance = props.fit.covariance || {requested:false,status:"not_requested"};
    var posterior = props.fit.posterior || {available:false,parameters:[]};
    var resultTabs = [{id:"parameters",label:"Parameters"}];
    if (covariance.requested) resultTabs.push({id:"covariance",label:"Covariance"});
    if (posterior.available) resultTabs.push({id:"posterior",label:"Posterior"});
    resultTabs.push({id:"run",label:"Run info"},{id:"report",label:"Report"});
    React.useEffect(function(){if(tab==="covariance"&&!covariance.requested)tabState[1]("parameters");},[covariance.requested]);
    React.useEffect(function(){if(tab==="posterior"&&!posterior.available)tabState[1]("parameters");},[posterior.available]);
    function toggle(section) { var next = Object.assign({}, sections[0]); next[section] = !next[section]; sections[1](next); }
    return e(Panel, { title: "Results", className: "lw-results-panel", bodyClass: "lw-results-body" },
      e(Tabs, { value: tab, onChange: tabState[1], items: resultTabs }),
      tab === "parameters" ? e("div", { className: "lw-results-tab" },
        props.result && props.result.kind === "comparison" ? e(React.Fragment, null,
          e("div", { className: "lw-fit-summary" }, e("strong", null, "Run comparison"), e("span", null, "Side-by-side estimates")),
          e(SimpleTable, { rows: props.result.parameters })) : props.fit.available ? e(React.Fragment, null,
          e("div", { className: "lw-fit-summary" }, e("strong", null, props.fit.method + " fit"), e("span", null, "OFV " + formatNumber(props.fit.objective)), e("span", null, props.fit.convergence === 0 ? "Converged" : "Code " + props.fit.convergence)),
          e(SimpleTable, { rows: parameters, columns: posterior.available?["name","value","posterior_sd","posterior_cv","median","lower_95","upper_95"]:covariance.status==="completed"?["name","value","se","rse"]:["name","value"] })) : e(Empty, { title: "No estimates", detail: "Open or run an estimation." })) : null,
      tab === "covariance" ? e("div", { className: "lw-results-tab lw-covariance-tab" }, covariance.status==="failed"?e("div",{className:"lw-destructive-note"},e("strong",null,"Covariance step failed"),e("p",null,value(covariance.error,"No covariance result was produced."))):e(React.Fragment,null,e("div",{className:"lw-fit-summary"},e("strong",null,(covariance.type||"covariance").toUpperCase()+" covariance"),e("span",null,"Condition "+formatNumber(covariance.condition)),e("span",null,"Regularization "+formatNumber(covariance.regularization))),e("h4",null,"Covariance matrix"),e(SimpleTable,{rows:covariance.covariance,empty:"No free parameters"}),e("h4",null,"Correlation matrix"),e(SimpleTable,{rows:covariance.correlation,empty:"No correlations"}))) : null,
      tab === "posterior" ? e("div", { className: "lw-results-tab lw-covariance-tab" }, e("div",{className:"lw-fit-summary"},e("strong",null,"Bayesian posterior uncertainty"),e("span",null,formatNumber(posterior.samples)+" saved samples"),e("span",null,"Outer acceptance "+formatNumber(posterior.outer_acceptance)),e("span",null,"ETA acceptance "+formatNumber(posterior.eta_acceptance))),e(SimpleTable,{rows:posterior.parameters,columns:["name","mean","posterior_sd","posterior_cv","median","lower_95","upper_95"]}),e("h4",null,"Posterior covariance"),e(SimpleTable,{rows:posterior.covariance,empty:"No posterior covariance"}),e("h4",null,"Posterior correlation"),e(SimpleTable,{rows:posterior.correlation,empty:"No posterior correlations"})) : null,
      tab === "run" ? e("div", { className: "lw-results-tab lw-run-info" }, props.result && props.result.kind === "comparison" ? e(SimpleTable, { rows: props.result.runs }) : props.fit.available ?
        e(SimpleTable, { rows: Object.keys(props.fit.run_info || {}).map(function (key) { return { Item: key, Value: props.fit.run_info[key] }; }), columns: ["Item","Value"] }) : e(Empty, { title: "No run information", detail: "Run an estimation first." })) : null,
      tab === "report" ? e("div", { className: "lw-results-tab lw-report-controls" },
        e("p", null, "PDF report with JSON manifest for interpretation and audit."),
        Object.keys(sections[0]).map(function (section) { return e("label", { key: section }, e("input", { type: "checkbox", checked: sections[0][section], onChange: function () { toggle(section); } }), section.replace(/_/g," ")); }),
        e(Field, { label: "Filename (no extension)" }, e("input", { value: reportName[0], onChange: function (event) { reportName[1](event.target.value); } })),
        e(Button, { className: "lw-button-primary", disabled: !props.fit.available, onClick: function () { emit(props, "report", { name: reportName[0], sections: Object.keys(sections[0]).filter(function (key) { return sections[0][key]; }) }); } }, "Generate PDF"),
        props.report && props.report.pdf ? e("div", { className: "lw-report-status" }, e("strong", null, "Report created"), e("span", null, props.report.pdf), props.report.json ? e("span", null, props.report.json) : null) : null) : null);
  }

  function ComparisonPlots(props) {
    var plots = props.plots || {};
    var labels = { gof:"Goodness-of-fit plots", vpc:"Visual predictive checks", vpc_categorical:"Categorical VPCs", vpc_tte:"Time-to-event VPCs", npde:"NPDE plots", npc:"NPC plots" };
    var kinds = ["gof","vpc","vpc_categorical","vpc_tte","npde","npc"].filter(function(kind){return list(plots[kind]).length === 2;});
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

  function HomePage(props) {
    var centerTab = React.useState("code"), tab = centerTab[0];
    var saveModal = React.useState(false), saveLabel = React.useState("Model revision");
    var comparisonModal = React.useState(false);
    var workspace = props.workspace || {};
    var diagnostics = props.diagnostics || {};
    var diagnosticAvailability = diagnostics.available || {};
    var centerTabs = [{id:"code",label:"Code"},{id:"gof",label:"GOF"}];
    if (diagnosticAvailability.vpc) centerTabs.push({id:"vpc",label:"VPC"});
    if (diagnosticAvailability.vpc_categorical) centerTabs.push({id:"vpc_categorical",label:"Cat VPC"});
    if (diagnosticAvailability.vpc_tte) centerTabs.push({id:"vpc_tte",label:"TTE VPC"});
    if (diagnosticAvailability.npde) centerTabs.push({id:"npde",label:"NPDE"});
    if (diagnosticAvailability.npc) centerTabs.push({id:"npc",label:"NPC"});
    if (diagnosticAvailability.bootstrap) centerTabs.push({id:"bootstrap",label:"Bootstrap"});
    if (diagnosticAvailability.profile) centerTabs.push({id:"profile",label:"Profile"});
    if (diagnosticAvailability.scm) centerTabs.push({id:"scm",label:"SCM"});
    React.useEffect(function(){if(!centerTabs.some(function(item){return item.id===centerTab[0];}))centerTab[1]("code");},[!!diagnosticAvailability.vpc,!!diagnosticAvailability.vpc_categorical,!!diagnosticAvailability.vpc_tte,!!diagnosticAvailability.npde,!!diagnosticAvailability.npc,!!diagnosticAvailability.bootstrap,!!diagnosticAvailability.profile,!!diagnosticAvailability.scm]);
    React.useEffect(function(){if(props.result&&props.result.kind==="comparison")comparisonModal[1](true);},[props.result&&props.result.comparison_id]);
    function closeComparison(){comparisonModal[1](false);emit(props,"comparison_close");}
    function selectCenterTab(next) {
      centerTab[1](next);
      if (next === "gof" && props.fit.available && !props.fit.gof_loaded) emit(props,"load_payload",{kind:"gof"});
      if (["vpc","npde","npc","vpc_categorical","vpc_tte","bootstrap","profile","scm"].indexOf(next) >= 0 && diagnosticAvailability[next] && !diagnostics[next]) emit(props,"load_payload",{kind:next});
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
          centerTabs.map(function(item){return e("div",{key:item.id,className:"lw-center-tab "+(tab===item.id?"":"lw-center-tab-hidden"),"aria-hidden":tab===item.id?"false":"true"},item.id==="code"?e(ModelEditor,props):e(CachedDiagnosticsPane,Object.assign({},props,{tab:item.id})));}))),
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

  function LibeRWorkbench(props) {
    var ribbon = React.useState("home"), theme = React.useState(function () { try { return window.localStorage.getItem("liberationDarkTheme") !== "0"; } catch (error) { return true; } });
    function toggleTheme() { var next = !theme[0]; theme[1](next); try { window.localStorage.setItem("liberationDarkTheme", next ? "1" : "0"); } catch (error) {} }
    function changeRibbon(next){ribbon[1](next);emit(props,"page_change",{page:next});}
    var ribbonItems = [{id:"home",label:"Home"},{id:"jobs",label:"Jobs"},{id:"data",label:"Data"}];
    return e("div", { className: "lw-legacy-shell " + (theme[0] ? "theme-dark" : "theme-light") },
      e("header", { className: "lw-app-header" }, e("div", null, props.server&&props.server.icon?e("img",{className:"lw-app-icon",src:props.server.icon,alt:""}):null, e("strong", null, "LibeRation"), e("span", {className:"lw-app-version"}, "v"+value(props.server&&props.server.package_version,"unknown")), e("span", null, "Population PK/PD modelling")),
        e("div", { className: "lw-app-header-right" }, e("label", { className: "lw-theme-toggle" }, e("span", null, theme[0] ? "Dark" : "Light"), e("input", { type: "checkbox", checked: theme[0], onChange: toggleTheme }), e("i", null)), e("span", { className: "lw-workspace-path" }, value(props.workspace && props.workspace.path, "No workspace selected")))),
      e("div", { className: "lw-ribbon" }, e(Tabs, { value: ribbon[0], onChange: changeRibbon, items: ribbonItems, className: "lw-ribbon-tabs" })),
      e(LogBanner, props),
      e("div", { className: "lw-page-host" }, ribbon[0] === "home" ? e(HomePage, props) : ribbon[0] === "jobs" ? e(JobsPage, props) : e(DataPage, props)),
      e("footer", { className: "lw-statusbar" }, e("span", null, e(StatusDot, { status: props.result && props.result.status === "error" ? "error" : "ready" }), value(props.result && props.result.message, "Engine ready")), e("span", null, value(props.server && props.server.mode, "local") + " | C++17 | Exact AD")));
  }

  reactR.reactWidget("liberWorkbench", "output", { LibeRWorkbench: LibeRWorkbench }, {});
}());
