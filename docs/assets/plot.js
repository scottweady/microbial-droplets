// Shared helpers for the interactive figures: SVG elements, axes, the
// stability coefficients of the paper, and play buttons for sliders.
(function () {
  "use strict";

  var NS = "http://www.w3.org/2000/svg";
  var uid = 0;

  function el(name, attrs, parent) {
    var e = document.createElementNS(NS, name);
    for (var a in attrs) e.setAttribute(a, attrs[a]);
    if (parent) parent.appendChild(e);
    return e;
  }

  function scale(d0, d1, r0, r1) {
    var f = function (v) { return r0 + (v - d0) / (d1 - d0) * (r1 - r0); };
    f.inv = function (p) { return d0 + (p - r0) / (r1 - r0) * (d1 - d0); };
    return f;
  }

  function path(xs, ys, sx, sy) {
    var s = "";
    for (var i = 0; i < xs.length; i++) {
      if (!isFinite(ys[i])) continue;
      s += (s ? "L" : "M") + sx(xs[i]).toFixed(1) + "," + sy(ys[i]).toFixed(1);
    }
    return s;
  }

  function fmt(v, d) {
    if (v === 0) return "0";
    var s = (+v.toFixed(d === undefined ? 2 : d)).toString();
    return s.replace("-", "−");
  }

  // Frame, ticks and tick labels. Returns scales and a clipped group for data.
  //   o = { W, H, M: {l,r,t,b}, x: [x0,x1], y: [y0,y1], xt, yt, xf, yf, logx, logy }
  function axes(svg, o) {
    while (svg.firstChild) svg.removeChild(svg.firstChild);
    var M = o.M, W = o.W, H = o.H;
    var tx = o.logx ? Math.log10 : function (v) { return v; };
    var ty = o.logy ? Math.log10 : function (v) { return v; };
    var sx0 = scale(tx(o.x[0]), tx(o.x[1]), M.l, W - M.r);
    var sy0 = scale(ty(o.y[0]), ty(o.y[1]), H - M.b, M.t);
    var sx = function (v) { return sx0(tx(v)); };
    var sy = function (v) { return sy0(ty(v)); };
    sx.inv = function (p) { var v = sx0.inv(p); return o.logx ? Math.pow(10, v) : v; };
    sy.inv = function (p) { var v = sy0.inv(p); return o.logy ? Math.pow(10, v) : v; };

    var id = "pl-clip-" + (++uid);
    var defs = el("defs", {}, svg);
    var cp = el("clipPath", { id: id }, defs);
    el("rect", { x: M.l, y: M.t, width: W - M.l - M.r, height: H - M.t - M.b }, cp);

    var ax = el("g", { "class": "rw-axes" }, svg);
    var xf = o.xf || function (v) { return fmt(v); };
    var yf = o.yf || function (v) { return fmt(v); };
    (o.xt || []).forEach(function (v) {
      var X = sx(v);
      el("line", { x1: X, x2: X, y1: H - M.b, y2: H - M.b - 5, "class": "rw-tick" }, ax);
      el("line", { x1: X, x2: X, y1: M.t, y2: M.t + 5, "class": "rw-tick" }, ax);
      el("text", { x: X, y: H - M.b + 17, "text-anchor": "middle" }, ax).textContent = xf(v);
    });
    (o.yt || []).forEach(function (v) {
      var Y = sy(v);
      el("line", { x1: M.l, x2: M.l + 5, y1: Y, y2: Y, "class": "rw-tick" }, ax);
      el("line", { x1: W - M.r, x2: W - M.r - 5, y1: Y, y2: Y, "class": "rw-tick" }, ax);
      el("text", { x: M.l - 7, y: Y + 4.5, "text-anchor": "end" }, ax).textContent = yf(v);
    });
    var g = el("g", { "clip-path": "url(#" + id + ")" }, svg);
    el("rect", { x: M.l, y: M.t, width: W - M.l - M.r, height: H - M.t - M.b, "class": "rw-frame" }, svg);
    return { sx: sx, sy: sy, g: g, ax: ax };
  }

  // Generalized eigenvalue lambda_m^m = Gamma(m+1/2) Gamma(1/2) / Gamma(m+1)
  // (recurrence from lambda_0^0 = pi), and lambda_1^0 = 4/pi.
  var LAM10 = 4 / Math.PI;
  function lamMM(m) {
    var v = Math.PI;
    for (var k = 1; k <= m; k++) v *= (k - 0.5) / k;
    return v;
  }
  // Stability coefficients (paper, Sec. 6) with beta*Ra = 1 for sigma_b.
  function sigmaG(m) { return 0.5 - m * LAM10 * lamMM(m) / 4; }
  function sigmaB(m) {
    if (m < 2) return 0;
    return (m * LAM10 * lamMM(m) / 4 - 9 / (2 * (m - 1) * (2 * m + 1))) / 96;
  }
  // Short-time critical beta*Ra for mode m >= 3, (1/2 - sigma_g) / sigma_b.
  function raCrit(m) { return (0.5 - sigmaG(m)) / sigmaB(m); }

  // Play button that sweeps a range input from its current value to max (or
  // from min if at max), calling onInput on each step.
  function player(btn, slider, onInput, stepsPerSec) {
    var timer = null;
    var rate = stepsPerSec || 30;
    function stop() { if (timer) { clearInterval(timer); timer = null; btn.textContent = "▶"; } }
    btn.addEventListener("click", function () {
      if (timer) { stop(); return; }
      var lo = +slider.min, hi = +slider.max, st = +slider.step || 1;
      if (+slider.value >= hi - st / 2) slider.value = lo;
      btn.textContent = "❚❚";
      timer = setInterval(function () {
        var v = +slider.value + st;
        if (v > hi + st / 2) { stop(); return; }
        slider.value = v;
        onInput();
      }, 1000 / rate);
    });
    slider.addEventListener("pointerdown", stop);
    return stop;
  }

  window.PL = {
    el: el, scale: scale, path: path, fmt: fmt, axes: axes,
    LAM10: LAM10, lamMM: lamMM, sigmaG: sigmaG, sigmaB: sigmaB, raCrit: raCrit,
    player: player
  };
})();
