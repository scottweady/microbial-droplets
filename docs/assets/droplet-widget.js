// Interactive Fig. 7: stability diagram in the (m, beta*Ra) plane (a), the
// perturbed droplet at time t (b), and the perturbation amplitude (c),
//   eps(t) = eps0 exp[(sigma_g^m - 1/2) t + beta*Ra sigma_b^m (e^t - 1)],
// which is exact for the linearized problem with R(t) = e^{t/2}. The edge is
// the image of the unit circle under f(zeta) = R (zeta + eps zeta^{m+1}),
// which is one-to-one only while eps (m+1) < 1: at eps (m+1) = 1 the edge
// forms cusps and would then fold over itself, so time stops there.
(function () {
  "use strict";

  var P = window.PL;
  var root = document.getElementById("drop-widget");
  if (!P || !root) return;

  var W = 420, H = 300, M = { l: 46, r: 14, t: 12, b: 30 };
  var EPS0 = 0.02, TMAX = 4, RAMAX = 200;

  var sM = root.querySelector("#dr-m"), sR = root.querySelector("#dr-ra"), sT = root.querySelector("#dr-t");
  var oM = root.querySelector("#dr-m-val"), oR = root.querySelector("#dr-ra-val"), oT = root.querySelector("#dr-t-val");
  var rd = {
    s: root.querySelector("#dr-s"), crit: root.querySelector("#dr-crit"),
    eps: root.querySelector("#dr-eps"), R: root.querySelector("#dr-R"),
    tstar: root.querySelector("#dr-tstar"), fold: root.querySelector("#dr-fold")
  };

  // (a) stability diagram
  var A = P.axes(root.querySelector("#dr-a"), { W: W, H: H, M: M, x: [2, 30], y: [0, RAMAX],
    xt: [2, 5, 10, 15, 20, 25, 30], yt: [0, 50, 96, 150, 200] });
  var cm = [], cr = [];
  for (var m = 3; m <= 30; m++) { cm.push(m); cr.push(P.raCrit(m)); }
  var poly = "M" + A.sx(2) + "," + A.sy(RAMAX) + "L" + A.sx(3) + "," + A.sy(RAMAX) +
    P.path(cm, cr, A.sx, A.sy).replace(/^M/, "L") + "L" + A.sx(30) + "," + A.sy(RAMAX) + "Z";
  P.el("path", { d: poly, "class": "dp-unstable" }, A.g);
  P.el("line", { x1: A.sx(2), x2: A.sx(30), y1: A.sy(96), y2: A.sy(96), "class": "rw-dash" }, A.g);
  P.el("path", { d: P.path(cm, cr, A.sx, A.sy), "class": "dr-crit" }, A.g);
  cm.forEach(function (m, k) { P.el("circle", { cx: A.sx(m), cy: A.sy(cr[k]), r: 3, "class": "dr-critdot" }, A.g); });
  P.el("text", { x: A.sx(16), y: A.sy(165), "text-anchor": "middle", "class": "sf-tag" }, A.g).textContent = "unstable";
  P.el("text", { x: A.sx(16), y: A.sy(45), "text-anchor": "middle", "class": "sf-tag" }, A.g).textContent = "stable";
  var mark = P.el("circle", { r: 6, "class": "rw-dot" }, A.g);
  var svgA = root.querySelector("#dr-a");

  function pick(ev) {
    var pt = svgA.createSVGPoint(); pt.x = ev.clientX; pt.y = ev.clientY;
    var p = pt.matrixTransform(svgA.getScreenCTM().inverse());
    sM.value = Math.round(Math.min(30, Math.max(2, A.sx.inv(p.x))));
    sR.value = Math.round(Math.min(RAMAX, Math.max(0, A.sy.inv(p.y))));
    update();
  }
  var dragging = false;
  svgA.addEventListener("pointerdown", function (e) { dragging = true; svgA.setPointerCapture(e.pointerId); pick(e); });
  svgA.addEventListener("pointermove", function (e) { if (dragging) pick(e); });
  svgA.addEventListener("pointerup", function () { dragging = false; });

  // (b) droplet at its true size; the view is fixed for each (m, beta*Ra)
  var svgB = root.querySelector("#dr-b");
  var c0 = P.el("circle", { cx: 150, cy: 150, "class": "dr-circle" }, svgB);
  var shape = P.el("path", { "class": "dr-shape" }, svgB);
  var c0lab = P.el("text", { x: 150, "text-anchor": "middle", "class": "sf-tag" }, svgB);
  c0lab.textContent = "R = 1";
  var scaleLab = P.el("text", { x: 296, y: 294, "text-anchor": "end", "class": "sf-tag" }, svgB);

  // (c) amplitude
  var C = P.axes(root.querySelector("#dr-c"), { W: W, H: H, M: M, x: [0, TMAX], y: [1e-3, 1e2], logy: true,
    xt: [0, 1, 2, 3, 4], yt: [1e-3, 1e-2, 1e-1, 1, 10, 100],
    yf: function (v) {
      var e = Math.round(Math.log10(v));
      if (e === 0) return "1";
      if (e === 1) return "10";
      return "10" + (e < 0 ? "⁻" : "") + "⁰¹²³".charAt(Math.abs(e));
    } });
  P.el("line", { x1: C.sx(0), x2: C.sx(TMAX), y1: C.sy(1), y2: C.sy(1), "class": "rw-dash" }, C.g);
  var foldLine = P.el("line", { x1: C.sx(0), x2: C.sx(TMAX), "class": "dr-foldline" }, C.g);
  var foldLab = P.el("text", { x: C.sx(0) + 6, "class": "sf-tag" }, C.g);
  var lin = P.el("path", { "class": "dr-lin" }, C.g);
  var amp = P.el("path", { "class": "rw-cur" }, C.g);
  var tline = P.el("line", { "class": "rw-guide" }, C.g);
  var tdot = P.el("circle", { r: 5, "class": "rw-dot" }, C.g);

  var lastKey = "", tEnd = TMAX, scale = 1;

  function update() {
    var m = +sM.value, ra = +sR.value;
    var g = P.sigmaG(m) - 0.5, b = P.sigmaB(m);
    var s = g + ra * b;
    var foldRatio = 1 / ((m + 1) * EPS0);              // eps (m+1) = 1
    function ratio(tt) { return Math.exp(g * tt + ra * b * (Math.exp(tt) - 1)); }

    var key = m + "," + ra;
    if (key !== lastKey) {
      lastKey = key;
      // time at which the edge first folds, if before TMAX
      tEnd = TMAX;
      var prev = 0;
      for (var i = 1; i <= 800; i++) {
        var ti = TMAX * i / 800;
        if (ratio(ti) >= foldRatio) {
          var lo = prev, hi = ti;
          for (var k = 0; k < 50; k++) { var mid = 0.5 * (lo + hi); if (ratio(mid) >= foldRatio) hi = mid; else lo = mid; }
          tEnd = lo;
          break;
        }
        prev = ti;
      }
      sT.max = tEnd.toFixed(3);
      if (+sT.value > tEnd) sT.value = tEnd;
      // view: the largest droplet over [0, tEnd] fills the panel
      var rMax = Math.exp(tEnd / 2) * (1 + EPS0 * Math.max(1, ratio(tEnd)));
      scale = 130 / rMax;
      c0.setAttribute("r", scale);
      c0lab.setAttribute("y", 150 - scale - 4);
      c0lab.style.display = scale > 14 ? "" : "none";

      var ts = [], es = [], ls = [];
      for (var j = 0; j <= 200; j++) {
        var tt = tEnd * j / 200; ts.push(tt);
        es.push(Math.min(Math.max(ratio(tt), 1e-4), 1e4));
      }
      for (var j2 = 0; j2 <= 200; j2++) ls.push(Math.min(Math.max(Math.exp(s * TMAX * j2 / 200), 1e-4), 1e4));
      amp.setAttribute("d", P.path(ts, es, C.sx, C.sy));
      lin.setAttribute("d", P.path(ls.map(function (v, j) { return TMAX * j / 200; }), ls, C.sx, C.sy));
      foldLine.setAttribute("y1", C.sy(foldRatio)); foldLine.setAttribute("y2", C.sy(foldRatio));
      foldLab.setAttribute("y", C.sy(foldRatio) - 6);
      foldLab.textContent = "edge folds, ε(m+1) = 1";
    }

    var t = +sT.value;
    oM.textContent = m; oR.textContent = ra; oT.textContent = t.toFixed(2);
    mark.setAttribute("cx", A.sx(m)); mark.setAttribute("cy", A.sy(ra));

    var e = ratio(t), ec = Math.min(Math.max(e, 1e-4), 1e4);
    tdot.setAttribute("cx", C.sx(t)); tdot.setAttribute("cy", C.sy(ec));
    tline.setAttribute("x1", C.sx(t)); tline.setAttribute("x2", C.sx(t));
    tline.setAttribute("y1", C.sy(ec)); tline.setAttribute("y2", H - M.b);

    var R = Math.exp(t / 2), eps = Math.min(EPS0 * e, 1 / (m + 1)), d = "";
    for (var n = 0; n <= 1440; n++) {
      var th = 2 * Math.PI * n / 1440;
      var rr = 1 + eps * Math.cos(m * th), tt2 = eps * Math.sin(m * th);
      var X = R * (rr * Math.cos(th) - tt2 * Math.sin(th)), Y = R * (rr * Math.sin(th) + tt2 * Math.cos(th));
      d += (n ? "L" : "M") + (150 + scale * X).toFixed(2) + "," + (150 - scale * Y).toFixed(2);
    }
    shape.setAttribute("d", d + "Z");
    scaleLab.textContent = "R(t) = " + P.fmt(R, 2);

    rd.s.textContent = P.fmt(s, 3) + (s > 0 ? " (unstable)" : " (stable)");
    rd.crit.textContent = m === 2 ? "none (stable for all βRa)" : P.fmt(P.raCrit(m), 1);
    rd.eps.textContent = P.fmt(e, e < 10 ? 3 : 1);
    rd.R.textContent = P.fmt(R, 2) + "  (βRa R² = " + P.fmt(ra * R * R, 0) + ")";
    if (b <= 0 || ra === 0) rd.tstar.textContent = "never";
    else {
      var ts2 = Math.log(-g / (ra * b));
      rd.tstar.textContent = ts2 <= 0 ? "immediately" : "t = " + P.fmt(ts2, 2);
    }
    rd.fold.textContent = tEnd < TMAX ? "t = " + P.fmt(tEnd, 2) : "not before t = " + TMAX;
  }

  [sM, sR].forEach(function (s) { s.addEventListener("input", update); });
  sT.addEventListener("input", update);
  P.player(root.querySelector("#dr-play"), sT, update, 30);
  update();
})();
