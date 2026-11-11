// Interactive Fig. 9: the matching constant kappa against the effective
// domain depth |H| (a), and the m -> infinity critical Rayleigh number
// Ra_{*,kappa} = 96 / (1 - 6 kappa) for beta = 1 (b). Data: assets/kappa-data.js.
(function () {
  "use strict";

  var D = window.KAPPA, P = window.PL;
  var root = document.getElementById("kappa-widget");
  if (!D || !P || !root) return;

  var W = 420, H = 300, M = { l: 46, r: 14, t: 12, b: 30 };
  var H0 = 7 / 8, H1 = 7;
  var hs = [], ks = [];
  D.H.forEach(function (h, i) { if (h >= H0 - 1e-9 && h <= H1 + 0.02) { hs.push(h); ks.push(D.kappa[i]); } });
  var ras = ks.map(function (k) { return 96 / (1 - 6 * k); });

  var slider = root.querySelector("#kp-h"), out = root.querySelector("#kp-h-val");
  var rdK = root.querySelector("#kp-k"), rdR = root.querySelector("#kp-ra");
  slider.max = hs.length - 1;
  var k0 = 0;
  hs.forEach(function (h, i) { if (Math.abs(h - 1) < Math.abs(hs[k0] - 1)) k0 = i; });
  slider.value = k0;

  var A = P.axes(root.querySelector("#kp-a"), { W: W, H: H, M: M, x: [H0, H1], y: [-2, 0.5],
    xt: [1, 2, 3, 4, 5, 6, 7], yt: [-2, -1.5, -1, -0.5, 0, 0.5] });
  P.el("line", { x1: A.sx(H0), x2: A.sx(H1), y1: A.sy(0), y2: A.sy(0), "class": "rw-dash" }, A.g);
  var asx = [], asy = [];
  for (var i = 0; i <= 100; i++) { var h = H0 + (H1 - H0) * i / 100; asx.push(h); asy.push(1 - Math.log(2 * h)); }
  P.el("path", { d: P.path(asx, asy, A.sx, A.sy), "class": "kp-asym" }, A.g);
  P.el("path", { d: P.path(hs, ks, A.sx, A.sy), "class": "dp-b" }, A.g);
  var dA = P.el("circle", { r: 5, "class": "rw-dot kp-dot" }, A.g);

  var B = P.axes(root.querySelector("#kp-b"), { W: W, H: H, M: M, x: [H0, H1], y: [0, 140],
    xt: [1, 2, 3, 4, 5, 6, 7], yt: [0, 20, 40, 60, 80, 96, 120, 140] });
  P.el("line", { x1: B.sx(H0), x2: B.sx(H1), y1: B.sy(96), y2: B.sy(96), "class": "rw-dash" }, B.g);
  P.el("path", { d: P.path(hs, ras, B.sx, B.sy), "class": "dp-b" }, B.g);
  var gB = P.el("line", { "class": "rw-guide kp-guide" }, B.g);
  var dB = P.el("circle", { r: 5, "class": "rw-dot kp-dot" }, B.g);

  function update() {
    var k = +slider.value;
    out.textContent = P.fmt(hs[k], 2);
    dA.setAttribute("cx", A.sx(hs[k])); dA.setAttribute("cy", A.sy(ks[k]));
    dB.setAttribute("cx", B.sx(hs[k])); dB.setAttribute("cy", B.sy(ras[k]));
    gB.setAttribute("x1", B.sx(hs[k])); gB.setAttribute("x2", B.sx(hs[k]));
    gB.setAttribute("y1", B.sy(ras[k])); gB.setAttribute("y2", H - M.b);
    rdK.textContent = P.fmt(ks[k], 3);
    rdR.textContent = P.fmt(ras[k], 1);
  }
  slider.addEventListener("input", update);
  P.player(root.querySelector("#kp-play"), slider, update, 30);
  update();
})();
