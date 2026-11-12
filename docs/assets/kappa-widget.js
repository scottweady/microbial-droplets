// Interactive Fig. 9: the matching constant kappa against the effective
// domain depth |H| (a), and the m -> infinity critical Rayleigh number
// Ra_{*,kappa} = 96 / (1 - 6 kappa) for beta = 1 (b), with a to-scale
// cross-section of the layer whose depth follows the slider. Data: assets/kappa-data.js.
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

  // Schematic: cross-section of the droplet (radius 1) on a layer of depth |H|,
  // drawn to scale. The surface sits at y = Y0 and one radius is S pixels.
  var sv = root.querySelector("#kp-s");
  var Y0 = 30, S = 32, XL = 24, XR = 560, CX = (XL + XR) / 2, XD = 584;
  var hZero = 1;   // depth where kappa = 0, by linear interpolation
  for (var j = 1; j < D.H.length; j++) {
    if (D.kappa[j - 1] > 0 !== D.kappa[j] > 0) {
      hZero = D.H[j - 1] + (D.H[j] - D.H[j - 1]) * D.kappa[j - 1] / (D.kappa[j - 1] - D.kappa[j]);
      break;
    }
  }
  var defs = P.el("defs", {}, sv);
  var clip = P.el("clipPath", { id: "kp-s-clip" }, defs);
  var clipR = P.el("rect", { x: XL, y: Y0, width: XR - XL }, clip);
  var mk = P.el("marker", { id: "kp-s-head", viewBox: "0 0 10 10", refX: 9, refY: 5, markerWidth: 7, markerHeight: 7, orient: "auto-start-reverse" }, defs);
  P.el("path", { d: "M0,1L10,5L0,9Z", "class": "kp-head" }, mk);
  var fluid = P.el("rect", { x: XL, y: Y0, width: XR - XL, "class": "kp-fluid" }, sv);
  var gp = P.el("g", { "clip-path": "url(#kp-s-clip)" }, sv);
  var seed = 7;
  function rnd() { seed = (seed * 16807) % 2147483647; return seed / 2147483647; }
  for (var q = 0; q < 320; q++) {
    P.el("circle", { cx: XL + 4 + rnd() * (XR - XL - 8), cy: Y0 + 6 + rnd() * (7 * S - 8), r: 3, "class": "kp-part" }, gp);
  }
  var refY = Y0 + S * hZero;
  var ref = P.el("g", {}, sv);
  P.el("line", { x1: XL, x2: XR, y1: refY, y2: refY, "class": "kp-ref" }, ref);
  P.el("text", { x: XL + 6, y: refY - 5, "class": "kp-lab" }, ref).textContent = "κ = 0";
  var axis = P.el("line", { x1: CX, x2: CX, y1: Y0 + 4, "class": "kp-axis" }, sv);
  var hatch = P.el("g", {}, sv);
  for (var xh = XL; xh < XR; xh += 10) P.el("line", { x1: xh + 8, x2: xh, y1: 0, y2: 8, "class": "kp-hatch" }, hatch);
  var wall = P.el("line", { x1: XL, x2: XR, "class": "kp-wall" }, sv);
  P.el("line", { x1: XL, x2: XR, y1: Y0, y2: Y0, "class": "kp-surf" }, sv);
  P.el("rect", { x: CX - S, y: Y0 - 5, width: 2 * S, height: 8, rx: 3, "class": "kp-drop" }, sv);
  P.el("text", { x: CX, y: Y0 - 11, "text-anchor": "middle", "class": "kp-lab" }, sv).textContent = "droplet, radius 1";
  var stag = P.el("circle", { cx: CX, r: 4, "class": "kp-stag" }, sv);
  var wLab = P.el("text", { x: CX + 9, "class": "kp-lab" }, sv);
  wLab.textContent = "w = 0";
  var dim = P.el("line", { x1: XD, x2: XD, y1: Y0 + 1, "marker-start": "url(#kp-s-head)", "marker-end": "url(#kp-s-head)", "class": "kp-dim" }, sv);
  var dimLab = P.el("text", { x: XD + 9, "class": "kp-math" }, sv);
  dimLab.textContent = "|H|";

  function drawSchematic(h) {
    var yb = Y0 + S * h;
    fluid.setAttribute("height", yb - Y0);
    clipR.setAttribute("height", yb - Y0);
    axis.setAttribute("y2", yb);
    wall.setAttribute("y1", yb); wall.setAttribute("y2", yb);
    hatch.setAttribute("transform", "translate(0," + (yb + 1) + ")");
    stag.setAttribute("cy", yb);
    wLab.setAttribute("y", yb - 7);
    dim.setAttribute("y2", yb - 1);
    dimLab.setAttribute("y", (Y0 + yb) / 2 + 5);
    ref.style.display = h > hZero + 0.1 ? "" : "none";   // only inside the fluid
  }

  function update() {
    var k = +slider.value;
    drawSchematic(hs[k]);
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
