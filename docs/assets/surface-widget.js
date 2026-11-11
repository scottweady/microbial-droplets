// Interactive Figs. 3 and 8: axisymmetric surface pressure (a) and radial
// surface velocity (b) for a chosen beta*Ra, by linear superposition of the
// growth (Ra = 0) and buoyancy solutions. Data: assets/surface-data.js.
(function () {
  "use strict";

  var D = window.SURFACE, P = window.PL;
  var root = document.getElementById("surface-widget");
  if (!D || !P || !root) return;

  var W = 420, H = 300, M = { l: 46, r: 14, t: 12, b: 30 };
  var slider = root.querySelector("#sf-ra"), out = root.querySelector("#sf-ra-val");
  var svgA = root.querySelector("#sf-a"), svgB = root.querySelector("#sf-b");
  var rdP0 = root.querySelector("#sf-p0"), rdEdge = root.querySelector("#sf-edge"),
      rdMax = root.querySelector("#sf-umax");

  // Pressure on the unit disk: p = lambda_1^0 w [1 - (beta Ra/96)(1 + 4/3 w^2)]
  var rIn = [];
  for (var i = 0; i <= 200; i++) rIn.push(1 - Math.pow(1 - i / 200, 1.5));
  function pressure(ra) {
    return rIn.map(function (r) {
      var w = Math.sqrt(Math.max(0, 1 - r * r));
      return P.LAM10 * w * (1 - ra / 96 * (1 + 4 / 3 * w * w));
    });
  }
  // Surface velocity: r/2 inside, u_g + beta Ra u_b outside
  function velocity(ra) {
    return D.ug.map(function (u, k) { return u + ra * D.ub[k]; });
  }

  var A = P.axes(svgA, { W: W, H: H, M: M, x: [0, 1], y: [-3.5, 1.5],
    xt: [0, 0.25, 0.5, 0.75, 1], yt: [-3, -2, -1, 0, 1] });
  var B = P.axes(svgB, { W: W, H: H, M: M, x: [0, 2], y: [-0.3, 0.7],
    xt: [0, 0.5, 1, 1.5, 2], yt: [-0.2, 0, 0.2, 0.4, 0.6] });

  // Static furniture
  P.el("line", { x1: A.sx(0), x2: A.sx(1), y1: A.sy(0), y2: A.sy(0), "class": "rw-dash" }, A.g);
  P.el("rect", { x: B.sx(0), y: M.t, width: B.sx(1) - B.sx(0), height: H - M.t - M.b, "class": "sf-drop" }, B.g);
  P.el("line", { x1: B.sx(0), x2: B.sx(2), y1: B.sy(0), y2: B.sy(0), "class": "rw-dash" }, B.g);
  P.el("text", { x: B.sx(0.5), y: M.t + 18, "text-anchor": "middle", "class": "sf-tag" }, B.g).textContent = "droplet";

  [0, 96].forEach(function (ra) {
    P.el("path", { d: P.path(rIn, pressure(ra), A.sx, A.sy), "class": "rw-ref" + (ra ? " sf-ref96" : "") }, A.g);
    P.el("path", { d: P.path(D.r, velocity(ra), B.sx, B.sy), "class": "rw-ref" + (ra ? " sf-ref96" : "") }, B.g);
  });
  P.el("line", { x1: B.sx(0), x2: B.sx(1), y1: B.sy(0), y2: B.sy(0.5), "class": "rw-line" }, B.g);

  var curA = P.el("path", { "class": "rw-cur" }, A.g);
  var dotA = P.el("circle", { r: 5, "class": "rw-dot" }, A.g);
  var curB = P.el("path", { "class": "rw-cur" }, B.g);
  var guideB = P.el("line", { "class": "rw-guide" }, B.g);
  var dotB = P.el("circle", { r: 5, "class": "rw-dot" }, B.g);

  function update() {
    var ra = +slider.value;
    out.textContent = ra;
    var p = pressure(ra);
    curA.setAttribute("d", P.path(rIn, p, A.sx, A.sy));
    dotA.setAttribute("cx", A.sx(0)); dotA.setAttribute("cy", A.sy(p[0]));

    var u = velocity(ra);
    curB.setAttribute("d", P.path(D.r, u, B.sx, B.sy));
    var kmax = -1;
    for (var k = 1; k < u.length - 1; k++) {
      if (u[k] > u[k - 1] && u[k] >= u[k + 1]) { kmax = k; break; }
    }
    if (kmax > 0) {
      dotB.style.display = guideB.style.display = "";
      dotB.setAttribute("cx", B.sx(D.r[kmax])); dotB.setAttribute("cy", B.sy(u[kmax]));
      guideB.setAttribute("x1", B.sx(D.r[kmax])); guideB.setAttribute("x2", B.sx(D.r[kmax]));
      guideB.setAttribute("y1", B.sy(u[kmax])); guideB.setAttribute("y2", H - M.b);
      rdMax.textContent = "r = " + D.r[kmax].toFixed(2) + ", u = " + P.fmt(u[kmax], 3);
    } else {
      dotB.style.display = guideB.style.display = "none";
      rdMax.textContent = "none (u decreases from 1/2)";
    }
    rdP0.textContent = P.fmt(p[0], 3);
    var c = P.LAM10 * (1 - ra / 96);
    rdEdge.textContent = P.fmt(c, 3) + (Math.abs(ra - 96) < 0.5 ? " (transition)" : c > 0 ? " (p > 0 at the edge)" : " (p < 0 at the edge)");
  }

  slider.addEventListener("input", update);
  P.player(root.querySelector("#sf-play"), slider, update, 40);
  update();
})();
