// Interactive Figs. 5 and 6: stability coefficients sigma_g^m and sigma_b^m
// (a), and the short-time growth rate of mode m for a chosen beta*Ra (b),
//   s_m = sigma_g^m - 1/2 + beta*Ra sigma_b^m.
// Everything is in closed form; see assets/plot.js.
(function () {
  "use strict";

  var P = window.PL;
  var root = document.getElementById("disp-widget");
  if (!P || !root) return;

  var W = 420, H = 300, M = { l: 46, r: 14, t: 12, b: 30 };
  var MMAX = 30, ms = [];
  for (var m = 0; m <= MMAX; m++) ms.push(m);
  var sg = ms.map(P.sigmaG), sb = ms.map(P.sigmaB);

  var slider = root.querySelector("#dp-ra"), out = root.querySelector("#dp-ra-val");
  var rdUnst = root.querySelector("#dp-unst"), rdS3 = root.querySelector("#dp-s3"),
      rdS30 = root.querySelector("#dp-s30");

  // (a) the two coefficients
  var A = P.axes(root.querySelector("#dp-a"), { W: W, H: H, M: M, x: [0, MMAX], y: [-3, 3.5],
    xt: [0, 5, 10, 15, 20, 25, 30], yt: [-3, -2, -1, 0, 1, 2, 3] });
  P.el("line", { x1: A.sx(0), x2: A.sx(MMAX), y1: A.sy(0), y2: A.sy(0), "class": "rw-dash" }, A.g);
  P.el("path", { d: P.path(ms, sg, A.sx, A.sy), "class": "dp-g" }, A.g);
  P.el("path", { d: P.path(ms, sb.map(function (v) { return 96 * v; }), A.sx, A.sy), "class": "dp-b" }, A.g);
  ms.forEach(function (m) {
    P.el("circle", { cx: A.sx(m), cy: A.sy(sg[m]), r: 3, "class": "dp-gdot" }, A.g);
    P.el("circle", { cx: A.sx(m), cy: A.sy(96 * sb[m]), r: 3, "class": "dp-bdot" }, A.g);
  });

  // (b) growth rate at the chosen beta*Ra
  var B = P.axes(root.querySelector("#dp-b"), { W: W, H: H, M: M, x: [0, MMAX], y: [-3.5, 3.5],
    xt: [0, 5, 10, 15, 20, 25, 30], yt: [-3, -2, -1, 0, 1, 2, 3] });
  P.el("rect", { x: M.l, y: B.sy(3.5), width: W - M.l - M.r, height: B.sy(0) - B.sy(3.5), "class": "dp-unstable" }, B.g);
  P.el("line", { x1: B.sx(0), x2: B.sx(MMAX), y1: B.sy(0), y2: B.sy(0), "class": "rw-dash" }, B.g);
  P.el("text", { x: W - M.r - 8, y: B.sy(3.5) + 18, "text-anchor": "end", "class": "sf-tag" }, B.g).textContent = "unstable";
  P.el("text", { x: W - M.r - 8, y: B.sy(-3.5) - 8, "text-anchor": "end", "class": "sf-tag" }, B.g).textContent = "stable";
  P.el("path", { d: P.path(ms, sg.map(function (v) { return v - 0.5; }), B.sx, B.sy), "class": "rw-ref" }, B.g);
  var cur = P.el("path", { "class": "rw-cur" }, B.g);
  var dots = ms.map(function () { return P.el("circle", { r: 3.6 }, B.g); });

  function update() {
    var ra = +slider.value;
    out.textContent = ra;
    var s = ms.map(function (m) { return sg[m] - 0.5 + ra * sb[m]; });
    cur.setAttribute("d", P.path(ms, s, B.sx, B.sy));
    var unst = [];
    dots.forEach(function (d, m) {
      d.setAttribute("cx", B.sx(m)); d.setAttribute("cy", B.sy(s[m]));
      var u = m >= 2 && s[m] > 1e-12;
      d.setAttribute("class", u ? "dp-on" : "dp-off");
      if (u) unst.push(m);
    });
    if (!unst.length) rdUnst.textContent = "none";
    else if (unst[unst.length - 1] === MMAX && unst.length === MMAX - unst[0] + 1)
      rdUnst.textContent = "m ≥ " + unst[0];
    else rdUnst.textContent = unst.join(", ");
    rdS3.textContent = P.fmt(s[3], 3);
    rdS30.textContent = P.fmt(s[30], 3);
  }
  slider.addEventListener("input", update);
  P.player(root.querySelector("#dp-play"), slider, update, 40);
  update();
})();
