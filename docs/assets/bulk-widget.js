// Interactive Figs. 2 and 4: cross-section of the bulk flow beneath the
// droplet for a chosen beta*Ra, coloured by speed or by nutrient
// concentration, with streamlines. Data: assets/bulk-data.js.
(function () {
  "use strict";

  var D = window.BULK, P = window.PL;
  var root = document.getElementById("bulk-widget");
  if (!D || !P || !root) return;

  var nx = D.x.length, nz = D.z.length;
  var x0 = D.x[0], x1 = D.x[nx - 1], z0 = D.z[0], z1 = D.z[nz - 1];
  var dx = (x1 - x0) / (nx - 1), dz = (z1 - z0) / (nz - 1);
  var jc = (nx - 1) / 2;                              // x = 0 column

  // beta*Ra values: 0, 10^(0 : 0.1 : 4), then infinity (buoyancy only)
  var RAS = [0];
  for (var q = 0; q <= 40; q++) RAS.push(Math.pow(10, q / 10));
  RAS.push(Infinity);

  var slider = root.querySelector("#bk-ra"), out = root.querySelector("#bk-ra-val");
  slider.max = RAS.length - 1;
  var canvas = root.querySelector(".bk-canvas"), ctx = canvas.getContext("2d");
  var bar = root.querySelector(".bk-bar"), barLabel = root.querySelector("#bk-bar-label");
  var modeBtns = root.querySelectorAll("[data-bk-mode]");
  var rdAxis = root.querySelector("#bk-axis"), rdStag = root.querySelector("#bk-stag");
  var svg = root.querySelector("#bk-axes");
  var mode = "speed";

  // Axes overlay (viewBox 450 x 240, plot area 400 x 200 at (40, 10))
  var AX = P.axes(svg, { W: 450, H: 240, M: { l: 40, r: 10, t: 10, b: 30 },
    x: [-2, 2], y: [-2, 0], xt: [-2, -1, 0, 1, 2], yt: [-2, -1, 0] });
  AX.g.parentNode.removeChild(AX.g);

  var U = new Float64Array(nx * nz), Wf = new Float64Array(nx * nz);
  var speedMax = 1;

  function field(ra) {
    var inf = !isFinite(ra);
    speedMax = 0;
    for (var k = 0; k < nx * nz; k++) {
      U[k] = inf ? D.ub[k] : D.ug[k] + ra * D.ub[k];
      Wf[k] = inf ? D.wb[k] : D.wg[k] + ra * D.wb[k];
      var s = Math.hypot(U[k], Wf[k]);
      if (s > speedMax) speedMax = s;
    }
  }

  function sample(A, x, z) {
    var fx = Math.min(Math.max((x - x0) / dx, 0), nx - 1.0001);
    var fz = Math.min(Math.max((z - z0) / dz, 0), nz - 1.0001);
    var j = Math.floor(fx), i = Math.floor(fz), a = fx - j, b = fz - i;
    var k = i * nx + j;
    return (1 - b) * ((1 - a) * A[k] + a * A[k + 1]) + b * ((1 - a) * A[k + nx] + a * A[k + nx + 1]);
  }

  function cmap(C, t) {
    t = Math.min(Math.max(t, 0), 1) * (C.length - 1);
    var i = Math.min(Math.floor(t), C.length - 2), a = t - i;
    return [0, 1, 2].map(function (c) { return 255 * ((1 - a) * C[i][c] + a * C[i + 1][c]); });
  }

  function drawBar() {
    var C = mode === "speed" ? D.ice : D.matter, stops = [];
    for (var k = 0; k <= 16; k++) {
      var c = cmap(C, 1 - k / 16);
      stops.push("rgb(" + c.map(Math.round).join(",") + ") " + (100 * k / 16).toFixed(2) + "%");
    }
    bar.style.background = "linear-gradient(to bottom, " + stops.join(", ") + ")";
    barLabel.innerHTML = mode === "speed" ? "|<i>u</i>|/|<i>u</i>|<sub>max</sub>" : "<i>c</i>";
  }

  // ── Evenly spaced streamlines (after Jobard & Lefer) ───────────────────────
  var DSEP = 0.13, DTEST = 0.06, H = 0.01, CELL = DTEST;
  var GX = Math.ceil(4 / CELL), GZ = Math.ceil(2 / CELL);

  function streamlines() {
    var grid = new Array(GX * GZ);
    var lines = [];
    function cellOf(x, z) {
      var cx = Math.floor((x + 2) / CELL), cz = Math.floor(-z / CELL);
      if (cx < 0 || cz < 0 || cx >= GX || cz >= GZ) return -1;
      return cz * GX + cx;
    }
    function near(x, z, d, id, step) {
      var cx = Math.floor((x + 2) / CELL), cz = Math.floor(-z / CELL), r = Math.ceil(d / CELL);
      for (var a = cx - r; a <= cx + r; a++) for (var b = cz - r; b <= cz + r; b++) {
        if (a < 0 || b < 0 || a >= GX || b >= GZ) continue;
        var L = grid[b * GX + a]; if (!L) continue;
        for (var n = 0; n < L.length; n++) {
          var p = L[n];
          if (p[2] === id && Math.abs(p[3] - step) < 3 * d / H) continue;
          if ((p[0] - x) * (p[0] - x) + (p[1] - z) * (p[1] - z) < d * d) return true;
        }
      }
      return false;
    }
    function add1(x, z, id, step) {
      var c = cellOf(x, z); if (c < 0) return;
      (grid[c] = grid[c] || []).push([x, z, id, step]);
    }
    function dir(x, z) {
      var u = sample(U, x, z), w = sample(Wf, x, z), s = Math.hypot(u, w);
      if (s < 1e-4 * speedMax) return null;
      return [u / s, w / s];
    }
    function trace(x, z, sgn, id) {
      var pts = [];
      for (var n = 0; n < 1200; n++) {
        var d1 = dir(x, z); if (!d1) break;
        var xm = x + 0.5 * sgn * H * d1[0], zm = z + 0.5 * sgn * H * d1[1];
        var d2 = dir(xm, zm); if (!d2) break;
        x += sgn * H * d2[0]; z += sgn * H * d2[1];
        if (x < 0 || x > 2 || z < -2 || z > 0) break;
        if (near(x, z, DTEST, id, sgn * (n + 1))) break;
        pts.push([x, z, sgn * (n + 1)]);
        add1(x, z, id, sgn * (n + 1));               // also stops closed orbits
      }
      return pts;
    }
    // Seeds in the right half; lines are mirrored, since u is odd in x.
    var seeds = [[0, -1.98], [0, -0.03]];
    for (var a = 0; a < 16; a++) for (var b = 0; b < 16; b++) {
      seeds.push([0.06 + a * 0.125 + (b % 2) * 0.06, -0.06 - b * 0.124]);
    }
    for (var s = 0; s < seeds.length; s++) {
      var sx = seeds[s][0], sz = seeds[s][1];
      if (near(sx, sz, DSEP, -1, 0)) continue;
      var id = s;
      add1(sx, sz, id, 0);
      var bwd = trace(sx, sz, -1, id).reverse(), fwd = trace(sx, sz, 1, id);
      var pts = bwd.concat([[sx, sz, 0]], fwd);
      if (pts.length >= 25) lines.push(pts);
    }
    return lines;
  }

  // ── Rendering ─────────────────────────────────────────────────────────────
  var CW = 800, CH = 400;
  canvas.width = CW; canvas.height = CH;
  var img = ctx.createImageData(CW, CH);
  function px(x) { return (x + 2) / 4 * CW; }
  function py(z) { return -z / 2 * CH; }

  function render() {
    var C = mode === "speed" ? D.ice : D.matter;
    var data = img.data;
    for (var r = 0; r < CH; r++) {
      var z = -(r + 0.5) / CH * 2;
      for (var c = 0; c < CW; c++) {
        var x = (c + 0.5) / CW * 4 - 2;
        var t = mode === "speed"
          ? Math.hypot(sample(U, x, z), sample(Wf, x, z)) / speedMax
          : sample(D.c, x, z);
        var col = cmap(C, t), k = 4 * (r * CW + c);
        data[k] = col[0]; data[k + 1] = col[1]; data[k + 2] = col[2]; data[k + 3] = 255;
      }
    }
    ctx.putImageData(img, 0, 0);

    var lines = streamlines();
    ctx.lineJoin = "round"; ctx.lineCap = "round";
    [[3.2, "rgba(0,0,0,0.35)"], [1.4, "rgba(255,255,255,0.92)"]].forEach(function (st) {
      ctx.lineWidth = st[0]; ctx.strokeStyle = st[1]; ctx.fillStyle = st[1];
      lines.forEach(function (L) {
        [1, -1].forEach(function (m) {
          ctx.beginPath();
          L.forEach(function (p, n) { ctx[n ? "lineTo" : "moveTo"](px(m * p[0]), py(p[1])); });
          ctx.stroke();
          // arrows every ~0.9 units of arc length
          for (var n = 45; n < L.length - 5; n += 90) {
            var a = L[n], b = L[n + 3];
            var ang = Math.atan2(py(b[1]) - py(a[1]), px(m * b[0]) - px(m * a[0]));
            ctx.save();
            ctx.translate(px(m * a[0]), py(a[1])); ctx.rotate(ang);
            ctx.beginPath(); ctx.moveTo(7, 0); ctx.lineTo(-4, 4.5); ctx.lineTo(-4, -4.5); ctx.closePath();
            ctx.fill();
            ctx.restore();
          }
        });
      });
    });

    // droplet
    ctx.fillStyle = "#cfe8c8"; ctx.strokeStyle = "#2f6b2a"; ctx.lineWidth = 1.5;
    ctx.fillRect(px(-1), 0, px(1) - px(-1), 7);
    ctx.strokeRect(px(-1), 0.75, px(1) - px(-1), 6.25);
  }

  function readouts(ra) {
    // vertical velocity on the axis x = 0, from the bottom up
    var w = [];
    for (var i = 0; i < nz; i++) w.push(Wf[i * nx + jc]);
    var stag = null;
    for (var i2 = nz - 2; i2 > 0; i2--) {
      if (w[i2] * w[i2 - 1] < 0) {
        stag = D.z[i2 - 1] + dz * w[i2 - 1] / (w[i2 - 1] - w[i2]);
        break;
      }
    }
    var wTop = w[nz - 8];
    rdAxis.textContent = wTop > 0 ? "rising" : "sinking";
    rdStag.textContent = stag === null ? "none" : "z = " + P.fmt(stag, 2);
  }

  function update() {
    var ra = RAS[+slider.value];
    out.textContent = !isFinite(ra) ? "∞" : ra < 10 ? P.fmt(ra, 1) : Math.round(ra);
    field(ra);
    render();
    readouts(ra);
  }

  modeBtns.forEach(function (b) {
    b.addEventListener("click", function () {
      mode = b.getAttribute("data-bk-mode");
      modeBtns.forEach(function (o) { o.setAttribute("aria-pressed", o === b ? "true" : "false"); });
      drawBar();
      render();
    });
  });

  var pending = false;
  slider.addEventListener("input", function () {
    if (pending) return;
    pending = true;
    requestAnimationFrame(function () { pending = false; update(); });
  });
  P.player(root.querySelector("#bk-play"), slider, update, 6);
  slider.value = 21;                                  // beta*Ra = 100
  drawBar();
  update();
})();
