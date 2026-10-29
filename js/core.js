/* core.js — 全部模块共用的数值方法、绘图、公式排版与界面组件 */
(function () {
'use strict';
const PI = Math.PI;
const PDE = window.PDE = { modules: {}, bench: {}, mode: 'class' };

/* ================= 格式化 ================= */
const $ = id => document.getElementById(id);
function fx(v, d = 3) {
  if (!isFinite(v)) return '—';
  const r = Math.abs(v) < 0.5 * 10 ** -d ? 0 : v;
  return r.toFixed(d);
}
const SUP = '⁰¹²³⁴⁵⁶⁷⁸⁹';
function sci(v, d = 1) {
  if (v === 0) return '0';
  if (!isFinite(v)) return '—';
  const e = Math.floor(Math.log10(Math.abs(v))), m = v / 10 ** e;
  if (e >= -2 && e <= 3) return v.toFixed(Math.max(0, d + 2 - e));
  return m.toFixed(d) + '×10' + (e < 0 ? '⁻' : '') + String(Math.abs(e)).replace(/\d/g, c => SUP[c]);
}
function stats(el, pairs) {
  el.innerHTML = pairs.map(([k, v]) => `<div class="stat"><span class="k">${k}</span><span class="v">${v}</span></div>`).join('');
}

/* ================= 数值方法 ================= */
// 复化 Simpson：u、v 长度为 NG+1，NG 为偶数
function dotS(u, v, h) {
  const n = u.length - 1;
  let s = u[0] * v[0] + u[n] * v[n];
  for (let i = 1; i < n; i++) s += (i & 1 ? 4 : 2) * u[i] * v[i];
  return s * h / 3;
}
function integ(u, h) {
  const n = u.length - 1;
  let s = u[0] + u[n];
  for (let i = 1; i < n; i++) s += (i & 1 ? 4 : 2) * u[i];
  return s * h / 3;
}
// 误差函数：|x|<2.5 用 Taylor 级数，否则用 erfc 连分式；相对精度约 1e-14
function erf(x) {
  const s = x < 0 ? -1 : 1; x = Math.abs(x);
  if (x < 2.5) {
    const x2 = x * x; let term = x, sum = x;
    for (let n = 1; n < 200; n++) {
      term *= -x2 / n; const add = term / (2 * n + 1); sum += add;
      if (Math.abs(add) < 1e-17 * Math.abs(sum)) break;
    }
    return s * 2 / Math.sqrt(PI) * sum;
  }
  let t = x;
  for (let n = 160; n >= 1; n--) t = x + (n / 2) / t;
  return s * (1 - Math.exp(-x * x) / (Math.sqrt(PI) * t));
}
// 经典四阶 Runge–Kutta，用于与解析解独立对照
function rk4(f, y0, t0, t1, n) {
  const hh = (t1 - t0) / n, ts = [t0], ys = [y0.slice()];
  let y = y0.slice(), t = t0;
  for (let i = 0; i < n; i++) {
    const k1 = f(t, y), k2 = f(t + hh / 2, y.map((v, j) => v + hh / 2 * k1[j]));
    const k3 = f(t + hh / 2, y.map((v, j) => v + hh / 2 * k2[j])), k4 = f(t + hh, y.map((v, j) => v + hh * k3[j]));
    y = y.map((v, j) => v + hh / 6 * (k1[j] + 2 * k2[j] + 2 * k3[j] + k4[j]));
    t += hh; ts.push(t); ys.push(y);
  }
  return { ts, ys };
}
function bisect(F, lo, hi, flo) {
  for (let it = 0; it < 90; it++) {
    const m = 0.5 * (lo + hi), fm = F(m);
    if (flo * fm <= 0) hi = m; else { lo = m; flo = fm; }
    if (hi - lo < 1e-15 * Math.max(1, hi)) break;
  }
  return 0.5 * (lo + hi);
}

/* ================= Sturm–Liouville 求解器 =================
   X'' + λX = 0, 0<x<L
   左端 a1 X(0) − b1 X'(0) = 0，右端 a2 X(L) + b2 X'(L) = 0
   D: (a,b)=(1,0)；N: (0,1)；R: (h,1)。h 可为负，此时可能出现负特征值。 */
function bcAB(t, hh) { return t === 'D' ? [1, 0] : t === 'N' ? [0, 1] : [hh, 1]; }
function solveSL(o) {
  const L = o.L, NG = o.NG || 2400, NM = o.NM || 60, h = L / NG;
  const [a1, b1] = bcAB(o.bcL, o.h1), [a2, b2] = bcAB(o.bcR, o.h2);
  const xs = new Float64Array(NG + 1);
  for (let i = 0; i <= NG; i++) xs[i] = L * i / NG;
  const rob = X => (o.bcL === 'R' ? o.h1 * X[0] * X[0] : 0) + (o.bcR === 'R' ? o.h2 * X[NG] * X[NG] : 0);
  const make = (kind, k, lam, fill) => {
    const X = new Float64Array(NG + 1), dX = new Float64Array(NG + 1);
    for (let i = 0; i <= NG; i++) fill(xs[i], X, dX, i);
    const n2 = dotS(X, X, h);
    return { kind, k, lam, X, dX, n2, ray: (dotS(dX, dX, h) + rob(X)) / n2 };
  };
  const modes = [];
  // λ = −μ² < 0：X = b1 μ cosh μx + a1 sinh μx
  const G = mu => (a2 * b1 + a1 * b2) * mu + (a1 * a2 + b1 * b2 * mu * mu) * Math.tanh(mu * L);
  const hmax = Math.abs(o.bcL === 'R' ? o.h1 : 0) + Math.abs(o.bcR === 'R' ? o.h2 : 0);
  if ((o.bcL === 'R' && o.h1 < 0) || (o.bcR === 'R' && o.h2 < 0)) {
    const muMax = 2 * hmax + 6 / L, n = 6000;
    let mPrev = 1e-7 / L, gPrev = G(mPrev);
    const neg = [];
    for (let j = 1; j <= n; j++) {
      const mu = muMax * (j + 0.37) / n, g = G(mu);
      if (gPrev * g < 0) neg.push(bisect(G, mPrev, mu, gPrev));
      mPrev = mu; gPrev = g;
    }
    neg.sort((p, q) => q - p).forEach(mu => {
      const s = Math.hypot(b1 * mu, a1), c = Math.cosh(mu * L);
      modes.push(make('neg', mu, -mu * mu, (x, X, dX, i) => {
        X[i] = (b1 * mu * Math.cosh(mu * x) + a1 * Math.sinh(mu * x)) / (s * c);
        dX[i] = (b1 * mu * mu * Math.sinh(mu * x) + a1 * mu * Math.cosh(mu * x)) / (s * c);
      }));
    });
  }
  // λ = 0：X = b1 + a1 x，条件 a2 b1 + a1 a2 L + a1 b2 = 0
  const z = a2 * b1 + a1 * a2 * L + a1 * b2;
  if (Math.abs(z) < 1e-10) {
    const s = Math.max(Math.abs(b1), Math.abs(b1 + a1 * L));
    modes.push(make('zero', 0, 0, (x, X, dX, i) => { X[i] = (b1 + a1 * x) / s; dX[i] = a1 / s; }));
  }
  // λ = k² > 0：X = b1 k cos kx + a1 sin kx
  const F = k => (a2 * b1 + a1 * b2) * k * Math.cos(k * L) + (a1 * a2 - b1 * b2 * k * k) * Math.sin(k * L);
  const dk = PI / (L * 40);
  let kPrev = 1e-7 / L, fPrev = F(kPrev);
  for (let j = 0; modes.length < NM && j < 200000; j++) {
    const k = (j + 0.37) * dk, f = F(k);
    if (fPrev * f < 0) {
      const kk = bisect(F, kPrev, k, fPrev), s = Math.hypot(b1 * kk, a1);
      modes.push(make('pos', kk, kk * kk, (x, X, dX, i) => {
        const c = Math.cos(kk * x), sn = Math.sin(kk * x);
        X[i] = (b1 * kk * c + a1 * sn) / s;
        dX[i] = (-b1 * kk * kk * sn + a1 * kk * c) / s;
      }));
    }
    kPrev = k; fPrev = f;
  }
  const nn = a1 === 0 && a2 === 0;
  return { L, NG, h, xs, modes, a1, b1, a2, b2, bcL: o.bcL, bcR: o.bcR, h1: o.h1, h2: o.h2,
    offset: nn ? 0 : 1, nNeg: modes.filter(m => m.kind === 'neg').length, F, G };
}
function project(sl, f) {
  return Float64Array.from(sl.modes, m => dotS(f, m.X, sl.h) / m.n2);
}
function maxAbs(a) { let m = 0; for (let i = 0; i < a.length; i++) { const v = Math.abs(a[i]); if (v > m) m = v; } return m || 1; }
const BCTEX = {
  L: { D: String.raw`X(0)=0`, N: String.raw`X'(0)=0`, R: String.raw`X'(0)-h_1X(0)=0` },
  R: { D: String.raw`X(L)=0`, N: String.raw`X'(L)=0`, R: String.raw`X'(L)+h_2X(L)=0` },
  uL: { D: String.raw`u(0,t)=0`, N: String.raw`u_x(0,t)=0`, R: String.raw`u_x(0,t)-h_1u(0,t)=0` },
  uR: { D: String.raw`u(L,t)=0`, N: String.raw`u_x(L,t)=0`, R: String.raw`u_x(L,t)+h_2u(L,t)=0` }
};

/* ================= 颜色与画布 ================= */
const C = {};
function readColors() {
  const cs = getComputedStyle(document.documentElement);
  ['ink', 'ink-2', 'ink-3', 'line', 'line-2', 'grid', 'accent', 'warm', 'surface', 'good', 'bad', 'violet', 'accent-soft'].forEach(k => C[k] = cs.getPropertyValue('--' + k).trim());
  C.mono = cs.getPropertyValue('--font-mono').trim();
}
function prep(cv) {
  const w = cv.clientWidth, h = cv.clientHeight;
  if (!w || !h) return null;
  const dpr = window.devicePixelRatio || 1;
  if (cv.width !== Math.round(w * dpr) || cv.height !== Math.round(h * dpr)) { cv.width = Math.round(w * dpr); cv.height = Math.round(h * dpr); }
  const ctx = cv.getContext('2d');
  ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
  ctx.clearRect(0, 0, w, h);
  return { ctx, w, h };
}
function niceTicks(min, max, count) {
  const span = max - min, raw = span / Math.max(1, count), mag = 10 ** Math.floor(Math.log10(raw)), r = raw / mag;
  const step = (r < 1.5 ? 1 : r < 3 ? 2 : r < 7 ? 5 : 10) * mag;
  const out = [];
  for (let v = Math.ceil(min / step - 1e-9) * step; v <= max + 1e-9 * span; v += step) out.push(Math.abs(v) < step * 1e-9 ? 0 : v);
  return out;
}
const tickLabel = v => String(+v.toPrecision(4));
// 建立坐标系：o = {xmin,xmax,ymin,ymax,xlabel,ylabel,xticks,yticks,equal,colorbar,padL}
function frame(cv, o) {
  const P = prep(cv); if (!P) return null;
  const { ctx, w, h } = P;
  const pad = { l: o.padL || 46, r: o.colorbar ? 58 : 14, t: 12, b: o.noX ? 12 : 30 };
  let pw = w - pad.l - pad.r, ph = h - pad.t - pad.b;
  if (o.equal) { // 等比例坐标
    const sx = pw / (o.xmax - o.xmin), sy = ph / (o.ymax - o.ymin), s = Math.min(sx, sy);
    const nw = s * (o.xmax - o.xmin), nh = s * (o.ymax - o.ymin);
    pad.l += (pw - nw) / 2; pad.t += (ph - nh) / 2; pw = nw; ph = nh;
  }
  const X = x => pad.l + (x - o.xmin) / (o.xmax - o.xmin) * pw;
  const Y = y => pad.t + (1 - (y - o.ymin) / (o.ymax - o.ymin)) * ph;
  ctx.font = '11px ' + (C.mono || 'monospace');
  ctx.lineWidth = 1;
  if (!o.noAxes) {
    const xt = o.xticks || niceTicks(o.xmin, o.xmax, Math.max(3, Math.floor(pw / 70))).map(v => ({ v, label: tickLabel(v) }));
    const yt = o.yticks || niceTicks(o.ymin, o.ymax, Math.max(3, Math.floor(ph / 42))).map(v => ({ v, label: tickLabel(v) }));
    ctx.strokeStyle = C.grid; ctx.fillStyle = C['ink-3'];
    ctx.textAlign = 'center'; ctx.textBaseline = 'top';
    xt.forEach(t => { const px = X(t.v); if (!o.noGrid) { ctx.beginPath(); ctx.moveTo(px, pad.t); ctx.lineTo(px, pad.t + ph); ctx.stroke(); } if (!o.noX) ctx.fillText(t.label, px, pad.t + ph + 6); });
    ctx.textAlign = 'right'; ctx.textBaseline = 'middle';
    yt.forEach(t => { const py = Y(t.v); if (!o.noGrid) { ctx.beginPath(); ctx.moveTo(pad.l, py); ctx.lineTo(pad.l + pw, py); ctx.stroke(); } ctx.fillText(t.label, pad.l - 6, py); });
    if (o.ymin < 0 && o.ymax > 0 && !o.noZero) { ctx.strokeStyle = C['ink-3']; ctx.beginPath(); ctx.moveTo(pad.l, Y(0)); ctx.lineTo(pad.l + pw, Y(0)); ctx.stroke(); }
    ctx.strokeStyle = C['line-2']; ctx.beginPath(); ctx.moveTo(pad.l, Math.round(pad.t + ph) + 0.5); ctx.lineTo(pad.l + pw, Math.round(pad.t + ph) + 0.5); ctx.stroke();
  }
  const f = { ctx, X, Y, pad, pw, ph, w, h, o };
  if (o.xlabel) label(f, o.xlabel, 'xr');
  if (o.ylabel) label(f, o.ylabel, 'yt');
  cv._inv = { x: px => o.xmin + (px - pad.l) / pw * (o.xmax - o.xmin), y: py => o.ymin + (1 - (py - pad.t) / ph) * (o.ymax - o.ymin), f };
  return f;
}
// 静态图层缓存：key 不变时直接贴回上次的位图，只重画随时间变化的覆盖层（播放动画时避免每帧重算热力图）
function layer(cv, key, o, draw) {
  const k = key + '|' + JSON.stringify(o) + '|' + cv.clientWidth + 'x' + cv.clientHeight + '|' + (window.devicePixelRatio || 1) + '|' + C.ink + C.surface + C.accent;
  const c = cv._layer;
  if (c && c.k === k) {
    const P = prep(cv); if (!P) return null;
    P.ctx.save(); P.ctx.setTransform(1, 0, 0, 1, 0, 0); P.ctx.drawImage(c.img, 0, 0); P.ctx.restore();
    cv._inv = c.inv; return c.f;
  }
  const f = frame(cv, o); if (!f) return null;
  if (draw) draw(f);
  const img = (c && c.img) || document.createElement('canvas'); img.width = cv.width; img.height = cv.height;
  img.getContext('2d').drawImage(cv, 0, 0);
  cv._layer = { k, img, f, inv: cv._inv };
  return f;
}
function label(f, s, where, color) {
  const { ctx, pad, pw, ph } = f;
  ctx.fillStyle = color || C['ink-2']; ctx.font = '11px ' + C.mono;
  if (where === 'xr') { ctx.textAlign = 'right'; ctx.textBaseline = 'bottom'; ctx.fillText(s, pad.l + pw - 4, pad.t + ph - 4); }
  else if (where === 'yt') { ctx.textAlign = 'left'; ctx.textBaseline = 'top'; ctx.fillText(s, pad.l + 5, pad.t + 4); }
}
function textAt(f, s, x, y, color, align = 'left', base = 'bottom') {
  if (!f) return;
  f.ctx.fillStyle = color || C['ink-2']; f.ctx.font = '11px ' + C.mono; f.ctx.textAlign = align; f.ctx.textBaseline = base;
  f.ctx.fillText(s, f.X(x), f.Y(y));
}
function clip(f) { f.ctx.save(); f.ctx.beginPath(); f.ctx.rect(f.pad.l, f.pad.t, f.pw, f.ph); f.ctx.clip(); }
function line(f, x, y, color, width = 1.6, dash = null, step = 1) {
  if (!f) return;
  const { ctx, X, Y } = f;
  clip(f);
  ctx.strokeStyle = color; ctx.lineWidth = width; ctx.setLineDash(dash || []); ctx.lineJoin = 'round'; ctx.lineCap = dash ? 'butt' : 'round';
  ctx.beginPath(); let pen = false;
  const n = x.length;
  for (let i = 0; i < n; i += step) {
    const yy = typeof y === 'function' ? y(x[i]) : y[i];
    if (!isFinite(yy)) { pen = false; continue; }
    const px = X(x[i]), py = Y(yy);
    if (pen) ctx.lineTo(px, py); else { ctx.moveTo(px, py); pen = true; }
  }
  if ((n - 1) % step) { const yy = typeof y === 'function' ? y(x[n - 1]) : y[n - 1]; if (isFinite(yy) && pen) ctx.lineTo(X(x[n - 1]), Y(yy)); }
  ctx.stroke(); ctx.setLineDash([]); ctx.restore();
}
function fnLine(f, xmin, xmax, fn, color, width, dash, n = 600) {
  const xs = new Float64Array(n + 1); for (let i = 0; i <= n; i++) xs[i] = xmin + (xmax - xmin) * i / n;
  line(f, xs, fn, color, width, dash);
}
function seg2(f, x0, y0, x1, y1, color, width = 1, dash = null) {
  if (!f) return;
  clip(f); const c = f.ctx; c.strokeStyle = color; c.lineWidth = width; c.setLineDash(dash || []);
  c.beginPath(); c.moveTo(f.X(x0), f.Y(y0)); c.lineTo(f.X(x1), f.Y(y1)); c.stroke(); c.setLineDash([]); c.restore();
}
function dot(f, x, y, color, r = 3.5, ring = false) {
  if (!f || !isFinite(y)) return;
  const c = f.ctx; c.beginPath(); c.arc(f.X(x), f.Y(y), r, 0, 2 * PI);
  if (ring) { c.strokeStyle = color; c.lineWidth = 1.6; c.stroke(); } else { c.fillStyle = color; c.fill(); }
}
function fillBetween(f, x, y0, y1, color, alpha = 0.15) {
  if (!f) return;
  clip(f); const c = f.ctx; c.globalAlpha = alpha; c.fillStyle = color; c.beginPath();
  for (let i = 0; i < x.length; i++) { const v = typeof y1 === 'function' ? y1(x[i]) : y1[i]; i ? c.lineTo(f.X(x[i]), f.Y(v)) : c.moveTo(f.X(x[i]), f.Y(v)); }
  for (let i = x.length - 1; i >= 0; i--) { const v = typeof y0 === 'function' ? y0(x[i]) : (y0.length ? y0[i] : y0); c.lineTo(f.X(x[i]), f.Y(v)); }
  c.closePath(); c.fill(); c.globalAlpha = 1; c.restore();
}
function bars(f, xs, ys, color, width, base = 0, outline = false) {
  if (!f) return;
  clip(f); const c = f.ctx;
  for (let i = 0; i < xs.length; i++) {
    const px = f.X(xs[i]) - width / 2, y0 = f.Y(base), y1 = f.Y(ys[i]);
    if (outline) { c.strokeStyle = color; c.lineWidth = 1; c.strokeRect(px + 0.5, Math.min(y0, y1), width - 1, Math.abs(y1 - y0)); }
    else { c.fillStyle = typeof color === 'function' ? color(i) : color; c.fillRect(px, Math.min(y0, y1), width, Math.abs(y1 - y0)); }
  }
  c.restore();
}
// 端点符号：D 夹持块，N 竖直滑杆 + 圆环，R 弹簧
function bcGlyphs(f, bcL, bcR) {
  if (!f) return;
  const { ctx, pad, pw, ph } = f;
  [[bcL, pad.l, -1], [bcR, pad.l + pw, 1]].forEach(([t, px, dir]) => {
    if (!t) return;
    const y0 = f.Y(0), cy = y0 >= pad.t && y0 <= pad.t + ph ? y0 : pad.t + ph / 2;
    ctx.strokeStyle = C['ink-2']; ctx.fillStyle = C['ink-2']; ctx.lineWidth = 1.5;
    if (t === 'D') ctx.fillRect(px + (dir < 0 ? 1 : -5), cy - 9, 4, 18);
    else if (t === 'N') { ctx.beginPath(); ctx.moveTo(px, pad.t + 4); ctx.lineTo(px, pad.t + ph - 4); ctx.setLineDash([2, 3]); ctx.stroke(); ctx.setLineDash([]); ctx.beginPath(); ctx.arc(px, cy, 4, 0, 2 * PI); ctx.stroke(); }
    else { ctx.beginPath(); const x1 = px - dir * 18; ctx.moveTo(px, cy); for (let i = 1; i <= 6; i++) ctx.lineTo(px + (x1 - px) * i / 6, cy + (i % 2 ? 4 : -4)); ctx.stroke(); }
  });
}
// 冷暖发散色图（在两种主题下都可读）
const CMAP = [[38, 84, 160], [110, 160, 210], [236, 238, 240], [240, 160, 110], [178, 44, 40]];
function cmap(t) {
  t = Math.max(0, Math.min(1, t)) * (CMAP.length - 1);
  const i = Math.min(CMAP.length - 2, Math.floor(t)), r = t - i, a = CMAP[i], b = CMAP[i + 1];
  return [a[0] + (b[0] - a[0]) * r, a[1] + (b[1] - a[1]) * r, a[2] + (b[2] - a[2]) * r];
}
// 热力图：data[j*nx+i]，j=0 为上边（y=ymax）；NaN 表示透明；bands>0 时画等值线
function heatmap(f, data, nx, ny, vmin, vmax, bands = 0) {
  if (!f) return;
  const off = document.createElement('canvas'); off.width = nx; off.height = ny;
  const oc = off.getContext('2d'), img = oc.createImageData(nx, ny), d = img.data;
  const span = vmax - vmin || 1, band = v => Math.floor((v - vmin) / span * bands);
  for (let j = 0; j < ny; j++) for (let i = 0; i < nx; i++) {
    const k = j * nx + i, v = data[k], p = 4 * k;
    if (!isFinite(v)) { d[p + 3] = 0; continue; }
    let [r, g, b] = cmap((v - vmin) / span);
    if (bands) {
      const bb = band(v);
      if ((i + 1 < nx && isFinite(data[k + 1]) && band(data[k + 1]) !== bb) || (j + 1 < ny && isFinite(data[k + nx]) && band(data[k + nx]) !== bb)) { r *= 0.55; g *= 0.55; b *= 0.55; }
    }
    d[p] = r; d[p + 1] = g; d[p + 2] = b; d[p + 3] = 255;
  }
  oc.putImageData(img, 0, 0);
  f.ctx.imageSmoothingEnabled = true;
  f.ctx.drawImage(off, f.pad.l, f.pad.t, f.pw, f.ph);
  if (f.o.colorbar) colorbar(f, vmin, vmax);
}
function colorbar(f, vmin, vmax) {
  const c = f.ctx, x = f.pad.l + f.pw + 12, w = 10, y0 = f.pad.t, hgt = f.ph;
  for (let i = 0; i < hgt; i++) { const [r, g, b] = cmap(1 - i / hgt); c.fillStyle = `rgb(${r | 0},${g | 0},${b | 0})`; c.fillRect(x, y0 + i, w, 1.5); }
  c.strokeStyle = C.line; c.strokeRect(x + 0.5, y0 + 0.5, w, hgt);
  c.fillStyle = C['ink-3']; c.font = '10px ' + C.mono; c.textAlign = 'left'; c.textBaseline = 'middle';
  [vmax, (vmin + vmax) / 2, vmin].forEach((v, k) => c.fillText(fx(v, 2), x + w + 4, y0 + k * hgt / 2));
}
// 3D 曲面（正交投影 + 画家算法）。pts[i][j] = [x,y,z,v]
function surface(cv, pts, view, vmin, vmax) {
  const P = prep(cv); if (!P) return;
  const { ctx, w, h } = P, ca = Math.cos(view.az), sa = Math.sin(view.az), ce = Math.cos(view.el), se = Math.sin(view.el);
  const sc = Math.min(w, h) * 0.38, cx = w / 2, cy = h / 2 + 10;
  const proj = ([x, y, z]) => { const x1 = x * ca - y * sa, y1 = x * sa + y * ca; return [cx + sc * x1, cy - sc * (y1 * se + z * ce), -y1 * ce + z * se]; };
  const quads = [];
  for (let i = 0; i + 1 < pts.length; i++) for (let j = 0; j + 1 < pts[i].length; j++) {
    const q = [pts[i][j], pts[i + 1][j], pts[i + 1][j + 1], pts[i][j + 1]], pp = q.map(proj);
    quads.push({ pp, d: (pp[0][2] + pp[1][2] + pp[2][2] + pp[3][2]) / 4, v: (q[0][3] + q[1][3] + q[2][3] + q[3][3]) / 4 });
  }
  quads.sort((a, b) => a.d - b.d);
  ctx.lineWidth = 0.5; ctx.strokeStyle = 'rgba(0,0,0,0.18)';
  quads.forEach(q => {
    const [r, g, b] = cmap((q.v - vmin) / (vmax - vmin || 1));
    ctx.fillStyle = `rgb(${r | 0},${g | 0},${b | 0})`;
    ctx.beginPath(); q.pp.forEach((p, k) => k ? ctx.lineTo(p[0], p[1]) : ctx.moveTo(p[0], p[1])); ctx.closePath(); ctx.fill(); ctx.stroke();
  });
  ctx.fillStyle = C['ink-3']; ctx.font = '11px ' + C.mono; ctx.textAlign = 'left'; ctx.textBaseline = 'top';
  ctx.fillText('拖动旋转', 8, 8);
}
function dragRotate(cv, view, redraw) {
  let last = null;
  cv.addEventListener('pointerdown', e => { cv.setPointerCapture(e.pointerId); last = [e.clientX, e.clientY]; });
  cv.addEventListener('pointermove', e => {
    if (!last || !cv.hasPointerCapture(e.pointerId)) return;
    view.az += (e.clientX - last[0]) * 0.01; view.el = Math.max(0.05, Math.min(1.5, view.el + (e.clientY - last[1]) * 0.01));
    last = [e.clientX, e.clientY]; redraw();
  });
  cv.addEventListener('pointerup', () => { last = null; });
}

/* ================= 公式排版 ================= */
// MathJax 异步加载：就绪前的排版请求先排队，就绪后统一处理
let mjChain = Promise.resolve();
const mjPending = new Set();
function typeset(el) {
  const MJ = window.MathJax;
  if (!MJ || !MJ.typesetPromise) { mjPending.add(el); return; }
  mjChain = mjChain.then(() => { if (!el.isConnected) return; MJ.typesetClear && MJ.typesetClear([el]); return MJ.typesetPromise([el]); }).catch(() => {});
}
window.__mjReady = () => { const w = document.getElementById('mj-warn'); if (w) w.hidden = true; const els = [...mjPending]; mjPending.clear(); els.forEach(typeset); };
if (window.MathJax && window.MathJax.typesetPromise) window.__mjReady();
function tex(el, s) { if (el._tex === s) return; el._tex = s; el.innerHTML = s; typeset(el); }

/* ================= 界面组件 ================= */
function segment(el, opts, value, onChange) {
  el.classList.add('seg');
  el.innerHTML = opts.map(o => `<button type="button" data-v="${o.v}">${o.label}</button>`).join('');
  const set = v => el.querySelectorAll('button').forEach(b => b.setAttribute('aria-pressed', String(b.dataset.v === String(v))));
  el.addEventListener('click', e => { const b = e.target.closest('button'); if (!b) return; set(b.dataset.v); onChange(b.dataset.v); });
  set(value);
  return { set };
}
// 滑块：返回 HTML，挂载后用 bindRange 绑定
function rangeHTML(id, labelHTML, min, max, step, value) {
  return `<div class="ctl" id="${id}-row"><label class="ctl-label" for="${id}">${labelHTML} <span class="num" id="${id}-v"></span></label><input type="range" id="${id}" min="${min}" max="${max}" step="${step}" value="${value}"></div>`;
}
// 滑块已走过的部分用强调色填充（WebKit 没有 range-progress 伪元素）
function fillRange(el) { const a = +el.min || 0, b = +el.max || 100; el.style.setProperty('--fill', ((+el.value - a) / (b - a || 1) * 100).toFixed(2) + '%'); }
function bindRange(id, onInput, fmt = v => v) {
  const el = $(id), out = $(id + '-v');
  const show = () => { out.textContent = fmt(+el.value); fillRange(el); };
  el.addEventListener('input', () => { show(); onInput(+el.value); });
  show();
  return { el, set(v) { el.value = v; show(); }, get: () => +el.value };
}
const ICON = {
  play: '<svg viewBox="0 0 12 12" aria-hidden="true"><path d="M3 1.6v8.8L10.4 6z" fill="currentColor"/></svg>',
  pause: '<svg viewBox="0 0 12 12" aria-hidden="true"><rect x="2.4" y="1.8" width="2.6" height="8.4" rx=".6" fill="currentColor"/><rect x="7" y="1.8" width="2.6" height="8.4" rx=".6" fill="currentColor"/></svg>',
  replay: '<svg viewBox="0 0 12 12" aria-hidden="true"><path d="M6 2a4 4 0 1 1-3.6 2.3" fill="none" stroke="currentColor" stroke-width="1.5" stroke-linecap="round"/><path d="M1.4 1.6l1.2 3 3-1.1z" fill="currentColor"/></svg>',
  step: '<svg viewBox="0 0 12 12" aria-hidden="true"><path d="M2.4 2v8L7.6 6z" fill="currentColor"/><rect x="8.4" y="2" width="1.6" height="8" rx=".5" fill="currentColor"/></svg>',
  reset: '<svg viewBox="0 0 12 12" aria-hidden="true"><rect x="2" y="2" width="1.6" height="8" rx=".5" fill="currentColor"/><path d="M10 2v8L4.8 6z" fill="currentColor"/></svg>'
};
// 时间播放器：滑块 + 播放/暂停/单帧/重置 + 速度
function player(el, o) {
  const id = o.id;
  el.innerHTML = `<div class="player"><div class="ctl"><label class="ctl-label" for="${id}">${o.label || '时间 \\(t\\)'} <span class="num" id="${id}-v"></span></label><input type="range" id="${id}" min="0" max="1000" step="1" value="0"></div>
  <div class="row"><button type="button" class="btn play" id="${id}-play"></button><button type="button" class="btn ghost icon" id="${id}-step" aria-label="前进一帧" title="前进一帧">${ICON.step}</button><button type="button" class="btn ghost icon" id="${id}-reset" aria-label="回到 t = 0" title="回到 t = 0">${ICON.reset}</button><div class="seg speed" id="${id}-speed" role="group" aria-label="播放速度"></div></div></div>`;
  const P = { s: 0, playing: false, raf: 0, last: 0, speed: 1 };
  segment($(id + '-speed'), [{ v: 0.5, label: '0.5×' }, { v: 1, label: '1×' }, { v: 2, label: '2×' }], 1, v => { P.speed = +v; });
  const sl = $(id), pb = $(id + '-play');
  const map = s => (o.map ? o.map(s) : s * o.tmax());
  P.t = () => map(P.s);
  const upd = () => { sl.value = Math.round(P.s * 1000); fillRange(sl); $(id + '-v').textContent = P.t().toFixed(o.digits || 3); const st = P.playing ? 'pause' : (P.s >= 1 ? 'replay' : 'play'); if (pb._st !== st) { pb._st = st; pb.innerHTML = ICON[st]; pb.setAttribute('aria-label', { play: '播放', pause: '暂停', replay: '重播' }[st]); pb.title = pb.getAttribute('aria-label'); } o.onChange(P.t()); };
  P.refresh = upd;
  P.set = s => { P.s = Math.max(0, Math.min(1, s)); upd(); };
  const tick = ts => {
    if (!P.playing) return;
    const dt = P.last ? (ts - P.last) / 1000 : 0; P.last = ts;
    P.s = Math.min(1, P.s + Math.min(dt, 0.1) * P.speed / (o.duration || 10));
    if (P.s >= 1) P.playing = false;
    upd();
    if (P.playing) P.raf = requestAnimationFrame(tick);
  };
  P.stop = () => { P.playing = false; cancelAnimationFrame(P.raf); };
  sl.addEventListener('input', () => { P.s = +sl.value / 1000; upd(); });
  pb.addEventListener('click', () => {
    if (P.playing) { P.stop(); upd(); return; }
    if (P.s >= 1) P.s = 0;
    P.playing = true; P.last = 0; P.raf = requestAnimationFrame(tick); upd();
  });
  $(id + '-step').addEventListener('click', () => { P.stop(); P.set(P.s + 1 / 200); });
  $(id + '-reset').addEventListener('click', () => { P.stop(); P.set(0); });
  return P;
}
// 理解检测
function quiz(box, items, scoreEl) {
  const answered = {};
  box.classList.add('quiz');
  box.innerHTML = items.map((it, i) => `<div class="q" data-i="${i}"><div class="q-head"><b class="mono q-no">Q${i + 1}</b><span class="pill info">${it.type || '概念'}</span></div><div>${it.q}</div><div class="q-opts">${it.o.map((o, j) => `<button type="button" data-j="${j}">${String.fromCharCode(65 + j)}. ${o}</button>`).join('')}</div><div class="q-exp" hidden></div></div>`).join('');
  typeset(box);
  const score = () => { const k = Object.keys(answered).length, r = Object.values(answered).filter(Boolean).length; scoreEl.textContent = k ? `已答 ${k}/${items.length}，正确 ${r}` : `共 ${items.length} 题，点击选项查看解析`; };
  box.addEventListener('click', e => {
    const b = e.target.closest('button'); if (!b) return;
    const qEl = b.closest('.q'), i = +qEl.dataset.i, j = +b.dataset.j, it = items[i];
    if (i in answered) return;
    answered[i] = j === it.a;
    qEl.querySelectorAll('.q-opts button').forEach((x, k) => { x.disabled = true; if (k === it.a) x.classList.add('right'); else if (k === j) x.classList.add('wrong'); });
    const ex = qEl.querySelector('.q-exp'); ex.hidden = false;
    tex(ex, (j === it.a ? '<span class="pill ok">正确</span> ' : '<span class="pill warn">再想想</span> ') + it.e);
    score();
  });
  score();
}
// 逐步推导：课堂模式逐步展开，探索模式全部展开
function stepper(box, steps) {
  box.classList.add('stepper');
  box.innerHTML = `<div class="steps">${steps.map((s, i) => `<div class="step" data-k="${i}"><h3>${s.t}</h3>${s.b}</div>`).join('')}</div>
    <div class="btn-row stepper-ctl cls-only"><button type="button" class="btn" data-a="next">下一步</button><button type="button" class="btn ghost" data-a="all">显示全部</button><button type="button" class="btn ghost" data-a="reset">从头开始</button><span class="small muted stepper-count"></span></div>`;
  let shown = 1;
  const els = box.querySelectorAll('.step'), cnt = box.querySelector('.stepper-count'), nb = box.querySelector('[data-a="next"]');
  const apply = () => {
    const all = PDE.mode === 'explore';
    els.forEach((e, i) => e.hidden = !all && i >= shown);
    cnt.textContent = `第 ${Math.min(shown, els.length)} / ${els.length} 步`;
    nb.disabled = shown >= els.length;
  };
  box.addEventListener('click', e => {
    const a = e.target.closest('[data-a]'); if (!a) return;
    if (a.dataset.a === 'next') shown++; else if (a.dataset.a === 'all') shown = els.length; else shown = 1;
    apply();
    if (a.dataset.a === 'next') els[Math.min(shown, els.length) - 1].scrollIntoView({ block: 'nearest', behavior: 'smooth' });
  });
  box._apply = apply;
  apply();
}
// 数值基准表
function benchTable(box, mod, rows) {
  PDE.bench[mod] = rows;
  const pass = rows.filter(r => r.pass).length;
  box.innerHTML = `<div class="bench-sum"><span class="pill ${pass === rows.length ? 'ok' : 'warn'}">${pass}/${rows.length} 通过</span><span class="small muted">每次打开模块或点击"重新运行"时实时计算，不是预存结果。</span><button type="button" class="btn ghost" data-rerun>重新运行</button></div>
  <div class="tbl-wrap"><table><thead><tr><th>基准案例</th><th>对照方法</th><th class="num">误差 / 指标</th><th class="num">容差</th><th>判定</th></tr></thead><tbody>${rows.map(r => `<tr><td>${r.name}</td><td class="small">${r.method}</td><td class="num">${typeof r.err === 'number' ? sci(r.err) : r.err}</td><td class="num">${r.tol === undefined ? '—' : sci(r.tol)}</td><td>${r.pass ? '<span class="pill ok">通过</span>' : '<span class="pill bad">未通过</span>'}</td></tr>`).join('')}</tbody></table></div>`;
  typeset(box);
}
function benchRow(name, method, err, tol) { return { name, method, err, tol, pass: isFinite(err) && err <= tol }; }

/* ================= 模块骨架 ================= */
const CHAIN = [['phys', '物理问题'], ['model', '数学模型'], ['cond', '定解条件'], ['idea', '求解思路'], ['deriv', '逐步推导'], ['viz', '可视化'], ['meaning', '物理意义'], ['quiz', '理解检测']];
const HERO = {
  '01': ['二阶线性方程的主部', String.raw`\begin{gathered}Au_{xx}+2Bu_{xy}+Cu_{yy}+\cdots=G\\[2pt] \Delta=B^2-AC\end{gathered}`, '\\(\\Delta>0\\) 双曲、\\(\\Delta=0\\) 抛物、\\(\\Delta<0\\) 椭圆；类型只由主部系数决定。'],
  '02': ['弦振动方程', String.raw`u_{tt}=a^2u_{xx}+f(x,t)`, '\\(a=\\sqrt{T/\\rho}\\) 是波速；无界弦的解是左右两列行波之和。'],
  '03': ['热传导方程', String.raw`u_t=D\,u_{xx}+\alpha u+f(x,t)`, '\\(D\\) 是热扩散率，\\(\\alpha u\\) 表示与温度成正比的产热或散热。'],
  '04': ['Sturm–Liouville 特征值问题', String.raw`X''+\lambda X=0,\qquad u=\sum_n T_n(t)\,X_n(x)`, '边界条件决定 \\(\\lambda_n\\) 与 \\(X_n\\)；解按特征函数展开。'],
  '05': ['极坐标下的 Laplace 方程', String.raw`u_{rr}+\frac1r u_r+\frac1{r^2}u_{\theta\theta}=0`, '圆盘内的调和函数完全由边界值决定，可以写成 Poisson 积分。'],
  '06': ['二阶常系数线性 ODE', String.raw`y''+p\,y'+q\,y=f(x)`, '分离变量后得到的每个方程都属于这一类。']
};
function modHead(id, num, title, lead) {
  const h = HERO[num];
  return `<div class="mod-head"><div class="eyebrow">模块 ${num} <span class="of">/ 06</span></div><h1>${title}</h1><p class="lead">${lead}</p>` +
    (h ? `<div class="hero-eq"><div class="label">${h[0]}</div><div class="eq">\\[ ${h[1]} \\]</div><div class="cap">${h[2]}</div></div>` : '') + `<p class="mode-hint"></p></div>
  <nav class="chain" aria-label="学习环节">${CHAIN.map(([k, l], i) => `<button type="button" data-go="${id}-${k}" data-stage="${i}"><span class="n">${String(i + 1).padStart(2, '0')}</span><span class="l">${l}</span></button>`).join('')}</nav>`;
}
function sec(id, key, title, inner, cls = '') {
  const i = CHAIN.findIndex(c => c[0] === key), tag = i >= 0 ? CHAIN[i][1] : key;
  return `<section class="sec ${cls}" id="${id}-${key}" data-stage="${i}"><header class="sec-head"><div class="k">${String(i + 1).padStart(2, '0')}<span>${tag}</span></div><h2>${title}</h2></header><div class="sec-body">${inner}</div></section>`;
}
function lab(tag, title, controls, instruments) {
  return `<div class="lab"><div class="lab-title"><span class="tag">${tag}</span><h3>${title}</h3></div><div class="lab-grid"><div class="instruments">${instruments}</div><div class="controls">${controls}</div></div></div>`;
}
function plotBlock(id, cap, legend = '', cls = '') {
  return `<div class="plot-block"><div class="plot-cap"><span>${cap}</span><span class="legend">${legend}</span></div><canvas id="${id}" class="${cls}"></canvas></div>`;
}
const lg = (color, text, dash) => `<span><i style="background:var(--${color})${dash ? ';opacity:.6' : ''}"></i>${text}</span>`;
function pitfalls(items) {
  return `<div class="pitfalls"><h3>易错点</h3><ul>${items.map(([a, b]) => `<li><div><b>${a}</b>${b}</div></li>`).join('')}</ul></div>`;
}
function register(id, def) { PDE.modules[id] = def; }

Object.assign(PDE, {
  PI, $, fx, sci, stats, dotS, integ, erf, rk4, bisect, solveSL, project, maxAbs, BCTEX, bcAB,
  C, readColors, prep, frame, layer, niceTicks, label, textAt, clip, line, fnLine, seg2, dot, fillBetween, bars, bcGlyphs,
  cmap, heatmap, colorbar, surface, dragRotate, typeset, tex,
  segment, rangeHTML, bindRange, fillRange, player, quiz, stepper, benchTable, benchRow,
  modHead, sec, lab, plotBlock, lg, pitfalls, register, CHAIN
});
})();
