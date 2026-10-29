/* m5-laplace.js — 模块 05：Laplace 方程、圆盘 Dirichlet 问题、Poisson 积分、最大值原理与 Neumann 问题 */
(function () {
'use strict';
const P = window.PDE, PI = Math.PI, TAU = 2 * PI, $ = P.$, fx = P.fx, sci = P.sci;
const ID = 'lp', NC = 400, MQ = 1024;
const S = { R: 1, f: 'step', N: 40, pr: 0.45, pt: 0.9, rho: 0.25, th0: 0.5, neu: 'cos', C: 0 };
const view = { az: -0.6, el: 0.75 };
const drawPts = new Float64Array(121);
const FNAME = { cos2: 'cos²θ', cos: 'cos θ', mix: 'cos 3θ + ½ sin 2θ', step: '上半圆 1、下半圆 0', tri: '三角波 |θ−π|/π', spot: '局部热点', draw: '手绘' };
let co = null; // {a: Float64Array, b, exact}

/* ---------- 边界数据与 Fourier 系数 ---------- */
function fB(th) {
  th = ((th % TAU) + TAU) % TAU;
  switch (S.f) {
    case 'cos2': return Math.cos(th) ** 2;
    case 'cos': return Math.cos(th);
    case 'mix': return Math.cos(3 * th) + 0.5 * Math.sin(2 * th);
    case 'step': return th === 0 || th === PI ? 0.5 : th < PI ? 1 : 0;
    case 'tri': return Math.abs(th - PI) / PI;
    case 'spot': { let d = th - PI / 3; d = Math.atan2(Math.sin(d), Math.cos(d)); return Math.exp(-d * d / (2 * 0.15 * 0.15)); }
    case 'draw': { const t = th / TAU * 120, i = Math.min(119, Math.floor(t)), r = t - i; return drawPts[i] * (1 - r) + drawPts[i + 1] * r; }
  }
  return 0;
}
function analyticCoef(type, n) { // 返回 [a_n, b_n]；n=0 时返回 [a_0, 0]
  switch (type) {
    case 'cos2': return n === 0 ? [1, 0] : n === 2 ? [0.5, 0] : [0, 0];
    case 'cos': return n === 1 ? [1, 0] : [0, 0];
    case 'mix': return n === 3 ? [1, 0] : n === 2 ? [0, 0.5] : [0, 0];
    case 'step': return n === 0 ? [1, 0] : [0, n % 2 ? 2 / (n * PI) : 0];
    case 'tri': return n === 0 ? [1, 0] : [n % 2 ? 4 / (PI * PI * n * n) : 0, 0];
  }
  return null;
}
function trapCoef(fn, n) {
  const a = new Float64Array(n + 1), b = new Float64Array(n + 1), fv = new Float64Array(MQ);
  for (let k = 0; k < MQ; k++) fv[k] = fn(TAU * k / MQ);
  for (let m = 0; m <= n; m++) {
    let sa = 0, sb = 0;
    for (let k = 0; k < MQ; k++) { const th = TAU * k / MQ * m; sa += fv[k] * Math.cos(th); sb += fv[k] * Math.sin(th); }
    a[m] = 2 * sa / MQ; b[m] = 2 * sb / MQ;
  }
  return { a, b };
}
function coefFor(type, fn, n = NC) {
  if (analyticCoef(type, 0)) {
    const a = new Float64Array(n + 1), b = new Float64Array(n + 1);
    for (let m = 0; m <= n; m++) [a[m], b[m]] = analyticCoef(type, m);
    return { a, b, exact: true };
  }
  return { ...trapCoef(fn, n), exact: false };
}
// u(ρ,θ) 的级数，ρ=r/R；用复数递推 z^n
function uSeries(c, N, x, y) { // x,y 已除以 R
  let s = c.a[0] / 2, zr = 1, zi = 0;
  for (let n = 1; n <= N; n++) { const t = zr * x - zi * y; zi = zr * y + zi * x; zr = t; s += c.a[n] * zr + c.b[n] * zi; }
  return s;
}
function poissonAt(fn, rho, th, M = 4096) { // R=1 归一化
  let s = 0;
  for (let k = 0; k < M; k++) { const ph = TAU * k / M; s += (1 - rho * rho) / (1 - 2 * rho * Math.cos(th - ph) + rho * rho) * fn(ph); }
  return s / M;
}

/* ---------- 渲染 ---------- */
let disk = null;
function computeDisk() {
  const n = 170, d = new Float64Array(n * n); let lo = Infinity, hi = -Infinity;
  for (let j = 0; j < n; j++) for (let i = 0; i < n; i++) {
    const x = -1 + 2 * (i + 0.5) / n, y = 1 - 2 * (j + 0.5) / n;
    if (x * x + y * y > 1) { d[j * n + i] = NaN; continue; }
    const v = uSeries(co, S.N, x, y); d[j * n + i] = v; lo = Math.min(lo, v); hi = Math.max(hi, v);
  }
  let flo = Infinity, fhi = -Infinity; for (let k = 0; k < 720; k++) { const v = fB(TAU * k / 720); flo = Math.min(flo, v); fhi = Math.max(fhi, v); }
  disk = { n, d, lo: Math.min(lo, flo), hi: Math.max(hi, fhi) };
}
let MX = null;
function computeMax() { // 400 项、ρ ≤ 0.95：截断误差 < 0.95^400 Σ|c_n| ≈ 1e-9
  const n = 120; let iMax = -Infinity, iMin = Infinity, pMax = null, pMin = null, strict = 0;
  const g = new Float64Array(n * n).fill(NaN);
  for (let j = 0; j < n; j++) for (let i = 0; i < n; i++) {
    const x = -1 + 2 * (i + 0.5) / n, y = 1 - 2 * (j + 0.5) / n;
    if (x * x + y * y > 0.95 * 0.95) continue;
    const v = uSeries(co, NC, x, y); g[j * n + i] = v;
    if (v > iMax) { iMax = v; pMax = [x, y]; } if (v < iMin) { iMin = v; pMin = [x, y]; }
  }
  for (let j = 1; j < n - 1; j++) for (let i = 1; i < n - 1; i++) { // 离散严格局部极值
    const v = g[j * n + i]; if (!isFinite(v)) continue;
    const nb = [g[j * n + i - 1], g[j * n + i + 1], g[(j - 1) * n + i], g[(j + 1) * n + i], g[(j - 1) * n + i - 1], g[(j - 1) * n + i + 1], g[(j + 1) * n + i - 1], g[(j + 1) * n + i + 1]];
    if (nb.some(q => !isFinite(q))) continue;
    if (nb.every(q => v > q + 1e-12) || nb.every(q => v < q - 1e-12)) strict++;
  }
  let bMax = -Infinity, bMin = Infinity, tMax = 0, tMin = 0;
  for (let k = 0; k < 4096; k++) { const th = TAU * k / 4096, v = fB(th); if (v > bMax) { bMax = v; tMax = th; } if (v < bMin) { bMin = v; tMin = th; } }
  MX = { iMax, iMin, pMax, pMin, bMax, bMin, tMax, tMin, strict, constant: bMax - bMin < 1e-12 };
}
function hessianAt(x, y) {
  const h = 1e-3, u = (a, b) => uSeries(co, NC, a, b);
  const uxx = (u(x + h, y) - 2 * u(x, y) + u(x - h, y)) / (h * h), uyy = (u(x, y + h) - 2 * u(x, y) + u(x, y - h)) / (h * h);
  const uxy = (u(x + h, y + h) - u(x + h, y - h) - u(x - h, y + h) + u(x - h, y - h)) / (4 * h * h);
  const m = (uxx + uyy) / 2, r = Math.hypot((uxx - uyy) / 2, uxy);
  return { uxx, uyy, uxy, l1: m + r, l2: m - r, lap: uxx + uyy };
}
function meanCircle(x0, y0, rr, K = 720) { let s = 0; for (let k = 0; k < K; k++) { const t = TAU * k / K; s += uSeries(co, NC, x0 + rr * Math.cos(t), y0 + rr * Math.sin(t)); } return s / K; }
const probeXY = () => [S.pr * Math.cos(S.pt), S.pr * Math.sin(S.pt)];

function render() {
  const R = String.raw, Rr = S.R;
  P.tex($('lp-prob'), R`\[ \begin{cases}\Delta u=u_{rr}+\frac1ru_r+\frac1{r^2}u_{\theta\theta}=0,& 0\le r\lt R=${fx(Rr, 2)},\\ u(R,\theta)=f(\theta)=\text{${FNAME[S.f]}}\end{cases}\quad u=\frac{a_0}2+\sum_{n=1}^{${S.N}}\Big(\frac rR\Big)^n(a_n\cos n\theta+b_n\sin n\theta) \]`);
  $('lp-coefSrc').innerHTML = co.exact ? '<span class="pill ok">解析系数</span> 该边界函数的 Fourier 系数有闭式，直接使用。' : '<span class="pill info">数值系数</span> 用 1024 点梯形公式计算；对光滑周期函数误差随点数指数下降。';
  // 圆盘热力图
  const f = P.frame($('lp-cvDisk'), { xmin: -Rr, xmax: Rr, ymin: -Rr, ymax: Rr, equal: true, colorbar: true, xlabel: 'x', ylabel: 'y' });
  const [px, py] = probeXY();
  if (f) {
    P.heatmap(f, disk.d, disk.n, disk.n, disk.lo, disk.hi, 14);
    const c = f.ctx, cx = f.X(0), cy = f.Y(0), rad = f.X(Rr) - cx;
    // Poisson 核权重画在边界环上
    const rho = S.pr, K = 240;
    for (let k = 0; k < K; k++) {
      const ph = TAU * k / K, w = (1 - rho * rho) / (1 - 2 * rho * Math.cos(S.pt - ph) + rho * rho);
      c.strokeStyle = P.C.violet; c.lineWidth = Math.min(10, 1 + 2.2 * w); c.globalAlpha = Math.min(1, 0.25 + 0.2 * w);
      c.beginPath(); c.arc(cx, cy, rad + 6, -ph, -ph - TAU / K, true); c.stroke();
    }
    c.globalAlpha = 1;
    if (MX) {
      P.dot(f, Rr * Math.cos(MX.tMax), Rr * Math.sin(MX.tMax), P.C.ink, 6, true);
      P.dot(f, Rr * Math.cos(MX.tMin), Rr * Math.sin(MX.tMin), P.C.ink, 6, true);
      if (MX.pMax) P.dot(f, MX.pMax[0] * Rr, MX.pMax[1] * Rr, P.C.bad, 3.5);
      if (MX.pMin) P.dot(f, MX.pMin[0] * Rr, MX.pMin[1] * Rr, P.C.accent, 3.5);
    }
    const mr = Math.min(S.rho, 1 - S.pr - 1e-3);
    c.strokeStyle = P.C.ink; c.lineWidth = 1.3; c.setLineDash([3, 3]); c.beginPath(); c.arc(f.X(px * Rr), f.Y(py * Rr), mr * rad, 0, TAU); c.stroke(); c.setLineDash([]);
    P.dot(f, px * Rr, py * Rr, P.C.ink, 4.5);
    P.seg2(f, -Rr * Math.cos(S.th0), -Rr * Math.sin(S.th0), Rr * Math.cos(S.th0), Rr * Math.sin(S.th0), P.C.ink, 1, [6, 4]);
  }
  // 3D 曲面
  const pts = [], nr = 22, nt = 72, span = (disk.hi - disk.lo) || 1;
  for (let i = 0; i <= nr; i++) { const rr = i / nr, row = []; for (let j = 0; j <= nt; j++) { const t = TAU * j / nt, x = rr * Math.cos(t), y = rr * Math.sin(t), v = i === nr ? fB(t) : uSeries(co, S.N, x, y); row.push([x, y, 0.9 * (v - (disk.lo + disk.hi) / 2) / span, v]); } pts.push(row); }
  P.surface($('lp-cv3d'), pts, view, disk.lo, disk.hi);
  // 径向截面
  const g = P.frame($('lp-cvRad'), { xmin: -Rr, xmax: Rr, ymin: disk.lo - 0.1 * span, ymax: disk.hi + 0.1 * span, xlabel: 's' });
  if (g) {
    const c0 = Math.cos(S.th0), s0 = Math.sin(S.th0);
    P.fnLine(g, -Rr, Rr, s => uSeries(co, S.N, s / Rr * c0, s / Rr * s0), P.C.accent, 2.2);
    P.dot(g, Rr, fB(S.th0), P.C.warm, 5, true); P.dot(g, -Rr, fB(S.th0 + PI), P.C.warm, 5, true);
    P.textAt(g, 'f(θ₀+π)', -Rr * 0.97, fB(S.th0 + PI), P.C.warm, 'left', 'bottom');
    P.textAt(g, 'f(θ₀)', Rr * 0.97, fB(S.th0), P.C.warm, 'right', 'bottom');
  }
  // 边界数据与部分和
  const b = P.frame($('lp-cvF'), { xmin: 0, xmax: TAU, ymin: S.f === 'draw' ? -1.25 : disk.lo - 0.15 * span, ymax: S.f === 'draw' ? 1.25 : disk.hi + 0.15 * span, xlabel: 'θ', xticks: [0, 1, 2, 3, 4].map(k => ({ v: k * PI / 2, label: ['0', 'π/2', 'π', '3π/2', '2π'][k] })) });
  if (b) {
    P.fnLine(b, 0, TAU, fB, P.C.warm, 2, null, 720);
    P.fnLine(b, 0, TAU, t => uSeries(co, S.N, Math.cos(t), Math.sin(t)), P.C.accent, 1.6, null, 720);
    P.seg2(b, S.pt, b.o.ymin, S.pt, b.o.ymax, P.C.violet, 1, [3, 3]);
    if (S.f === 'draw') P.textAt(b, '按住拖动来绘制 f(θ)', 0.2, b.o.ymax, P.C['ink-3'], 'left', 'top');
  }
  // Poisson 核
  const rho = S.pr, k = P.frame($('lp-cvK'), { xmin: 0, xmax: TAU, ymin: 0, ymax: Math.min(40, (1 + rho) / (1 - rho)) * 1.1, xlabel: 'φ', xticks: [0, 1, 2, 3, 4].map(q => ({ v: q * PI / 2, label: ['0', 'π/2', 'π', '3π/2', '2π'][q] })) });
  if (k) {
    const xs = []; for (let i = 0; i <= 720; i++) xs.push(TAU * i / 720);
    const wk = ph => (1 - rho * rho) / (1 - 2 * rho * Math.cos(S.pt - ph) + rho * rho);
    P.fillBetween(k, xs, 0, wk, P.C.violet, 0.18); P.line(k, xs, wk, P.C.violet, 2);
    P.seg2(k, 0, 1, TAU, 1, P.C['ink-3'], 1, [4, 4]);
  }
  const uS = uSeries(co, NC, px, py), uP = poissonAt(fB, S.pr, S.pt);
  P.stats($('lp-pStats'), [
    ['选点 (r/R, θ)', `(${fx(S.pr, 3)}, ${fx(S.pt, 3)})`],
    ['级数值（400 项）', fx(uS, 8)],
    ['Poisson 积分（4096 点）', fx(uP, 8)],
    ['两者之差', sci(Math.abs(uS - uP))],
    ['核的峰值 (R+r)/(R−r)', fx((1 + rho) / (1 - rho), 3)],
    ['核的最小值 (R−r)/(R+r)', fx((1 - rho) / (1 + rho), 4)]
  ]);
  // 最大值原理
  if (MX) {
    const H = hessianAt(px, py), mr = Math.min(S.rho, 1 - S.pr - 1e-3), mv = meanCircle(px, py, mr);
    const okM = MX.iMax <= MX.bMax + 1e-9, okm = MX.iMin >= MX.bMin - 1e-9;
    P.stats($('lp-mStats'), [
      ['边界最大值 max f', fx(MX.bMax, 6)],
      ['内部最大值（r ≤ 0.95R）', fx(MX.iMax, 6)],
      ['最大值判定', okM ? (MX.constant ? '<span class="pill info">常数解</span>' : '<span class="pill ok">内部 &lt; 边界（严格）</span>') : '<span class="pill bad">违反</span>'],
      ['边界最小值 min f', fx(MX.bMin, 6)],
      ['内部最小值', fx(MX.iMin, 6)],
      ['最小值判定', okm ? '<span class="pill ok">内部 &gt; 边界</span>' : '<span class="pill bad">违反</span>'],
      ['网格上的严格局部极值点', MX.strict === 0 ? '<span class="pill ok">0 个</span>' : MX.strict + ' 个'],
      ['选点处 Hessian 特征值', `${fx(H.l1, 4)}, ${fx(H.l2, 4)}`],
      ['Δu = 迹（差分）', sci(H.lap)],
      [`圆周平均（半径 ${fx(mr, 2)}R）`, fx(mv, 8)],
      ['与中心值之差', sci(Math.abs(mv - uS))]
    ]);
  }
  renderNeumann();
}

/* ---------- Neumann 问题 ---------- */
const GNAME = { cos: 'cos θ', sin2: 'sin 2θ', sign: 'sgn(cos θ)', bad: 'cos θ + 0.5' };
function gN(th) { switch (S.neu) { case 'cos': return Math.cos(th); case 'sin2': return Math.sin(2 * th); case 'sign': { const c = Math.cos(th); return Math.abs(c) < 1e-14 ? 0 : Math.sign(c); } case 'bad': return Math.cos(th) + 0.5; } return 0; }
function neuCoef(n = 200) {
  const a = new Float64Array(n + 1), b = new Float64Array(n + 1);
  if (S.neu === 'cos' || S.neu === 'bad') { a[1] = 1; a[0] = S.neu === 'bad' ? 1 : 0; }
  else if (S.neu === 'sin2') b[2] = 1;
  else for (let m = 1; m <= n; m += 2) a[m] = 4 / (m * PI) * (((m - 1) / 2) % 2 ? -1 : 1);
  return { a, b };
}
function uNeu(c, x, y, N) { // R=1：u = C + Σ (1/n) ρ^n (a_n cos + b_n sin)
  let s = S.C, zr = 1, zi = 0;
  for (let n = 1; n <= N; n++) { const t = zr * x - zi * y; zi = zr * y + zi * x; zr = t; s += (c.a[n] * zr + c.b[n] * zi) / n; }
  return s;
}
function renderNeumann() {
  const c = neuCoef(), N = 120, flux = c.a[0] * PI; // ∫ g dθ = π a_0
  const R = String.raw;
  P.tex($('lp-nProb'), R`\[ \Delta u=0\ (r\lt1),\quad \frac{\partial u}{\partial r}\Big|_{r=1}=g(\theta)=${GNAME[S.neu]},\qquad \int_0^{2\pi}g\,d\theta=${fx(flux, 4)}\ ${Math.abs(flux) > 1e-12 ? R`\neq0\ \Rightarrow\ \text{无解}` : R`=0\ \Rightarrow\ u=C+\sum_{n\ge1}\frac{r^n}{n}(a_n\cos n\theta+b_n\sin n\theta)`} \]`);
  const n = 140, d = new Float64Array(n * n); let lo = Infinity, hi = -Infinity;
  for (let j = 0; j < n; j++) for (let i = 0; i < n; i++) { const x = -1 + 2 * (i + 0.5) / n, y = 1 - 2 * (j + 0.5) / n; if (x * x + y * y > 1) { d[j * n + i] = NaN; continue; } const v = uNeu(c, x, y, N); d[j * n + i] = v; lo = Math.min(lo, v); hi = Math.max(hi, v); }
  const f = P.frame($('lp-cvN'), { xmin: -1, xmax: 1, ymin: -1, ymax: 1, equal: true, colorbar: true, xlabel: 'x', ylabel: 'y' });
  if (f) {
    P.heatmap(f, d, n, n, Math.min(lo, -1 + S.C), Math.max(hi, 1 + S.C), 14);
    for (let k = 0; k < 24; k++) { const t = TAU * k / 24, gv = gN(t), x0 = Math.cos(t), y0 = Math.sin(t); P.seg2(f, x0, y0, x0 * (1 + 0.18 * gv), y0 * (1 + 0.18 * gv), gv >= 0 ? P.C.bad : P.C.accent, 2); }
  }
  let e = 0; const h = 1e-5;
  for (let k = 0; k < 64; k++) { const t = TAU * (k + 0.5) / 64; const dr = (uNeu(c, Math.cos(t), Math.sin(t), N) - uNeu(c, (1 - h) * Math.cos(t), (1 - h) * Math.sin(t), N)) / h; e = Math.max(e, Math.abs(dr - (gN(t) - c.a[0] / 2))); }
  P.stats($('lp-nStats'), [
    ['净热流 ∫g dθ', fx(flux, 6)],
    ['相容性', Math.abs(flux) < 1e-12 ? '<span class="pill ok">满足</span>' : '<span class="pill bad">不满足：稳态不存在</span>'],
    ['图中实际求解的 g', Math.abs(flux) < 1e-12 ? 'g' : 'g − 平均值（仅作对照）'],
    ['加法常数 C', fx(S.C, 2)],
    ['边界 ∂u/∂r 与 g 的偏差', S.neu === 'sign' ? '间断 g，级数在跳跃处收敛慢' : sci(e)],
    ['u(0) = C（平均值性质）', fx(uNeu(c, 0, 0, N), 4)]
  ]);
}

/* ---------- 基准 ---------- */
function runBench() {
  const rows = [], keep = S.f;
  const err = (type, exact) => { S.f = type; const c = coefFor(type, fB); let e = 0; for (let k = 0; k < 400; k++) { const r = Math.sqrt((k % 20) / 19) * 0.999, t = TAU * Math.floor(k / 20) / 20, x = r * Math.cos(t), y = r * Math.sin(t); e = Math.max(e, Math.abs(uSeries(c, NC, x, y) - exact(x, y))); } return e; };
  rows.push(P.benchRow('f = cos θ：u = (r/R) cos θ', '解析解', err('cos', x => x), 1e-14));
  rows.push(P.benchRow('f = cos²θ：u = ½ + (x² − y²)/2', '解析解', err('cos2', (x, y) => 0.5 + (x * x - y * y) / 2), 1e-14));
  S.f = 'mix'; let c = coefFor('mix', fB), e = 0;
  [[0.3, 1], [0.6, 2.5], [0.8, 4]].forEach(([r, t]) => { e = Math.max(e, Math.abs(poissonAt(fB, r, t) - uSeries(c, NC, r * Math.cos(t), r * Math.sin(t)))); });
  rows.push(P.benchRow('Poisson 积分 vs 级数（三角多项式数据）', '4096 点梯形公式', e, 1e-12));
  S.f = 'step'; c = coefFor('step', fB); e = 0;
  [[0.3, 0.4], [0.5, 2], [0.7, 4]].forEach(([r, t]) => { e = Math.max(e, Math.abs(poissonAt(fB, r, t) - uSeries(c, NC, r * Math.cos(t), r * Math.sin(t)))); });
  rows.push(P.benchRow('Poisson 积分 vs 级数（阶跃数据）', '跳跃点取平均的梯形公式，误差 O(h²)', e, 1e-5));
  const tc = trapCoef(fB, 60); e = 0; for (let m = 1; m <= 60; m++) e = Math.max(e, Math.abs(tc.b[m] - analyticCoef('step', m)[1]), Math.abs(tc.a[m]));
  rows.push(P.benchRow('阶跃数据：梯形系数 vs 解析 b_n = 2/(nπ)', '1024 点梯形公式（n ≤ 60）', e, 2e-4));
  co = c;
  const mv = meanCircle(0.4, 0.2, 0.35), uc = uSeries(c, NC, 0.4, 0.2);
  rows.push(P.benchRow('平均值性质：圆周平均 = 圆心值', '720 点梯形公式', Math.abs(mv - uc), 1e-12));
  e = 0; [[0.2, 0.1], [-0.3, 0.4], [0.1, -0.5]].forEach(([x, y]) => { e = Math.max(e, Math.abs(hessianAt(x, y).lap)); });
  rows.push(P.benchRow('Δu ≈ 0：五点差分残差（阶跃数据，400 项）', '差分 h = 10⁻³，误差 O(h²)', e, 1e-4));
  let okAll = true;['step', 'spot', 'tri'].forEach(t => { S.f = t; co = coefFor(t, fB); computeMax(); okAll = okAll && MX.iMax <= MX.bMax + 1e-9 && MX.iMin >= MX.bMin - 1e-9 && MX.strict === 0; });
  rows.push({ name: '最大值原理：3 种边界数据，内部无超出、无严格局部极值', method: 'r ≤ 0.95R 网格，400 项级数', err: okAll ? 0 : 1, tol: 0, pass: okAll });
  const keepN = S.neu; S.neu = 'cos'; const nc = neuCoef(); e = 0;
  for (let k = 0; k < 32; k++) { const t = TAU * k / 32, h = 1e-6; e = Math.max(e, Math.abs((uNeu(nc, Math.cos(t), Math.sin(t), 10) - uNeu(nc, (1 - h) * Math.cos(t), (1 - h) * Math.sin(t), 10)) / h - Math.cos(t))); }
  rows.push(P.benchRow('Neumann g = cos θ：∂u/∂r|_{r=1} = g', '单侧差分 h = 10⁻⁶', e, 2e-6));
  S.neu = 'bad'; const bc = neuCoef();
  rows.push({ name: 'Neumann 相容性：g = cos θ + 0.5 被识别为无解', method: `∫g dθ = ${fx(bc.a[0] * PI, 4)}`, err: 0, tol: 0, pass: Math.abs(bc.a[0] * PI - PI) < 1e-14 });
  S.neu = keepN; S.f = keep; co = coefFor(S.f, fB);
  P.benchTable($('lp-bench'), 'laplace', rows);
}

/* ---------- 文字 ---------- */
const QUIZ = [
  { type: '极坐标公式', q: '极坐标下的 Laplace 算子是？', o: [String.raw`\(u_{rr}+u_{\theta\theta}\)`, String.raw`\(u_{rr}+\frac1ru_r+\frac1{r^2}u_{\theta\theta}\)`, String.raw`\(u_{rr}+\frac1{r^2}u_r+\frac1ru_{\theta\theta}\)`, String.raw`\(u_{rr}+ru_r+r^2u_{\theta\theta}\)`], a: 1, e: String.raw`由链式法则得到。常见错误是漏掉 \(\frac1ru_r\) 一项，或把 \(r\) 的幂次写反。` },
  { type: '推导辨析', q: String.raw`圆盘内部问题中，为什么舍去 \(r^{-n}\) 和 \(\ln r\)？`, o: ['它们不满足 Laplace 方程', '它们在 r=0 处无界，而解在圆心应有界', '它们不满足边界条件', '为了简化计算'], a: 1, e: String.raw`它们都满足径向 Euler 方程，但在 \(r=0\) 处趋于无穷。圆环或圆外问题中这些项必须保留。` },
  { type: '特征值', q: String.raw`角向方程 \(\Theta''+\lambda\Theta=0\) 中 \(\lambda\) 只能取 \(n^2\) 的原因是？`, o: ['Θ 必须以 2π 为周期（周期边界条件）', 'Θ 在 θ=0 处为零', '物理上要求', '径向方程要求'], a: 0, e: String.raw`\(\Theta(\theta+2\pi)=\Theta(\theta)\) 是隐含的边界条件，它给出 \(\lambda=n^2\)，每个 \(n\ge1\) 有两个特征函数 \(\cos n\theta,\sin n\theta\)。` },
  { type: '最大值原理', q: '圆盘内非常数的调和函数能在内部取到局部极大值吗？', o: ['可以，只要边界数据合适', '不能：强最大值原理', '只能在圆心取到', '只有当 f 不连续时可以'], a: 1, e: String.raw`若在内部点 \(P\) 取最大值，由平均值性质，\(u\) 在以 \(P\) 为心的小圆盘上恒等于 \(u(P)\)，再由连通性推出 \(u\) 为常数。Hessian 的迹为 0，所以不可能两个特征值都为负。` },
  { type: 'Poisson 积分', q: '选点越靠近边界，Poisson 核的形状如何变化？', o: ['越来越平', '在最近的边界点附近越来越集中', '变成负值', '与位置无关'], a: 1, e: String.raw`核的最大值 \((R+r)/(R-r)\to\infty\)，最小值 \((R-r)/(R+r)\to0\)，总积分保持为 1：核趋于 δ 函数，所以 \(u\to f\)。` },
  { type: 'Neumann 问题', q: String.raw`圆盘的 Neumann 问题 \(\partial_ru|_{r=R}=g\) 有解的必要条件是？`, o: ['g 处处为正', String.raw`\(\int_0^{2\pi}g\,d\theta=0\)`, 'g 连续', 'g 的最大值等于最小值'], a: 1, e: String.raw`散度定理：\(\int_D\Delta u=\oint\partial_nu=R\int g\,d\theta\)，左边为 0。物理上：稳态时净热流必须为零。解只确定到加一个常数。` },
  { type: '物理现象', q: '圆盘中心的温度等于什么？', o: ['边界温度的最大值', '边界温度的平均值', '0', '边界最高温与最低温的平均'], a: 1, e: String.raw`\(u(0)=\frac{a_0}2=\frac1{2\pi}\int_0^{2\pi}f\,d\theta\)：平均值性质，也是 Poisson 核在 \(r=0\) 时等于常数 \(\frac1{2\pi}\) 的结果。` }
];
const STEPS = [
  { t: '写出问题', b: String.raw`<p>圆盘 \(r\lt R\) 内的稳态温度，边界温度为 \(f(\theta)\)：</p><div class="formula-box">\[ u_{rr}+\frac1ru_r+\frac1{r^2}u_{\theta\theta}=0,\qquad u(R,\theta)=f(\theta). \]</div><p>隐含条件：\(u\) 关于 \(\theta\) 以 \(2\pi\) 为周期，在圆心有界。</p>` },
  { t: '分离变量', b: String.raw`<div class="formula-box">\[ u=\mathcal R(r)\Theta(\theta):\quad \frac{r^2\mathcal R''+r\mathcal R'}{\mathcal R}=-\frac{\Theta''}{\Theta}=\lambda . \]</div>` },
  { t: '角向特征值问题（周期条件）', b: String.raw`<div class="formula-box">\[ \Theta''+\lambda\Theta=0,\ \Theta(\theta+2\pi)=\Theta(\theta)\ \Rightarrow\ \lambda=n^2,\ \ \Theta_n=a_n\cos n\theta+b_n\sin n\theta\ (n=0,1,2,\dots) \]</div>` },
  { t: '径向 Euler 方程', b: String.raw`<div class="formula-box">\[ r^2\mathcal R''+r\mathcal R'-n^2\mathcal R=0\ \Rightarrow\ \mathcal R=\begin{cases}c_1+c_2\ln r,&n=0\\ c_1r^n+c_2r^{-n},&n\ge1\end{cases} \]</div><p>在圆心有界 ⇒ \(c_2=0\)。这是模块 06 中 Euler 方程的直接应用。</p>` },
  { t: '叠加并匹配边界', b: String.raw`<div class="formula-box">\[ u=\frac{a_0}2+\sum_{n=1}^\infty\Big(\frac rR\Big)^n(a_n\cos n\theta+b_n\sin n\theta),\quad a_n=\frac1\pi\int_0^{2\pi}f\cos n\theta\,d\theta,\ b_n=\frac1\pi\int_0^{2\pi}f\sin n\theta\,d\theta . \]</div>` },
  { t: '例：f = cos²θ', b: String.raw`<p>\(\cos^2\theta=\tfrac12+\tfrac12\cos2\theta\)，所以 \(a_0=1,\ a_2=\tfrac12\)，其余为 0：</p><div class="formula-box">\[ u=\frac12+\frac12\Big(\frac rR\Big)^2\cos2\theta=\frac12+\frac{x^2-y^2}{2R^2}. \]</div><p>可直接验证 \(\Delta u=0\)。中心温度 \(\tfrac12\) 正是边界平均温度。</p>` },
  { t: '求和得到 Poisson 积分', b: String.raw`<p>把系数公式代回级数，交换求和与积分，对几何级数求和：</p><div class="formula-box">\[ u(r,\theta)=\frac1{2\pi}\int_0^{2\pi}\frac{R^2-r^2}{R^2-2Rr\cos(\theta-\varphi)+r^2}\,f(\varphi)\,d\varphi . \]</div><p>核为正、积分为 1，因此 \(u\) 是边界值的加权平均。最大值原理和平均值性质都从这里看得出来。</p>` }
];
function html() {
  const R = String.raw;
  return P.modHead(ID, '05', 'Laplace 方程、调和函数与圆盘问题', R`\(\Delta u=0\) 描述没有源的稳态：稳态温度、静电势、不可压缩无旋流的势函数。本模块在圆盘上用极坐标分离变量求解，引出 Poisson 积分，并用可交互的图像检验平均值性质和最大值原理。`) +
  P.sec(ID, 'phys', '稳态问题', R`<div class="cols">
    <div class="card"><h3>稳态温度</h3><p class="small muted">热方程 \(u_t=D\Delta u\) 在 \(t\to\infty\) 时 \(u_t\to0\)，得 \(\Delta u=0\)。边界温度固定为 \(f(\theta)\) 的金属圆盘，内部最终温度满足 Laplace 方程。</p></div>
    <div class="card"><h3>静电势</h3><p class="small muted">无电荷区域中 \(\nabla\cdot E=0\)，\(E=-\nabla u\) ⇒ \(\Delta u=0\)。导体边界上给定电势就是 Dirichlet 条件。</p></div>
    <div class="card"><h3>几何意义</h3><p class="small muted">\(\Delta u(P)\) 正比于"\(P\) 附近小圆周上 \(u\) 的平均值 − \(u(P)\)"。\(\Delta u=0\) 说明每点的值恰好等于周围的平均值：没有局部的"凸起"或"凹陷"。</p></div></div>`) +
  P.sec(ID, 'model', '调和函数与 Laplace 算子', R`<div class="formula-box">\[ \Delta u=u_{xx}+u_{yy}=u_{rr}+\frac1ru_r+\frac1{r^2}u_{\theta\theta}=0 \]</div><p class="prose small muted">在区域内二阶连续可微且满足 \(\Delta u=0\) 的函数称为调和函数。例：\(1,\ x,\ y,\ x^2-y^2,\ xy,\ r^n\cos n\theta,\ \ln r\)（\(r\neq0\)）。调和函数的实部、虚部来自解析函数：\(\operatorname{Re}z^n=r^n\cos n\theta\)。</p>`) +
  P.sec(ID, 'cond', '边界条件', R`<div class="tbl-wrap"><table><thead><tr><th>问题</th><th>边界条件</th><th>解的存在唯一性</th></tr></thead><tbody>
    <tr><td>Dirichlet</td><td>\(u|_{r=R}=f\)</td><td>\(f\) 连续时存在唯一</td></tr>
    <tr><td>Neumann</td><td>\(\partial_nu|_{r=R}=g\)</td><td>需 \(\oint g=0\)；解差一个常数</td></tr>
    <tr><td>Robin</td><td>\(\partial_nu+hu=g\)，\(h\gt0\)</td><td>存在唯一</td></tr></tbody></table></div>
    <p class="prose small muted">圆盘问题还有两个"隐含"条件：关于 \(\theta\) 的 \(2\pi\) 周期性和圆心处有界。它们不写在题目里，却决定了特征值和解的形式。</p>`) +
  P.sec(ID, 'idea', '求解思路', R`<div class="prose small muted"><p>区域是圆盘 ⇒ 用极坐标，边界 \(r=R\) 变成坐标线。分离变量后角向是周期特征值问题（得 Fourier 级数），径向是 Euler 方程（得 \(r^n\)）。边界数据的 Fourier 系数乘以 \((r/R)^n\) 就是内部的解：越往里，高频成分衰减得越厉害，所以内部的解比边界数据光滑。</p></div>`) +
  P.sec(ID, 'deriv', '圆盘 Dirichlet 问题的推导', '<div id="lp-steps"></div>') +
  P.sec(ID, 'viz', '交互实验',
    `<h3 class="sub">B–D　圆盘 Dirichlet 问题、Poisson 积分与最大值原理</h3>` +
    P.lab('DISK', '点击圆盘选点；拖动 3D 图旋转',
      `<div class="ctl"><label class="ctl-label" for="lp-f">边界函数 \\(f(\\theta)\\)</label><select id="lp-f">${Object.entries(FNAME).map(([k, v]) => `<option value="${k}">${v}</option>`).join('')}</select></div>
       ${P.rangeHTML('lp-R', '半径 \\(R\\)', 0.5, 3, 0.1, 1)}
       ${P.rangeHTML('lp-N', '截断阶数 \\(N\\)', 1, 120, 1, 40)}
       ${P.rangeHTML('lp-pr', '选点半径 \\(r/R\\)', 0, 0.98, 0.01, 0.45)}
       ${P.rangeHTML('lp-pt', '选点角度 \\(\\theta\\)', 0, 6.28, 0.01, 0.9)}
       ${P.rangeHTML('lp-rho', '平均值圆半径 \\(\\rho/R\\)', 0.02, 0.9, 0.01, 0.25)}
       ${P.rangeHTML('lp-th0', '截面方向 \\(\\theta_0\\)', 0, 3.14, 0.01, 0.5)}
       <div class="tasks"><b>试一试</b><ul><li>把选点拖向边界：紫色核在最近的边界点处集中，级数与 Poisson 积分仍一致。</li><li>阶跃数据：边界上的部分和有 Gibbs 振荡，内部却很光滑。</li><li>任选数据：红点（内部最大）总在最大边界值之下；Hessian 两特征值异号。</li><li>选"手绘"后在 \\(f(\\theta)\\) 图上拖动，自己设计边界温度。</li></ul></div>`,
      `<div class="problem" id="lp-prob"></div>
       <p class="small" id="lp-coefSrc"></p>
       <div class="two">
         <div class="plot-block"><div class="plot-cap"><span>圆盘内的解 \\(u(x,y)\\)，环上紫色粗细 = Poisson 核权重</span></div><canvas id="lp-cvDisk" class="square clickable"></canvas></div>
         <div class="plot-block"><div class="plot-cap"><span>三维曲面（边界一圈为 f）</span></div><canvas id="lp-cv3d" class="square"></canvas></div>
       </div>
       <div class="two">
         ${P.plotBlock('lp-cvF', '边界数据 \\(f(\\theta)\\) 与 \\(r=R\\) 处的部分和', P.lg('warm', 'f(θ)') + P.lg('accent', 'S_N(θ)') + P.lg('violet', '选点角度'), 'short')}
         ${P.plotBlock('lp-cvRad', '沿直径的截面 \\(u(s\\cos\\theta_0, s\\sin\\theta_0)\\)', '', 'short')}
       </div>
       ${P.plotBlock('lp-cvK', 'C　Poisson 核 \\(2\\pi P(r,\\theta;\\varphi)\\)：边界各点对选点的权重（虚线为平均值 1）', '', 'short')}
       <div class="stats" id="lp-pStats"></div>
       <h3 class="sub">D　最大值原理</h3>
       <p class="small muted">圆盘图中：黑色空心圈 = 边界上的最大、最小值点；红点 / 蓝点 = 内部（r ≤ 0.95R）的最大 / 最小值点；虚线圆 = 平均值性质检验用的圆。</p>
       <div class="stats" id="lp-mStats"></div>`) +
    `<div class="note"><b>最大值原理（前提与结论）</b>　设 \\(\\Omega\\) 有界、\\(u\\in C^2(\\Omega)\\cap C(\\bar\\Omega)\\) 且 \\(\\Delta u=0\\)。<br>弱最大值原理：\\(\\max_{\\bar\\Omega}u=\\max_{\\partial\\Omega}u\\)。证明：令 \\(v=u+\\varepsilon|x|^2\\)，\\(\\Delta v=4\\varepsilon\\gt0\\)；若 \\(v\\) 在内部取最大值则该点 \\(\\Delta v\\le0\\)，矛盾；再令 \\(\\varepsilon\\to0\\)。<br>强最大值原理：若 \\(\\Omega\\) 连通且 \\(u\\) 在内部某点取到最大值，则 \\(u\\) 为常数（由平均值性质推出）。<br>图像只是验证，不能代替证明；网格检查也只能说明在采样点上没有反例。</div>` +
    `<h3 class="sub">E　Neumann 问题</h3>` +
    P.lab('NEUMANN', '\\(\\partial u/\\partial r=g(\\theta)\\)：法向导数 = 流入的热流密度',
      `<div class="ctl"><label class="ctl-label" for="lp-g">边界热流 \\(g(\\theta)\\)</label><select id="lp-g">${Object.entries(GNAME).map(([k, v]) => `<option value="${k}">${v}${k === 'bad' ? '（不相容）' : ''}</option>`).join('')}</select></div>
       ${P.rangeHTML('lp-C', '加法常数 \\(C\\)', -1, 1, 0.05, 0)}
       <div class="tasks"><b>试一试</b><ul><li>改变 \\(C\\)：边界导数不变，所以它们都是解。</li><li>选 cos θ + 0.5：净热流不为零，持续流入的热量无法达到稳态。</li></ul></div>`,
      `<div class="problem" id="lp-nProb"></div>
       <div class="plot-block"><div class="plot-cap"><span>解（箭头：红色为热量流出 g&gt;0 方向的法向导数，蓝色相反）</span></div><canvas id="lp-cvN" class="square"></canvas></div>
       <div class="stats" id="lp-nStats"></div>
       <p class="small muted">例题：\\(g=\\cos\\theta\\) 时 \\(u=C+r\\cos\\theta=C+x\\)。验证：\\(\\partial_ru|_{r=1}=\\cos\\theta\\)，且 \\(\\int_0^{2\\pi}\\cos\\theta\\,d\\theta=0\\)。</p>`)) +
  P.sec(ID, 'meaning', '物理意义', R`<div class="cols">
    <div class="card"><h3>平均值性质</h3><p class="small muted">稳态温度在任一点等于周围任一圆周上的平均温度：热量没有理由在某点聚集或亏空。</p></div>
    <div class="card"><h3>最大值原理</h3><p class="small muted">没有内部热源时，最高温和最低温都在边界上。若内部出现局部热点，热量会向四周流走，就不是稳态。</p></div>
    <div class="card"><h3>Poisson 核</h3><p class="small muted">点 \(P\) 的温度是边界温度的加权平均，离 \(P\) 越近的边界点权重越大。</p></div>
    <div class="card"><h3>\((r/R)^n\) 的光滑化</h3><p class="small muted">边界上的高频起伏在内部按 \((r/R)^n\) 衰减，所以即使边界数据间断，内部的解也是无穷光滑的。</p></div></div>` +
    P.pitfalls([
      ['极坐标公式：', String.raw`\(\Delta u=u_{rr}+\frac1ru_r+\frac1{r^2}u_{\theta\theta}\)，不要漏 \(\frac1ru_r\)。`],
      ['隐含条件：', '周期性给出 λ = n²；圆心有界舍去 r⁻ⁿ 和 ln r。圆环问题两者都要保留。'],
      ['n = 0 项：', String.raw`常数项是 \(\frac{a_0}2\)，\(a_0=\frac1\pi\int f\)。漏掉 1/2 会使中心温度加倍。`],
      ['最大值原理的前提：', String.raw`需要 \(u\) 在闭区域上连续、区域有界；对 \(\Delta u=f\neq0\) 不直接成立（\(f\le0\) 时关于最大值的部分仍成立，即上调和函数）。`],
      ['Neumann 问题：', '先检查 ∮g = 0；解不唯一，差一个常数。']
    ]) + `<div class="bench exp-only"><h3 class="sub">数值基准</h3><div id="lp-bench"></div></div>`) +
  P.sec(ID, 'quiz', '检查你的理解 <span class="score" id="lp-score"></span>', '<div id="lp-quiz"></div>');
}

/* ---------- 初始化 ---------- */
function recompute() { co = coefFor(S.f, fB); computeDisk(); computeMax(); render(); }
P.register('laplace', {
  init(root) {
    root.innerHTML = html();
    P.stepper($('lp-steps'), STEPS);
    $('lp-f').value = S.f;
    $('lp-f').addEventListener('change', e => {
      if (e.target.value === 'draw') for (let i = 0; i <= 120; i++) drawPts[i] = Math.max(-1.2, Math.min(1.2, fB(TAU * i / 120)));
      S.f = e.target.value; recompute();
    });
    const f2 = v => v.toFixed(2);
    P.bindRange('lp-R', v => { S.R = v; render(); }, f2);
    P.bindRange('lp-N', v => { S.N = v; computeDisk(); render(); });
    P.bindRange('lp-pr', v => { S.pr = v; render(); }, f2);
    P.bindRange('lp-pt', v => { S.pt = v; render(); }, f2);
    P.bindRange('lp-rho', v => { S.rho = v; render(); }, f2);
    P.bindRange('lp-th0', v => { S.th0 = v; render(); }, f2);
    $('lp-g').addEventListener('change', e => { S.neu = e.target.value; renderNeumann(); });
    P.bindRange('lp-C', v => { S.C = v; renderNeumann(); }, f2);
    const cv = $('lp-cvDisk');
    cv.addEventListener('pointerdown', e => {
      if (!cv._inv) return; const r = cv.getBoundingClientRect(), x = cv._inv.x(e.clientX - r.left) / S.R, y = cv._inv.y(e.clientY - r.top) / S.R, rr = Math.hypot(x, y);
      if (rr >= 1) return;
      S.pr = Math.min(0.98, rr); S.pt = ((Math.atan2(y, x) % TAU) + TAU) % TAU;
      $('lp-pr').value = S.pr; $('lp-pr-v').textContent = S.pr.toFixed(2); $('lp-pt').value = S.pt; $('lp-pt-v').textContent = S.pt.toFixed(2);
      render();
    });
    P.dragRotate($('lp-cv3d'), view, () => render());
    const fc = $('lp-cvF'); let last = null, pend = false;
    const paint = e => {
      const r = fc.getBoundingClientRect(), th = Math.max(0, Math.min(TAU, fc._inv.x(e.clientX - r.left))), y = Math.max(-1.2, Math.min(1.2, fc._inv.y(e.clientY - r.top))), i = Math.round(th / TAU * 120);
      if (last) { const [i0, y0] = last, d = Math.abs(i - i0); for (let k = 0; k <= d; k++) { const t = d ? k / d : 1; drawPts[Math.round(i0 + (i - i0) * t)] = y0 + (y - y0) * t; } } else drawPts[i] = y;
      drawPts[120] = drawPts[0] = (i === 0 || i === 120) ? y : drawPts[0];
      last = [i, y];
      if (!pend) { pend = true; requestAnimationFrame(() => { pend = false; recompute(); }); }
    };
    fc.addEventListener('pointerdown', e => { if (S.f !== 'draw' || !fc._inv) return; fc.setPointerCapture(e.pointerId); last = null; paint(e); });
    fc.addEventListener('pointermove', e => { if (S.f === 'draw' && fc.hasPointerCapture(e.pointerId)) paint(e); });
    fc.addEventListener('pointerup', () => { last = null; });
    P.quiz($('lp-quiz'), QUIZ, $('lp-score'));
    runBench();
    $('lp-bench').addEventListener('click', e => { if (e.target.closest('[data-rerun]')) { runBench(); recompute(); } });
    co = coefFor(S.f, fB); computeDisk(); computeMax();
  },
  render
});
})();
