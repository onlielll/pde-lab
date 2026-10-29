/* m3-heat.js — 模块 03：热传导方程 u_t = D u_xx + αu + f(x,t)，三类边界、非齐次问题与极值原理 */
(function () {
'use strict';
const P = window.PDE, PI = Math.PI, $ = P.$, fx = P.fx, sci = P.sci;
const ID = 'ht', NM = 60;
const S = { bcL: 'D', bcR: 'D', h1: 2, h2: 2, T0: 0, T1: 0, D: 1, alpha: 0, src: 'none', q: 10, Om: 30, init: 'gauss', m: 1, N: 40, L: 1 };
const PRE = {
  basic: { label: '高斯初温，两端 0°', set: { bcL: 'D', bcR: 'D', T0: 0, T1: 0, alpha: 0, src: 'none', init: 'gauss' } },
  pw: { label: '分段常数初温，两端绝热（热量守恒）', set: { bcL: 'N', bcR: 'N', alpha: 0, src: 'none', init: 'pw' } },
  robin: { label: '两端对流换热（牛顿冷却）', set: { bcL: 'R', bcR: 'R', h1: 2, h2: 2, alpha: 0, src: 'none', init: 'const' } },
  mode: { label: '单一正弦模态', set: { bcL: 'D', bcR: 'D', T0: 0, T1: 0, alpha: 0, src: 'none', init: 'sine', m: 2 } },
  nonhom: { label: '非齐次边界：左端 1°、右端 0.2°', set: { bcL: 'D', bcR: 'D', T0: 1, T1: 0.2, alpha: 0, src: 'none', init: 'zero' } },
  source: { label: '恒定均匀热源 f = q', set: { bcL: 'D', bcR: 'D', T0: 0, T1: 0, alpha: 0, src: 'const', q: 10, init: 'zero' } },
  growth: { label: '反应项 α = 12 > π²（失稳）', set: { bcL: 'D', bcR: 'D', T0: 0, T1: 0, alpha: 12, src: 'none', init: 'sine', m: 1 } },
  periodic: { label: '周期热源 q sin Ωt', set: { bcL: 'D', bcR: 'D', T0: 0, T1: 0, alpha: 0, src: 'periodic', q: 10, Om: 30, init: 'zero' } }
};
let sl, M, player;

/* ---------- 数学核心 ---------- */
function initF(xi) {
  switch (S.init) {
    case 'gauss': return Math.exp(-((xi - 0.4) ** 2) / (2 * 0.08 * 0.08));
    case 'pw': return xi < 0.3 ? 1 : xi > 0.3 && xi < 0.6 ? 0.2 : xi === 0.3 ? 0.6 : xi === 0.6 ? 0.45 : 0.7;
    case 'sine': return Math.sin(S.m * PI * xi);
    case 'tri': return 1 - Math.abs(2 * xi - 1);
    case 'const': return 1;
    default: return 0;
  }
}
function srcF(xi) {
  if (S.src === 'const' || S.src === 'periodic') return S.q;
  if (S.src === 'local') return S.q * Math.exp(-((xi - 0.5) ** 2) / (2 * 0.05 * 0.05));
  return 0;
}
// 边界提升 w = p + s x：满足非齐次边界条件（D 端取给定温度，N/R 端齐次）
function lift() {
  const { a1, b1, a2, b2 } = sl, g0 = S.bcL === 'D' ? S.T0 : 0, g1 = S.bcR === 'D' ? S.T1 : 0;
  if (!g0 && !g1) return { p: 0, s: 0, ok: true };
  const det = a1 * (a2 * S.L + b2) + b1 * a2;
  if (Math.abs(det) < 1e-12) return { p: 0, s: 0, ok: false };
  return { p: (g0 * (a2 * S.L + b2) + b1 * g1) / det, s: (a1 * g1 - a2 * g0) / det, ok: true };
}
function update() {
  sl = P.solveSL({ L: S.L, bcL: S.bcL, bcR: S.bcR, h1: S.h1, h2: S.h2, NM, NG: 1200 });
  const NG = sl.NG, phi = new Float64Array(NG + 1), w = new Float64Array(NG + 1), s = new Float64Array(NG + 1);
  const lf = lift();
  for (let i = 0; i <= NG; i++) { const xi = i / NG; phi[i] = initF(xi); w[i] = lf.p + lf.s * sl.xs[i]; s[i] = srcF(xi); }
  const v0 = phi.map((v, i) => v - w[i]);
  M = { phi, w, s, lf, V0: P.project(sl, v0), Sn: P.project(sl, s), Wn: P.project(sl, w), sig: sl.modes.map(m => S.alpha - S.D * m.lam) };
  M.tmax = tmaxCalc();
}
function I1(sg, t) { return Math.abs(sg) < 1e-12 ? t : Math.expm1(sg * t) / sg; }
function Isin(sg, Om, t) { return (-sg * Math.sin(Om * t) - Om * Math.cos(Om * t) + Om * Math.exp(sg * t)) / (sg * sg + Om * Om); }
function vn(j, t) {
  const sg = M.sig[j];
  let v = Math.exp(sg * t) * M.V0[j] + S.alpha * M.Wn[j] * I1(sg, t);
  if (S.src === 'periodic') v += M.Sn[j] * Isin(sg, S.Om, t); else if (S.src !== 'none') v += M.Sn[j] * I1(sg, t);
  return v;
}
function uAt(t, N = S.N, deriv = false) {
  const NG = sl.NG, out = new Float64Array(NG + 1);
  for (let i = 0; i <= NG; i++) out[i] = deriv ? M.lf.s : M.w[i];
  for (let j = 0; j < Math.min(N, sl.modes.length); j++) {
    const c = vn(j, t); if (c === 0) continue;
    const X = deriv ? sl.modes[j].dX : sl.modes[j].X;
    for (let i = 0; i <= NG; i++) out[i] += c * X[i];
  }
  return out;
}
function steady(N = S.N) {
  if (S.src === 'periodic') return null;
  const out = new Float64Array(sl.NG + 1); for (let i = 0; i <= sl.NG; i++) out[i] = M.w[i];
  for (let j = 0; j < Math.min(N, sl.modes.length); j++) {
    const sg = M.sig[j], F = (S.src !== 'none' ? M.Sn[j] : 0) + S.alpha * M.Wn[j];
    let c;
    if (sg < -1e-12) c = -F / sg;
    else if (Math.abs(sg) <= 1e-12 && Math.abs(F) < 1e-12) c = M.V0[j];
    else return null;
    const X = sl.modes[j].X; for (let i = 0; i <= sl.NG; i++) out[i] += c * X[i];
  }
  return out;
}
function tmaxCalc() {
  const s1 = M.sig.find(s => Math.abs(s) > 1e-9) ?? -1;
  let tm = 1.5 / Math.abs(s1);
  if (S.src === 'periodic') tm = Math.max(tm, 6 * 2 * PI / S.Om);
  return Math.min(tm, 5);
}
// 总热量收支：dQ/dt = D[u_x]_0^L + α∫u + ∫f
function heatBudget(t) {
  const u = uAt(t), ux = uAt(t, S.N, true), h = sl.h, NG = sl.NG;
  const tf = S.src === 'periodic' ? Math.sin(S.Om * t) : 1;
  return { Q: P.integ(u, h), J0: -S.D * ux[0], JL: -S.D * ux[NG], src: S.src === 'none' ? 0 : tf * P.integ(M.s, h), react: S.alpha * P.integ(u, h), u, ux };
}

/* ---------- 渲染：主实验 ---------- */
const BCU = {
  L: { D: () => `u(0,t)=${fx(S.T0, 2)}`, N: () => 'u_x(0,t)=0', R: () => 'u_x(0,t)-h_1u(0,t)=0' },
  R: { D: () => `u(L,t)=${fx(S.T1, 2)}`, N: () => 'u_x(L,t)=0', R: () => 'u_x(L,t)+h_2u(L,t)=0' }
};
function srcTex() { return { none: '0', const: `${fx(S.q, 1)}`, local: `${fx(S.q, 1)}e^{-(x-0.5)^2/0.005}`, periodic: `${fx(S.q, 1)}\\sin(${fx(S.Om, 0)}t)` }[S.src]; }
function render() {
  const t = player.t(), R = String.raw;
  P.tex($('ht-prob'), R`\[ \begin{cases} u_t=${fx(S.D, 2)}\,u_{xx}${S.alpha ? (S.alpha > 0 ? '+' : '') + fx(S.alpha, 1) + 'u' : ''}+f,\quad f=${srcTex()},\quad 0\lt x\lt ${fx(S.L, 2)},\\ ${BCU.L[S.bcL]()},\qquad ${BCU.R[S.bcR]()}\end{cases}\qquad u=\underbrace{w(x)}_{\text{边界提升}}+\sum_n v_n(t)X_n(x),\ \ v_n'=\sigma_nv_n+F_n(t),\ \sigma_n=\alpha-D\lambda_n \]`);
  const lfNote = M.lf.ok ? (M.lf.s || M.lf.p ? `\\(w(x)=${fx(M.lf.p, 4)}${M.lf.s >= 0 ? '+' : ''}${fx(M.lf.s, 4)}x\\) 满足非齐次边界，\\(v=u-w\\) 的边界条件齐次。${S.alpha ? '注意 \\(\\alpha w\\) 成为 \\(v\\) 方程的附加源项。' : ''}` : '边界条件齐次，无需提升。') : '<span class="pill bad">该边界组合下提升函数不存在（零特征值情形），已按齐次边界计算</span>';
  P.tex($('ht-lift'), lfNote);
  const u = uAt(t), st = steady();
  let lo = Math.min(0, ...M.phi), hi = Math.max(...M.phi, S.T0, S.T1);
  if (st) { lo = Math.min(lo, ...st); hi = Math.max(hi, ...st); }
  for (let i = 0; i < u.length; i++) { if (u[i] > hi) hi = u[i]; if (u[i] < lo) lo = u[i]; }
  const pad = 0.12 * (hi - lo || 1);
  const f = P.frame($('ht-cvU'), { xmin: 0, xmax: S.L, ymin: lo - pad, ymax: hi + pad, xlabel: 'x' });
  if (f) {
    P.line(f, sl.xs, M.phi, P.C.warm, 1.2, [4, 4], 2);
    if (st) P.line(f, sl.xs, st, P.C.good, 1.6, [2, 3], 2);
    P.line(f, sl.xs, u, P.C.accent, 2.4, null, 2);
    P.bcGlyphs(f, S.bcL, S.bcR);
  }
  const b = heatBudget(t), sg1 = M.sig[0];
  P.stats($('ht-stats'), [
    ['t', fx(t, 4)],
    ['总热量 Q(t) = ∫u dx', fx(b.Q, 5)],
    ['左端流入 −Du_x(0)', fx(b.J0, 4)],
    ['右端流出 −Du_x(L)', fx(b.JL, 4)],
    ['最慢模态 σ₁ = α − Dλ₁', fx(sg1, 4) + (sg1 > 1e-9 ? ' <span class="pill bad">增长</span>' : Math.abs(sg1) <= 1e-9 ? ' <span class="pill info">中性</span>' : ' <span class="pill ok">衰减</span>')],
    ['弛豫时间 1/|σ₁|', Math.abs(sg1) > 1e-9 ? fx(1 / Math.abs(sg1), 4) : '∞'],
    ['稳态', S.src === 'periodic' ? '周期强迫：趋于周期解' : st ? '存在（绿色虚线）' : '<span class="pill warn">不存在</span>'],
    ['热量守恒', S.bcL === 'N' && S.bcR === 'N' && !S.alpha && S.src === 'none' ? '<span class="pill ok">守恒（绝热、无源）</span>' : '不守恒']
  ]);
  if (S.bcL === 'D' && S.bcR === 'D') {
    const Lc = S.alpha > 0 ? PI * Math.sqrt(S.D / S.alpha) : Infinity;
    $('ht-crit').innerHTML = S.alpha > 0 ? `两端 D 时，σ₁ = α − Dπ²/L² &gt; 0 ⇔ L &gt; L_c = π√(D/α) = <span class="mono">${fx(Lc, 4)}</span>（当前 L = ${fx(S.L, 2)}）。反应项产生的热量超过从两端散失的热量时，温度指数增长。` : 'α ≤ 0 时所有 σ_n &lt; 0，解衰减到稳态。';
  } else $('ht-crit').textContent = '临界长度的闭式只针对两端 D 给出；其他边界看 σ₁ 的符号即可。';
  // 模态衰减
  const nb = Math.min(20, sl.modes.length), amp0 = Array.from(M.V0.slice(0, nb)).map(Math.abs), ampt = Array.from({ length: nb }, (_, j) => Math.abs(vn(j, t)));
  const cm = Math.max(...amp0, ...ampt.map(v => Math.min(v, 1e6))) || 1;
  const g = P.frame($('ht-cvM'), { xmin: sl.offset - 0.8, xmax: nb + sl.offset - 0.2, ymin: 0, ymax: cm * 1.12, xlabel: 'n', xticks: P.niceTicks(sl.offset, nb - 1 + sl.offset, 5).filter(v => v % 1 === 0).map(v => ({ v, label: v + '' })) });
  if (g) {
    const bw = Math.max(3, g.pw / nb * 0.55), xn = amp0.map((_, j) => j + sl.offset);
    P.bars(g, xn, amp0, P.C['ink-3'], bw, 0, true);
    P.bars(g, xn, ampt.map(v => Math.min(v, cm * 1.12)), j => j < S.N ? P.C.accent : P.C['ink-3'], bw);
  }
  renderQ(t); renderMax(t);
}
let Qhist = null;
function renderQ(t) {
  if (!Qhist) {
    const n = 160, ts = [], Q = [], bal = []; const tm = M.tmax;
    for (let i = 0; i <= n; i++) { const tt = tm * (i / n) ** 2; const b = heatBudget(tt); ts.push(tt); Q.push(b.Q); bal.push(b.J0 - b.JL + b.src + b.react); }
    Qhist = { ts, Q, bal };
  }
  const lo = Math.min(...Qhist.Q), hi = Math.max(...Qhist.Q), pad = 0.15 * (hi - lo || Math.abs(hi) || 1);
  const f = P.frame($('ht-cvQ'), { xmin: 0, xmax: M.tmax, ymin: lo - pad, ymax: hi + pad, xlabel: 't' });
  if (!f) return;
  P.line(f, Qhist.ts, Qhist.Q, P.C.accent, 2);
  P.seg2(f, t, f.o.ymin, t, f.o.ymax, P.C.ink, 1, [3, 3]);
}

/* ---------- 极值原理 ---------- */
let MX = null, mxId = 0;
function computeMax() {
  const nx = 160, nt = 120, tm = M.tmax, step = sl.NG / nx;
  const lamN = sl.modes[NM - 1].lam, t1 = Math.max(tm / nt * 0.5, 14 / (S.D * lamN)); // 截断误差 < e^{-14}
  const grid = [], ts = [];
  for (let j = 0; j <= nt; j++) {
    const t = tm * j / nt; ts.push(t);
    const row = j === 0 ? M.phi : uAt(t, NM);
    const r = new Float64Array(nx + 1); for (let i = 0; i <= nx; i++) r[i] = row[Math.round(i * step)];
    if (j > 0) { if (S.bcL === 'D') r[0] = S.T0; if (S.bcR === 'D') r[nx] = S.T1; }
    grid.push(r);
  }
  let bMax = -Infinity, bMin = Infinity, iMax = -Infinity, iMin = Infinity, pMax = null, pMin = null, ipMax = null, ipMin = null;
  for (let j = 0; j <= nt; j++) for (let i = 0; i <= nx; i++) {
    const v = grid[j][i], onB = j === 0 || i === 0 || i === nx;
    if (onB) { if (v > bMax) { bMax = v; pMax = [i, j]; } if (v < bMin) { bMin = v; pMin = [i, j]; } }
    else if (ts[j] >= t1) { if (v > iMax) { iMax = v; ipMax = [i, j]; } if (v < iMin) { iMin = v; ipMin = [i, j]; } }
  }
  const scale = Math.max(Math.abs(bMax), Math.abs(bMin), 1e-12), tol = 1e-7 * scale;
  MX = { id: ++mxId, grid, nx, nt, tm, t1, bMax, bMin, iMax, iMin, pMax, pMin, ipMax, ipMin, okMax: iMax <= bMax + tol, okMin: iMin >= bMin - tol };
}
function renderMax(t) {
  if (!MX) computeMax();
  const { grid, nx, nt, tm } = MX;
  // 极值图与 t 无关：缓存为静态图层，播放时只重画时间线
  const f = P.layer($('ht-cvXT'), 'mx' + MX.id, { xmin: 0, xmax: S.L, ymin: 0, ymax: tm, xlabel: 'x', ylabel: 't', colorbar: true }, f => {
    const d = new Float64Array((nx + 1) * (nt + 1)); let lo = Infinity, hi = -Infinity;
    for (let j = 0; j <= nt; j++) for (let i = 0; i <= nx; i++) { const v = grid[nt - j][i]; d[j * (nx + 1) + i] = v; lo = Math.min(lo, v); hi = Math.max(hi, v); }
    P.heatmap(f, d, nx + 1, nt + 1, lo, hi, 12);
    const X = i => i / nx * S.L, T = j => j / nt * tm;
    f.ctx.lineWidth = 3; f.ctx.strokeStyle = P.C.good; f.ctx.beginPath();
    f.ctx.moveTo(f.X(0), f.Y(tm)); f.ctx.lineTo(f.X(0), f.Y(0)); f.ctx.lineTo(f.X(S.L), f.Y(0)); f.ctx.lineTo(f.X(S.L), f.Y(tm)); f.ctx.stroke();
    if (MX.pMax) P.dot(f, X(MX.pMax[0]), T(MX.pMax[1]), P.C.ink, 6, true);
    if (MX.ipMax) P.dot(f, X(MX.ipMax[0]), T(MX.ipMax[1]), MX.okMax ? P.C['ink-2'] : P.C.bad, 4);
    if (MX.pMin) P.dot(f, X(MX.pMin[0]), T(MX.pMin[1]), P.C.violet, 6, true);
  });
  if (f) P.seg2(f, 0, t, S.L, t, P.C.ink, 1, [3, 3]);
  const hyp = !S.alpha && S.src === 'none';
  P.stats($('ht-maxStats'), [
    ['定理前提（α=0, f=0）', hyp ? '<span class="pill ok">满足</span>' : '<span class="pill warn">不满足</span>'],
    ['抛物边界上的最大值', fx(MX.bMax, 5)],
    ['内部（t ≥ t₁）最大值', fx(MX.iMax, 5)],
    ['最大值判定', MX.okMax ? '<span class="pill ok">内部 ≤ 边界</span>' : `<span class="pill bad">内部超出 ${sci(MX.iMax - MX.bMax)}</span>`],
    ['抛物边界上的最小值', fx(MX.bMin, 5)],
    ['最小值判定', MX.okMin ? '<span class="pill ok">内部 ≥ 边界</span>' : `<span class="pill bad">内部低于 ${sci(MX.bMin - MX.iMin)}</span>`],
    ['排除的初始薄层 t₁', sci(MX.t1)]
  ]);
  $('ht-maxNote').innerHTML = hyp
    ? (MX.okMax && MX.okMin ? '前提成立，数值结果与定理一致：最大值、最小值都落在抛物边界（绿色折线：底边 t=0 与两条侧边）上。' : '前提成立但数值检查未通过：请减小 N 截断带来的影响。')
    : (MX.okMax && MX.okMin ? '前提不成立，但这组数据恰好仍满足结论。定理不保证这一点，不能由图像反推定理在此情形成立。' : '前提不成立，结论确实失效：热源或反应项在区域内部产生热量，最高温出现在内部。');
}

/* ---------- 基准 ---------- */
function runBench() {
  const keep = { ...S }, rows = [];
  const set = o => { Object.assign(S, { bcL: 'D', bcR: 'D', h1: 2, h2: 2, T0: 0, T1: 0, D: 1, alpha: 0, src: 'none', q: 10, Om: 30, init: 'gauss', m: 1, N: 60, L: 1 }, o); update(); };
  set({ init: 'sine', m: 1 });
  let e = 0; [0.01, 0.05, 0.2].forEach(t => { const u = uAt(t); for (let i = 0; i <= sl.NG; i += 10) e = Math.max(e, Math.abs(u[i] - Math.exp(-PI * PI * t) * Math.sin(PI * sl.xs[i]))); });
  rows.push(P.benchRow('单模态 sin πx：u = e^{−π²Dt} sin πx', '解析解', e, 1e-12));
  // 窄高斯 vs 无界热核解
  set({ init: 'gauss' }); const s0 = 0.05; M.phi = M.phi.map((_, i) => Math.exp(-((sl.xs[i] - 0.5) ** 2) / (2 * s0 * s0))); M.V0 = P.project(sl, M.phi);
  e = 0; const tg = 0.002, sv = s0 * s0 + 2 * tg, u = uAt(tg);
  for (let i = 0; i <= sl.NG; i += 4) e = Math.max(e, Math.abs(u[i] - s0 / Math.sqrt(sv) * Math.exp(-((sl.xs[i] - 0.5) ** 2) / (2 * sv))));
  rows.push(P.benchRow('窄高斯短时间：有界杆 vs 无界热核解', '无界解；边界影响 ~e^{−27}，误差主要来自系数的 Simpson 积分 O(h⁴k⁴)', e, 1e-8));
  set({ bcL: 'N', bcR: 'N', init: 'pw' });
  const Q0 = P.integ(M.phi, sl.h); e = 0; [0.001, 0.05, 0.5, 2].forEach(t => { e = Math.max(e, Math.abs(P.integ(uAt(t), sl.h) - Q0)); });
  rows.push(P.benchRow('N–N 绝热：总热量 Q(t) 守恒', '∫u dx 与 ∫φ dx 比较', e, 1e-12));
  const ss = uAt(50);
  rows.push(P.benchRow('N–N：t→∞ 趋于初始平均温度', '平均值 ∫φ/L', Math.abs(Math.max(...ss) - Q0) + Math.abs(Math.min(...ss) - Q0), 1e-10));
  set({ src: 'const', q: 10, init: 'zero' });
  let st = steady(60); e = 0; for (let i = 0; i <= sl.NG; i++) { const x = sl.xs[i]; e = Math.max(e, Math.abs(st[i] - 10 * x * (1 - x) / 2)); }
  let tail = 0; for (let n = 61; n < 400000; n += 2) tail += 10 * 4 / (PI ** 3 * n ** 3);
  rows.push(P.benchRow('恒定热源稳态 u_s = q x(L−x)/(2D)', '解析稳态；容差 = 截断尾项上界 Σ_{n>60}|c_n|', e, tail));
  set({ T0: 1, T1: 0.2, init: 'zero' });
  const ul = uAt(30); e = 0; for (let i = 0; i <= sl.NG; i++) e = Math.max(e, Math.abs(ul[i] - (1 - 0.8 * sl.xs[i])));
  rows.push(P.benchRow('非齐次边界：长时间趋于线性稳态 1 − 0.8x', '解析稳态', e, 1e-12));
  set({ alpha: 12, init: 'sine', m: 1 });
  const ug = uAt(0.5); e = Math.abs(ug[600] - Math.exp((12 - PI * PI) * 0.5));
  rows.push(P.benchRow('反应项 α=12：sin πx 按 e^{(α−π²)t} 增长', '解析解', e, 1e-12));
  set({ bcL: 'R', bcR: 'D', h1: 3, T1: 0.5, src: 'local', q: 8, init: 'gauss' });
  const tb = 0.03, dtb = 1e-5, b = heatBudget(tb), dQ = (P.integ(uAt(tb + dtb), sl.h) - P.integ(uAt(tb - dtb), sl.h)) / (2 * dtb);
  rows.push(P.benchRow('热量收支 dQ/dt = 边界热流 + 源项（R–D，局部热源）', '时间差分 vs 边界通量公式', Math.abs(dQ - (b.J0 - b.JL + b.src)) / Math.abs(dQ), 1e-3));
  set({ init: 'pw' }); computeMax();
  rows.push({ name: '极值原理（α=0, f=0，分段常数初值）', method: `内部最大 ${fx(MX.iMax, 6)} ≤ 边界最大 ${fx(MX.bMax, 6)}`, err: Math.max(0, MX.iMax - MX.bMax), tol: 1e-7, pass: MX.okMax && MX.okMin });
  set({ src: 'const', q: 30, init: 'zero' }); computeMax();
  rows.push({ name: '反例：有热源时最大值出现在内部', method: `内部最大 ${fx(MX.iMax, 4)} > 边界最大 ${fx(MX.bMax, 4)}`, err: MX.okMax ? 'bad' : 0, tol: 0, pass: !MX.okMax });
  Object.assign(S, keep); update(); Qhist = null; MX = null;
  P.benchTable($('ht-bench'), 'heat', rows);
}

/* ---------- 文字 ---------- */
const QUIZ = [
  { type: '预测演化', q: '两端保持 0°，初温为高 1 的方波。哪个现象最先发生？', o: ['方波整体平移', '尖角立即被抹平，高频模态最先消失', '温度先升高后降低', '端点处出现温度跳跃并保持'], a: 1, e: String.raw`第 \(n\) 个模态按 \(e^{-D(n\pi/L)^2t}\) 衰减，高频衰减得最快，所以任意 \(t\gt0\) 时解都是光滑的。` },
  { type: '边界条件', q: '两端绝热、无热源时，以下哪个量守恒？', o: ['最高温度', '温度的 L² 范数', '总热量 ∫u dx', '中点温度'], a: 2, e: String.raw`\(\frac{d}{dt}\int u\,dx=D[u_x]_0^L=0\)。数学上对应 N–N 的零特征值模态 \(X_0=1\) 不衰减。` },
  { type: '非齐次边界', q: String.raw`\(u(0,t)=1,\ u(L,t)=0.2\)。分离变量前应先做什么？`, o: ['直接设 u = XT', '减去满足边界的稳态 w = 1 − 0.8x/L', '对初始条件求导', '把边界值加到初始条件上'], a: 1, e: String.raw`\(v=u-w\) 满足齐次边界，\(v(x,0)=\varphi-w\)。若方程含 \(\alpha u\)，\(v\) 的方程会多出源项 \(\alpha w\)，本页的计算包含这一项。` },
  { type: '反应项', q: String.raw`两端 0°，\(u_t=Du_{xx}+\alpha u\)，\(\alpha\gt0\)。温度何时会指数增长？`, o: ['只要 α > 0', String.raw`当 \(\alpha\gt D\pi^2/L^2\)`, '当 α > D', '永远不会'], a: 1, e: String.raw`最慢模态的增长率 \(\sigma_1=\alpha-D\pi^2/L^2\)。杆越长，边界散热越弱；\(L\gt\pi\sqrt{D/\alpha}\) 时失稳（核反应堆临界尺寸的一维模型）。` },
  { type: '极值原理', q: '对 u_t = Du_xx 在 [0,L]×[0,T] 上，最大值原理断言什么？', o: ['最大值在 t = T 时取到', '最大值在抛物边界（t=0 或 x=0、x=L）上取到', '最大值在区域中心取到', '最大值随时间增长'], a: 1, e: String.raw`弱最大值原理：\(\max_{\bar Q_T}u=\max_{\Gamma_T}u\)，\(\Gamma_T\) 为底边和两条侧边，不含顶边 \(t=T\)。前提是无源、无反应项（或 \(f\le0\)）。` },
  { type: '物理解释', q: '高阶模态为什么衰减更快？', o: ['因为它们的初始系数更小', '因为空间变化越剧烈，温度梯度越大，热流把起伏抹平得越快', '因为边界条件只作用于高阶模态', '因为数值误差'], a: 1, e: String.raw`衰减率 \(D\lambda_n\propto n^2\)：波长减半，梯度加倍，同时需要抹平的距离减半，衰减快 4 倍。` }
];
const STEPS = [
  { t: '建立模型', b: String.raw`<p>长 \(L\) 的均匀细杆，侧面绝热，两端保持 0°，初温处处为 \(u_0\)：</p><div class="formula-box">\[ u_t=Du_{xx},\quad u(0,t)=u(L,t)=0,\quad u(x,0)=u_0 . \]</div>` },
  { t: '分离变量', b: String.raw`<div class="formula-box">\[ u=X(x)T(t):\quad \frac{T'}{DT}=\frac{X''}{X}=-\lambda,\qquad X''+\lambda X=0,\ X(0)=X(L)=0 . \]</div>` },
  { t: '特征值问题（模块 04）', b: String.raw`<div class="formula-box">\[ \lambda_n=\Big(\frac{n\pi}{L}\Big)^2,\quad X_n=\sin\frac{n\pi x}{L},\quad T_n=e^{-D\lambda_nt}. \]</div><p>时间部分是一阶 ODE，所以出现指数而不是三角函数（模块 06）。</p>` },
  { t: '叠加并确定系数', b: String.raw`<div class="formula-box">\[ b_n=\frac2L\int_0^Lu_0\sin\frac{n\pi x}{L}dx=\frac{2u_0}{n\pi}\big(1-(-1)^n\big)=\begin{cases}\dfrac{4u_0}{n\pi},&n\text{ 奇}\\0,&n\text{ 偶}\end{cases} \]</div>` },
  { t: '写出解', b: String.raw`<div class="formula-box">\[ u(x,t)=\frac{4u_0}{\pi}\sum_{k=0}^{\infty}\frac{1}{2k+1}\,e^{-D(2k+1)^2\pi^2t/L^2}\sin\frac{(2k+1)\pi x}{L}. \]</div><p>初值与边界在角点 \((0,0)\) 不相容（\(u_0\neq0\)），但对任意 \(t\gt0\) 级数快速收敛。</p>` },
  { t: '长期行为', b: String.raw`<p>\(t\) 较大时只剩第一项：\(u\approx\frac{4u_0}{\pi}e^{-D\pi^2t/L^2}\sin\frac{\pi x}{L}\)。中点温度降到一半所需时间约 \(t_{1/2}\approx\frac{L^2}{D\pi^2}\ln\frac{8}{\pi}\)。杆长加倍，冷却时间变为 4 倍。</p>` },
  { t: '推广', b: String.raw`<p>边界非齐次：先减去稳态 \(w\)。有源项 \(f\)：把 \(f\) 也按 \(X_n\) 展开，每个模态解一阶线性 ODE \(v_n'=-D\lambda_nv_n+f_n(t)\)，用常数变易法（即 Duhamel 原理）得：</p><div class="formula-box">\[ v_n(t)=e^{-D\lambda_nt}v_n(0)+\int_0^te^{-D\lambda_n(t-\tau)}f_n(\tau)\,d\tau . \]</div>` }
];
function html() {
  const R = String.raw;
  return P.modHead(ID, '03', '热传导方程与扩散问题', R`研究 \(u_t=Du_{xx}+\alpha u+f(x,t)\)。重点是：模态为什么衰减、边界条件如何决定热量收支和稳态、非齐次问题如何分解，以及最大值原理在什么条件下成立。`) +
  P.sec(ID, 'phys', '热量守恒 + Fourier 定律', R`<div class="cols">
    <div class="card"><h3>Fourier 热传导定律</h3><p class="small muted">热流密度与温度梯度成正比、方向相反：\(q=-\kappa u_x\)。热量从高温流向低温。</p></div>
    <div class="card"><h3>能量守恒</h3><p class="small muted">对小段 \([x,x+dx]\)：\(c\rho u_t\,dx=-[q(x+dx)-q(x)]+(\text{热源})\,dx\)，代入 Fourier 定律得 \(u_t=Du_{xx}+f\)，\(D=\kappa/(c\rho)\)。</p></div>
    <div class="card"><h3>反应项 αu</h3><p class="small muted">产热率与温度成正比（\(\alpha\gt0\)，如放热反应、中子增殖）或侧面按温度散热（\(\alpha\lt0\)，侧面不绝热时的牛顿冷却）。</p></div></div>`) +
  P.sec(ID, 'model', '数学模型', R`<div class="formula-box">\[ u_t=D\,\Delta u+\alpha u+f(x,t)\ \xrightarrow{\ \text{一维}\ }\ u_t=Du_{xx}+\alpha u+f(x,t),\qquad 0\lt x\lt L,\ t\gt0 . \]</div><p class="prose small muted">抛物型方程：时间一阶，空间二阶。只需一个初始条件 \(u(x,0)=\varphi(x)\)。本模块以一维为主；高维问题在矩形、圆盘上同样可用分离变量。</p>`) +
  P.sec(ID, 'cond', '三类边界条件的热学含义', R`<div class="tbl-wrap"><table><thead><tr><th>类型</th><th>形式</th><th>物理含义</th><th>对总热量的影响</th></tr></thead><tbody>
    <tr><td>Dirichlet</td><td>\(u(0,t)=T_0\)</td><td>端点接恒温热库</td><td>热量可以自由进出</td></tr>
    <tr><td>Neumann</td><td>\(u_x=0\)（或给定热流）</td><td>绝热（或给定热流）</td><td>两端绝热且无源时 \(Q\) 守恒</td></tr>
    <tr><td>Robin</td><td>\(-\kappa\partial_nu=H(u-u_{\text{env}})\)</td><td>对流换热（牛顿冷却）</td><td>热量按 \(hu\) 的速率流出</td></tr></tbody></table></div>`) +
  P.sec(ID, 'idea', '求解思路：稳态 + 瞬态，再按模态分解', R`<div class="prose small muted"><p>(1) 用线性函数 \(w(x)\) 吸收非齐次边界值，\(v=u-w\) 的边界齐次。(2) 把 \(v\) 和源项都按模块 04 的特征函数展开：\(v=\sum v_n(t)X_n\)。(3) 每个系数满足一阶 ODE \(v_n'=\sigma_nv_n+F_n(t)\)，\(\sigma_n=\alpha-D\lambda_n\)。(4) 恒定、简谐的 \(F_n\) 都有闭式解，本页逐模态精确计算，唯一的近似是截断阶数 \(N\)。</p></div>`) +
  P.sec(ID, 'deriv', '典型例题：两端 0°、初温均匀的杆', '<div id="ht-steps"></div>') +
  P.sec(ID, 'viz', '交互实验',
    `<h3 class="sub">A–D　扩散、模态衰减、边界条件与非齐次问题</h3>` +
    P.lab('HEAT', '\\(u_t=Du_{xx}+\\alpha u+f\\)，逐模态精确求解',
      `<div class="ctl"><label class="ctl-label" for="ht-pre">典型情形</label><select id="ht-pre">${Object.entries(PRE).map(([k, p]) => `<option value="${k}">${p.label}</option>`).join('')}</select></div>
       <div class="ctl"><label class="ctl-label" for="ht-init">初始温度 \\(\\varphi\\)</label><select id="ht-init"><option value="gauss">高斯</option><option value="pw">分段常数</option><option value="sine">正弦模态 sin mπx/L</option><option value="tri">三角形</option><option value="const">常数 1</option><option value="zero">0</option></select></div>
       ${P.rangeHTML('ht-m', '模态数 \\(m\\)', 1, 6, 1, 1)}
       <div class="ctl"><div class="ctl-label">左端</div><div id="ht-bcL"></div></div>
       <div class="ctl"><div class="ctl-label">右端</div><div id="ht-bcR"></div></div>
       ${P.rangeHTML('ht-h1', '左端换热系数 \\(h_1\\)', 0, 10, 0.1, 2)}
       ${P.rangeHTML('ht-h2', '右端换热系数 \\(h_2\\)', 0, 10, 0.1, 2)}
       ${P.rangeHTML('ht-T0', '左端温度 \\(T_0\\)', -1, 1, 0.05, 0)}
       ${P.rangeHTML('ht-T1', '右端温度 \\(T_1\\)', -1, 1, 0.05, 0)}
       ${P.rangeHTML('ht-D', '扩散系数 \\(D\\)', 0.1, 2, 0.05, 1)}
       ${P.rangeHTML('ht-alpha', '反应系数 \\(\\alpha\\)', -20, 20, 0.5, 0)}
       <div class="ctl"><label class="ctl-label" for="ht-src">热源 \\(f\\)</label><select id="ht-src"><option value="none">无</option><option value="const">均匀恒定 q</option><option value="local">中点局部恒定</option><option value="periodic">均匀周期 q sin Ωt</option></select></div>
       ${P.rangeHTML('ht-q', '热源强度 \\(q\\)', -30, 30, 0.5, 10)}
       ${P.rangeHTML('ht-Om', '热源频率 \\(\\Omega\\)', 2, 80, 1, 30)}
       ${P.rangeHTML('ht-N', '截断阶数 \\(N\\)', 1, 60, 1, 40)}
       <div id="ht-player"></div>`,
      `<div class="problem" id="ht-prob"></div>
       <p class="small" id="ht-lift"></p>
       ${P.plotBlock('ht-cvU', '温度分布 \\(u(x,t)\\)', P.lg('accent', 'u(x,t)') + P.lg('warm', 'φ(x)') + P.lg('good', '稳态'))}
       <div class="stats" id="ht-stats"></div>
       <p class="small muted" id="ht-crit"></p>
       <div class="two">
         ${P.plotBlock('ht-cvM', '各模态 \\(|v_n(t)|\\)（空心为 \\(t=0\\)）：高阶先消失', '', 'short')}
         ${P.plotBlock('ht-cvQ', '总热量 \\(Q(t)=\\int_0^Lu\\,dx\\)', '', 'short')}
       </div>`) +
    `<h3 class="sub">E　极值原理</h3>` +
    P.lab('MAX', '在 x–t 矩形上检查最大值与最小值的位置（与上方实验共用参数）',
      `<div class="small muted">绿色折线为抛物边界 \\(\\Gamma_T\\)：底边 \\(t=0\\) 和两条侧边。黑色圆圈为整个闭矩形上的最大值位置，紫色圆圈为最小值位置，实心点为内部最大值（红色表示超过边界最大值）。</div>
       <div class="btn-row"><button type="button" class="btn ghost" data-pre="pw">纯扩散</button><button type="button" class="btn ghost" data-pre="source">有热源</button><button type="button" class="btn ghost" data-pre="growth">α &gt; 0</button></div>
       <div class="formula-box small">\\[ u_t-Du_{xx}\\le0\\ \\text{于}\\ Q_T\\ \\Longrightarrow\\ \\max_{\\bar Q_T}u=\\max_{\\Gamma_T}u \\]</div>
       <p class="small muted">证明思路：令 \\(v=u+\\varepsilon x^2\\)，则 \\(v_t-Dv_{xx}\\lt0\\)。若 \\(v\\) 在内部或顶边取最大值，那里 \\(v_t\\ge0\\)、\\(v_{xx}\\le0\\)，矛盾。再令 \\(\\varepsilon\\to0\\)。</p>`,
      `${P.plotBlock('ht-cvXT', 'x–t 平面上的温度（等值线间隔为色阶的 1/12）', '', 'tall')}
       <div class="stats" id="ht-maxStats"></div>
       <p class="small" id="ht-maxNote"></p>
       <p class="small muted exp-only">数值说明：使用全部 60 个模态；t=0 行用精确初值。截断级数在 t→0 时有 Gibbs 振荡，因此内部检查从 t₁（截断项衰减到 e⁻¹⁴ 以下）开始。</p>`)) +
  P.sec(ID, 'meaning', '物理意义', R`<div class="cols">
    <div class="card"><h3>平滑化与不可逆</h3><p class="small muted">任意 \(t\gt0\) 时解无穷光滑；初值的细节（高频）按 \(e^{-Dn^2\pi^2t/L^2}\) 丢失，无法从 \(u(\cdot,t)\) 反推（逆热传导问题不适定）。</p></div>
    <div class="card"><h3>无限传播速度</h3><p class="small muted">局部加热后，远处温度立刻变为正（虽然极小）。这是抛物型方程的数学性质，物理上由于该影响指数小而可接受。</p></div>
    <div class="card"><h3>趋于平衡</h3><p class="small muted">\(\sigma_n\lt0\) 时各瞬态模态衰减，解趋于稳态 \(Du_s''+\alpha u_s+f=0\)。弛豫时间由最慢模态 \(1/|\sigma_1|\) 决定。</p></div>
    <div class="card"><h3>Fourier 系数的含义</h3><p class="small muted">\(v_n(0)\) 是初温中"第 \(n\) 种温度起伏"的含量；高阶起伏空间尺度小，散热快。</p></div></div>` +
    P.pitfalls([
      ['边界非齐次：', '不能直接分离变量；先减稳态（或任一满足边界的函数）。'],
      ['零模态：', String.raw`两端绝热时 \(\lambda_0=0\) 的常数模态不衰减，平均温度守恒；漏掉它会得到"温度趋于 0"的错误结论。`],
      ['反应项与提升：', String.raw`有 \(\alpha u\) 时，\(v=u-w\) 的方程多出 \(\alpha w\) 项。`],
      ['最大值原理的前提：', String.raw`需要 \(u_t-Du_{xx}\le0\)。有正热源或 \(\alpha\gt0\) 时结论可能失效（本页"有热源"按钮给出反例）。`],
      ['时间因子：', '热方程时间一阶 ⇒ 指数衰减；不要照搬波动方程的 cos ωt。']
    ]) + `<div class="bench exp-only"><h3 class="sub">数值基准</h3><div id="ht-bench"></div></div>`) +
  P.sec(ID, 'quiz', '检查你的理解 <span class="score" id="ht-score"></span>', '<div id="ht-quiz"></div>');
}

/* ---------- 初始化 ---------- */
const R_ = {};
const BCO = [{ v: 'D', label: 'D 恒温' }, { v: 'N', label: 'N 绝热' }, { v: 'R', label: 'R 换热' }];
let segL, segR;
function syncRows() {
  $('ht-h1-row').hidden = S.bcL !== 'R'; $('ht-h2-row').hidden = S.bcR !== 'R';
  $('ht-T0-row').hidden = S.bcL !== 'D'; $('ht-T1-row').hidden = S.bcR !== 'D';
  $('ht-m-row').hidden = S.init !== 'sine';
  $('ht-q-row').hidden = S.src === 'none'; $('ht-Om-row').hidden = S.src !== 'periodic';
}
function all() { update(); syncRows(); Qhist = null; MX = null; player.refresh(); }
function applyPreset(k) {
  Object.assign(S, PRE[k].set);
  $('ht-pre').value = k; $('ht-init').value = S.init; $('ht-src').value = S.src;
  segL.set(S.bcL); segR.set(S.bcR);
  ['h1', 'h2', 'T0', 'T1', 'D', 'alpha', 'q', 'Om', 'm'].forEach(key => R_[key].set(S[key]));
  player.stop && player.stop(); player.s = 0;
  all();
}
P.register('heat', {
  init(root) {
    root.innerHTML = html();
    P.stepper($('ht-steps'), STEPS);
    segL = P.segment($('ht-bcL'), BCO, S.bcL, v => { S.bcL = v; all(); });
    segR = P.segment($('ht-bcR'), BCO, S.bcR, v => { S.bcR = v; all(); });
    const f2 = v => v.toFixed(2);
    R_.h1 = P.bindRange('ht-h1', v => { S.h1 = v; all(); }, f2);
    R_.h2 = P.bindRange('ht-h2', v => { S.h2 = v; all(); }, f2);
    R_.T0 = P.bindRange('ht-T0', v => { S.T0 = v; all(); }, f2);
    R_.T1 = P.bindRange('ht-T1', v => { S.T1 = v; all(); }, f2);
    R_.D = P.bindRange('ht-D', v => { S.D = v; all(); }, f2);
    R_.alpha = P.bindRange('ht-alpha', v => { S.alpha = v; all(); }, v => v.toFixed(1));
    R_.q = P.bindRange('ht-q', v => { S.q = v; all(); }, v => v.toFixed(1));
    R_.Om = P.bindRange('ht-Om', v => { S.Om = v; all(); }, v => v.toFixed(0));
    R_.m = P.bindRange('ht-m', v => { S.m = v; all(); }, v => v.toFixed(0));
    R_.N = P.bindRange('ht-N', v => { S.N = v; Qhist = null; player.refresh(); });
    $('ht-pre').addEventListener('change', e => applyPreset(e.target.value));
    $('ht-init').addEventListener('change', e => { S.init = e.target.value; all(); });
    $('ht-src').addEventListener('change', e => { S.src = e.target.value; all(); });
    root.querySelectorAll('[data-pre]').forEach(b => b.addEventListener('click', () => applyPreset(b.dataset.pre)));
    update();
    player = P.player($('ht-player'), { id: 'ht-t', tmax: () => M.tmax, map: s => M.tmax * s * s, digits: 4, duration: 10, onChange: render });
    syncRows();
    P.quiz($('ht-quiz'), QUIZ, $('ht-score'));
    runBench();
    $('ht-bench').addEventListener('click', e => { if (e.target.closest('[data-rerun]')) runBench(); });
  },
  render() { player.refresh(); }
});
})();
