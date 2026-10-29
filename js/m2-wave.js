/* m2-wave.js — 模块 02：一维波动方程（达朗贝尔公式、反射法、驻波与能量、Duhamel 原理） */
(function () {
'use strict';
const P = window.PDE, PI = Math.PI, $ = P.$, fx = P.fx, sci = P.sci;
const ID = 'wv', SQ = Math.sqrt(PI / 2), SIG = 0.5;

/* ================= 初始数据（解析定义，含 ψ 的原函数） ================= */
const PHI = {
  tri: x => Math.max(0, 1 - Math.abs(x)),
  gauss: x => Math.exp(-x * x / (2 * SIG * SIG)),
  sine: x => Math.abs(x) <= 1 ? Math.sin(PI * x) : 0,
  none: () => 0
};
const PSI = { // [ψ, Ψ(x)=∫_0^x ψ]
  none: [() => 0, () => 0],
  box: [x => Math.abs(x) <= 1 ? 1 : 0, x => Math.max(-1, Math.min(1, x))],
  gauss: [x => Math.exp(-x * x / (2 * SIG * SIG)), x => SIG * SQ * P.erf(x / (SIG * Math.SQRT2))],
  tri: [x => Math.max(0, 1 - Math.abs(x)), x => { const s = Math.sign(x), z = Math.min(1, Math.abs(x)); return s * (z - z * z / 2); }]
};

/* ================= A. 达朗贝尔公式 ================= */
const SA = { phi: 'tri', psi: 'none', A: 1, V: 0.6, a: 1, x0: 2.5, t0: 1.5, view: 'lr' };
const XA = 10;
function dal(x, t, s = SA) {
  const phi = PHI[s.phi], Psi = PSI[s.psi][1], a = s.a;
  return 0.5 * s.A * (phi(x - a * t) + phi(x + a * t)) + s.V * (Psi(x + a * t) - Psi(x - a * t)) / (2 * a);
}
const tmaxA = () => 7 / SA.a;
let plA;
function renderA() {
  const t = plA.t(), a = SA.a, phi = PHI[SA.phi], Psi = PSI[SA.psi][1];
  const R = String.raw;
  P.tex($('wv-aProb'), R`\[ \begin{cases} u_{tt}=a^2u_{xx}, & -\infty\lt x\lt\infty,\ t\gt0,\\ u(x,0)=\varphi(x),\ u_t(x,0)=\psi(x) \end{cases}\quad a=${fx(a, 2)},\ \ u=\underbrace{\tfrac12[\varphi(x-at)+\varphi(x+at)]}_{\text{初始位移}}+\underbrace{\tfrac1{2a}\int_{x-at}^{x+at}\psi(s)\,ds}_{\text{初始速度}} \]`);
  const f = P.frame($('wv-cvA'), { xmin: -XA, xmax: XA, ymin: -1.3, ymax: 1.6, xlabel: 'x' });
  if (f) {
    const N = 800, xs = new Float64Array(N + 1); for (let i = 0; i <= N; i++) xs[i] = -XA + 2 * XA * i / N;
    P.line(f, xs, x => SA.A * phi(x), P.C['ink-3'], 1, [3, 3]);
    if (SA.view === 'lr') {
      P.line(f, xs, x => 0.5 * SA.A * phi(x - a * t) - SA.V * Psi(x - a * t) / (2 * a), P.C.warm, 1.5);
      P.line(f, xs, x => 0.5 * SA.A * phi(x + a * t) + SA.V * Psi(x + a * t) / (2 * a), P.C.violet, 1.5);
    } else {
      P.line(f, xs, x => 0.5 * SA.A * (phi(x - a * t) + phi(x + a * t)), P.C.warm, 1.5);
      P.line(f, xs, x => SA.V * (Psi(x + a * t) - Psi(x - a * t)) / (2 * a), P.C.violet, 1.5);
    }
    P.line(f, xs, x => dal(x, t), P.C.accent, 2.4);
  }
  $('wv-aLeg').innerHTML = P.lg('accent', 'u(x,t)') + (SA.view === 'lr' ? P.lg('warm', '右行波 F(x−at)') + P.lg('violet', '左行波 G(x+at)') : P.lg('warm', '位移贡献') + P.lg('violet', '速度贡献')) + P.lg('ink-3', 'φ(x)');
  // x–t 平面
  // 热力图与 t 无关：缓存为静态图层，播放时每帧只重画时间线和选点
  const tm = tmaxA(), g = P.layer($('wv-cvXT'), ['A', SA.phi, SA.psi, SA.A, SA.V, SA.a].join('|'), { xmin: -XA, xmax: XA, ymin: 0, ymax: tm, xlabel: 'x', ylabel: 't', colorbar: true }, g => {
    const nx = 240, nt = 120, d = new Float64Array(nx * nt); let lo = Infinity, hi = -Infinity;
    for (let j = 0; j < nt; j++) for (let i = 0; i < nx; i++) { const v = dal(-XA + 2 * XA * (i + 0.5) / nx, tm * (1 - (j + 0.5) / nt)); d[j * nx + i] = v; lo = Math.min(lo, v); hi = Math.max(hi, v); }
    const m = Math.max(Math.abs(lo), Math.abs(hi), 1e-9);
    P.heatmap(g, d, nx, nt, -m, m);
    const sup = SA.phi === 'gauss' || SA.psi === 'gauss' ? 1.5 : 1;
    [-sup, sup].forEach(x => { P.seg2(g, x, 0, x - a * tm, tm, P.C['ink-2'], 1, [4, 4]); P.seg2(g, x, 0, x + a * tm, tm, P.C['ink-2'], 1, [4, 4]); });
  });
  if (g) {
    P.seg2(g, -XA, t, XA, t, P.C.ink, 1.2);
    const { x0, t0 } = SA;
    P.seg2(g, x0, t0, x0 - a * t0, 0, P.C.warm, 2); P.seg2(g, x0, t0, x0 + a * t0, 0, P.C.violet, 2);
    P.seg2(g, x0 - a * t0, 0, x0 + a * t0, 0, P.C.accent, 5);
    P.dot(g, x0, t0, P.C.ink, 4.5);
  }
  const lo = SA.x0 - a * SA.t0, hi = SA.x0 + a * SA.t0;
  P.stats($('wv-aStats'), [
    ['当前 t', fx(t, 3)],
    ['选点 (x₀, t₀)', `(${fx(SA.x0, 2)}, ${fx(SA.t0, 2)})`],
    ['依赖区间', `[${fx(lo, 2)}, ${fx(hi, 2)}]`],
    ['u(x₀, t₀)', fx(dal(SA.x0, SA.t0), 5)],
    ['位移贡献', fx(0.5 * SA.A * (phi(lo) + phi(hi)), 5)],
    ['速度贡献', fx(SA.V * (Psi(hi) - Psi(lo)) / (2 * a), 5)],
    ['信号最早到达 x₀ 的时间', SA.phi === 'gauss' || SA.psi === 'gauss' ? '高斯无紧支集' : fx(Math.max(0, (Math.abs(SA.x0) - 1) / a), 3)]
  ]);
}

/* ================= B. 半无界弦：反射法 ================= */
const SB = { end: 'D', shape: 'tri', x0: 5, mode: 'left', a: 1 };
const XB = 10;
function bump(x) { const z = x - SB.x0; return SB.shape === 'tri' ? Math.max(0, 1 - Math.abs(z)) : Math.exp(-z * z / (2 * 0.35 * 0.35)); }
// ψ = aφ'（纯左行）时 Ψ(x)=∫_0^x ψ = a[φ(x)−φ(0)]
function PsiB(x) { return SB.mode === 'left' ? SB.a * (bump(x) - bump(0)) : 0; }
function extPhi(x) { return x >= 0 ? bump(x) : (SB.end === 'D' ? -1 : 1) * bump(-x); }
function extPsi(x) { return x >= 0 ? PsiB(x) : (SB.end === 'D' ? 1 : -1) * PsiB(-x); } // φ 奇 ⇒ ψ 奇 ⇒ Ψ 偶
function uB(x, t) { const a = SB.a; return 0.5 * (extPhi(x - a * t) + extPhi(x + a * t)) + (extPsi(x + a * t) - extPsi(x - a * t)) / (2 * a); }
function incB(x, t) {
  const a = SB.a, p0 = y => y >= 0 ? bump(y) : 0, Q0 = y => y >= 0 ? PsiB(y) : 0;
  return 0.5 * (p0(x - a * t) + p0(x + a * t)) + (Q0(x + a * t) - Q0(x - a * t)) / (2 * a);
}
const tmaxB = () => (SB.x0 + 6) / SB.a;
let plB;
function renderB() {
  const t = plB.t(), a = SB.a, R = String.raw;
  P.tex($('wv-bProb'), R`\[ \begin{cases} u_{tt}=a^2u_{xx},\ x\gt0,\\ ${SB.end === 'D' ? 'u(0,t)=0\\ (\\text{固定端})' : 'u_x(0,t)=0\\ (\\text{自由端})'},\\ u(x,0)=\varphi(x),\ u_t(x,0)=${SB.mode === 'left' ? 'a\\varphi\'(x)' : '0'} \end{cases}\quad\Rightarrow\quad \Phi(x)=\begin{cases}\varphi(x),&x\ge0\\ ${SB.end === 'D' ? '-' : ''}\varphi(-x),&x\lt0\end{cases}\ (\text{${SB.end === 'D' ? '奇' : '偶'}延拓}) \]`);
  const f = P.frame($('wv-cvB'), { xmin: -XB, xmax: XB, ymin: -1.4, ymax: 2.2, xlabel: 'x' });
  if (f) {
    f.ctx.fillStyle = P.C.grid; f.ctx.globalAlpha = 0.6; f.ctx.fillRect(f.pad.l, f.pad.t, f.X(0) - f.pad.l, f.ph); f.ctx.globalAlpha = 1;
    P.textAt(f, '虚拟区域（延拓）', -XB + 0.3, 2.05, P.C['ink-3'], 'left', 'top');
    P.textAt(f, '物理区域 x ≥ 0', 0.3, 2.05, P.C['ink-3'], 'left', 'top');
    const n = 800, xl = [], xr = [];
    for (let i = 0; i <= n; i++) { xl.push(-XB * (1 - i / n)); xr.push(XB * i / n); }
    P.line(f, xl, x => uB(x, t), P.C['ink-3'], 1.6, [4, 3]);
    P.line(f, xr, x => incB(x, t), P.C.violet, 1.3);
    P.line(f, xr, x => uB(x, t) - incB(x, t), P.C.warm, 1.3);
    P.line(f, xr, x => uB(x, t), P.C.accent, 2.4);
    P.seg2(f, 0, -1.4, 0, 2.2, P.C.ink, 1.5);
    P.dot(f, 0, uB(0, t), P.C.ink, 4, SB.end === 'N');
  }
  const g = P.layer($('wv-cvBxt'), 'B|' + JSON.stringify(SB), { xmin: 0, xmax: XB, ymin: 0, ymax: tmaxB(), xlabel: 'x', ylabel: 't', colorbar: true }, g => {
    const nx = 200, nt = 120, tm = tmaxB(), d = new Float64Array(nx * nt);
    for (let j = 0; j < nt; j++) for (let i = 0; i < nx; i++) d[j * nx + i] = uB(XB * (i + 0.5) / nx, tm * (1 - (j + 0.5) / nt));
    P.heatmap(g, d, nx, nt, -1.2, 1.2);
    const tr = SB.x0 / a;
    P.seg2(g, SB.x0, 0, 0, tr, P.C.ink, 1.4, [5, 4]); P.seg2(g, 0, tr, XB, tr + XB / a, P.C.ink, 1.4, [5, 4]);
    if (SB.mode === 'rest') P.seg2(g, SB.x0, 0, XB, (XB - SB.x0) / a, P.C.ink, 1.4, [5, 4]);
  });
  if (g) {
    P.seg2(g, 0, t, XB, t, P.C.ink, 1.2);
  }
  const h = 1e-5, ux0 = (uB(h, t) - uB(-h, t)) / (2 * h);
  P.stats($('wv-bStats'), [
    ['当前 t', fx(t, 3)],
    ['端点 u(0,t)', sci(uB(0, t))],
    ['端点 u_x(0,t)（差分）', sci(ux0)],
    ['脉冲中心到达端点', 't = ' + fx(SB.x0 / a, 3)],
    ['反射波', SB.end === 'D' ? '<span class="pill warn">倒相（相位 π）</span>' : '<span class="pill ok">同相</span>'],
    ['端点附近叠加', SB.end === 'D' ? '入射 + 反射抵消，u(0)=0' : '斜率抵消，端点振幅加倍']
  ]);
}

/* ================= C. 有界弦：驻波、模态叠加与能量 ================= */
const SC = { bcL: 'D', bcR: 'D', h: 2, L: 1, a: 1, phi: 'tri', p: 0.3, psi: 'none', xh: 0.7, N: 30, sel: 0 };
let slC, AC, BC, plC, energy = null;
function phiC(xi) {
  switch (SC.phi) {
    case 'tri': return xi < SC.p ? xi / SC.p : (1 - xi) / (1 - SC.p);
    case 'gauss': return Math.exp(-((xi - 0.5) ** 2) / (2 * 0.06 * 0.06));
    case 'mode': return null;
    default: return 0;
  }
}
function psiC(xi) {
  switch (SC.psi) {
    case 'box': return Math.abs(xi - SC.xh) <= 0.05 ? 2 : 0;
    case 'gauss': return 2 * Math.exp(-((xi - SC.xh) ** 2) / (2 * 0.04 * 0.04));
    case 'uniform': return 0.5;
    default: return 0;
  }
}
function updateC() {
  slC = P.solveSL({ L: SC.L, bcL: SC.bcL, bcR: SC.bcR, h1: SC.h, h2: SC.h, NM: 60, NG: 1200 });
  const NG = slC.NG, f = new Float64Array(NG + 1), g = new Float64Array(NG + 1);
  if (SC.phi === 'mode') { const X = slC.modes[Math.min(2, slC.modes.length - 1)].X, m = P.maxAbs(X); for (let i = 0; i <= NG; i++) f[i] = 0.8 * X[i] / m; }
  else for (let i = 0; i <= NG; i++) f[i] = 0.8 * phiC(i / NG);
  for (let i = 0; i <= NG; i++) g[i] = psiC(i / NG) * SC.a;
  AC = P.project(slC, f); const gp = P.project(slC, g);
  BC = Float64Array.from(slC.modes, (m, j) => m.kind === 'pos' ? gp[j] / (SC.a * m.k) : gp[j]);
  slC.f = f; slC.g = g;
  computeEnergy();
}
const omega = m => SC.a * m.k;
function uC(t, deriv) { // deriv: 0=u, 1=u_t, 2=u_x
  const NG = slC.NG, out = new Float64Array(NG + 1), n = Math.min(SC.N, slC.modes.length);
  for (let j = 0; j < n; j++) {
    const m = slC.modes[j], w = omega(m);
    let T;
    if (m.kind === 'pos') T = deriv === 1 ? w * (-AC[j] * Math.sin(w * t) + BC[j] * Math.cos(w * t)) : AC[j] * Math.cos(w * t) + BC[j] * Math.sin(w * t);
    else T = deriv === 1 ? BC[j] : AC[j] + BC[j] * t;
    if (T === 0) continue;
    const X = deriv === 2 ? m.dX : m.X;
    for (let i = 0; i <= NG; i++) out[i] += T * X[i];
  }
  return out;
}
const tmaxC = () => 2 * 2 * PI / omega(slC.modes.find(m => m.kind === 'pos'));
function energyAt(t) {
  const ut = uC(t, 1), ux = uC(t, 2), u = SC.bcL === 'R' || SC.bcR === 'R' ? uC(t, 0) : null, h = slC.h, a2 = SC.a * SC.a;
  const K = 0.5 * P.dotS(ut, ut, h);
  let V = 0.5 * a2 * P.dotS(ux, ux, h);
  if (u) V += 0.5 * a2 * ((SC.bcL === 'R' ? SC.h * u[0] * u[0] : 0) + (SC.bcR === 'R' ? SC.h * u[slC.NG] ** 2 : 0));
  return [K, V];
}
function computeEnergy() {
  const tm = tmaxC(), n = 120, ts = [], K = [], V = [], E = [];
  for (let i = 0; i <= n; i++) { const t = tm * i / n, [k, v] = energyAt(t); ts.push(t); K.push(k); V.push(v); E.push(k + v); }
  let modal = 0; const n0 = Math.min(SC.N, slC.modes.length);
  for (let j = 0; j < n0; j++) { const m = slC.modes[j], w = omega(m); modal += m.kind === 'pos' ? 0.5 * m.n2 * w * w * (AC[j] ** 2 + BC[j] ** 2) : 0.5 * m.n2 * BC[j] ** 2; }
  const drift = Math.max(...E.map(e => Math.abs(e - E[0]))) / (E[0] || 1);
  energy = { ts, K, V, E, modal, drift };
}
function renderC() {
  const t = plC.t(), R = String.raw, a = SC.a;
  P.tex($('wv-cProb'), R`\[ \begin{cases} u_{tt}=a^2u_{xx},\ 0\lt x\lt L,\\ ${P.BCTEX.uL[SC.bcL]},\ \ ${P.BCTEX.uR[SC.bcR]},\\ u(x,0)=\varphi(x),\ u_t(x,0)=\psi(x)\end{cases}\qquad u=\sum_n\big(A_n\cos\omega_nt+B_n\sin\omega_nt\big)X_n(x),\ \ \omega_n=ak_n,\ \ A_n=\frac{\langle\varphi,X_n\rangle}{\|X_n\|^2},\ B_n=\frac{\langle\psi,X_n\rangle}{\omega_n\|X_n\|^2} \]`);
  const m = slC.modes[SC.sel], w = omega(m), X = m.X, mx = P.maxAbs(X);
  const f = P.frame($('wv-cvMode'), { xmin: 0, xmax: SC.L, ymin: -1.2, ymax: 1.2, xlabel: 'x' });
  if (f) {
    P.fillBetween(f, slC.xs, X.map(v => -Math.abs(v) / mx), X.map(v => Math.abs(v) / mx), P.C.accent, 0.12);
    P.line(f, slC.xs, X.map(v => v / mx * (m.kind === 'pos' ? Math.cos(w * t) : 1)), P.C.accent, 2.4, null, 2);
    P.bcGlyphs(f, SC.bcL, SC.bcR);
    const NG = slC.NG;
    for (let i = 0; i <= NG; i++) {
      if (i > 0 && X[i - 1] * X[i] < 0 || Math.abs(X[i] / mx) < 1e-9) P.dot(f, slC.xs[i], 0, P.C.ink, 3.5);
      if (i > 0 && i < NG && Math.abs(X[i]) >= Math.abs(X[i - 1]) && Math.abs(X[i]) > Math.abs(X[i + 1]) || ((i === 0 || i === NG) && Math.abs(X[i] / mx) > 0.999)) { P.dot(f, slC.xs[i], X[i] / mx, P.C.warm, 4, true); }
    }
  }
  const pos = slC.modes.filter(q => q.kind === 'pos'), w1 = omega(pos[0]);
  const ratios = pos.slice(1, 4).map(q => omega(q) / w1), harm = ratios.every(r => Math.abs(r - Math.round(r)) < 1e-6);
  P.stats($('wv-cModeStats'), [
    ['模态 n', SC.sel + slC.offset],
    ['k_n', fx(m.k, 5)],
    ['ω_n = a k_n', fx(w, 5)],
    ['频率 f_n = ω_n/2π', fx(w / (2 * PI), 5)],
    ['波长 λ_n = 2π/k_n', m.k ? fx(2 * PI / m.k, 4) : '∞'],
    ['f_n λ_n（应等于 a）', m.k ? fx(w / (2 * PI) * 2 * PI / m.k, 6) : '—']
  ]);
  const u = uC(t, 0), g = P.frame($('wv-cvC'), { xmin: 0, xmax: SC.L, ymin: -1.2, ymax: 1.2, xlabel: 'x' });
  if (g) { P.line(g, slC.xs, slC.f, P.C.warm, 1.2, [4, 4], 2); P.line(g, slC.xs, u, P.C.accent, 2.2); P.bcGlyphs(g, SC.bcL, SC.bcR); }
  const nb = Math.min(20, slC.modes.length), cm = Math.max(P.maxAbs(AC.slice(0, nb)), P.maxAbs(BC.slice(0, nb)));
  const b = P.frame($('wv-cvCoef'), { xmin: slC.offset - 0.8, xmax: nb + slC.offset - 0.2, ymin: -cm * 1.1, ymax: cm * 1.1, xlabel: 'n', xticks: P.niceTicks(slC.offset, nb - 1 + slC.offset, 5).filter(v => v % 1 === 0).map(v => ({ v, label: v + '' })) });
  if (b) {
    const bw = Math.max(2, b.pw / nb * 0.3), xsN = Array.from({ length: nb }, (_, j) => j + slC.offset);
    P.bars(b, xsN.map(x => x - 0.17), Array.from(AC.slice(0, nb)), j => j < SC.N ? P.C.accent : P.C['ink-3'], bw);
    P.bars(b, xsN.map(x => x + 0.17), Array.from(BC.slice(0, nb)), j => j < SC.N ? P.C.warm : P.C['ink-3'], bw);
  }
  const e = P.frame($('wv-cvE'), { xmin: 0, xmax: energy.ts[energy.ts.length - 1], ymin: 0, ymax: Math.max(...energy.E) * 1.15 || 1, xlabel: 't' });
  if (e) {
    P.line(e, energy.ts, energy.K, P.C.warm, 1.6); P.line(e, energy.ts, energy.V, P.C.violet, 1.6); P.line(e, energy.ts, energy.E, P.C.accent, 2.4);
    P.seg2(e, t, 0, t, e.o.ymax, P.C.ink, 1, [3, 3]);
  }
  const [k, v] = energyAt(t);
  P.stats($('wv-cStats'), [
    ['当前 t', fx(t, 3)],
    ['周期性', harm ? `<span class="pill ok">ω_n/ω₁ 为整数，周期 ${fx(2 * PI / w1, 4)}</span>` : `<span class="pill warn">ω₂/ω₁ = ${fx(ratios[0], 4)}，非周期</span>`],
    ['动能 K(t)', fx(k, 5)],
    ['势能 V(t)', fx(v, 5)],
    ['总能量 E = K + V', fx(k + v, 6)],
    ['E 的最大相对漂移', sci(energy.drift)],
    ['模态能量和 ½Σ‖X_n‖²ω_n²(A_n²+B_n²)', fx(energy.modal, 6)]
  ]);
}

/* ================= D. Duhamel 原理 ================= */
const SD = { force: 'harm', prof: 'local', xf: 0.3, Om: 0.9, tau0: 0.3, M: 8, a: 1, N: 40 };
let slD, gD, plD, histD = null;
function updateD() {
  slD = slD || P.solveSL({ L: 1, bcL: 'D', bcR: 'D', NM: 60, NG: 1200 });
  const NG = slD.NG, g = new Float64Array(NG + 1);
  for (let i = 0; i <= NG; i++) { const x = i / NG; g[i] = SD.prof === 'uniform' ? 1 : 4 * Math.exp(-((x - SD.xf) ** 2) / (2 * 0.03 * 0.03)); }
  gD = P.project(slD, g); slD.g = g;
  histD = null;
}
const wD = j => SD.a * slD.modes[j].k, OmD = () => SD.Om * wD(0);
function forcedHarm(w, Om, t) { // ü + ω²u = sin Ωt, u(0)=u'(0)=0
  if (Math.abs(w - Om) < 1e-6 * w) return (Math.sin(w * t) - w * t * Math.cos(w * t)) / (2 * w * w);
  return (Math.sin(Om * t) - (Om / w) * Math.sin(w * t)) / (w * w - Om * Om);
}
function unD(j, t) {
  const w = wD(j), c = gD[j];
  if (SD.force === 'const') return c * (1 - Math.cos(w * t)) / (w * w);
  if (SD.force === 'harm') return c * forcedHarm(w, OmD(), t);
  return t > SD.tau0 ? c * Math.sin(w * (t - SD.tau0)) / w : 0;
}
const fnD = (j, tau) => SD.force === 'const' ? gD[j] : gD[j] * Math.sin(OmD() * tau);
function sumModes(coefFn) { const NG = slD.NG, out = new Float64Array(NG + 1); for (let j = 0; j < SD.N; j++) { const c = coefFn(j); if (!c) continue; const X = slD.modes[j].X; for (let i = 0; i <= NG; i++) out[i] += c * X[i]; } return out; }
const tmaxD = () => (SD.force === 'harm' ? 8 : 3) * 2 * PI / wD(0);
const probeD = () => SD.prof === 'uniform' ? 0.5 : SD.xf;
function renderD() {
  const t = plD.t(), R = String.raw;
  const ftex = SD.force === 'const' ? R`F\,g(x)\,H(t)` : SD.force === 'harm' ? R`F\,g(x)\sin\Omega t` : R`I\,g(x)\,\delta(t-\tau_0)`;
  P.tex($('wv-dProb'), R`\[ \begin{cases} u_{tt}-a^2u_{xx}=f(x,t)=${ftex},\\ u(0,t)=u(L,t)=0,\ \ u(x,0)=u_t(x,0)=0\end{cases}\qquad u(x,t)=\int_0^t w(x,t;\tau)\,d\tau \]`);
  const u = sumModes(j => unD(j, t));
  let slices = [], approx = null;
  if (SD.force !== 'impulse' && t > 0) {
    const M = SD.M, dt = t / M; approx = new Float64Array(slD.NG + 1);
    for (let k = 0; k < M; k++) {
      const tau = (k + 0.5) * dt, w = sumModes(j => dt * fnD(j, tau) * Math.sin(wD(j) * (t - tau)) / wD(j));
      for (let i = 0; i < w.length; i++) approx[i] += w[i];
      if (k % Math.max(1, Math.ceil(M / 8)) === 0) slices.push(w);
    }
  }
  const static_ = SD.force === 'const' ? sumModes(j => gD[j] / wD(j) ** 2) : null;
  let ym = Math.max(P.maxAbs(u), static_ ? 2 * P.maxAbs(static_) : 0, 1e-3);
  if (!histD) {
    const tm = tmaxD(), n = 400, xp = probeD(), Xp = slD.modes.map(m => { const i = Math.round(xp * slD.NG); return m.X[i]; });
    histD = { ts: [], us: [] };
    for (let i = 0; i <= n; i++) { const tt = tm * i / n; let s = 0; for (let j = 0; j < SD.N; j++) s += unD(j, tt) * Xp[j]; histD.ts.push(tt); histD.us.push(s); }
  }
  ym = Math.max(ym, P.maxAbs(histD.us) * 0.6);
  const f = P.frame($('wv-cvD'), { xmin: 0, xmax: 1, ymin: -ym * 1.15, ymax: ym * 1.15, xlabel: 'x' });
  if (f) {
    slices.forEach(w => P.line(f, slD.xs, w, P.C['ink-3'], 1, null, 4));
    if (static_) P.line(f, slD.xs, static_, P.C.warm, 1.4, [5, 4], 4);
    if (approx) P.line(f, slD.xs, approx, P.C.violet, 1.6, [3, 3], 2);
    P.line(f, slD.xs, u, P.C.accent, 2.4, null, 2);
    if (SD.prof === 'local') P.seg2(f, SD.xf, -ym * 1.15, SD.xf, ym * 1.15, P.C.warm, 1, [2, 4]);
    P.bcGlyphs(f, 'D', 'D');
  }
  const hh = P.frame($('wv-cvDh'), { xmin: 0, xmax: histD.ts[histD.ts.length - 1], ymin: -P.maxAbs(histD.us) * 1.1, ymax: P.maxAbs(histD.us) * 1.1, xlabel: 't' });
  if (hh) { P.line(hh, histD.ts, histD.us, P.C.accent, 1.6); P.seg2(hh, t, hh.o.ymin, t, hh.o.ymax, P.C.ink, 1, [3, 3]); }
  let err = NaN; if (approx) { err = 0; for (let i = 0; i < u.length; i++) err = Math.max(err, Math.abs(approx[i] - u[i])); err /= (P.maxAbs(u) || 1); }
  let nearest = 0; for (let j = 0; j < 8; j++) if (Math.abs(wD(j) - OmD()) < Math.abs(wD(nearest) - OmD())) nearest = j;
  const reson = SD.force === 'harm' && Math.abs(wD(nearest) - OmD()) < 0.02 * wD(nearest) && Math.abs(gD[nearest]) > 1e-8;
  P.stats($('wv-dStats'), [
    ['当前 t', fx(t, 3)],
    ['Ω / ω₁', SD.force === 'harm' ? fx(SD.Om, 3) : '—'],
    ['最近的固有频率', SD.force === 'harm' ? `ω_${nearest + 1} = ${fx(wD(nearest), 3)}` : '—'],
    ['共振状态', SD.force !== 'harm' ? '—' : reson ? '<span class="pill warn">共振：振幅随 t 线性增长</span>' : Math.abs(gD[nearest]) < 1e-8 ? '<span class="pill info">该模态未被激发（g_n=0）</span>' : '<span class="pill ok">非共振：振幅有界</span>'],
    ['黎曼和相对误差（M 片）', isFinite(err) ? sci(err) : '—'],
    [`探针 x=${fx(probeD(), 2)} 处 u`, fx(histD.us[Math.round(t / histD.ts[histD.ts.length - 1] * 400)] ?? 0, 5)]
  ]);
}

/* ================= 基准 ================= */
function runBench() {
  const rows = [], s = { phi: 'gauss', psi: 'gauss', A: 1, V: 0.7, a: 1.3 };
  let e = 0; for (let x = -6; x <= 6; x += 0.05) e = Math.max(e, Math.abs(dal(x, 0, s) - PHI.gauss(x)));
  rows.push(P.benchRow('达朗贝尔：初始位移 u(x,0)=φ', '直接代入 t=0', e, 1e-15));
  e = 0; const h = 1e-5; for (let x = -4; x <= 4; x += 0.1) e = Math.max(e, Math.abs((dal(x, h, s) - dal(x, -h, s)) / (2 * h) - 0.7 * PSI.gauss[0](x)));
  rows.push(P.benchRow('达朗贝尔：初始速度 u_t(x,0)=ψ', '中心差分 h=10⁻⁵', e, 1e-8));
  e = 0; const hd = 2e-3;
  for (let x = -3; x <= 3; x += 0.25) for (let t = 0.3; t <= 2; t += 0.4) {
    const utt = (dal(x, t + hd, s) - 2 * dal(x, t, s) + dal(x, t - hd, s)) / hd ** 2, uxx = (dal(x + hd, t, s) - 2 * dal(x, t, s) + dal(x - hd, t, s)) / hd ** 2;
    e = Math.max(e, Math.abs(utt - 1.69 * uxx));
  }
  rows.push(P.benchRow('达朗贝尔：方程残差 u_tt − a²u_xx', '二阶差分 h=2×10⁻³（截断误差 O(h²)）', e, 1e-4));
  const sb = { phi: 'tri', psi: 'box', A: 1, V: 1, a: 1 }; e = 0;
  for (let x = 3.0001; x < 10; x += 0.1) e = Math.max(e, Math.abs(dal(x, 2, sb)), Math.abs(dal(-x, 2, sb)));
  rows.push(P.benchRow('有限传播速度：t=2 时 |x|>1+at 处 u≡0', '三角形 φ + 方波 ψ', e, 0));
  rows.push(P.benchRow('误差函数 erf(1)、erf(3)（用于高斯 ψ 的原函数）', '参考值 0.8427007929497149, 0.9999779095030014', Math.abs(P.erf(1) - 0.8427007929497149) + Math.abs(P.erf(3) - 0.9999779095030014), 1e-14));
  // 无界弦能量守恒（高斯数据，[−30,30] 上积分）
  const En = t => { const N = 6000, hx = 60 / N, ut = new Float64Array(N + 1), ux = new Float64Array(N + 1), d = 1e-5;
    for (let i = 0; i <= N; i++) { const x = -30 + i * hx; ut[i] = (dal(x, t + d, s) - dal(x, t - d, s)) / (2 * d); ux[i] = (dal(x + d, t, s) - dal(x - d, t, s)) / (2 * d); }
    return 0.5 * (P.dotS(ut, ut, hx) + 1.69 * P.dotS(ux, ux, hx)); };
  const E0 = En(0), E3 = En(3);
  rows.push(P.benchRow('无界弦能量守恒 E(3) vs E(0)', '差分 + Simpson 积分', Math.abs(E3 - E0) / E0, 1e-7));
  // 反射法
  const keep = { ...SB }; Object.assign(SB, { end: 'D', shape: 'gauss', x0: 5, mode: 'left', a: 1 });
  e = 0; for (let t = 0; t <= 12; t += 0.05) e = Math.max(e, Math.abs(uB(0, t)));
  rows.push(P.benchRow('固定端：u(0,t) ≡ 0（奇延拓）', '直接计算', e, 1e-15));
  e = 0; for (let x = 0; x <= 10; x += 0.05) e = Math.max(e, Math.abs(uB(x, 10) + bump(10 - x)));
  rows.push(P.benchRow('固定端：t=2x₀/a 后 u = −φ(2x₀−x)（倒相）', '解析反射公式', e, 1e-12));
  SB.end = 'N'; e = 0; for (let t = 0; t <= 12; t += 0.05) e = Math.max(e, Math.abs((uB(1e-5, t) - uB(-1e-5, t)) / 2e-5));
  rows.push(P.benchRow('自由端：u_x(0,t) ≡ 0（偶延拓）', '中心差分', e, 1e-6));
  Object.assign(SB, keep);
  // 有界弦：级数 vs 周期奇延拓的达朗贝尔解
  const keepC = { ...SC }; Object.assign(SC, { bcL: 'D', bcR: 'D', L: 1, a: 1, phi: 'tri', p: 0.3, psi: 'none' });
  updateC();
  const ext = x => { let y = ((x % 2) + 2) % 2; const sgn = y > 1 ? -1 : 1; if (y > 1) y = 2 - y; return sgn * 0.8 * (y < 0.3 ? y / 0.3 : (1 - y) / 0.7); };
  const errN = N => { SC.N = N; const u = uC(0.37, 0); let m = 0; for (let i = 0; i <= slC.NG; i++) { const x = slC.xs[i]; m = Math.max(m, Math.abs(u[i] - 0.5 * (ext(x - 0.37) + ext(x + 0.37)))); } return m; };
  const e15 = errN(15), e30 = errN(30), e60 = errN(60);
  rows.push({ name: '有界弦 D–D：级数解 vs 达朗贝尔周期延拓', method: `N=15/30/60 误差 ${sci(e15)} / ${sci(e30)} / ${sci(e60)}`, err: e60, tol: 5e-3, pass: e60 < e30 && e30 < e15 && e60 < 5e-3 });
  let we = 0; slC.modes.slice(0, 20).forEach((m, j) => { we = Math.max(we, Math.abs(m.k - (j + 1) * PI) / ((j + 1) * PI)); });
  rows.push(P.benchRow('有界弦 D–D：ω_n = nπa/L（前 20 个）', '闭式解', we, 1e-12));
  Object.assign(SC, { psi: 'gauss', N: 40 }); updateC();
  rows.push(P.benchRow('能量守恒：E(t) 的最大相对漂移（位移+速度激励）', '直接对 u_t、u_x 积分', energy.drift, 1e-8));
  rows.push(P.benchRow('能量 = 模态能量之和 ½Σ‖X_n‖²ω_n²(A_n²+B_n²)', '正交性推出的恒等式', Math.abs(energy.E[0] - energy.modal) / energy.modal, 1e-8));
  Object.assign(SC, keepC); updateC();
  // Duhamel
  const keepD = { ...SD }; Object.assign(SD, { prof: 'uniform', a: 1, N: 60 }); updateD();
  const st = sumModes(j => gD[j] / wD(j) ** 2); e = 0; for (let i = 0; i <= slD.NG; i++) { const x = slD.xs[i]; e = Math.max(e, Math.abs(st[i] - x * (1 - x) / 2)); }
  let tail = 0; for (let n = 61; n < 400000; n += 2) tail += 4 / (PI ** 3 * n ** 3);
  rows.push(P.benchRow('恒力作用的平均位置 = 静挠度 x(L−x)/(2a²)', `解析静力解；容差 = 截断尾项上界 Σ_{n>60}|c_n|`, e, tail));
  const quad = (w, Om, t) => { const n = 4000, hq = t / n, v = new Float64Array(n + 1); for (let i = 0; i <= n; i++) { const tau = i * hq; v[i] = Math.sin(Om * tau) * Math.sin(w * (t - tau)) / w; } return P.integ(v, hq); };
  const w1 = PI;
  rows.push(P.benchRow('Duhamel 闭式（非共振 Ω=0.9ω₁）vs 直接积分', 'Simpson 4000 段', Math.abs(forcedHarm(w1, 0.9 * w1, 3.3) - quad(w1, 0.9 * w1, 3.3)) / Math.abs(quad(w1, 0.9 * w1, 3.3)), 1e-10));
  rows.push(P.benchRow('Duhamel 闭式（共振 Ω=ω₁）vs 直接积分', 'Simpson 4000 段', Math.abs(forcedHarm(w1, w1, 3.3) - quad(w1, w1, 3.3)) / Math.abs(quad(w1, w1, 3.3)), 1e-10));
  Object.assign(SD, { force: 'harm', Om: 0.9, prof: 'local' }); updateD();
  const rerr = M => { const t = 2.5, dt = t / M, ap = new Float64Array(slD.NG + 1), ex = sumModes(j => unD(j, t));
    for (let k = 0; k < M; k++) { const tau = (k + 0.5) * dt, w = sumModes(j => dt * fnD(j, tau) * Math.sin(wD(j) * (t - tau)) / wD(j)); for (let i = 0; i < w.length; i++) ap[i] += w[i]; }
    let m = 0; for (let i = 0; i < ap.length; i++) m = Math.max(m, Math.abs(ap[i] - ex[i])); return m; };
  const r1 = rerr(20) / rerr(40);
  rows.push({ name: 'Duhamel 黎曼和（中点）收敛阶', method: `误差比 e(20)/e(40) = ${fx(r1, 2)}，二阶方法应约为 4`, err: Math.abs(r1 - 4), tol: 0.6, pass: Math.abs(r1 - 4) < 0.6 });
  Object.assign(SD, keepD); updateD();
  P.benchTable($('wv-bench'), 'wave', rows);
}

/* ================= 文字内容 ================= */
const QUIZ = [
  { type: '预测演化', q: String.raw`无界弦 \(\varphi\) 为高 1 的三角形，\(\psi=0\)。足够长时间后看到什么？`, o: ['一个高 1 的三角形向右移动', '两个高 1/2 的三角形分别向左右移动', '两个高 1 的三角形', '三角形原地衰减'], a: 1, e: String.raw`\(u=\tfrac12[\varphi(x-at)+\varphi(x+at)]\)，初始形状一分为二，各带一半高度。实验 A 可直接验证。` },
  { type: '依赖区间', q: String.raw`点 \((x_0,t_0)\) 处的解只依赖于哪段初始数据？`, o: ['整条直线', String.raw`\([x_0-at_0,\ x_0+at_0]\)`, String.raw`\([0,\ x_0]\)`, String.raw`只依赖 \(x_0\) 一点`], a: 1, e: String.raw`由达朗贝尔公式，\(u(x_0,t_0)\) 只用到 \(\varphi(x_0\pm at_0)\) 和 \(\psi\) 在 \([x_0-at_0,x_0+at_0]\) 上的积分。这是有限传播速度的数学表达。` },
  { type: '边界与延拓', q: '半无界弦的固定端 u(0,t)=0，应对初始数据作什么延拓？反射波怎样？', o: ['偶延拓，同相反射', '奇延拓，倒相反射', '周期延拓，不反射', '奇延拓，同相反射'], a: 1, e: String.raw`奇延拓保证 \(\Phi(-x)=-\Phi(x)\)，于是 \(u(0,t)=0\) 自动成立；镜像脉冲符号相反，反射波倒相。自由端 \(u_x=0\) 用偶延拓，同相反射。` },
  { type: '识别推导错误', q: String.raw`某同学写出 \(u=\varphi(x-at)+\varphi(x+at)\)（\(\psi=0\)）。错在哪里？`, o: ['应为减号', '缺少系数 1/2：t=0 时会得到 2φ', '应写成 φ(x−t)', '没有错误'], a: 1, e: String.raw`代入 \(t=0\) 检验初始条件：\(u(x,0)=2\varphi(x)\neq\varphi(x)\)。检验初值是发现这类错误最快的方法。` },
  { type: '预测演化', q: String.raw`\(\varphi=0\)，\(\psi\) 是 \([-1,1]\) 上的常数 \(V\)。很长时间后 \(x=0\) 附近 \(u\) 等于？`, o: ['0', 'V', 'V/a', '随时间线性增长'], a: 2, e: String.raw`\(u=\tfrac1{2a}\int_{x-at}^{x+at}\psi\). 当 \([-1,1]\subset[x-at,x+at]\) 时积分为 \(2V\)，\(u=V/a\)：一个高为 \(V/a\) 的平台向两侧扩展。` },
  { type: 'Duhamel 原理', q: String.raw`Duhamel 原理中，时刻 \(\tau\) 的外力冲量 \(f(x,\tau)\,d\tau\) 产生的"子问题"初始条件是？`, o: ['位移为 f dτ，速度为 0', '位移为 0，速度为 f dτ', '位移、速度都为 f dτ', '边界条件变为 f'], a: 1, e: String.raw`冲量改变的是动量（速度），不能瞬间改变位置。所以 \(w(x,\tau;\tau)=0\)，\(w_t(x,\tau;\tau)=f(x,\tau)\)，\(u=\int_0^tw\,d\tau\)。` },
  { type: '驻波周期', q: '长 L 的弦一端固定、一端自由（D–N），自由振动的周期是？', o: ['L/a', '2L/a', '4L/a', '不是周期运动'], a: 2, e: String.raw`\(k_n=(n-\tfrac12)\pi/L\)，\(\omega_n=(2n-1)\omega_1\) 都是 \(\omega_1=\pi a/(2L)\) 的奇数倍，周期 \(2\pi/\omega_1=4L/a\)。` }
];
const STEPS = [
  { t: '引入特征坐标', b: String.raw`<p>波动方程是双曲型，两族特征线 \(x\pm at=c\)。令</p><div class="formula-box">\[ \xi=x-at,\quad \eta=x+at\;\Longrightarrow\; u_{tt}-a^2u_{xx}=-4a^2u_{\xi\eta}=0 . \]</div>` },
  { t: '积分得通解', b: String.raw`<div class="formula-box">\[ u_{\xi\eta}=0\;\Rightarrow\;u=F(\xi)+G(\eta)=F(x-at)+G(x+at). \]</div><p>\(F(x-at)\) 形状不变地以速度 \(a\) 向右平移，\(G(x+at)\) 向左平移。</p>` },
  { t: '代入初始条件', b: String.raw`<div class="formula-box">\[ F(x)+G(x)=\varphi(x),\qquad -aF'(x)+aG'(x)=\psi(x)\;\Rightarrow\;G(x)-F(x)=\frac1a\int_0^x\psi(s)\,ds+C . \]</div>` },
  { t: '解出 F、G 并合成', b: String.raw`<div class="formula-box">\[ u(x,t)=\frac{\varphi(x-at)+\varphi(x+at)}2+\frac1{2a}\int_{x-at}^{x+at}\psi(s)\,ds \]</div><p>常数 \(C\) 在合成时抵消。这是达朗贝尔公式，适用于无界弦、\(\varphi\in C^2,\ \psi\in C^1\)（对分段光滑数据给出弱解）。</p>` },
  { t: '例题：三角形初始位移、静止释放', b: String.raw`<p>\(\varphi(x)=\max(0,1-|x|)\)，\(\psi=0\)：</p><div class="formula-box">\[ u=\tfrac12\varphi(x-at)+\tfrac12\varphi(x+at). \]</div><p>\(t\gt1/a\) 后两个半高三角形分开，各自以速度 \(a\) 远离。尖角不会被抹平（与热方程对比）。</p>` },
  { t: '依赖区间与影响区域', b: String.raw`<p>\(u(x_0,t_0)\) 只依赖初始数据在 \([x_0-at_0,x_0+at_0]\) 上的值；反过来，初始点 \(x_1\) 只影响 \(|x-x_1|\le at\) 的区域。</p>` },
  { t: '推广：半无界弦与有界弦', b: String.raw`<p>固定端 \(u(0,t)=0\)：把 \(\varphi,\psi\) 奇延拓到 \(x\lt0\) 再用达朗贝尔公式；自由端 \(u_x(0,t)=0\)：偶延拓。两端固定的有界弦：奇延拓并以 \(2L\) 为周期延拓，它与分离变量的级数解相同（数值基准中有对照）。</p>` }
];
function html() {
  const R = String.raw;
  return P.modHead(ID, '02', '一维波动方程', R`从达朗贝尔公式出发理解行波、特征线和有限传播速度；用反射法处理半无界弦；用分离变量法处理有界弦的驻波和能量；用 Duhamel 原理处理受迫振动。`) +
  P.sec(ID, 'phys', '弦为什么满足波动方程', R`<div class="cols">
    <div class="card"><h3>受力分析</h3><p class="small muted">取弦元 \([x,x+dx]\)，张力 \(T\) 沿切线。小振幅下竖直合力约为 \(T[u_x(x+dx)-u_x(x)]\approx Tu_{xx}dx\)，由牛顿第二定律 \(\rho\,dx\,u_{tt}=Tu_{xx}\,dx+f\,dx\)。</p></div>
    <div class="card"><h3>适用条件</h3><p class="small muted">小振幅（\(|u_x|\ll1\)）、张力近似不变、忽略重力和阻尼、弦完全柔软。超出这些条件，方程会变成非线性或带阻尼项。</p></div>
    <div class="card"><h3>其他例子</h3><p class="small muted">细杆的纵振动（\(a^2=E/\rho\)）、管中声波（\(a\) 为声速）、传输线电压。</p></div></div>`) +
  P.sec(ID, 'model', '数学模型', R`<div class="formula-box">\[ u_{tt}=a^2u_{xx}+F(x,t),\qquad a=\sqrt{T/\rho},\quad F=f/\rho . \]</div><p class="prose small muted">\(F=0\) 为齐次方程（自由振动），\(F\neq0\) 为受迫振动。方程线性，解可以叠加，这是 Duhamel 原理的基础。</p>`) +
  P.sec(ID, 'cond', '定解条件', R`<div class="tbl-wrap"><table><thead><tr><th>问题</th><th>区域</th><th>需要的条件</th><th>方法</th></tr></thead><tbody>
    <tr><td>初值问题（Cauchy 问题）</td><td>\(-\infty\lt x\lt\infty\)</td><td>\(u(x,0)=\varphi,\ u_t(x,0)=\psi\)</td><td>达朗贝尔公式</td></tr>
    <tr><td>半无界问题</td><td>\(x\gt0\)</td><td>初值 + \(x=0\) 处一个边界条件</td><td>奇 / 偶延拓（反射法）</td></tr>
    <tr><td>初边值问题</td><td>\(0\lt x\lt L\)</td><td>初值 + 两端边界条件</td><td>分离变量（模块 04 的特征函数）</td></tr>
    <tr><td>受迫振动</td><td>同上</td><td>零初值 + 源项 \(f\)</td><td>Duhamel 原理 / 特征函数展开</td></tr></tbody></table></div>
    <p class="prose small muted">二阶时间导数需要两个初始条件：位置和速度。只给位置会使问题不定。</p>`) +
  P.sec(ID, 'idea', '求解思路', R`<div class="prose small muted"><p>沿特征线 \(x-at=c\) 和 \(x+at=c\) 换元，方程变成 \(u_{\xi\eta}=0\)，可以直接积分两次。初始条件决定两个任意函数。边界的影响用"镜像"处理：在虚拟区域放一个对称或反对称的镜像波，使边界条件自动满足。</p></div>`) +
  P.sec(ID, 'deriv', '推导达朗贝尔公式', '<div id="wv-steps"></div>') +
  P.sec(ID, 'viz', '交互实验',
    `<h3 class="sub">A　达朗贝尔公式：行波分解与特征线</h3>` +
    P.lab('D\'ALEMBERT', '无界弦的初值问题',
      `<div class="ctl"><label class="ctl-label" for="wv-phi">初始位移 \\(\\varphi\\)</label><select id="wv-phi"><option value="tri">三角形</option><option value="gauss">高斯</option><option value="sine">一个周期的正弦</option><option value="none">0</option></select></div>
       ${P.rangeHTML('wv-A', '位移幅值', 0, 1.5, 0.05, 1)}
       <div class="ctl"><label class="ctl-label" for="wv-psi">初始速度 \\(\\psi\\)</label><select id="wv-psi"><option value="none">0（静止释放）</option><option value="box">方波（锤击）</option><option value="gauss">高斯</option><option value="tri">三角形</option></select></div>
       ${P.rangeHTML('wv-V', '速度幅值 \\(V\\)', -1.5, 1.5, 0.05, 0.6)}
       ${P.rangeHTML('wv-a', '波速 \\(a\\)', 0.5, 2, 0.05, 1)}
       <div class="ctl"><div class="ctl-label">分解方式</div><div id="wv-view"></div></div>
       <div id="wv-plA"></div>
       <div class="tasks"><b>试一试</b><ul><li>在 x–t 图上点击任意一点，看它的依赖区间。</li><li>只给方波速度：远处出现高为 \\(V/a\\) 的平台。</li><li>φ 为三角形、ψ 为同形三角形且 V=a：只剩左行波。</li></ul></div>`,
      `<div class="problem" id="wv-aProb"></div>
       <div class="plot-block"><div class="plot-cap"><span>\\(u(x,t)\\) 及其分解</span><span class="legend" id="wv-aLeg"></span></div><canvas id="wv-cvA"></canvas></div>
       ${P.plotBlock('wv-cvXT', 'x–t 平面上的解：虚线为初始支集的影响区域边界，彩色线为选点的两条特征线，粗线为依赖区间（点击图像选点）', '', 'tall clickable')}
       <div class="stats" id="wv-aStats"></div>`) +
    `<h3 class="sub">B　半无界弦：反射法</h3>` +
    P.lab('REFLECT', '端点处的反射与延拓',
      `<div class="ctl"><div class="ctl-label">端点类型</div><div id="wv-end"></div></div>
       <div class="ctl"><div class="ctl-label">初始状态</div><div id="wv-bmode"></div></div>
       <div class="ctl"><div class="ctl-label">脉冲形状</div><div id="wv-shape"></div></div>
       ${P.rangeHTML('wv-x0', '脉冲中心 \\(x_0\\)', 2, 8, 0.1, 5)}
       <div id="wv-plB"></div>
       <div class="tasks"><b>试一试</b><ul><li>固定端：脉冲到达端点时与倒相的镜像相遇，端点始终不动。</li><li>自由端：到达时端点位移加倍，反射脉冲同相。</li><li>选"静止释放"：一半直接向右走，另一半先向左、反射后再向右。</li></ul></div>`,
      `<div class="problem" id="wv-bProb"></div>
       ${P.plotBlock('wv-cvB', '实线：物理区域的解；虚线：延拓后的虚拟区域', P.lg('accent', '总位移') + P.lg('violet', '入射波') + P.lg('warm', '反射波（镜像贡献）'))}
       ${P.plotBlock('wv-cvBxt', 'x–t 平面：虚线为脉冲中心的特征线及其反射', '', 'short')}
       <div class="stats" id="wv-bStats"></div>`) +
    `<h3 class="sub">C　有界弦：驻波、模态叠加与能量守恒</h3>` +
    P.lab('STANDING', '与模块 04 共用特征值求解器',
      `<div class="ctl"><div class="ctl-label">左端</div><div id="wv-cL"></div></div>
       <div class="ctl"><div class="ctl-label">右端</div><div id="wv-cR"></div></div>
       ${P.rangeHTML('wv-ch', 'Robin 系数 \\(h\\)（两端相同）', 0, 10, 0.1, 2)}
       ${P.rangeHTML('wv-cLen', '弦长 \\(L\\)', 0.5, 2, 0.05, 1)}
       ${P.rangeHTML('wv-ca', '波速 \\(a\\)', 0.5, 2, 0.05, 1)}
       <div class="ctl"><label class="ctl-label" for="wv-cphi">初始位移 \\(\\varphi\\)</label><select id="wv-cphi"><option value="tri">拨弦（三角形）</option><option value="gauss">窄高斯</option><option value="mode">第 3 个特征函数</option><option value="none">0</option></select></div>
       ${P.rangeHTML('wv-cp', '拨弦位置 \\(p\\)', 0.05, 0.95, 0.05, 0.3)}
       <div class="ctl"><label class="ctl-label" for="wv-cpsi">初始速度 \\(\\psi\\)</label><select id="wv-cpsi"><option value="none">0</option><option value="box">锤击（窄方波）</option><option value="gauss">锤击（高斯）</option><option value="uniform">均匀速度</option></select></div>
       ${P.rangeHTML('wv-cxh', '锤击位置', 0.05, 0.95, 0.05, 0.7)}
       ${P.rangeHTML('wv-cN', '截断阶数 \\(N\\)', 1, 60, 1, 30)}
       <div class="ctl"><div class="ctl-label">单独查看的正则模态</div><div id="wv-csel"></div></div>
       <div id="wv-plC"></div>`,
      `<div class="problem" id="wv-cProb"></div>
       ${P.plotBlock('wv-cvMode', '正则模态 \\(X_n(x)\\cos\\omega_nt\\)：阴影为包络，黑点为节点，空心圈为腹点', '')}
       <div class="stats" id="wv-cModeStats"></div>
       ${P.plotBlock('wv-cvC', '多模态叠加 \\(u(x,t)\\)', P.lg('accent', 'u(x,t)') + P.lg('warm', 'φ(x)'))}
       ${P.plotBlock('wv-cvCoef', '系数 \\(A_n\\)（位移）与 \\(B_n\\)（速度）', P.lg('accent', 'A_n') + P.lg('warm', 'B_n'), 'short')}
       ${P.plotBlock('wv-cvE', '能量随时间变化（ρ=1）', P.lg('warm', '动能 K') + P.lg('violet', '势能 V（含弹簧）') + P.lg('accent', '总能量 E'), 'short')}
       <div class="stats" id="wv-cStats"></div>
       <p class="small muted">能量 \\(E=\\tfrac12\\int_0^L(u_t^2+a^2u_x^2)\\,dx\\)（Robin 端另加弹簧势能 \\(\\tfrac12a^2hu^2\\)）。守恒的条件：无阻尼、无外力、边界条件为 D、N 或 \\(h\\ge0\\) 的 R（边界不做功）。</p>`) +
    `<h3 class="sub">D　Duhamel 原理：受迫振动</h3>` +
    P.lab('DUHAMEL', '两端固定的弦，零初始条件，外力 \\(f(x,t)\\)',
      `<div class="ctl"><div class="ctl-label">外力的时间形式</div><div id="wv-dforce"></div></div>
       <div class="ctl"><div class="ctl-label">外力的空间分布 \\(g(x)\\)</div><div id="wv-dprof"></div></div>
       ${P.rangeHTML('wv-dxf', '作用位置 \\(x_f\\)', 0.1, 0.9, 0.05, 0.3)}
       ${P.rangeHTML('wv-dOm', '驱动频率 \\(\\Omega/\\omega_1\\)', 0.1, 3.5, 0.01, 0.9)}
       ${P.rangeHTML('wv-dtau', '冲量时刻 \\(\\tau_0\\)', 0, 1.5, 0.05, 0.3)}
       ${P.rangeHTML('wv-dM', '黎曼和切片数 \\(M\\)', 2, 60, 1, 8)}
       <div id="wv-plD"></div>
       <div class="tasks"><b>试一试</b><ul><li>把 \\(\\Omega/\\omega_1\\) 拉到 1.00 和 3.00：探针振幅线性增长（共振）。</li><li>均匀力 + \\(\\Omega/\\omega_1=2\\)：不共振，因为 \\(g_2=0\\)（偶模态不被均匀力激发）。</li><li>增大 \\(M\\)：紫色黎曼和逼近精确解。</li></ul></div>`,
      `<div class="problem" id="wv-dProb"></div>
       ${P.plotBlock('wv-cvD', '灰线为各时刻冲量产生的响应 \\(w(x,t;\\tau_k)\\Delta\\tau\\)，它们的和即 Duhamel 积分', P.lg('accent', '精确解（模态闭式）') + P.lg('violet', 'M 片黎曼和') + P.lg('ink-3', '单片响应') + P.lg('warm', '静挠度'))}
       ${P.plotBlock('wv-cvDh', '探针点的位移历史 \\(u(x_p,t)\\)', '', 'short')}
       <div class="stats" id="wv-dStats"></div>
       <div class="formula-box small">\\(w\\) 满足 \\(w_{tt}=a^2w_{xx}\\)，\\(w|_{t=\\tau}=0\\)，\\(w_t|_{t=\\tau}=f(x,\\tau)\\)。按模态：\\(u_n(t)=\\int_0^t f_n(\\tau)\\dfrac{\\sin\\omega_n(t-\\tau)}{\\omega_n}d\\tau\\)，恒力、简谐力和冲量都有闭式；共振 \\(\\Omega=\\omega_n\\) 时 \\(u_n=\\dfrac{\\sin\\omega_nt-\\omega_nt\\cos\\omega_nt}{2\\omega_n^2}g_n\\)。</div>`)) +
  P.sec(ID, 'meaning', '物理意义', R`<div class="cols">
    <div class="card"><h3>有限传播速度</h3><p class="small muted">扰动以速度 \(a\) 传播，依赖区间之外的初始数据不影响 \((x_0,t_0)\)。这是双曲型方程与热方程的根本区别。</p></div>
    <div class="card"><h3>间断沿特征线传播</h3><p class="small muted">三角形的尖角、方波的跳跃都原样保留并沿 \(x\pm at=c\) 移动，不会被平滑。</p></div>
    <div class="card"><h3>反射的相位</h3><p class="small muted">固定端产生倒相反射（端点力迫使位移为零）；自由端同相反射（端点斜率为零，位移可以加倍）。</p></div>
    <div class="card"><h3>驻波与音色</h3><p class="small muted">有界弦只能以离散频率 \(\omega_n\) 振动；系数 \(A_n,B_n\) 决定各泛音的强弱。拨弦位置在某阶节点上时，该阶不被激发。</p></div></div>` +
    P.pitfalls([
      ['达朗贝尔公式的适用范围：', '只适用于无界弦；有边界时必须先延拓。'],
      ['系数：', String.raw`位移项有 \(\tfrac12\)，速度项是 \(\tfrac1{2a}\int\)。用 \(t=0\) 代入检验两个初始条件。`],
      ['延拓的奇偶：', String.raw`固定端（D）奇延拓，自由端（N）偶延拓；\(\psi\) 与 \(\varphi\) 用同样的奇偶性。`],
      ['有界弦的周期：', String.raw`D–D、N–N 周期 \(2L/a\)，D–N 周期 \(4L/a\)，含 R 时一般不是周期运动。`],
      ['Duhamel 原理：', '冲量给的是初速度而不是初位移；子问题从时刻 τ 开始，所以是 sin ω(t−τ)。'],
      ['能量守恒：', '需要无阻尼、无外力、边界不做功；有外力时能量变化等于外力做的功。']
    ]) + `<div class="bench exp-only"><h3 class="sub">数值基准</h3><div id="wv-bench"></div></div>`) +
  P.sec(ID, 'quiz', '检查你的理解 <span class="score" id="wv-score"></span>', '<div id="wv-quiz"></div>');
}

/* ================= 初始化 ================= */
const BCO = [{ v: 'D', label: 'D 固定' }, { v: 'N', label: 'N 自由' }, { v: 'R', label: 'R 弹性' }];
function buildSelC() { P.segment($('wv-csel'), Array.from({ length: 6 }, (_, j) => ({ v: j, label: j + slC.offset })), SC.sel, v => { SC.sel = +v; renderC(); }); }
function syncRowsC() { $('wv-ch-row').hidden = SC.bcL !== 'R' && SC.bcR !== 'R'; $('wv-cp-row').hidden = SC.phi !== 'tri'; $('wv-cxh-row').hidden = SC.psi === 'none' || SC.psi === 'uniform'; }
function syncRowsD() { $('wv-dxf-row').hidden = SD.prof !== 'local'; $('wv-dOm-row').hidden = SD.force !== 'harm'; $('wv-dtau-row').hidden = SD.force !== 'impulse'; $('wv-dM-row').hidden = SD.force === 'impulse'; }
P.register('wave', {
  init(root) {
    root.innerHTML = html();
    P.stepper($('wv-steps'), STEPS);
    // A
    $('wv-phi').addEventListener('change', e => { SA.phi = e.target.value; renderA(); });
    $('wv-psi').addEventListener('change', e => { SA.psi = e.target.value; renderA(); });
    P.bindRange('wv-A', v => { SA.A = v; renderA(); }, v => v.toFixed(2));
    P.bindRange('wv-V', v => { SA.V = v; renderA(); }, v => v.toFixed(2));
    P.bindRange('wv-a', v => { SA.a = v; plA.refresh(); }, v => v.toFixed(2));
    P.segment($('wv-view'), [{ v: 'lr', label: '左右行波' }, { v: 'pv', label: '位移 / 速度贡献' }], SA.view, v => { SA.view = v; renderA(); });
    plA = P.player($('wv-plA'), { id: 'wv-tA', tmax: tmaxA, duration: 10, onChange: renderA });
    const cv = $('wv-cvXT');
    cv.addEventListener('pointerdown', e => { if (!cv._inv) return; const r = cv.getBoundingClientRect(); SA.x0 = Math.max(-XA, Math.min(XA, cv._inv.x(e.clientX - r.left))); SA.t0 = Math.max(0, Math.min(tmaxA(), cv._inv.y(e.clientY - r.top))); renderA(); });
    // B
    P.segment($('wv-end'), [{ v: 'D', label: '固定端 u=0' }, { v: 'N', label: '自由端 u_x=0' }], SB.end, v => { SB.end = v; renderB(); });
    P.segment($('wv-bmode'), [{ v: 'left', label: '向端点行进' }, { v: 'rest', label: '静止释放' }], SB.mode, v => { SB.mode = v; renderB(); });
    P.segment($('wv-shape'), [{ v: 'tri', label: '三角形' }, { v: 'gauss', label: '高斯' }], SB.shape, v => { SB.shape = v; renderB(); });
    P.bindRange('wv-x0', v => { SB.x0 = v; plB.refresh(); }, v => v.toFixed(1));
    plB = P.player($('wv-plB'), { id: 'wv-tB', tmax: tmaxB, duration: 10, onChange: renderB });
    // C
    const allC = () => { updateC(); buildSelC(); syncRowsC(); plC.refresh(); };
    P.segment($('wv-cL'), BCO, SC.bcL, v => { SC.bcL = v; allC(); });
    P.segment($('wv-cR'), BCO, SC.bcR, v => { SC.bcR = v; allC(); });
    P.bindRange('wv-ch', v => { SC.h = v; allC(); }, v => v.toFixed(1));
    P.bindRange('wv-cLen', v => { SC.L = v; allC(); }, v => v.toFixed(2));
    P.bindRange('wv-ca', v => { SC.a = v; allC(); }, v => v.toFixed(2));
    $('wv-cphi').addEventListener('change', e => { SC.phi = e.target.value; allC(); });
    $('wv-cpsi').addEventListener('change', e => { SC.psi = e.target.value; allC(); });
    P.bindRange('wv-cp', v => { SC.p = v; allC(); }, v => v.toFixed(2));
    P.bindRange('wv-cxh', v => { SC.xh = v; allC(); }, v => v.toFixed(2));
    P.bindRange('wv-cN', v => { SC.N = v; computeEnergy(); plC.refresh(); });
    updateC();
    plC = P.player($('wv-plC'), { id: 'wv-tC', tmax: tmaxC, duration: 12, onChange: renderC });
    buildSelC(); syncRowsC();
    // D
    const allD = () => { updateD(); syncRowsD(); plD.refresh(); };
    P.segment($('wv-dforce'), [{ v: 'harm', label: '简谐 sin Ωt' }, { v: 'const', label: '恒力' }, { v: 'impulse', label: '冲量 δ(t−τ₀)' }], SD.force, v => { SD.force = v; allD(); });
    P.segment($('wv-dprof'), [{ v: 'local', label: '局部（窄高斯）' }, { v: 'uniform', label: '均匀' }], SD.prof, v => { SD.prof = v; allD(); });
    P.bindRange('wv-dxf', v => { SD.xf = v; allD(); }, v => v.toFixed(2));
    P.bindRange('wv-dOm', v => { SD.Om = v; histD = null; plD.refresh(); }, v => v.toFixed(2));
    P.bindRange('wv-dtau', v => { SD.tau0 = v; histD = null; plD.refresh(); }, v => v.toFixed(2));
    P.bindRange('wv-dM', v => { SD.M = v; plD.refresh(); });
    updateD();
    plD = P.player($('wv-plD'), { id: 'wv-tD', tmax: tmaxD, duration: 14, onChange: renderD });
    syncRowsD();
    P.quiz($('wv-quiz'), QUIZ, $('wv-score'));
    runBench();
    $('wv-bench').addEventListener('click', e => { if (e.target.closest('[data-rerun]')) runBench(); });
  },
  render() { plA.refresh(); plB.refresh(); plC.refresh(); plD.refresh(); }
});
})();
