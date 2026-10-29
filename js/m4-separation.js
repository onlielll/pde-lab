/* m4-separation.js — 模块 04：分离变量法与特征值问题 */
(function () {
'use strict';
const P = window.PDE, PI = Math.PI, $ = P.$, fx = P.fx, sci = P.sci;
const ID = 'sv', NM = 60;
const S = { bcL: 'D', bcR: 'D', L: 1, h1: 2, h2: 2, sel: 0, fType: 'tri', p: 0.3, N: 10, eq: 'heat', a: 1 };
const drawPts = new Float64Array(121);
let sl = null, fv = null, coef = null, player = null, yRange = [-0.2, 1.2];

/* ---------- 计算 ---------- */
function updateEigen() { sl = P.solveSL({ L: S.L, bcL: S.bcL, bcR: S.bcR, h1: S.h1, h2: S.h2, NM, NG: 2400 }); if (S.sel >= sl.modes.length) S.sel = 0; }
function fShape(xi) {
  switch (S.fType) {
    case 'tri': return xi < S.p ? xi / S.p : (1 - xi) / (1 - S.p);
    case 'para': return 4 * xi * (1 - xi);
    case 'step': return xi < 0.5 ? 1 : xi > 0.5 ? 0 : 0.5;
    case 'const': return 1;
    case 'gauss': return Math.exp(-((xi - 0.5) ** 2) / (2 * 0.07 * 0.07));
    case 'draw': { const t = xi * 120, i = Math.min(119, Math.floor(t)), r = t - i; return drawPts[i] * (1 - r) + drawPts[i + 1] * r; }
  }
  return 0;
}
function updateCoeffs() {
  const NG = sl.NG; fv = new Float64Array(NG + 1);
  if (S.fType === 'mode') { const M = sl.modes[1].X, mx = P.maxAbs(M); for (let i = 0; i <= NG; i++) fv[i] = M[i] / mx; }
  else for (let i = 0; i <= NG; i++) fv[i] = fShape(i / NG);
  coef = P.project(sl, fv);
}
function analyticDD(n, type = S.fType, p = S.p) {
  switch (type) {
    case 'para': return n % 2 ? 32 / (n * PI) ** 3 : 0;
    case 'tri': return 2 * Math.sin(n * PI * p) / ((n * PI) ** 2 * p * (1 - p));
    case 'step': return 2 / (n * PI) * (1 - Math.cos(n * PI / 2));
    case 'const': return n % 2 ? 4 / (n * PI) : 0;
  }
  return null;
}
function timeFactor(m, t) {
  if (S.eq === 'heat') return Math.exp(-S.a * S.a * m.lam * t);
  if (m.kind === 'pos') return Math.cos(S.a * m.k * t);
  if (m.kind === 'zero') return 1;
  return Math.cosh(S.a * m.k * t);
}
function series(Nt, t) {
  const NG = sl.NG, out = new Float64Array(NG + 1), n = Math.min(Nt, sl.modes.length);
  for (let j = 0; j < n; j++) {
    const c = t == null ? coef[j] : coef[j] * timeFactor(sl.modes[j], t);
    if (c === 0) continue;
    const X = sl.modes[j].X;
    for (let i = 0; i <= NG; i++) out[i] += c * X[i];
  }
  return out;
}
const nLabel = j => j + sl.offset;
const firstPos = () => sl.modes.find(m => m.kind === 'pos');
function tMax() {
  if (S.eq === 'heat') { const m = sl.modes.find(m => Math.abs(m.lam) > 1e-9); return 1.5 / (S.a * S.a * Math.abs(m.lam)); }
  return 2 * (2 * PI / (S.a * firstPos().k));
}

/* ---------- 公式 ---------- */
function eqTex() {
  const m = {
    DD: String.raw`\sin kL=0\;\Rightarrow\;k_n=\frac{n\pi}{L},\ n=1,2,\dots`,
    NN: String.raw`\sin kL=0\;\Rightarrow\;k_n=\frac{n\pi}{L},\ n=0,1,2,\dots\quad(\lambda_0=0,\ X_0=1)`,
    DN: String.raw`\cos kL=0\;\Rightarrow\;k_n=\frac{(n-\frac12)\pi}{L},\ n=1,2,\dots`,
    ND: String.raw`\cos kL=0\;\Rightarrow\;k_n=\frac{(n-\frac12)\pi}{L},\ n=1,2,\dots`,
    DR: String.raw`\tan kL=-\frac{k}{h_2}`, RD: String.raw`\tan kL=-\frac{k}{h_1}`,
    NR: String.raw`\tan kL=\frac{h_2}{k}`, RN: String.raw`\tan kL=\frac{h_1}{k}`,
    RR: String.raw`\tan kL=\frac{(h_1+h_2)\,k}{k^2-h_1h_2}`
  };
  const c = S.bcL + S.bcR, rob = c.includes('R');
  return m[c] + (rob ? String.raw`\quad\text{（超越方程，数值求根）}` : '');
}
function efTex() { return { D: String.raw`X_n=\sin k_nx`, N: String.raw`X_n=\cos k_nx`, R: String.raw`X_n=k_n\cos k_nx+h_1\sin k_nx` }[S.bcL]; }

/* ---------- 实验 A ---------- */
function renderSystem() {
  P.tex($('sv-sys'), String.raw`\[ \begin{cases} X''+\lambda X=0,\quad 0\lt x\lt L,\\ ${P.BCTEX.L[S.bcL]},\qquad ${P.BCTEX.R[S.bcR]} \end{cases}\qquad \lambda=k^2\gt 0\text{ 时 } ${efTex()}\ (\text{差常数因子}) \]`);
  P.tex($('sv-eq'), String.raw`\[ ${eqTex()} \]`);
  const neg = sl.nNeg, z = sl.modes.some(m => m.kind === 'zero');
  let note = '';
  if (neg) note = `<span class="pill warn">出现 ${neg} 个负特征值</span> Robin 系数为负，Rayleigh 商中的边界项 \\(h X^2\\lt 0\\) 可以让 \\(\\lambda\\lt0\\)。对应 \\(X=\\mu\\cosh\\mu x+h_1\\sinh\\mu x\\) 型双曲函数，由 \\(\\lambda=-\\mu^2\\) 的方程 \\(G(\\mu)=0\\) 数值求出。热传导中该模态按 \\(e^{a^2\\mu^2t}\\) 增长。`;
  else if (z) note = `<span class="pill info">存在零特征值</span> \\(\\lambda=0\\) 对应线性函数 \\(X=b_1+a_1x\\)${S.bcL === 'N' && S.bcR === 'N' ? '，N–N 时就是常数模态' : '，这是 Robin 系数恰好满足 \\(a_2b_1+a_1a_2L+a_1b_2=0\\) 的特殊情形'}。`;
  else note = '<span class="pill ok">全部特征值为正</span> 在当前边界条件下 \\(\\lambda\\le0\\) 只有零解。';
  P.tex($('sv-signnote'), note);
}
function renderEigen() {
  const f = P.frame($('sv-cvEig'), { xmin: 0, xmax: S.L, ymin: -1.15, ymax: 1.15, xlabel: 'x' });
  if (!f) return;
  const show = Math.min(6, sl.modes.length);
  for (let j = 0; j < show; j++) if (j !== S.sel) { const X = sl.modes[j].X, mx = P.maxAbs(X); P.line(f, sl.xs, X.map(v => v / mx), P.C['ink-3'], 1, [3, 3], 4); }
  const m = sl.modes[S.sel], Xs = m.X, mx = P.maxAbs(Xs);
  P.line(f, sl.xs, Xs.map(v => v / mx), P.C.accent, 2.4, null, 2);
  P.bcGlyphs(f, S.bcL, S.bcR);
  let zeros = 0;
  for (let i = 2; i < sl.NG - 1; i++) if (Xs[i] * Xs[i + 1] < 0) { zeros++; P.dot(f, sl.xs[i], 0, P.C.accent, 3); }
  const fp = firstPos();
  P.stats($('sv-eigStats'), [
    ['模态编号 n', nLabel(S.sel) + (m.kind === 'neg' ? ' <span class="pill warn">λ&lt;0</span>' : m.kind === 'zero' ? ' <span class="pill info">λ=0</span>' : '')],
    [m.kind === 'neg' ? 'μ_n（λ=−μ²）' : 'k_n', fx(m.k, 6)],
    ['λ_n', fx(m.lam, 5)],
    ['Rayleigh 商', fx(m.ray, 5)],
    ['内部零点数', zeros + (zeros === S.sel ? ' <span class="pill ok">= 序号−1 ✓</span>' : '')],
    ['基频 ω₁/a（首个正 k）', fx(fp.k, 4)],
    ['负特征值个数', sl.nNeg]
  ]);
  renderTan(); renderTable(); renderSystem();
}
function renderTan() {
  const L = S.L, qmax = 6, Ym = 6, a1 = sl.a1, b1 = sl.b1, a2 = sl.a2, b2 = sl.b2;
  const cotForm = a1 * a2 === 0 && b1 * b2 === 0;
  const xt = []; for (let q = 0; q <= qmax; q++) xt.push({ v: q, label: q === 0 ? '0' : (q === 1 ? '' : q) + 'π' });
  const f = P.frame($('sv-cvTan'), { xmin: 0, xmax: qmax, ymin: -Ym, ymax: Ym, xticks: xt, xlabel: 'kL' });
  if (!f) return;
  const n = 1500, q = new Float64Array(n + 1), y1 = new Float64Array(n + 1), y2 = new Float64Array(n + 1);
  for (let i = 0; i <= n; i++) {
    const qq = qmax * i / n, k = qq * PI / L; q[i] = qq;
    if (cotForm) { const c = Math.cos(qq * PI) / Math.sin(qq * PI); y1[i] = Math.abs(c) > 40 ? NaN : c; y2[i] = 0; }
    else {
      const t = Math.tan(qq * PI); y1[i] = Math.abs(t) > 40 ? NaN : t;
      const r = -(a2 * b1 + a1 * b2) * k / (a1 * a2 - b1 * b2 * k * k);
      y2[i] = Math.abs(r) > 40 || !isFinite(r) ? NaN : r;
    }
  }
  P.line(f, q, y1, P.C['ink-2'], 1.3);
  P.line(f, q, y2, P.C.warm, 1.8);
  sl.modes.forEach((m, j) => {
    if (m.kind !== 'pos') return;
    const qq = m.k * L / PI; if (qq > qmax) return;
    const y = cotForm ? 0 : Math.tan(m.k * L);
    if (Math.abs(y) <= Ym) P.dot(f, qq, y, j === S.sel ? P.C.accent : P.C.ink, j === S.sel ? 5 : 3.5);
  });
  $('sv-tanLegend').innerHTML = cotForm ? P.lg('ink-2', 'cot kL') + P.lg('warm', '0') : P.lg('ink-2', 'tan kL') + P.lg('warm', '右端函数') + P.lg('ink', '交点 = 正特征值');
}
function renderTable() {
  const rows = Math.min(8, sl.modes.length);
  let h = '<thead><tr><th>n</th><th>类型</th><th class="num">kL/π 或 μL</th><th class="num">λ_n</th><th class="num">Rayleigh 商</th><th class="num">相对差</th></tr></thead><tbody>';
  for (let j = 0; j < rows; j++) {
    const m = sl.modes[j], rel = Math.abs(m.ray - m.lam) / Math.max(1, Math.abs(m.lam));
    h += `<tr class="${j === S.sel ? 'hl' : ''}"><td>${nLabel(j)}</td><td>${m.kind === 'neg' ? 'λ&lt;0' : m.kind === 'zero' ? 'λ=0' : 'λ&gt;0'}</td><td class="num">${m.kind === 'pos' ? fx(m.k * S.L / PI, 5) : m.kind === 'neg' ? fx(m.k * S.L, 5) : '0'}</td><td class="num">${fx(m.lam, 4)}</td><td class="num">${fx(m.ray, 4)}</td><td class="num">${sci(rel)}</td></tr>`;
  }
  $('sv-eigTable').innerHTML = h + '</tbody>';
  const g = Math.min(6, sl.modes.length);
  let gh = '<tbody>', mo = 0;
  for (let a = 0; a < g; a++) {
    gh += '<tr>';
    for (let b = 0; b < g; b++) {
      const v = P.dotS(sl.modes[a].X, sl.modes[b].X, sl.h) / Math.sqrt(sl.modes[a].n2 * sl.modes[b].n2);
      if (a !== b) mo = Math.max(mo, Math.abs(v));
      gh += a === b ? `<td class="diag">${v.toFixed(3)}</td>` : `<td>${Math.abs(v) < 5e-4 ? '≈0' : v.toFixed(3)}</td>`;
    }
    gh += '</tr>';
  }
  $('sv-gram').innerHTML = gh + '</tbody>';
  $('sv-gramNote').innerHTML = `最大非对角元 <span class="mono">${sci(mo)}</span>：在 Simpson 积分误差范围内为 0，特征函数两两正交。`;
}

/* ---------- 实验 B ---------- */
function computeYRange() {
  if (S.fType === 'draw') { yRange = [-1.25, 1.25]; return; }
  let lo = 0, hi = 0; for (let i = 0; i < fv.length; i++) { lo = Math.min(lo, fv[i]); hi = Math.max(hi, fv[i]); }
  const pad = 0.22 * (hi - lo || 1); yRange = [lo - pad, hi + pad];
}
function compat() {
  const NG = sl.NG, hh = sl.h, f0 = fv[0], fL = fv[NG];
  const d0 = (-3 * fv[0] + 4 * fv[1] - fv[2]) / (2 * hh), dL = (3 * fv[NG] - 4 * fv[NG - 1] + fv[NG - 2]) / (2 * hh);
  const res = [
    S.bcL === 'D' ? Math.abs(f0) : S.bcL === 'N' ? Math.abs(d0) : Math.abs(d0 - S.h1 * f0),
    S.bcR === 'D' ? Math.abs(fL) : S.bcR === 'N' ? Math.abs(dL) : Math.abs(dL + S.h2 * fL)
  ];
  const bad = res.map(r => r > 1e-3);
  if (!bad[0] && !bad[1]) return '<span class="pill ok">f 满足两端边界条件</span>';
  const which = [bad[0] ? '左端' : '', bad[1] ? '右端' : ''].filter(Boolean).join('、');
  const dBad = (bad[0] && S.bcL === 'D') || (bad[1] && S.bcR === 'D');
  return dBad ? `<span class="pill warn">${which}不相容</span> D 端要求 f=0：级数在端点收敛到 0，附近出现 Gibbs 过冲，非一致收敛`
    : `<span class="pill info">${which}导数条件不满足</span> 级数仍一致收敛到 f（f 连续），只是系数衰减变慢`;
}
function renderFourier() {
  computeYRange();
  const f = P.frame($('sv-cvF'), { xmin: 0, xmax: S.L, ymin: yRange[0], ymax: yRange[1], xlabel: 'x' });
  const SN = series(S.N, null);
  if (f) {
    P.line(f, sl.xs, fv, P.C.warm, 2, null, 2);
    P.line(f, sl.xs, SN, P.C.accent, 1.8);
    P.bcGlyphs(f, S.bcL, S.bcR);
    if (S.fType === 'draw') { f.ctx.fillStyle = P.C['ink-3']; f.ctx.textAlign = 'left'; f.ctx.textBaseline = 'top'; f.ctx.fillText('按住拖动来绘制 f(x)', f.pad.l + 8, f.pad.t + 6); }
  }
  const NG = sl.NG, e = new Float64Array(NG + 1);
  let maxErr = 0, fmax = -Infinity, smax = -Infinity;
  for (let i = 0; i <= NG; i++) { e[i] = SN[i] - fv[i]; maxErr = Math.max(maxErr, Math.abs(e[i])); fmax = Math.max(fmax, fv[i]); smax = Math.max(smax, SN[i]); }
  const l2 = Math.sqrt(P.dotS(e, e, sl.h) / (P.dotS(fv, fv, sl.h) || 1));
  P.stats($('sv-fStats'), [
    ['截断阶数 N', S.N],
    ['相对 L² 误差', sci(l2)],
    ['最大误差 max|S_N−f|', fx(maxErr, 4)],
    ['过冲 max S_N − max f', fx(Math.max(0, smax - fmax), 4)],
    ['f(0) / S_N(0)', fx(fv[0], 3) + ' / ' + fx(SN[0], 3)],
    ['f(L) / S_N(L)', fx(fv[NG], 3) + ' / ' + fx(SN[NG], 3)],
    ['系数衰减', decayFit()]
  ]);
  P.tex($('sv-compat'), compat());
  renderCoef();
  const el = $('sv-fCheck');
  if (S.bcL === 'D' && S.bcR === 'D' && analyticDD(1) !== null) {
    let md = 0; for (let n = 1; n <= NM; n++) md = Math.max(md, Math.abs(coef[n - 1] - analyticDD(n)));
    P.tex(el, `<span class="pill ok">解析校验</span> 两端 D 时 \\(c_n\\) 即正弦级数系数 \\(b_n=\\tfrac{2}{L}\\int_0^L f\\sin\\tfrac{n\\pi x}{L}dx\\)。前 60 项数值与闭式的最大偏差 <span class="mono">${sci(md)}</span>。`);
  } else P.tex(el, '<span class="pill info">说明</span> 系数 \\(c_n=\\langle f,X_n\\rangle/\\langle X_n,X_n\\rangle\\) 由 2400 段 Simpson 积分计算；阶跃和尖角落在网格节点上，保持高精度。两端 D 且选三角形、抛物线、阶跃、常数时可与闭式对比。');
}
function decayFit() {
  const mx = P.maxAbs(coef), pts = [];
  for (let j = 4; j < coef.length; j++) { const c = Math.abs(coef[j]); if (c > 1e-10 * mx && c > 1e-13) pts.push([Math.log(j + 1), Math.log(c)]); }
  if (pts.length < 6) return S.fType === 'mode' ? '只有单个模态' : '快于任意幂次';
  const n = pts.length, sx = pts.reduce((s, p) => s + p[0], 0) / n, sy = pts.reduce((s, p) => s + p[1], 0) / n;
  let num = 0, den = 0; pts.forEach(([x, y]) => { num += (x - sx) * (y - sy); den += (x - sx) ** 2; });
  const p = -num / den;
  return p > 4.5 ? '快于任意幂次' : '≈ n^(−' + p.toFixed(2) + ')';
}
function renderCoef() {
  const yt = [-10, -8, -6, -4, -2, 0].map(v => ({ v, label: v === 0 ? '1' : '1e' + v }));
  const f = P.frame($('sv-cvC'), { xmin: 0, xmax: NM + 1, ymin: -11, ymax: 0.8, yticks: yt, xlabel: 'n' });
  if (!f) return;
  const { ctx, X, Y } = f, mx = P.maxAbs(coef);
  [1, 2, 3].forEach(pw => {
    const xx = [], yy = []; for (let k = 1; k <= NM; k++) { xx.push(k); yy.push(Math.log10(mx / k ** pw)); }
    P.line(f, xx, yy, P.C['ink-3'], 1, [4, 4]);
    const yl = Y(Math.log10(mx / NM ** pw));
    if (yl > f.pad.t + 8 && yl < f.pad.t + f.ph) { ctx.fillStyle = P.C['ink-3']; ctx.textAlign = 'left'; ctx.textBaseline = 'bottom'; ctx.fillText('n⁻' + '¹²³'[pw - 1], X(NM - 5), yl - 2); }
  });
  P.clip(f);
  for (let j = 0; j < coef.length; j++) {
    const c = Math.abs(coef[j]), nn = nLabel(j) || 0.5; if (c < 1e-11) continue;
    const y = Math.log10(c), on = j < S.N;
    ctx.strokeStyle = on ? P.C.accent : P.C['ink-3']; ctx.lineWidth = on ? 2 : 1;
    ctx.beginPath(); ctx.moveTo(X(nn), Y(-11)); ctx.lineTo(X(nn), Y(y)); ctx.stroke();
    ctx.fillStyle = ctx.strokeStyle; ctx.beginPath(); ctx.arc(X(nn), Y(y), on ? 2.6 : 2, 0, 2 * PI); ctx.fill();
  }
  ctx.restore();
}

/* ---------- 实验 C ---------- */
function renderTime() {
  const t = player.t();
  const u = series(S.N, t);
  let umax = 0; for (let i = 0; i < u.length; i++) umax = Math.max(umax, Math.abs(u[i]));
  const lim = S.eq === 'heat' && sl.nNeg ? [Math.min(yRange[0], -umax * 1.1), Math.max(yRange[1], umax * 1.1)] : yRange;
  const f = P.frame($('sv-cvU'), { xmin: 0, xmax: S.L, ymin: lim[0], ymax: lim[1], xlabel: 'x' });
  if (f) { P.line(f, sl.xs, fv, P.C.warm, 1.2, [4, 4], 3); P.line(f, sl.xs, u, P.C.accent, 2.2); P.bcGlyphs(f, S.bcL, S.bcR); }
  const fp = firstPos();
  let st;
  if (S.eq === 'heat') {
    const m1 = sl.modes[0], tau = 1 / (S.a * S.a * Math.abs(fp.lam)), nxt = sl.modes[Math.min(S.N, sl.modes.length - 1)];
    const lim0 = sl.nNeg ? '指数增长（λ₁&lt;0）' : sl.modes[0].kind === 'zero' ? '常数模态 ' + fx(coef[0], 4) : '0';
    st = [
      ['t / τ（τ=1/(a²λ_+)）', fx(t / tau, 3)],
      ['∫₀ᴸ u dx（总热量）', fx(P.integ(u, sl.h), 5)],
      ['max |u|', fx(umax, 4)],
      ['首个截断项衰减因子', t === 0 ? '1（t=0 截断误差最大）' : sci(Math.exp(-S.a * S.a * nxt.lam * t))],
      ['最慢模态增长/衰减率', fx(-S.a * S.a * m1.lam, 4)],
      ['t → ∞', lim0]
    ];
  } else {
    const pos = sl.modes.filter(m => m.kind === 'pos'), T1 = 2 * PI / (S.a * fp.k);
    const r2 = pos[1].k / fp.k, r3 = pos[2].k / fp.k;
    const harmonic = Math.abs(r2 - Math.round(r2)) < 1e-6 && Math.abs(r3 - Math.round(r3)) < 1e-6;
    st = [
      ['t / T₁', fx(t / T1, 3)],
      ['T₁ = 2π/(a k₁)', fx(T1, 4)],
      ['k₂/k₁, k₃/k₁', fx(r2, 4) + ', ' + fx(r3, 4)],
      ['频率关系', harmonic ? '<span class="pill ok">整数倍：周期</span>' : '<span class="pill warn">非整数倍：不周期</span>'],
      ['max |u|', fx(umax, 4)],
      ['非振荡模态', sl.nNeg ? sl.nNeg + ' 个按 cosh 增长' : sl.modes[0].kind === 'zero' ? '常数模态（不振动）' : '无']
    ];
  }
  P.stats($('sv-tStats'), st);
  const nb = Math.min(20, sl.modes.length), cmax = Math.max(...Array.from(coef).slice(0, nb).map(Math.abs)) || 1;
  const ts = Array.from({ length: nb }, (_, j) => Math.abs(coef[j] * timeFactor(sl.modes[j], t)));
  const ymax = cmax * 1.12;
  const g = P.frame($('sv-cvM'), { xmin: sl.offset - 0.8, xmax: nb + sl.offset - 0.2, ymin: 0, ymax, xlabel: 'n', xticks: P.niceTicks(sl.offset, nb - 1 + sl.offset, 5).filter(v => v % 1 === 0).map(v => ({ v, label: v + '' })) });
  if (g) {
    const bw = Math.max(3, g.pw / nb * 0.55), xsN = ts.map((_, j) => nLabel(j));
    P.bars(g, xsN, Array.from(coef).slice(0, nb).map(Math.abs), P.C['ink-3'], bw, 0, true);
    P.bars(g, xsN, ts.map(v => Math.min(v, ymax)), j => j < S.N ? P.C.accent : P.C['ink-3'], bw);
  }
  P.tex($('sv-tFormula'), S.eq === 'heat'
    ? String.raw`\(T_n(t)=e^{-a^2\lambda_n t}\)：第 \(n\) 项衰减时间 \(1/(a^2\lambda_n)\propto 1/n^2\)。时间轴前段加密。`
    : String.raw`\(T_n=\cos(ak_nt)\)（\(\lambda_n\gt0\)）；\(\lambda=0\) 时 \(T=1\)；\(\lambda=-\mu^2\) 时 \(T=\cosh(a\mu t)\)。`);
}

/* ---------- 基准 ---------- */
function runBench() {
  const rows = [];
  let s = P.solveSL({ L: 1, bcL: 'D', bcR: 'D', NM: 20, NG: 2400 }), e = 0;
  s.modes.forEach((m, j) => { e = Math.max(e, Math.abs(m.lam - ((j + 1) * PI) ** 2) / ((j + 1) * PI) ** 2); });
  rows.push(P.benchRow('D–D：λ_n = (nπ)²，前 20 个', '闭式解', e, 1e-12));
  s = P.solveSL({ L: 1, bcL: 'D', bcR: 'R', h2: 1, NM: 4, NG: 2400 });
  rows.push(P.benchRow('D–R, h=1：tan k = −k 首根', '文献值 2.028757838110434', Math.abs(s.modes[0].k - 2.028757838110434), 1e-12));
  s = P.solveSL({ L: 1, bcL: 'R', bcR: 'R', h1: 1, h2: 1, NM: 4, NG: 2400 });
  rows.push(P.benchRow('R–R, h=1：tan k = 2k/(k²−1) 首根', '文献值 1.306542374188806', Math.abs(s.modes[0].k - 1.306542374188806), 1e-12));
  s = P.solveSL({ L: 1, bcL: 'N', bcR: 'N', NM: 6, NG: 2400 });
  rows.push(P.benchRow('N–N：λ₀ = 0 且 X₀ 为常数', '检查 λ₀ 与 X₀ 的变化量', Math.abs(s.modes[0].lam) + (P.maxAbs(s.modes[0].X) - Math.min(...s.modes[0].X)), 1e-14));
  // D–R, h2 = −2：负特征值 μ = 2 tanh μ，用独立 Newton 迭代对照
  s = P.solveSL({ L: 1, bcL: 'D', bcR: 'R', h2: -2, NM: 4, NG: 2400 });
  let mu = 2; for (let i = 0; i < 50; i++) mu -= (mu - 2 * Math.tanh(mu)) / (1 - 2 / Math.cosh(mu) ** 2);
  rows.push(P.benchRow('D–R, h=−2：负特征值 μ = 2 tanh μ', 'Newton 迭代 μ=' + mu.toFixed(10), s.modes[0].kind === 'neg' ? Math.abs(s.modes[0].k - mu) : Infinity, 1e-11));
  rows.push(P.benchRow('上例负模态的 Rayleigh 商（含负边界项）', 'Rayleigh 商 vs λ = −μ²', Math.abs(s.modes[0].ray - s.modes[0].lam) / Math.abs(s.modes[0].lam), 1e-8));
  s = P.solveSL({ L: 1, bcL: 'R', bcR: 'R', h1: -0.5, h2: 1, NM: 4, NG: 2400 });
  const zm = s.modes.find(m => m.kind === 'zero');
  rows.push(P.benchRow('R–R, h₁=−½, h₂=1：零特征值 X=1−x/2', '条件 h₁+h₂+h₁h₂L=0', zm ? Math.abs(zm.ray) : Infinity, 1e-9));
  s = P.solveSL({ L: 1.7, bcL: 'R', bcR: 'N', h1: 3, NM: 8, NG: 2400 });
  let go = 0; for (let a = 0; a < 8; a++) for (let b = 0; b < a; b++) go = Math.max(go, Math.abs(P.dotS(s.modes[a].X, s.modes[b].X, s.h)) / Math.sqrt(s.modes[a].n2 * s.modes[b].n2));
  rows.push(P.benchRow('R–N 正交性：前 8 个模态的 Gram 矩阵', '最大非对角元', go, 1e-9));
  s = P.solveSL({ L: 1, bcL: 'D', bcR: 'D', NM: 60, NG: 2400 });
  const pv = new Float64Array(2401); for (let i = 0; i <= 2400; i++) pv[i] = 4 * (i / 2400) * (1 - i / 2400);
  const cp = P.project(s, pv); let ce = 0; cp.forEach((c, j) => { ce = Math.max(ce, Math.abs(c - analyticDD(j + 1, 'para'))); });
  rows.push(P.benchRow('抛物线 4x(1−x) 的正弦系数（60 项）', '闭式 32/(nπ)³（n 奇）', ce, 1e-10));
  const one = new Float64Array(2401).fill(1), c1 = P.project(s, one), Sn = new Float64Array(2401);
  c1.forEach((c, j) => { const X = s.modes[j].X; for (let i = 0; i <= 2400; i++) Sn[i] += c * X[i]; });
  rows.push(P.benchRow('Gibbs：f≡1 (D–D)，N=60 的过冲', '理论极限 0.17898（跳跃 2 × 8.949%）', Math.abs((Math.max(...Sn) - 1) - 0.178979744), 5e-3));
  P.benchTable($('sv-bench'), 'separation', rows);
}

/* ---------- 题目 ---------- */
const QUIZ = [
  { type: '推导辨析', q: String.raw`两端 Dirichlet 时，为什么不需要考虑 \(\lambda\lt0\)？`, o: ['因为负数没有平方根', String.raw`代入边界条件后只剩零解 \(X\equiv0\)`, '因为物理上时间不能倒流', '因为分离常数规定必须为正'], a: 1,
    e: String.raw`\(\lambda=-\mu^2\) 时 \(X=A\cosh\mu x+B\sinh\mu x\)，两端为零迫使 \(A=B=0\)。更一般地，Rayleigh 商 \(\lambda=\big(\int X'^2+\text{边界项}\big)/\int X^2\)；边界项非负时 \(\lambda\ge0\)。` },
  { type: '收敛与相容', q: String.raw`两端 D，把 \(f\equiv1\) 展开成 \(\sum b_n\sin\frac{n\pi x}{L}\)。部分和在 \(x=0\) 处的值是？`, o: ['1', '0', '1/2', '随 N 增大而发散'], a: 1,
    e: String.raw`每一项 \(\sin 0=0\)，\(S_N(0)=0\)。\(f\) 与边界条件不相容，端点附近出现 Gibbs 过冲（约为奇延拓跳跃 2 的 8.95%，即 0.179），过冲不随 \(N\) 消失，只是变窄。` },
  { type: '长期行为', q: '两端绝热（N–N）的杆，t → ∞ 时温度趋于？', o: ['0', '初始温度的平均值', '初始温度的最大值', '取决于 a 的大小'], a: 1,
    e: String.raw`N–N 有 \(\lambda_0=0\)、\(X_0=1\)，\(T_0\equiv1\) 不衰减，\(c_0=\frac1L\int_0^Lf\,dx\)。绝热意味着总热量守恒。` },
  { type: '边界条件', q: String.raw`右端 Robin 条件 \(X'(L)+hX(L)=0\)，当 \(h\to\infty\) 和 \(h\to0\) 时分别退化为？`, o: ['都退化为 Dirichlet', 'h→∞ 为 Neumann，h→0 为 Dirichlet', 'h→∞ 为 Dirichlet，h→0 为 Neumann', '都退化为 Neumann'], a: 2,
    e: String.raw`除以 \(h\) 得 \(X'(L)/h+X(L)=0\)，\(h\to\infty\) 时就是 \(X(L)=0\)；\(h=0\) 时就是 \(X'(L)=0\)。` },
  { type: '非齐次边界', q: String.raw`\(u_t=a^2u_{xx}\)，\(u(0,t)=100\)，\(u(L,t)=0\)，能否直接分离变量？`, o: ['可以，直接设 u = X T', '不能，边界非齐次；先减去稳态解 w = 100(1 − x/L)', '不能，热传导方程不可分离', '可以，但只能求数值解'], a: 1,
    e: String.raw`若 \(u=XT\)，则 \(X(0)T(t)=100\) 要对所有 \(t\) 成立，与 \(T\) 随时间变化矛盾。令 \(u=w(x)+v\)，\(w''=0\) 满足原边界，\(v\) 的边界齐次。模块 03 有交互演示。` },
  { type: '系数与光滑性', q: String.raw`两端 D 时，抛物线系数按 \(n^{-3}\)、三角形按 \(n^{-2}\)、阶跃按 \(n^{-1}\) 衰减。这反映了什么？`, o: ['区间长度不同', 'f 的奇延拓的光滑程度', '初始条件的最大值不同', '数值误差的大小'], a: 1,
    e: String.raw`分部积分一次多一个 \(1/n\) 因子，直到遇到不连续为止。跳跃 → \(n^{-1}\)，尖角 → \(n^{-2}\)，导数连续而二阶导跳跃 → \(n^{-3}\)。` },
  { type: '特征值符号', q: String.raw`左 D 右 R，若 \(h_2=-2\)、\(L=1\)，最小特征值的符号是？`, o: ['正', '零', '负', '无法判断'], a: 2,
    e: String.raw`\(h_2\lt0\) 时 Rayleigh 商分子中出现 \(h_2X(L)^2\lt0\)。当 \(h_2L\lt-1\) 时出现负特征值，\(\lambda=-\mu^2\)，\(\mu=2\tanh\mu\Rightarrow\mu\approx1.915\)。物理上对应边界向杆内"放热"，温度会指数增长。` }
];

/* ---------- 页面 ---------- */
function html() {
  const R = String.raw;
  return P.modHead(ID, '04', '分离变量法与特征值问题', '有界区间上的弦振动和杆的热传导，都可以拆成"空间形状 × 时间因子"的叠加。本模块从三类边界条件出发求特征值和特征函数，再用它们展开初始条件、推进时间。实验 A–C 的计算核心也被模块 02、03 共用。') +
  P.sec(ID, 'phys', '两个典型的有界问题', `<div class="cols">
    <div class="card"><h3>弦的横振动</h3><p class="muted small">长 \\(L\\) 的弦，张力 \\(T\\)，线密度 \\(\\rho\\)。小振幅下横向位移 \\(u(x,t)\\) 满足波动方程。端点可以固定、套在光滑竖杆上自由滑动，或接弹簧。</p></div>
    <div class="card"><h3>细杆的热传导</h3><p class="muted small">长 \\(L\\) 的均匀细杆，侧面绝热，温度 \\(u(x,t)\\)（相对环境温度）满足热传导方程。端点可以恒温、绝热，或与空气对流换热。</p></div>
    <div class="card"><h3>共同点</h3><p class="muted small">方程线性齐次，区域是 \\([0,L]\\times[0,\\infty)\\)，两端边界条件线性齐次。这正是分离变量法能直接使用的条件。</p></div></div>`) +
  P.sec(ID, 'model', '数学模型', R`<div class="cols">
    <div class="card"><h3>热传导方程（抛物型）</h3><div class="formula-box">\[ u_t = a^2 u_{xx},\qquad a^2=\frac{\kappa}{c\rho} \]</div><p class="muted small">\(a^2\) 是热扩散率（钢约 \(1.2\times10^{-5}\,\mathrm{m^2/s}\)）。</p></div>
    <div class="card"><h3>波动方程（双曲型）</h3><div class="formula-box">\[ u_{tt} = a^2 u_{xx},\qquad a=\sqrt{T/\rho} \]</div><p class="muted small">\(a\) 是波速。实验中 \(a\)、\(L\) 均为无量纲量。</p></div></div>`) +
  P.sec(ID, 'cond', '三类边界条件与初始条件', R`<div class="tbl-wrap"><table>
    <thead><tr><th>类型</th><th>左端 \(x=0\)</th><th>右端 \(x=L\)</th><th>弦振动</th><th>杆热传导</th></tr></thead><tbody>
    <tr><td>第一类 Dirichlet (D)</td><td>\(u=0\)</td><td>\(u=0\)</td><td>端点固定</td><td>端点保持环境温度</td></tr>
    <tr><td>第二类 Neumann (N)</td><td>\(u_x=0\)</td><td>\(u_x=0\)</td><td>端点在光滑竖杆上自由滑动</td><td>端点绝热</td></tr>
    <tr><td>第三类 Robin (R)</td><td>\(u_x-h_1u=0\)</td><td>\(u_x+h_2u=0\)</td><td>端点接弹簧，\(h=k_s/T\)</td><td>牛顿冷却，\(h=H/\kappa\)</td></tr></tbody></table></div>
    <p class="prose muted small">Robin 条件两端符号不同，因为外法向不同：统一写成 \(-\kappa\,\partial_n u = H u\)。物理上 \(h\ge0\)（热量流向较冷的环境）；实验允许 \(h\lt0\) 以观察负特征值。初始条件：热传导给 \(u(x,0)=f(x)\)；波动再给 \(u_t(x,0)=g(x)\)，本模块取 \(g=0\)。</p>`) +
  P.sec(ID, 'idea', '求解思路', R`<div class="prose small muted"><p>线性齐次方程 + 齐次边界条件 ⇒ 满足方程和边界条件的解组成线性空间。先找最简单的"驻波型"解 \(X(x)T(t)\)，再叠加它们去匹配初始条件。边界条件全部落到 \(X\) 上，得到特征值问题；初始条件只用来确定叠加系数。</p></div>`) +
  P.sec(ID, 'deriv', '分离变量法的逐步推导', `<div id="sv-steps"></div>
    <div class="note"><b>适用条件</b>　(1) 方程线性齐次且可分离（这里是常系数）；(2) 两端边界条件线性齐次；(3) 区域是"区间 × 时间"的乘积区域。边界非齐次时先减去稳态解；有源项时用特征函数展开法或 Duhamel 原理（模块 02、03）。</div>`) +
  P.sec(ID, 'viz', '实验：特征值、展开与演化',
    `<h3 class="sub">实验 A　特征值与特征函数</h3>` +
    P.lab('EIGEN', '改变边界条件，观察 \\(\\lambda_n\\) 与 \\(X_n\\)',
      `<div class="ctl"><div class="ctl-label">左端 \\(x=0\\)</div><div id="sv-bcL"></div></div>
       <div class="ctl"><div class="ctl-label">右端 \\(x=L\\)</div><div id="sv-bcR"></div></div>
       ${P.rangeHTML('sv-h1', '左端 Robin 系数 \\(h_1\\)', -5, 20, 0.1, 2)}
       ${P.rangeHTML('sv-h2', '右端 Robin 系数 \\(h_2\\)', -5, 20, 0.1, 2)}
       ${P.rangeHTML('sv-L', '区间长度 \\(L\\)', 0.5, 3, 0.05, 1)}
       <div class="ctl"><div class="ctl-label">查看第几个特征函数</div><div id="sv-sel"></div></div>
       <div class="tasks"><b>试一试</b><ul>
         <li>左 D 右 R，把 \\(h_2\\) 从 0 拉到 20：\\(k_1L/\\pi\\) 从 0.5 走向 1。</li>
         <li>把 \\(h_2\\) 拉到 −1.5：\\(h_2L\\lt-1\\) 时出现负特征值，X₁ 变成双曲函数形状。</li>
         <li>R–R 取 \\(h_1=-0.5,\\ h_2=1\\)：出现零特征值，\\(X=1-x/2\\)。</li>
       </ul></div>`,
      `<div class="problem" id="sv-sys"></div>
       <p class="small" id="sv-signnote"></p>
       ${P.plotBlock('sv-cvEig', '特征函数 \\(X_n(x)\\)（按最大值归一化显示）', P.lg('accent', '选中模态') + P.lg('ink-3', '前 6 个模态'))}
       <div class="stats" id="sv-eigStats"></div>
       <div class="plot-block"><div class="plot-cap"><span>正特征值方程的图解</span><span class="legend" id="sv-tanLegend"></span></div><div class="formula-box" id="sv-eq"></div><canvas id="sv-cvTan" class="short"></canvas></div>
       <div class="two">
         <div class="plot-block"><div class="plot-cap"><span>前 8 个特征值与 Rayleigh 商校验</span></div><div class="tbl-wrap"><table id="sv-eigTable"></table></div></div>
         <div class="plot-block"><div class="plot-cap"><span>归一化 Gram 矩阵 \\(\\langle X_m,X_n\\rangle/\\|X_m\\|\\|X_n\\|\\)</span></div><div class="tbl-wrap"><table class="gram" id="sv-gram"></table></div><p class="small muted" id="sv-gramNote"></p></div>
       </div>
       <div class="small muted exp-only"><p>数值方法：把边界条件统一写成 \\(a_1X(0)-b_1X'(0)=0\\)、\\(a_2X(L)+b_2X'(L)=0\\)。\\(\\lambda=k^2\\gt0\\) 时以步长 \\(\\pi/(40L)\\) 扫描下式的变号区间再二分；\\(\\lambda=-\\mu^2\\lt0\\) 时扫描 \\(G\\)：</p>
       <div class="formula-box">\\[ \\begin{aligned} F(k)&=(a_2b_1+a_1b_2)k\\cos kL+(a_1a_2-b_1b_2k^2)\\sin kL,\\\\ G(\\mu)&=(a_2b_1+a_1b_2)\\mu+(a_1a_2+b_1b_2\\mu^2)\\tanh\\mu L . \\end{aligned} \\]</div>
       <p>\\(\\lambda=0\\) 当且仅当 \\(a_2b_1+a_1a_2L+a_1b_2=0\\)。积分用 2400 段复化 Simpson 公式。</p></div>`) +
    `<h3 class="sub">实验 B　用特征函数展开初始条件</h3>` +
    P.lab('EXPAND', '\\(f(x)\\approx S_N(x)=\\sum_{n}^{N}c_nX_n(x)\\)，边界条件沿用实验 A',
      `<div class="ctl"><label class="ctl-label" for="sv-f">初始条件 \\(f(x)\\)</label><select id="sv-f">
         <option value="tri">三角形（拨弦）</option><option value="para">抛物线 4ξ(1−ξ)</option><option value="step">阶跃：左半段为 1</option>
         <option value="const">常数 1</option><option value="gauss">高斯脉冲</option><option value="mode">恰好是第 2 个特征函数</option><option value="draw">手绘（在图上拖动）</option></select></div>
       ${P.rangeHTML('sv-p', '拨弦位置 \\(p=x_0/L\\)', 0.05, 0.95, 0.05, 0.3)}
       ${P.rangeHTML('sv-N', '截断阶数 \\(N\\)', 1, 60, 1, 10)}
       <div class="tasks"><b>试一试</b><ul>
         <li>常数 1 + 两端 D：端点处 \\(S_N=0\\neq f\\)，把 \\(N\\) 加到 60，过冲不消失（Gibbs 现象），只是变窄。</li>
         <li>同一个常数 1 换成两端 N：只需要 1 项。</li>
         <li>比较三角形、抛物线、阶跃的系数衰减指数 \\(p\\)：约为 2、3、1。</li></ul></div>`,
      `${P.plotBlock('sv-cvF', '初始条件与部分和', P.lg('warm', 'f(x)') + P.lg('accent', 'S_N(x)'))}
       <p class="small" id="sv-compat"></p>
       <div class="stats" id="sv-fStats"></div>
       ${P.plotBlock('sv-cvC', '系数 \\(|c_n|\\)（对数坐标），虚线为参考衰减 \\(n^{-1},n^{-2},n^{-3}\\)', P.lg('accent', '已保留') + P.lg('ink-3', '被截断'), 'short')}
       <p class="small muted" id="sv-fCheck"></p>`) +
    `<div class="note"><b>收敛定理</b>　\\(f\\in L^2\\)：\\(S_N\\to f\\) 均方收敛。\\(f\\) 分段光滑：每点收敛到 \\(\\tfrac12[f(x^-)+f(x^+)]\\)，D 端收敛到 0。\\(f\\) 连续、分段光滑且满足 D 端条件：一致收敛，无 Gibbs 过冲。N、R 端的导数条件不影响一致收敛，只影响收敛速度。</div>` +
    `<h3 class="sub">实验 C　每个模态各自随时间演化</h3>` +
    P.lab('EVOLVE', '\\(u(x,t)=\\sum_{n}^{N}c_n\\,T_n(t)\\,X_n(x)\\)，\\(N\\)、\\(f\\) 与边界条件沿用实验 A、B',
      `<div class="ctl"><div class="ctl-label">方程</div><div id="sv-eqSeg"></div></div>
       ${P.rangeHTML('sv-a', '\\(a\\)（扩散率平方根 / 波速）', 0.2, 2, 0.05, 1)}
       <div id="sv-player"></div>
       <div class="formula-box small" id="sv-tFormula"></div>
       <div class="tasks"><b>试一试</b><ul>
         <li>热传导 + 阶跃：极短时间后高频模态全部消失，曲线立刻变光滑。</li>
         <li>热传导两端 N：温度趋于初始平均值，\\(\\int u\\,dx\\) 不变。</li>
         <li>波动两端 D：\\(k_n/k_1\\) 是整数，运动周期为 \\(2L/a\\)；换成 R–R 不再周期。</li></ul></div>`,
      `${P.plotBlock('sv-cvU', '\\(u(x,t)\\)', P.lg('accent', 'u(x,t)') + P.lg('warm', 'f(x)'))}
       <div class="stats" id="sv-tStats"></div>
       ${P.plotBlock('sv-cvM', '各模态振幅 \\(|c_nT_n(t)|\\)（空心为 \\(t=0\\) 的 \\(|c_n|\\)）', '', 'short')}`)) +
  P.sec(ID, 'meaning', '每个数学对象对应什么', R`<div class="tbl-wrap"><table>
    <thead><tr><th>数学对象</th><th>弦振动</th><th>杆热传导</th></tr></thead><tbody>
    <tr><td>特征函数 \(X_n\)</td><td>第 \(n\) 阶振型（驻波形状）</td><td>第 \(n\) 个温度分布模式</td></tr>
    <tr><td>特征值 \(\lambda_n=k_n^2\)</td><td>角频率 \(\omega_n=a k_n\)</td><td>衰减率 \(a^2\lambda_n\)</td></tr>
    <tr><td>最小正特征值</td><td>基频，决定音高</td><td>最慢衰减，决定冷却多久</td></tr>
    <tr><td>\(\lambda_0=0\)（N–N）</td><td>整体平移，不振动</td><td>平均温度，总热量守恒</td></tr>
    <tr><td>\(\lambda\lt0\)（\(h\lt0\)）</td><td>负刚度弹簧，失稳</td><td>边界向杆内供热，温度指数增长</td></tr>
    <tr><td>零点数 = 序号 − 1</td><td>驻波节点</td><td>温度符号变化处</td></tr>
    <tr><td>正交性</td><td>各阶振动能量互不交换</td><td>各模式独立衰减</td></tr>
    <tr><td>系数 \(c_n\)</td><td>初始形状激发第 \(n\) 阶的程度（音色）</td><td>初始温度中第 \(n\) 模式的含量</td></tr></tbody></table></div>` +
    P.pitfalls([
      ['分离常数的符号：', String.raw`写成 \(-\lambda\) 是为了让 \(X''+\lambda X=0\) 的特征值为正；写成 \(+\lambda\) 也对，但结论中符号全部反过来。`],
      ['零特征值不能漏：', String.raw`N–N 的 \(\lambda_0=0\)、\(X_0=1\) 是真实模态；D–D 时 \(\lambda=0\) 只有平凡解 \(X\equiv0\)，不是特征值。`],
      ['三类边界条件：', String.raw`D 约束函数值，N 约束导数（通量），R 是两者的线性组合；左右端 Robin 的符号因外法向不同而相反。`],
      ['系数公式的分母：', String.raw`\(c_n=\int fX_n/\int X_n^2\)。只有 D–D 正弦级数才有 \(\int X_n^2=L/2\)；N–N 的 \(n=0\) 项分母是 \(L\)。`],
      ['相容性：', String.raw`\(f\) 不满足 D 端条件时，级数在端点收敛到 0 而不是 \(f\)，并出现 Gibbs 现象。`]
    ]) +
    `<div class="bench exp-only"><h3 class="sub">数值基准</h3><div id="sv-bench"></div></div>`) +
  P.sec(ID, 'quiz', '检查你的理解 <span class="score" id="sv-score"></span>', '<div id="sv-quiz"></div>');
}
const STEPS = [
  { t: '假设变量分离形式的解', b: String.raw`<p>寻找非零解 \(u(x,t)=X(x)T(t)\)，代入热传导方程：</p><div class="formula-box">\[ X\,T' = a^2X''\,T \;\Longrightarrow\; \frac{T'}{a^2T}=\frac{X''}{X}=-\lambda . \]</div><p>左边只依赖 \(t\)，右边只依赖 \(x\)，两边只能等于同一常数。波动方程同理得 \(T''/(a^2T)=X''/X=-\lambda\)。</p>` },
  { t: '边界条件落到 X 上：特征值问题', b: String.raw`<p>\(u(0,t)=X(0)T(t)=0\) 对所有 \(t\) 成立且 \(T\not\equiv0\)，所以 \(X(0)=0\)。得到 Sturm–Liouville 问题：</p><div class="formula-box">\[ X''+\lambda X=0,\quad 0\lt x\lt L,\qquad \text{加两端边界条件}. \]</div><p>只有特定的 \(\lambda\) 才有非零解，这一步是整个方法的核心。</p>` },
  { t: '分情况讨论 λ：以 D–D 为例', b: String.raw`<div class="formula-box">\[ \begin{aligned} \lambda=-\mu^2\lt0:&\quad X=A\cosh\mu x+B\sinh\mu x,\; X(0)=0\Rightarrow A=0,\; B\sinh\mu L=0\Rightarrow B=0;\\ \lambda=0:&\quad X=Ax+B\Rightarrow A=B=0;\\ \lambda=k^2\gt0:&\quad X=A\cos kx+B\sin kx,\; A=0,\; B\sin kL=0\Rightarrow k_n=\tfrac{n\pi}{L}. \end{aligned} \]</div>` },
  { t: 'Rayleigh 商：λ 的符号由边界项决定', b: String.raw`<p>方程乘 \(X\) 后分部积分，代入边界条件：</p><div class="formula-box">\[ \lambda=\frac{\int_0^L X'^2\,dx + h_1X(0)^2 + h_2X(L)^2}{\int_0^L X^2\,dx} \]</div><p>D、N 端边界项为零，R 端贡献 \(hX^2\)。若 \(h_1,h_2\ge0\)，则 \(\lambda\ge0\)，等号只在 N–N（或 \(h=0\)）时由常数函数取到。若某个 \(h\lt0\)，分子可以为负，于是可能出现负特征值，必须单独讨论 \(\lambda\lt0\) 的情形。</p>` },
  { t: '正交性', b: String.raw`<div class="formula-box">\[ (\lambda_m-\lambda_n)\int_0^L X_mX_n\,dx=\Big[X_m'X_n-X_mX_n'\Big]_0^L=0 . \]</div><p>三类齐次边界条件（任意符号的 \(h\)）都使边界项为零，所以不同特征值的特征函数正交。</p>` },
  { t: '叠加与系数', b: String.raw`<div class="formula-box">\[ u(x,t)=\sum_n c_n\,T_n(t)\,X_n(x),\qquad c_n=\frac{\int_0^L f\,X_n\,dx}{\int_0^L X_n^2\,dx}, \]</div><div class="formula-box">\[ \text{热传导: } T_n=e^{-a^2\lambda_n t},\qquad \text{波动 } (g=0)\text{: } T_n=\cos\!\big(a\sqrt{\lambda_n}\,t\big). \]</div><p>Sturm–Liouville 定理保证 \(\{X_n\}\) 在 \(L^2(0,L)\) 中完备，级数按均方收敛。</p>` },
  { t: '物理意义回顾', b: String.raw`<p>\(X_n\) 是振型或温度模式，\(\lambda_n\) 决定频率或衰减率。高阶模态 \(\lambda_n\sim n^2\)：热传导中衰减得更快，波动中振得更快。</p>` }
];

/* ---------- 初始化 ---------- */
function allBC() { updateEigen(); updateCoeffs(); buildSel(); syncRows(); renderEigen(); renderFourier(); player.refresh(); }
function allF() { updateCoeffs(); syncRows(); renderFourier(); player.refresh(); }
let selSeg = null;
function buildSel() {
  const opts = Array.from({ length: 6 }, (_, j) => ({ v: j, label: nLabel(j) }));
  selSeg = P.segment($('sv-sel'), opts, S.sel, v => { S.sel = +v; renderEigen(); });
}
function syncRows() {
  $('sv-h1-row').hidden = S.bcL !== 'R'; $('sv-h2-row').hidden = S.bcR !== 'R'; $('sv-p-row').hidden = S.fType !== 'tri';
}
const BC_OPTS = [{ v: 'D', label: 'D 固定' }, { v: 'N', label: 'N 自由' }, { v: 'R', label: 'R 弹性' }];
P.register('separation', {
  init(root) {
    root.innerHTML = html();
    P.stepper($('sv-steps'), STEPS);
    P.segment($('sv-bcL'), BC_OPTS, S.bcL, v => { S.bcL = v; allBC(); });
    P.segment($('sv-bcR'), BC_OPTS, S.bcR, v => { S.bcR = v; allBC(); });
    P.segment($('sv-eqSeg'), [{ v: 'heat', label: '热传导' }, { v: 'wave', label: '波动 (g=0)' }], S.eq, v => { S.eq = v; player.set(0); });
    P.bindRange('sv-h1', v => { S.h1 = v; allBC(); }, v => v.toFixed(1));
    P.bindRange('sv-h2', v => { S.h2 = v; allBC(); }, v => v.toFixed(1));
    P.bindRange('sv-L', v => { S.L = v; allBC(); }, v => v.toFixed(2));
    P.bindRange('sv-p', v => { S.p = v; allF(); }, v => v.toFixed(2));
    P.bindRange('sv-N', v => { S.N = v; renderFourier(); player.refresh(); });
    P.bindRange('sv-a', v => { S.a = v; player.refresh(); }, v => v.toFixed(2));
    $('sv-f').addEventListener('change', e => {
      const v = e.target.value;
      if (v === 'draw') for (let i = 0; i <= 120; i++) drawPts[i] = Math.max(-1.2, Math.min(1.2, fv[Math.round(i * sl.NG / 120)]));
      S.fType = v; allF();
    });
    // 手绘
    const cv = $('sv-cvF'); let last = null, pending = false;
    const at = e => { const r = cv.getBoundingClientRect(), inv = cv._inv; return [Math.round(Math.max(0, Math.min(1, inv.x(e.clientX - r.left) / S.L)) * 120), Math.max(-1.2, Math.min(1.2, inv.y(e.clientY - r.top)))]; };
    const paint = e => {
      const [i, y] = at(e);
      if (last) { const [i0, y0] = last, d = Math.abs(i - i0); for (let k = 0; k <= d; k++) { const t = d ? k / d : 1; drawPts[Math.round(i0 + (i - i0) * t)] = y0 + (y - y0) * t; } } else drawPts[i] = y;
      last = [i, y];
      if (!pending) { pending = true; requestAnimationFrame(() => { pending = false; allF(); }); }
    };
    cv.addEventListener('pointerdown', e => { if (S.fType !== 'draw' || !cv._inv) return; cv.setPointerCapture(e.pointerId); last = null; paint(e); });
    cv.addEventListener('pointermove', e => { if (S.fType === 'draw' && cv.hasPointerCapture(e.pointerId)) paint(e); });
    cv.addEventListener('pointerup', () => { last = null; });
    player = P.player($('sv-player'), { id: 'sv-t', tmax: tMax, map: s => S.eq === 'heat' ? tMax() * s * s : tMax() * s, digits: 4, duration: 9, onChange: renderTime });
    updateEigen(); updateCoeffs(); buildSel(); syncRows();
    P.quiz($('sv-quiz'), QUIZ, $('sv-score'));
    runBench();
    $('sv-bench').addEventListener('click', e => { if (e.target.closest('[data-rerun]')) runBench(); });
  },
  render() { renderEigen(); renderFourier(); player.refresh(); }
});
})();
