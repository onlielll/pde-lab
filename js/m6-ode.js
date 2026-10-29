/* m6-ode.js — 模块 06：二阶常系数 ODE、待定系数法与常数变易法、与 PDE 分离变量的联系 */
(function () {
'use strict';
const P = window.PDE, PI = Math.PI, $ = P.$, fx = P.fx, sci = P.sci;
const ID = 'od', TOL = 1e-9, T = 10;
const S = { p: 1, q: 9, y0: 1, v0: 0, g: 'sin', A: 1, beta: -0.5, w: 2, lam: 12, Ls: 1, sbc: 'D' };
const PRE = {
  under: { label: '欠阻尼 p=1, q=9', p: 1, q: 9 }, crit: { label: '临界阻尼 p=4, q=4', p: 4, q: 4 },
  over: { label: '过阻尼 p=5, q=4', p: 5, q: 4 }, free: { label: '无阻尼 p=0, q=4', p: 0, q: 4 },
  saddle: { label: '不稳定 p=0, q=−4', p: 0, q: -4 }, ex: { label: '例题 p=2, q=5', p: 2, q: 5 }
};

/* ---------- 齐次方程 y'' + p y' + q y = 0 ---------- */
function basis(p, q) {
  const D = p * p - 4 * q;
  if (D > TOL) {
    const s = Math.sqrt(D), r1 = (-p + s) / 2, r2 = (-p - s) / 2;
    return { kind: 'real', D, r1, r2, y1: t => [Math.exp(r1 * t), r1 * Math.exp(r1 * t)], y2: t => [Math.exp(r2 * t), r2 * Math.exp(r2 * t)] };
  }
  if (D >= -TOL) {
    const r = -p / 2;
    return { kind: 'double', D, r1: r, r2: r, y1: t => [Math.exp(r * t), r * Math.exp(r * t)], y2: t => [t * Math.exp(r * t), (1 + r * t) * Math.exp(r * t)] };
  }
  const a = -p / 2, b = Math.sqrt(-D) / 2;
  return { kind: 'complex', D, a, b, y1: t => { const e = Math.exp(a * t); return [e * Math.cos(b * t), e * (a * Math.cos(b * t) - b * Math.sin(b * t))]; }, y2: t => { const e = Math.exp(a * t); return [e * Math.sin(b * t), e * (a * Math.sin(b * t) + b * Math.cos(b * t))]; } };
}
function constants(B, y0, v0) {
  const [a1, b1] = B.y1(0), [a2, b2] = B.y2(0), det = a1 * b2 - a2 * b1;
  return [(y0 * b2 - a2 * v0) / det, (a1 * v0 - b1 * y0) / det];
}
function solH(B, C1, C2, t) { const u = B.y1(t), v = B.y2(t); return [C1 * u[0] + C2 * v[0], C1 * u[1] + C2 * v[1]]; }
function classify(p, q, B) {
  if (q < 0) return ['一正一负实根：不稳定（鞍点）', 'bad'];
  if (Math.abs(q) < TOL) return ['有零根：解含常数项', 'info'];
  if (Math.abs(p) < TOL) return ['纯虚根：无阻尼简谐振动', 'info'];
  if (p < 0) return ['负阻尼：振幅增长', 'bad'];
  return B.kind === 'complex' ? ['共轭复根：欠阻尼振动', 'ok'] : B.kind === 'double' ? ['重根：临界阻尼', 'warn'] : ['两个负实根：过阻尼', 'warn'];
}

/* ---------- 非齐次：g(t) ---------- */
function gFn(t) {
  switch (S.g) {
    case 'const': return S.A;
    case 'exp': return S.A * Math.exp(S.beta * t);
    case 'sin': return S.A * Math.sin(S.w * t);
    case 'poly': return S.A * t;
    case 'cauchy': return S.A / (1 + t * t);
  }
  return 0;
}
// 待定系数法的特解及其导数；不适用时返回 null
function ucParticular(p, q) {
  const A = S.A, z = v => Math.abs(v) < TOL;
  switch (S.g) {
    case 'const':
      if (!z(q)) return { tex: String.raw`y_p=\frac{A}{q}`, f: () => [A / q, 0] };
      if (!z(p)) return { tex: String.raw`y_p=\frac{A}{p}\,t\quad(q=0\text{，常数是特征根，乘 }t)`, f: t => [A * t / p, A / p] };
      return { tex: String.raw`y_p=\frac A2t^2`, f: t => [A * t * t / 2, A * t] };
    case 'exp': {
      const b = S.beta, Pb = b * b + p * b + q, dP = 2 * b + p;
      if (!z(Pb)) return { tex: String.raw`y_p=\frac{A\,e^{\beta t}}{\beta^2+p\beta+q}`, f: t => { const e = A * Math.exp(b * t) / Pb; return [e, b * e]; } };
      if (!z(dP)) return { tex: String.raw`y_p=\frac{A\,t\,e^{\beta t}}{2\beta+p}\quad(\beta\text{ 是单根，乘 }t)`, f: t => { const e = A * Math.exp(b * t) / dP; return [t * e, (1 + b * t) * e]; } };
      return { tex: String.raw`y_p=\frac A2t^2e^{\beta t}\quad(\beta\text{ 是重根，乘 }t^2)`, f: t => { const e = A * Math.exp(b * t) / 2; return [t * t * e, (2 * t + b * t * t) * e]; } };
    }
    case 'sin': {
      const w = S.w, a = q - w * w, c = p * w, det = a * a + c * c;
      if (!z(det)) { const Cc = -c * A / det, Dd = a * A / det; return { tex: String.raw`y_p=C\cos\omega t+D\sin\omega t,\ \ C=${fx(Cc, 4)},\ D=${fx(Dd, 4)}`, f: t => [Cc * Math.cos(w * t) + Dd * Math.sin(w * t), w * (-Cc * Math.sin(w * t) + Dd * Math.cos(w * t))] }; }
      if (z(w)) return { tex: 'y_p=0', f: () => [0, 0] };
      return { tex: String.raw`y_p=-\frac{A\,t\cos\omega t}{2\omega}\quad(\text{共振：}i\omega\text{ 是特征根})`, f: t => [-A * t * Math.cos(w * t) / (2 * w), -A * (Math.cos(w * t) - w * t * Math.sin(w * t)) / (2 * w)] };
    }
    case 'poly':
      if (!z(q)) return { tex: String.raw`y_p=\frac Aq\Big(t-\frac pq\Big)`, f: t => [A * (t / q - p / (q * q)), A / q] };
      if (!z(p)) return { tex: String.raw`y_p=\frac{A}{2p}t^2-\frac{A}{p^2}t`, f: t => [A * t * t / (2 * p) - A * t / (p * p), A * t / p - A / (p * p)] };
      return { tex: String.raw`y_p=\frac A6t^3`, f: t => [A * t ** 3 / 6, A * t * t / 2] };
  }
  return null;
}
// 常数变易法：累积 Simpson 积分求 u1, u2（零初值特解）
function vop(B, n = 2000) {
  const h = T / n, i1 = new Float64Array(n + 1), i2 = new Float64Array(n + 1);
  for (let k = 0; k <= n; k++) { const t = k * h, a = B.y1(t), b = B.y2(t), W = a[0] * b[1] - a[1] * b[0], g = gFn(t); i1[k] = -b[0] * g / W; i2[k] = a[0] * g / W; }
  const ts = [], yp = [], W = [];
  let U1 = 0, U2 = 0;
  for (let k = 0; k <= n; k += 2) {
    if (k > 0) { U1 += h / 3 * (i1[k - 2] + 4 * i1[k - 1] + i1[k]); U2 += h / 3 * (i2[k - 2] + 4 * i2[k - 1] + i2[k]); }
    const t = k * h, a = B.y1(t), b = B.y2(t);
    ts.push(t); yp.push(U1 * a[0] + U2 * b[0]); W.push(a[0] * b[1] - a[1] * b[0]);
  }
  return { ts, yp, W };
}
function rkSolve(p, q, withG, n = 2000) {
  return P.rk4((t, y) => [y[1], -p * y[1] - q * y[0] + (withG ? gFn(t) : 0)], [S.y0, S.v0], 0, T, n);
}

/* ---------- 渲染 ---------- */
const num = v => (+v.toFixed(4)).toString();
function rootsTex(B) {
  if (B.kind === 'real') return String.raw`r_1=${num(B.r1)},\ r_2=${num(B.r2)}`;
  if (B.kind === 'double') return String.raw`r_1=r_2=${num(B.r1)}`;
  return String.raw`r=${num(B.a)}\pm ${num(B.b)}i`;
}
function basisTex(B) {
  if (B.kind === 'real') return String.raw`y_1=e^{${num(B.r1)}t},\ y_2=e^{${num(B.r2)}t}`;
  if (B.kind === 'double') return String.raw`y_1=e^{${num(B.r1)}t},\ y_2=t\,e^{${num(B.r1)}t}`;
  return String.raw`y_1=e^{${num(B.a)}t}\cos ${num(B.b)}t,\ y_2=e^{${num(B.a)}t}\sin ${num(B.b)}t`;
}
function renderA() {
  const { p, q } = S, B = basis(p, q), [C1, C2] = constants(B, S.y0, S.v0), R = String.raw;
  const [cl, pill] = classify(p, q, B);
  P.tex($('od-aEq'), R`\[ y''+${num(p)}y'+${num(q)}y=0,\ \ y(0)=${num(S.y0)},\ y'(0)=${num(S.v0)}\quad\Rightarrow\quad r^2+${num(p)}r+${num(q)}=0,\ \ \Delta=p^2-4q=${num(B.D)} \]`);
  P.tex($('od-aSol'), R`\[ ${rootsTex(B)};\qquad ${basisTex(B)};\qquad y=C_1y_1+C_2y_2,\ \ C_1=${num(C1)},\ C_2=${num(C2)} \]`);
  const rk = rkSolve(p, q, false), ex = rk.ts.map(t => solH(B, C1, C2, t)[0]);
  let err = 0; rk.ts.forEach((t, i) => { err = Math.max(err, Math.abs(rk.ys[i][0] - ex[i])); });
  const ym = Math.min(8, Math.max(1.2, ...ex.map(Math.abs)) * 1.1);
  const f = P.frame($('od-cvA'), { xmin: 0, xmax: T, ymin: -ym, ymax: ym, xlabel: 't' });
  if (f) {
    P.line(f, rk.ts, ex, P.C.accent, 2.2);
    for (let i = 0; i < rk.ts.length; i += 50) P.dot(f, rk.ts[i], rk.ys[i][0], P.C.warm, 2.5);
    if (B.kind === 'complex' && p > 0) { P.line(f, rk.ts, t => Math.hypot(C1, C2) * Math.exp(B.a * t), P.C['ink-3'], 1, [4, 4]); P.line(f, rk.ts, t => -Math.hypot(C1, C2) * Math.exp(B.a * t), P.C['ink-3'], 1, [4, 4]); }
  }
  const vs = rk.ts.map(t => solH(B, C1, C2, t)[1]), vm = Math.min(20, Math.max(1, ...vs.map(Math.abs)) * 1.1);
  const g = P.frame($('od-cvPh'), { xmin: -ym, xmax: ym, ymin: -vm, ymax: vm, xlabel: 'y', ylabel: "y'" });
  if (g) { P.line(g, ex, vs, P.C.accent, 1.8); P.dot(g, S.y0, S.v0, P.C.warm, 4.5); }
  const k = P.frame($('od-cvRoots'), { xmin: -6, xmax: 3, ymin: -5.5, ymax: 5.5, equal: true, xlabel: 'Re r', ylabel: 'Im r' });
  if (k) {
    P.seg2(k, 0, -5.5, 0, 5.5, P.C['ink-3'], 1);
    if (B.kind === 'complex') { P.dot(k, B.a, B.b, P.C.accent, 5); P.dot(k, B.a, -B.b, P.C.accent, 5); }
    else { P.dot(k, B.r1, 0, P.C.accent, 5); P.dot(k, B.r2, 0, B.kind === 'double' ? P.C.warm : P.C.accent, B.kind === 'double' ? 8 : 5, B.kind === 'double'); }
  }
  P.stats($('od-aStats'), [
    ['类型', `<span class="pill ${pill}">${cl}</span>`],
    ['阻尼比 ζ = p/(2√q)', q > 0 ? fx(p / (2 * Math.sqrt(q)), 4) : '—'],
    ['固有角频率 ω₀ = √q', q > 0 ? fx(Math.sqrt(q), 4) : '—'],
    ['阻尼振动频率 ω_d', B.kind === 'complex' ? fx(B.b, 4) : '—（不振荡）'],
    ['RK4（h=0.005）与解析解最大差', sci(err)],
    ['重根判据', `|Δ| &lt; 10⁻⁹ 视为重根`]
  ]);
}
const GNAME = { const: '常数 A', exp: 'A e^{βt}', sin: 'A sin ωt', poly: 'A t', cauchy: 'A/(1+t²)' };
function renderB() {
  const { p, q } = S, B = basis(p, q), R = String.raw, uc = ucParticular(p, q), vp = vop(B);
  const [Ch1, Ch2] = constants(B, S.y0, S.v0);
  const yVoP = vp.ts.map((t, i) => solH(B, Ch1, Ch2, t)[0] + vp.yp[i]);
  let yUC = null, cU = null;
  if (uc) { const [a0, b0] = uc.f(0); cU = constants(B, S.y0 - a0, S.v0 - b0); yUC = vp.ts.map(t => solH(B, cU[0], cU[1], t)[0] + uc.f(t)[0]); }
  const rk = rkSolve(p, q, true), yRK = vp.ts.map((_, i) => rk.ys[2 * i][0]);
  const gtex = { const: num(S.A), exp: `${num(S.A)}e^{${num(S.beta)}t}`, sin: `${num(S.A)}\\sin ${num(S.w)}t`, poly: `${num(S.A)}t`, cauchy: `\\dfrac{${num(S.A)}}{1+t^2}` }[S.g];
  P.tex($('od-bEq'), R`\[ y''+${num(p)}y'+${num(q)}y=g(t)=${gtex},\qquad y(0)=${num(S.y0)},\ y'(0)=${num(S.v0)} \]`);
  const W0 = vp.W[0];
  const steps = [
    R`\textbf{1. 齐次通解：}\ ${basisTex(B)}`,
    R`\textbf{2. Wronskian：}\ W(t)=y_1y_2'-y_1'y_2=${num(W0)}\,e^{${num(-p)}t}\neq0\ \ (\text{Abel 公式 } W=W(0)e^{-pt})`,
    R`\textbf{3. 常数变易：}\ y_p=u_1y_1+u_2y_2,\ \ u_1'=-\frac{y_2\,g}{W},\ \ u_2'=\frac{y_1\,g}{W}`,
    R`\textbf{4. 积分：}\ y_p(t)=\int_0^t\frac{y_1(s)y_2(t)-y_1(t)y_2(s)}{W(s)}\,g(s)\,ds\ \ (\text{零初值特解})`,
    R`\textbf{5. 通解：}\ y=C_1y_1+C_2y_2+y_p,\ \ \text{由初值得 }C_1=${num(Ch1)},\ C_2=${num(Ch2)}`
  ];
  P.tex($('od-bVoP'), steps.map(s => `\\[ ${s} \\]`).join(''));
  P.tex($('od-bUC'), uc ? R`\[ \text{待定系数法：}\ ${uc.tex} \]` : R`<span class="pill warn">待定系数法不适用</span> \(g=A/(1+t^2)\) 不是多项式、指数、正余弦及其乘积，猜不出特解的形式；常数变易法仍然适用（积分没有初等原函数时用数值积分）。`);
  const all = [...yVoP, ...yRK].map(Math.abs), ym = Math.min(20, Math.max(1, ...all) * 1.1);
  const f = P.frame($('od-cvB'), { xmin: 0, xmax: T, ymin: -ym, ymax: ym, xlabel: 't' });
  if (f) {
    P.line(f, vp.ts, vp.ts.map(gFn), P.C['ink-3'], 1, [4, 4], 4);
    P.line(f, vp.ts, vp.yp, P.C.violet, 1.4);
    if (yUC) P.line(f, vp.ts, yUC, P.C.warm, 3.5);
    P.line(f, vp.ts, yVoP, P.C.accent, 1.8);
    for (let i = 0; i < vp.ts.length; i += 25) P.dot(f, vp.ts[i], yRK[i], P.C.ink, 2.2);
  }
  let eV = 0, eU = 0, eW = 0;
  yVoP.forEach((v, i) => { eV = Math.max(eV, Math.abs(v - yRK[i]) / Math.max(1, Math.abs(yRK[i]))); if (yUC) eU = Math.max(eU, Math.abs(v - yUC[i]) / Math.max(1, Math.abs(yUC[i]))); });
  vp.W.forEach((w, i) => { eW = Math.max(eW, Math.abs(w - W0 * Math.exp(-p * vp.ts[i])) / Math.abs(W0 * Math.exp(-p * vp.ts[i]))); });
  P.stats($('od-bStats'), [
    ['常数变易 vs 待定系数', yUC ? sci(eU) : '待定系数法不适用'],
    ['常数变易 vs RK4', sci(eV)],
    ['Abel 公式相对误差', sci(eW)],
    ['共振', S.g === 'sin' && Math.abs(p) < TOL && Math.abs(q - S.w * S.w) < TOL ? '<span class="pill warn">ω² = q：振幅线性增长</span>' : S.g === 'exp' && Math.abs(S.beta ** 2 + p * S.beta + q) < TOL ? '<span class="pill warn">β 是特征根</span>' : '无']
  ]);
}
function shootF(lam, x) { // X(0)=0, X'(0)=1
  if (lam > 1e-12) { const k = Math.sqrt(lam); return [Math.sin(k * x) / k, Math.cos(k * x)]; }
  if (lam < -1e-12) { const m = Math.sqrt(-lam); return [Math.sinh(m * x) / m, Math.cosh(m * x)]; }
  return [x, 1];
}
const shootVal = lam => shootF(lam, S.Ls)[S.sbc === 'D' ? 0 : 1];
function shootRoots(n = 5) {
  const out = []; let a = -30, fa = shootVal(a);
  for (let i = 1; i <= 40000 && out.length < n; i++) { const b = -30 + i * 0.01 + 0.0037, fb = shootVal(b); if (fa * fb < 0) out.push(P.bisect(shootVal, a, b, fa)); a = b; fa = fb; }
  return out;
}
function renderC() {
  const R = String.raw, lam = S.lam, L = S.Ls;
  const form = lam > 0 ? R`X=\frac{\sin\sqrt\lambda\,x}{\sqrt\lambda}\ (\text{振荡})` : lam < 0 ? R`X=\frac{\sinh\sqrt{-\lambda}\,x}{\sqrt{-\lambda}}\ (\text{单调增长})` : R`X=x\ (\text{线性})`;
  P.tex($('od-cEq'), R`\[ X''+\lambda X=0,\ X(0)=0,\ X'(0)=1\ \Rightarrow\ ${form};\quad \text{要求 } ${S.sbc === 'D' ? 'X(L)=0' : "X'(L)=0"} \]`);
  const f = P.frame($('od-cvX'), { xmin: 0, xmax: L, ymin: -1.2, ymax: 1.2, xlabel: 'x' });
  if (f) {
    const mx = Math.max(...Array.from({ length: 201 }, (_, i) => Math.abs(shootF(lam, L * i / 200)[0]))) || 1;
    P.fnLine(f, 0, L, x => shootF(lam, x)[0] / mx, P.C.accent, 2.2);
    P.dot(f, L, shootF(lam, L)[0] / mx, Math.abs(shootVal(lam)) < 0.02 * (S.sbc === 'D' ? mx : 1) ? P.C.good : P.C.bad, 5);
    P.bcGlyphs(f, 'D', S.sbc);
  }
  const roots = shootRoots(5), lmax = Math.max(60, roots[4] * 1.08 || 60);
  const g = P.frame($('od-cvS'), { xmin: -30, xmax: lmax, ymin: -1.3, ymax: 1.3, xlabel: 'λ' });
  if (g) {
    P.fnLine(g, -30, lmax, l => { const v = shootVal(l); return Math.sign(v) * Math.min(1.25, Math.abs(v) * (S.sbc === 'D' ? 3 : 1)); }, P.C.accent, 1.8, null, 1500);
    roots.forEach(r => P.dot(g, r, 0, P.C.good, 4.5));
    P.seg2(g, lam, -1.3, lam, 1.3, P.C.warm, 1.4, [4, 3]);
  }
  const exact = n => S.sbc === 'D' ? (n * PI / L) ** 2 : ((n - 0.5) * PI / L) ** 2;
  $('od-cTable').innerHTML = '<thead><tr><th>n</th><th class="num">打靶求得 λ_n</th><th class="num">解析值</th><th class="num">相对误差</th></tr></thead><tbody>' + roots.map((r, i) => `<tr><td>${i + 1}</td><td class="num">${fx(r, 8)}</td><td class="num">${fx(exact(i + 1), 8)}</td><td class="num">${sci(Math.abs(r - exact(i + 1)) / exact(i + 1))}</td></tr>`).join('') + '</tbody>';
}
function render() { renderA(); renderB(); renderC(); }

/* ---------- 基准 ---------- */
function runBench() {
  const keep = { ...S }, rows = [];
  const homErr = (p, q) => { S.p = p; S.q = q; S.y0 = 1; S.v0 = -0.5; const B = basis(p, q), [c1, c2] = constants(B, 1, -0.5), rk = rkSolve(p, q, false); let e = 0; rk.ts.forEach((t, i) => { e = Math.max(e, Math.abs(rk.ys[i][0] - solH(B, c1, c2, t)[0])); }); return e; };
  rows.push(P.benchRow('欠阻尼 p=1, q=9：解析解 vs RK4', 'RK4，h = 0.005', homErr(1, 9), 1e-8));
  rows.push(P.benchRow('临界阻尼 p=4, q=4（重根 t e^{rt}）', 'RK4，h = 0.005', homErr(4, 4), 1e-8));
  rows.push(P.benchRow('过阻尼 p=5, q=4', 'RK4，h = 0.005', homErr(5, 4), 1e-8));
  const cmp = () => { const B = basis(S.p, S.q), uc = ucParticular(S.p, S.q), vp = vop(B), [c1, c2] = constants(B, S.y0, S.v0), [a0, b0] = uc.f(0), cu = constants(B, S.y0 - a0, S.v0 - b0); let e = 0; vp.ts.forEach((t, i) => { e = Math.max(e, Math.abs(solH(B, c1, c2, t)[0] + vp.yp[i] - solH(B, cu[0], cu[1], t)[0] - uc.f(t)[0])); }); return e; };
  Object.assign(S, { p: 0.5, q: 4, g: 'sin', A: 1, w: 1.5, y0: 1, v0: 0 });
  rows.push(P.benchRow('常数变易法 vs 待定系数法（sin 强迫）', '累积 Simpson，h = 0.005', cmp(), 1e-9));
  Object.assign(S, { p: 0, q: 4, w: 2 });
  rows.push(P.benchRow('共振 ω² = q：常数变易 vs −At cos ωt/(2ω)', '累积 Simpson', cmp(), 1e-9));
  Object.assign(S, { p: 2, q: 1, g: 'exp', beta: -1 });
  rows.push(P.benchRow('β 为重根：y_p = A t² e^{βt}/2', '常数变易法对照', cmp(), 1e-9));
  Object.assign(S, { p: 1, q: 9, g: 'cauchy', A: 1 });
  const B = basis(1, 9), vp = vop(B), rk = rkSolve(1, 9, true), [c1, c2] = constants(B, S.y0, S.v0); let e = 0;
  vp.ts.forEach((t, i) => { e = Math.max(e, Math.abs(solH(B, c1, c2, t)[0] + vp.yp[i] - rk.ys[2 * i][0])); });
  rows.push(P.benchRow('g = 1/(1+t²)（只能常数变易）vs RK4', '两种独立数值方法', e, 1e-8));
  e = 0; vp.ts.forEach((t, i) => { e = Math.max(e, Math.abs(vp.W[i] - vp.W[0] * Math.exp(-t))); });
  rows.push(P.benchRow('Abel 公式 W(t) = W(0)e^{−pt}', '由 y₁、y₂ 直接计算', e, 1e-13));
  Object.assign(S, { Ls: 1, sbc: 'D' }); let r = shootRoots(5); e = 0; r.forEach((v, i) => { e = Math.max(e, Math.abs(v - ((i + 1) * PI) ** 2) / ((i + 1) * PI) ** 2); });
  rows.push(P.benchRow('打靶法：X(L)=0 的根 = (nπ/L)²', '二分法', e, 1e-12));
  S.sbc = 'N'; r = shootRoots(5); e = 0; r.forEach((v, i) => { e = Math.max(e, Math.abs(v - ((i + 0.5) * PI) ** 2) / ((i + 0.5) * PI) ** 2); });
  rows.push(P.benchRow("打靶法：X'(L)=0 的根 = ((n−½)π/L)²", '二分法', e, 1e-12));
  Object.assign(S, keep);
  P.benchTable($('od-bench'), 'ode', rows);
}

/* ---------- 文字 ---------- */
const QUIZ = [
  { type: '特征方程', q: String.raw`\(y''+4y'+4y=0\) 的通解是？`, o: [String.raw`\(C_1e^{-2t}+C_2e^{-2t}\)`, String.raw`\((C_1+C_2t)e^{-2t}\)`, String.raw`\(C_1\cos2t+C_2\sin2t\)`, String.raw`\(C_1e^{2t}+C_2te^{2t}\)`], a: 1, e: String.raw`\(r^2+4r+4=(r+2)^2\)，重根 \(r=-2\)。\(e^{-2t}\) 与 \(e^{-2t}\) 线性相关，第二个解要乘 \(t\)。选项 A 只有一个独立解。` },
  { type: '物理解释', q: '质量—弹簧—阻尼系统 y″ + py′ + qy = 0，什么情况下回到平衡位置最快且不振荡？', o: ['p = 0', 'p² < 4q', 'p² = 4q（临界阻尼）', 'p² 远大于 4q'], a: 2, e: String.raw`临界阻尼时衰减率为 \(p/2=\sqrt q\)。过阻尼时慢根 \(r=(-p+\sqrt{p^2-4q})/2\) 接近 0，回零反而更慢。汽车减震器按接近临界阻尼设计。` },
  { type: '待定系数法', q: String.raw`\(y''+4y=\sin2t\) 的特解应设为？`, o: [String.raw`\(A\sin2t\)`, String.raw`\(A\cos2t+B\sin2t\)`, String.raw`\(t(A\cos2t+B\sin2t)\)`, String.raw`\(At^2\sin2t\)`], a: 2, e: String.raw`\(\pm2i\) 是特征根（共振），普通试探式代入后为 0，必须乘 \(t\)。结果 \(y_p=-\tfrac t4\cos2t\)，振幅线性增长。` },
  { type: '常数变易法', q: '常数变易法中，为什么可以附加条件 u₁′y₁ + u₂′y₂ = 0？', o: ['因为它总是自动成立', '两个未知函数只需一个方程确定，附加条件是自由选择，用来避免出现 u″', '因为 Wronskian 为零', '因为初始条件要求'], a: 1, e: String.raw`\(y_p=u_1y_1+u_2y_2\) 有两个未知函数，代入方程只得一个条件。附加 \(u_1'y_1+u_2'y_2=0\) 使 \(y_p'\) 中不含 \(u'\)，再代入方程得 \(u_1'y_1'+u_2'y_2'=g\)。两式联立，系数行列式就是 \(W\)。` },
  { type: '与 PDE 的联系', q: '分离变量后，热方程的时间部分和波动方程的时间部分分别是？', o: ['都是指数函数', '指数衰减；三角函数振荡', '三角函数；指数函数', '都是三角函数'], a: 1, e: String.raw`热方程 \(T'+a^2\lambda T=0\) 一阶 ⇒ \(e^{-a^2\lambda t}\)；波动方程 \(T''+a^2\lambda T=0\) 二阶、特征根 \(\pm ia\sqrt\lambda\) ⇒ 正弦余弦。` },
  { type: '边界条件', q: String.raw`\(X''+\lambda X=0\)，\(X(0)=X(L)=0\)。为什么 \(\lambda=-4\) 不是特征值？`, o: ['因为 λ 必须是整数', String.raw`通解 \(A\cosh2x+B\sinh2x\) 只有 \(A=B=0\) 满足边界条件`, '因为 −4 太小', '因为对应 ODE 无解'], a: 1, e: String.raw`ODE 本身总有解，边界条件才是筛选：\(X(0)=0\Rightarrow A=0\)，\(B\sinh2L=0\Rightarrow B=0\)。实验 C 的打靶函数在 \(\lambda\lt0\) 时没有零点。` }
];
const STEPS = [
  { t: '写特征方程', b: String.raw`<p>例题：\(y''+2y'+5y=0\)，\(y(0)=1\)，\(y'(0)=0\)。设 \(y=e^{rt}\)：</p><div class="formula-box">\[ r^2+2r+5=0,\qquad \Delta=4-20=-16\lt0 . \]</div>` },
  { t: '求根并写通解', b: String.raw`<div class="formula-box">\[ r=-1\pm2i\ \Rightarrow\ y=e^{-t}(C_1\cos2t+C_2\sin2t). \]</div><p>复根 \(\alpha\pm i\beta\) 对应 \(e^{\alpha t}\cos\beta t,\ e^{\alpha t}\sin\beta t\)（Euler 公式取实部和虚部）。</p>` },
  { t: '由初值定常数', b: String.raw`<div class="formula-box">\[ y(0)=C_1=1,\quad y'(0)=-C_1+2C_2=0\Rightarrow C_2=\tfrac12,\qquad y=e^{-t}\Big(\cos2t+\tfrac12\sin2t\Big). \]</div>` },
  { t: '物理解释', b: String.raw`<p>欠阻尼振动：阻尼比 \(\zeta=p/(2\sqrt q)=1/\sqrt5\approx0.447\)，振幅按 \(e^{-t}\) 衰减，振动角频率 \(\omega_d=2\)，比无阻尼的 \(\sqrt5\) 略小。</p>` },
  { t: '加上外力：待定系数法', b: String.raw`<p>\(y''+2y'+5y=\sin t\)。\(\pm i\) 不是特征根，设 \(y_p=C\cos t+D\sin t\)：</p><div class="formula-box">\[ 4C+2D=0,\quad -2C+4D=1\ \Rightarrow\ C=-\tfrac1{10},\ D=\tfrac15 . \]</div>` },
  { t: '同一问题用常数变易法', b: String.raw`<div class="formula-box">\[ W=y_1y_2'-y_1'y_2=2e^{-2t},\quad y_p(t)=\int_0^t\frac{e^{-(t-s)}\sin2(t-s)}{2}\sin s\,ds . \]</div><p>它与待定系数法的特解相差一个齐次解，满足初值后两者给出同一个 \(y\)。常数变易法适用于任意连续的 \(g\)。</p>` },
  { t: '与 PDE 的联系', b: String.raw`<p>分离变量得到的 \(X''+\lambda X=0\)、\(T'+a^2\lambda T=0\)、\(T''+a^2\lambda T=0\)、\(r^2R''+rR'-n^2R=0\) 都是本模块的方程。有源项时，每个模态的 \(T_n\) 满足非齐次 ODE，用常数变易法就得到 Duhamel 公式。</p>` }
];
function html() {
  const R = String.raw;
  return P.modHead(ID, '06', '常微分方程与求解工具箱', R`分离变量把偏微分方程拆成若干常微分方程。本模块整理二阶常系数方程的解法：特征方程、三种根的情形、待定系数法、常数变易法，以及边界条件如何从 ODE 的解中筛选出特征值。`) +
  P.sec(ID, 'phys', '振动系统', R`<div class="cols">
    <div class="card"><h3>质量—弹簧—阻尼</h3><p class="small muted">\(m\ddot y+c\dot y+ky=F(t)\)，除以 \(m\)：\(y''+py'+qy=g\)，\(p=c/m\)，\(q=k/m\)。</p></div>
    <div class="card"><h3>RLC 电路</h3><p class="small muted">\(L\ddot Q+R\dot Q+Q/C=E(t)\)，与机械振动是同一个方程。</p></div>
    <div class="card"><h3>PDE 的模态</h3><p class="small muted">弦的第 \(n\) 阶模态就是一个频率为 \(\omega_n\) 的振子；受迫振动时每个模态独立地响应外力。</p></div></div>`) +
  P.sec(ID, 'model', '二阶线性常系数方程', R`<div class="formula-box">\[ y''+py'+qy=g(t),\qquad \text{齐次通解 } y_h=C_1y_1+C_2y_2,\qquad \text{通解 } y=y_h+y_p . \]</div><p class="prose small muted">线性 ⇒ 叠加原理成立。\(y_1,y_2\) 线性无关当且仅当 Wronskian \(W=y_1y_2'-y_1'y_2\neq0\)。</p>`) +
  P.sec(ID, 'cond', '初值问题与边值问题', R`<div class="tbl-wrap"><table><thead><tr><th>类型</th><th>条件</th><th>解的情况</th><th>在 PDE 中的角色</th></tr></thead><tbody>
    <tr><td>初值问题</td><td>\(y(0),\ y'(0)\)</td><td>系数连续时存在唯一</td><td>时间部分 \(T(t)\)</td></tr>
    <tr><td>边值问题</td><td>\(X(0),\ X(L)\) 等</td><td>可能无解、唯一或无穷多解</td><td>空间部分：只有特征值 λ 才有非零解</td></tr></tbody></table></div>`) +
  P.sec(ID, 'idea', '求解思路', R`<div class="prose small muted"><p>齐次方程：代入 \(y=e^{rt}\) 化为代数方程。非齐次方程：右端为多项式、指数、正余弦时用待定系数法（注意共振时乘 \(t\)）；一般的 \(g\) 用常数变易法，它给出积分形式的特解，也就是 Green 函数（Duhamel 原理）。</p></div>`) +
  P.sec(ID, 'deriv', '典型例题：y″ + 2y′ + 5y = 0 及其受迫振动', '<div id="od-steps"></div>') +
  P.sec(ID, 'viz', '交互实验',
    `<h3 class="sub">A　齐次方程：三种根与阻尼</h3>` +
    P.lab('HOMOG', '\\(y\'\'+py\'+qy=0\\)',
      `<div class="ctl"><label class="ctl-label" for="od-pre">预设</label><select id="od-pre">${Object.entries(PRE).map(([k, v]) => `<option value="${k}">${v.label}</option>`).join('')}</select></div>
       ${P.rangeHTML('od-p', '阻尼 \\(p\\)', -1, 6, 0.1, 1)}
       ${P.rangeHTML('od-q', '刚度 \\(q\\)', -4, 25, 0.1, 9)}
       ${P.rangeHTML('od-y0', '\\(y(0)\\)', -2, 2, 0.1, 1)}
       ${P.rangeHTML('od-v0', '\\(y\'(0)\\)', -4, 4, 0.1, 0)}
       <div class="tasks"><b>试一试</b><ul><li>固定 \\(q=4\\)，把 \\(p\\) 从 0 增加到 6：两个复根在 \\(p=4\\) 处汇合成重根，再分裂成两个负实根。</li><li>比较 p=4 和 p=6 哪个先回到 0 附近。</li></ul></div>`,
      `<div class="problem" id="od-aEq"></div>
       <div class="formula-box" id="od-aSol"></div>
       ${P.plotBlock('od-cvA', '解 \\(y(t)\\)：实线为解析解，点为 RK4 独立数值解，虚线为包络 \\(\\pm Ce^{\\alpha t}\\)', P.lg('accent', '解析解') + P.lg('warm', 'RK4'))}
       <div class="two">${P.plotBlock('od-cvPh', '相图 \\((y, y\')\\)', '', 'short')}${P.plotBlock('od-cvRoots', '特征根在复平面上的位置', '', 'short')}</div>
       <div class="stats" id="od-aStats"></div>`) +
    `<h3 class="sub">B　非齐次方程：待定系数法与常数变易法</h3>` +
    P.lab('FORCED', '\\(y\'\'+py\'+qy=g(t)\\)，\\(p,q\\) 与初值沿用实验 A',
      `<div class="ctl"><label class="ctl-label" for="od-g">右端 \\(g(t)\\)</label><select id="od-g">${Object.entries(GNAME).map(([k, v]) => `<option value="${k}">${v}</option>`).join('')}</select></div>
       ${P.rangeHTML('od-A', '幅值 \\(A\\)', -3, 3, 0.1, 1)}
       ${P.rangeHTML('od-beta', '指数 \\(\\beta\\)', -3, 1, 0.1, -0.5)}
       ${P.rangeHTML('od-w', '频率 \\(\\omega\\)', 0, 6, 0.05, 2)}
       <div class="tasks"><b>试一试</b><ul><li>选"无阻尼 p=0, q=4"，再选 sin 且 ω=2：特解变为 \\(-\\tfrac{At}{4}\\cos2t\\)，振幅线性增长。</li><li>选临界阻尼 p=4, q=4，再选指数且 β=−2：β 是重根，特解乘 \\(t^2\\)。</li><li>选 A/(1+t²)：只能用常数变易法。</li></ul></div>`,
      `<div class="problem" id="od-bEq"></div>
       <div class="formula-box small" id="od-bVoP"></div>
       <div class="formula-box small" id="od-bUC"></div>
       ${P.plotBlock('od-cvB', '三种方法得到的同一个解', P.lg('warm', '待定系数法（粗）') + P.lg('accent', '常数变易法') + P.lg('ink', 'RK4（点）') + P.lg('violet', '零初值特解 y_p') + P.lg('ink-3', 'g(t)'))}
       <div class="stats" id="od-bStats"></div>`) +
    `<h3 class="sub">C　与 PDE 的联系：边界条件如何选出特征值</h3>` +
    P.lab('SHOOT', '打靶法：从左端出发，看能否"打中"右端条件',
      `${P.rangeHTML('od-lam', '\\(\\lambda\\)', -30, 120, 0.1, 12)}
       ${P.rangeHTML('od-Ls', '区间长度 \\(L\\)', 0.5, 2, 0.05, 1)}
       <div class="ctl"><div class="ctl-label">右端条件</div><div id="od-sbc"></div></div>
       <div class="tasks"><b>试一试</b><ul><li>把 λ 拖到约 9.87：X(L) 恰好为 0，这就是 \\(\\lambda_1=\\pi^2\\)。</li><li>λ &lt; 0 时解是 sinh，单调增长，永远打不中 0。</li></ul></div>`,
      `<div class="problem" id="od-cEq"></div>
       ${P.plotBlock('od-cvX', '当前 λ 的解 \\(X(x)\\)（右端点绿色表示满足右端条件）', '')}
       ${P.plotBlock('od-cvS', '打靶函数（右端值，经缩放）随 λ 的变化：零点就是特征值；橙色虚线为当前 λ', '', 'short')}
       <div class="tbl-wrap"><table id="od-cTable"></table></div>`) +
    `<div class="tbl-wrap"><table><thead><tr><th>PDE 中的 ODE</th><th>方程</th><th>解的形式</th><th>原因</th></tr></thead><tbody>
      <tr><td>波动方程的时间部分</td><td>\\(T''+a^2\\lambda T=0\\)</td><td>\\(\\cos a\\sqrt\\lambda t,\\ \\sin a\\sqrt\\lambda t\\)</td><td>特征根纯虚 \\(\\pm ia\\sqrt\\lambda\\)</td></tr>
      <tr><td>热方程的时间部分</td><td>\\(T'+a^2\\lambda T=0\\)</td><td>\\(e^{-a^2\\lambda t}\\)</td><td>一阶方程，单个负实根</td></tr>
      <tr><td>空间部分</td><td>\\(X''+\\lambda X=0\\)</td><td>\\(\\lambda\\gt0\\)：sin、cos；\\(\\lambda=0\\)：\\(A+Bx\\)；\\(\\lambda\\lt0\\)：sinh、cosh</td><td>根的符号由 λ 决定；边界条件筛选 λ</td></tr>
      <tr><td>矩形区域 Laplace 方程</td><td>\\(Y''-\\lambda Y=0\\)</td><td>\\(\\sinh,\\cosh\\)</td><td>另一方向的分离常数符号相反</td></tr>
      <tr><td>圆盘的径向部分</td><td>\\(r^2R''+rR'-n^2R=0\\)</td><td>\\(r^n,\\ r^{-n}\\)（\\(n=0\\)：\\(1,\\ln r\\)）</td><td>Euler 方程，设 \\(R=r^m\\) 得 \\(m^2=n^2\\)</td></tr></tbody></table></div>`) +
  P.sec(ID, 'meaning', '物理意义', R`<div class="cols">
    <div class="card"><h3>特征根 = 系统的"指纹"</h3><p class="small muted">实部决定衰减或增长率，虚部决定振荡频率。所有根实部为负 ⇔ 系统渐近稳定。</p></div>
    <div class="card"><h3>共振</h3><p class="small muted">外力频率等于固有频率时，每个周期外力都做正功，无阻尼时振幅线性增长。有阻尼时振幅有界，但在 \(\omega\approx\omega_0\) 附近出现峰值。</p></div>
    <div class="card"><h3>Green 函数</h3><p class="small muted">常数变易法的积分核 \(G(t,s)\) 是"在时刻 \(s\) 敲一下"产生的响应。总响应是各时刻冲击响应的叠加：这就是 Duhamel 原理。</p></div></div>` +
    P.pitfalls([
      ['重根：', String.raw`第二个解是 \(te^{rt}\)，不是再写一次 \(e^{rt}\)。`],
      ['复根：', String.raw`\(\alpha\pm i\beta\) 对应 \(e^{\alpha t}\cos\beta t\)、\(e^{\alpha t}\sin\beta t\)，不要漏 \(e^{\alpha t}\)。`],
      ['待定系数的共振：', '右端的指数或频率与特征根重合时，试探式乘 t（重根时乘 t²）。'],
      ['常数变易的公式：', String.raw`\(u_1'=-y_2g/W\)，\(u_2'=y_1g/W\)。公式要求方程首项系数为 1。`],
      ['初值条件的使用顺序：', '先写出完整通解 y_h + y_p，再代入初值；只对 y_h 代入初值是常见错误。'],
      ['边值问题：', '与初值问题不同，边值问题可能无非零解；特征值正是使它有非零解的参数。']
    ]) + `<div class="bench exp-only"><h3 class="sub">数值基准</h3><div id="od-bench"></div></div>`) +
  P.sec(ID, 'quiz', '检查你的理解 <span class="score" id="od-score"></span>', '<div id="od-quiz"></div>');
}

/* ---------- 初始化 ---------- */
P.register('ode', {
  init(root) {
    root.innerHTML = html();
    P.stepper($('od-steps'), STEPS);
    const f1 = v => v.toFixed(1), f2 = v => v.toFixed(2);
    const rp = P.bindRange('od-p', v => { S.p = v; render(); }, f1), rq = P.bindRange('od-q', v => { S.q = v; render(); }, f1);
    P.bindRange('od-y0', v => { S.y0 = v; render(); }, f1);
    P.bindRange('od-v0', v => { S.v0 = v; render(); }, f1);
    $('od-pre').addEventListener('change', e => { const p = PRE[e.target.value]; S.p = p.p; S.q = p.q; rp.set(p.p); rq.set(p.q); render(); });
    const rows = () => { $('od-beta-row').hidden = S.g !== 'exp'; $('od-w-row').hidden = S.g !== 'sin'; };
    $('od-g').value = S.g;
    $('od-g').addEventListener('change', e => { S.g = e.target.value; rows(); renderB(); });
    P.bindRange('od-A', v => { S.A = v; renderB(); }, f1);
    P.bindRange('od-beta', v => { S.beta = v; renderB(); }, f1);
    P.bindRange('od-w', v => { S.w = v; renderB(); }, f2);
    P.bindRange('od-lam', v => { S.lam = v; renderC(); }, v => v.toFixed(1));
    P.bindRange('od-Ls', v => { S.Ls = v; renderC(); }, f2);
    P.segment($('od-sbc'), [{ v: 'D', label: "X(L)=0" }, { v: 'N', label: "X'(L)=0" }], S.sbc, v => { S.sbc = v; renderC(); });
    rows();
    P.quiz($('od-quiz'), QUIZ, $('od-score'));
    runBench();
    $('od-bench').addEventListener('click', e => { if (e.target.closest('[data-rerun]')) runBench(); });
  },
  render
});
})();
