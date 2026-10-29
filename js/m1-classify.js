/* m1-classify.js — 模块 01：二阶线性 PDE 的分类、特征线与标准型 */
(function () {
'use strict';
const P = window.PDE, $ = P.$, fx = P.fx, sci = P.sci;
const ID = 'cl', TOL = 1e-12;
const S = { A: 4, B: 2, C: 1, sample: true };
const PRESETS = [
  { k: 'ex', label: '例题 4u_xx+4u_xy+u_yy', A: 4, B: 2, C: 1 },
  { k: 'wave', label: '波动 u_xx − u_yy', A: 1, B: 0, C: -1 },
  { k: 'lap', label: 'Laplace u_xx + u_yy', A: 1, B: 0, C: 1 },
  { k: 'hyp', label: 'u_xx + 4u_xy + 3u_yy', A: 1, B: 2, C: 3 },
  { k: 'mixed', label: 'u_xy = 0', A: 0, B: 0.5, C: 0 },
  { k: 'heat', label: '热传导型 u_xx（y 为时间）', A: 1, B: 0, C: 0 },
  { k: 'ell', label: 'u_xx + 2u_xy + 5u_yy', A: 1, B: 1, C: 5 }
];

/* ---------- 数学核心 ---------- */
const Q = (A, B, C, a, b) => A * a * a + 2 * B * a * b + C * b * b;
const Bil = (A, B, C, a1, b1, a2, b2) => A * a1 * a2 + B * (a1 * b2 + a2 * b1) + C * b1 * b2;
function classify(A, B, C) {
  const D = B * B - A * C;
  if (Math.abs(A) < TOL && Math.abs(B) < TOL && Math.abs(C) < TOL) return { D, type: 'none' };
  return { D, type: Math.abs(D) < TOL ? 'para' : D > 0 ? 'hyp' : 'ell' };
}
function normCov([a, b]) { // 让第一个非零分量为 1，得到教材式的 ξ = x + βy
  if (Math.abs(a) > TOL) return [1, b / a];
  return [0, 1];
}
// 特征坐标：返回 {xi:[α,β], eta:[α,β], kind}
function charCoords(A, B, C) {
  const { D, type } = classify(A, B, C);
  if (type === 'hyp') {
    const s = Math.sqrt(D);
    let c1, c2;
    if (Math.abs(A) > TOL) { c1 = [(-B + s) / A, 1]; c2 = [(-B - s) / A, 1]; }
    else if (Math.abs(C) > TOL) { c1 = [1, (-B + s) / C]; c2 = [1, (-B - s) / C]; }
    else { c1 = [1, 0]; c2 = [0, 1]; }
    return { type, xi: normCov(c1), eta: normCov(c2), real: 2 };
  }
  if (type === 'para') {
    const c = Math.abs(A) > TOL ? normCov([-B / A, 1]) : [1, 0];
    const eta = Math.abs(c[1]) > TOL ? [1, 0] : [0, 1]; // 任取与 ξ 无关的 η
    return { type, xi: c, eta, real: 1 };
  }
  if (type === 'ell') { const s = Math.sqrt(-D); return { type, xi: [-B / A, 1], eta: [s / A, 0], real: 0 }; }
  return { type, xi: [1, 0], eta: [0, 1], real: 0 };
}
function transformed(A, B, C, xi, eta) {
  return { A: Q(A, B, C, xi[0], xi[1]), B: Bil(A, B, C, xi[0], xi[1], eta[0], eta[1]), C: Q(A, B, C, eta[0], eta[1]), J: xi[0] * eta[1] - xi[1] * eta[0] };
}
function eig2(A, B, C) {
  const m = (A + C) / 2, r = Math.hypot((A - C) / 2, B), l1 = m + r, l2 = m - r;
  let v = Math.abs(B) > TOL ? [B, l1 - A] : (A >= C ? [1, 0] : [0, 1]);
  const n = Math.hypot(v[0], v[1]); v = [v[0] / n, v[1] / n];
  return { l1, l2, e1: v, e2: [-v[1], v[0]] };
}
// 样例解：在标准坐标里写出主部方程的解，再回到 (x,y)
function sampleU(ch, x, y) {
  const xi = ch.xi[0] * x + ch.xi[1] * y, eta = ch.eta[0] * x + ch.eta[1] * y;
  if (ch.type === 'hyp') return Math.exp(-xi * xi) + 0.6 * Math.sin(1.4 * eta);
  if (ch.type === 'para') return Math.exp(-xi * xi) * (1 + 0.25 * eta) ;
  if (ch.type === 'ell') return Math.exp(0.45 * xi) * Math.cos(0.45 * eta);
  return 0;
}
function residual(A, B, C, ch) {
  const h = 1e-3; let mr = 0, sc = 0;
  for (let i = 0; i < 9; i++) for (let j = 0; j < 9; j++) {
    const x = -2 + i * 0.5, y = -2 + j * 0.5, u = (a, b) => sampleU(ch, a, b);
    const uxx = (u(x + h, y) - 2 * u(x, y) + u(x - h, y)) / (h * h), uyy = (u(x, y + h) - 2 * u(x, y) + u(x, y - h)) / (h * h);
    const uxy = (u(x + h, y + h) - u(x + h, y - h) - u(x - h, y + h) + u(x - h, y - h)) / (4 * h * h);
    mr = Math.max(mr, Math.abs(A * uxx + 2 * B * uxy + C * uyy));
    sc = Math.max(sc, Math.abs(A * uxx) + Math.abs(2 * B * uxy) + Math.abs(C * uyy));
  }
  return mr / (sc || 1);
}

/* ---------- TeX ---------- */
const num = v => { const s = (+v.toFixed(4)).toString(); return s; };
function poly(terms) { // terms: [[coef, symbol]]
  let s = '';
  terms.forEach(([c, sym]) => {
    if (Math.abs(c) < 1e-10) return;
    const a = Math.abs(c), sign = c < 0 ? '-' : (s ? '+' : '');
    s += sign + (Math.abs(a - 1) < 1e-10 && sym ? '' : num(a)) + sym;
  });
  return s || '0';
}
const lin = ([a, b]) => poly([[a, 'x'], [b, 'y']]);
const TYPE_NAME = { hyp: '双曲型', para: '抛物型', ell: '椭圆型', none: '非二阶方程' };

/* ---------- 渲染 ---------- */
function render() {
  const { A, B, C } = S, cls = classify(A, B, C), ch = charCoords(A, B, C), T = transformed(A, B, C, ch.xi, ch.eta), E = eig2(A, B, C);
  const R = String.raw;
  P.tex($('cl-eq'), R`\[ ${poly([[A, 'u_{xx}'], [2 * B, 'u_{xy}'], [C, 'u_{yy}']])}+Du_x+Eu_y+Fu=G \qquad (A,\;B,\;C)=(${num(A)},\;${num(B)},\;${num(C)}) \]`);
  P.tex($('cl-disc'), R`\[ \Delta=B^2-AC=(${num(B)})^2-(${num(A)})(${num(C)})=${num(cls.D)}\;\Rightarrow\;\text{${TYPE_NAME[cls.type]}} \qquad M=\begin{pmatrix}${num(A)}&${num(B)}\\${num(B)}&${num(C)}\end{pmatrix},\ \lambda_{1,2}=${num(E.l1)},\ ${num(E.l2)} \]`);
  let chTex = '';
  if (cls.type === 'hyp') chTex = R`\[ A\,dy^2-2B\,dx\,dy+C\,dx^2=0\ \text{有两族实特征线：}\ ${lin(ch.xi)}=c_1,\quad ${lin(ch.eta)}=c_2 \]`;
  else if (cls.type === 'para') chTex = R`\[ \text{判别式为零，只有一族实特征线：}\ ${lin(ch.xi)}=c \]`;
  else if (cls.type === 'ell') chTex = R`\[ \text{特征方程的根是共轭复数，没有实特征线；取复特征 } \varphi=\xi+i\eta \text{ 的实部、虚部：}\ \xi=${lin(ch.xi)},\ \eta=${lin(ch.eta)} \]`;
  else chTex = R`\[ A=B=C=0\text{：没有二阶项，不属于二阶方程的分类范围} \]`;
  P.tex($('cl-char'), chTex);
  let can = '';
  const jac = R`J=\dfrac{\partial(\xi,\eta)}{\partial(x,y)}=${num(T.J)}\neq0`;
  if (cls.type === 'hyp') can = R`\[ \xi=${lin(ch.xi)},\ \eta=${lin(ch.eta)},\ ${jac}:\quad \bar A=${num(T.A)},\ \bar B=${num(T.B)},\ \bar C=${num(T.C)}\ \Rightarrow\ u_{\xi\eta}=\text{低阶项} \]`;
  else if (cls.type === 'para') can = R`\[ \xi=${lin(ch.xi)},\ \eta=${lin(ch.eta)}\ (\text{任取无关}),\ ${jac}:\quad \bar A=${num(T.A)},\ \bar B=${num(T.B)},\ \bar C=${num(T.C)}\ \Rightarrow\ u_{\eta\eta}=\text{低阶项} \]`;
  else if (cls.type === 'ell') can = R`\[ ${jac}:\quad \bar A=${num(T.A)},\ \bar B=${num(T.B)},\ \bar C=${num(T.C)}\ \Rightarrow\ u_{\xi\xi}+u_{\eta\eta}=\text{低阶项} \]`;
  P.tex($('cl-can'), can);
  const rot = cls.type === 'none' ? '' : R`\[ \text{主轴旋转 } s=${lin(E.e1)},\ t=${lin(E.e2)}:\ \ ${poly([[E.l1, 'u_{ss}'], [E.l2, 'u_{tt}']])}+\cdots \]`;
  P.tex($('cl-rot'), rot + (cls.type === 'hyp' ? '<p class="small muted">旋转后得到双曲型的第二标准型 \\(u_{ss}-u_{tt}\\)（伸缩后），但 \\(s\\)、\\(t\\) 不是特征坐标：混合导数被消去了，\\(u_{ss}\\)、\\(u_{tt}\\) 却仍在。再令 \\(\\xi=s+t,\\ \\eta=s-t\\) 才回到特征坐标下的 \\(u_{\\xi\\eta}\\)。</p>' : cls.type === 'para' ? '<p class="small muted">旋转后有一个特征值为 0：零特征值对应的主轴方向恰好垂直于特征线。抛物型只有这一种化简。</p>' : cls.type === 'ell' ? '<p class="small muted">两个特征值同号，伸缩后得到 \\(u_{\\xi\\xi}+u_{\\eta\\eta}\\)。旋转只是化简的一种方式，与"实特征线"无关。</p>' : ''));
  const res = cls.type === 'none' ? NaN : residual(A, B, C, ch);
  P.stats($('cl-stats'), [
    ['Δ = B² − AC', fx(cls.D, 4)],
    ['类型', `<span class="pill ${cls.type === 'hyp' ? 'warn' : cls.type === 'para' ? 'info' : cls.type === 'ell' ? 'ok' : 'bad'}">${TYPE_NAME[cls.type]}</span>`],
    ['实特征线族数', cls.type === 'none' ? '—' : ch.real],
    ['特征线斜率 dy/dx', cls.type === 'hyp' ? slope(ch.xi) + ', ' + slope(ch.eta) : cls.type === 'para' ? slope(ch.xi) : '无（复）'],
    ['变换后 (Ā, B̄, C̄)', cls.type === 'none' ? '—' : `${fx(T.A, 3)}, ${fx(T.B, 3)}, ${fx(T.C, 3)}`],
    ['样例解主部残差（相对）', isFinite(res) ? sci(res) : '—']
  ]);
  drawXY(cls, ch); drawSymbol(cls, ch); drawXiEta(cls, ch);
}
function slope([a, b]) { return Math.abs(b) < TOL ? '∞（竖直）' : fx(-a / b, 3); }
function drawXY(cls, ch) {
  const f = P.frame($('cl-cvXY'), { xmin: -3, xmax: 3, ymin: -3, ymax: 3, equal: true, xlabel: 'x', ylabel: 'y' });
  if (!f) return;
  if (S.sample && cls.type !== 'none') {
    const n = 140, d = new Float64Array(n * n); let lo = Infinity, hi = -Infinity;
    for (let j = 0; j < n; j++) for (let i = 0; i < n; i++) { const v = sampleU(ch, -3 + 6 * (i + 0.5) / n, 3 - 6 * (j + 0.5) / n); d[j * n + i] = v; lo = Math.min(lo, v); hi = Math.max(hi, v); }
    f.ctx.globalAlpha = 0.55; P.heatmap(f, d, n, n, lo, hi); f.ctx.globalAlpha = 1;
  }
  const fam = (cov, color, dash) => {
    for (let c = -6; c <= 6; c += 0.75) {
      const [a, b] = cov;
      if (Math.abs(b) > TOL) P.seg2(f, -3, (c + 3 * a) / b, 3, (c - 3 * a) / b, color, 1.6, dash);
      else P.seg2(f, c / a, -3, c / a, 3, color, 1.6, dash);
    }
  };
  if (cls.type === 'hyp') { fam(ch.xi, P.C.accent); fam(ch.eta, P.C.warm); }
  else if (cls.type === 'para') { fam(ch.xi, P.C.accent); fam(ch.eta, P.C['ink-3'], [4, 4]); }
  else if (cls.type === 'ell') { fam(ch.xi, P.C['ink-3'], [4, 4]); fam(ch.eta, P.C['ink-3'], [4, 4]); }
  $('cl-xyLegend').innerHTML = cls.type === 'hyp' ? P.lg('accent', 'ξ = 常数（特征线）') + P.lg('warm', 'η = 常数（特征线）')
    : cls.type === 'para' ? P.lg('accent', 'ξ = 常数（特征线）') + P.lg('ink-3', 'η = 常数（非特征线）')
    : cls.type === 'ell' ? P.lg('ink-3', 'ξ、η 坐标线（不是特征线）') : '';
}
function drawSymbol(cls, ch) {
  const f = P.frame($('cl-cvQ'), { xmin: -2.2, xmax: 2.2, ymin: -2.2, ymax: 2.2, equal: true, xlabel: 'ξ₁', ylabel: 'ξ₂' });
  if (!f) return;
  const { A, B, C } = S, n = 720, pos = [[], []], neg = [[], []];
  for (let i = 0; i <= n; i++) {
    const th = 2 * Math.PI * i / n, c = Math.cos(th), s = Math.sin(th), q = Q(A, B, C, c, s);
    const r = Math.abs(q) > 1e-9 ? 1 / Math.sqrt(Math.abs(q)) : Infinity, ok = r < 2.6;
    pos[0].push(c * r); pos[1].push(q > 0 && ok ? s * r : NaN);
    neg[0].push(c * r); neg[1].push(q < 0 && ok ? s * r : NaN);
  }
  P.line(f, pos[0], pos[1], P.C.accent, 2);
  P.line(f, neg[0], neg[1], P.C.warm, 2);
  if (cls.type === 'hyp' || cls.type === 'para') [ch.xi, ch.eta].slice(0, ch.real).forEach(([a, b]) => { const m = Math.hypot(a, b); P.seg2(f, -2.2 * a / m, -2.2 * b / m, 2.2 * a / m, 2.2 * b / m, P.C.violet, 1.4, [5, 4]); });
}
function drawXiEta(cls, ch) {
  const f = P.frame($('cl-cvXE'), { xmin: -3, xmax: 3, ymin: -3, ymax: 3, equal: true, xlabel: 'ξ', ylabel: 'η' });
  if (!f || cls.type === 'none') return;
  const n = 140, d = new Float64Array(n * n); let lo = Infinity, hi = -Infinity;
  const canon = { type: ch.type, xi: [1, 0], eta: [0, 1] };
  for (let j = 0; j < n; j++) for (let i = 0; i < n; i++) { const v = sampleU(canon, -3 + 6 * (i + 0.5) / n, 3 - 6 * (j + 0.5) / n); d[j * n + i] = v; lo = Math.min(lo, v); hi = Math.max(hi, v); }
  f.ctx.globalAlpha = 0.55; P.heatmap(f, d, n, n, lo, hi); f.ctx.globalAlpha = 1;
  for (let c = -3; c <= 3; c += 0.75) {
    P.seg2(f, c, -3, c, 3, cls.type === 'ell' ? P.C['ink-3'] : P.C.accent, 1.4, cls.type === 'ell' ? [4, 4] : null);
    P.seg2(f, -3, c, 3, c, cls.type === 'hyp' ? P.C.warm : P.C['ink-3'], 1.4, cls.type === 'hyp' ? null : [4, 4]);
  }
  const msg = { hyp: '样例解 u = F(ξ) + G(η)：沿特征线分别为常数', para: '样例解 u = F(ξ) + ηG(ξ)：u_ηη = 0', ell: '样例解 u = e^{0.45ξ}cos(0.45η)：u_ξξ + u_ηη = 0' }[cls.type];
  P.textAt(f, msg, -2.9, -2.9, P.C.ink, 'left', 'bottom');
}

/* ---------- 基准 ---------- */
function runBench() {
  const rows = [];
  let ch = charCoords(4, 2, 1), T = transformed(4, 2, 1, ch.xi, ch.eta);
  rows.push(P.benchRow('例题 4u_xx+4u_xy+u_yy：Δ = 0', '精确算术', Math.abs(classify(4, 2, 1).D), 0));
  rows.push(P.benchRow('例题：ξ = x − 2y, η = x 后 Ā = B̄ = 0, C̄ = 4', '代入变换公式', Math.abs(T.A) + Math.abs(T.B) + Math.abs(T.C - 4) + Math.abs(ch.xi[1] + 2), 1e-14));
  ch = charCoords(1, 0, -1);
  rows.push(P.benchRow('u_xx − u_yy：特征线斜率 ±1', '解析 dy/dx = ±1', Math.abs(Math.abs(ch.xi[0] / ch.xi[1]) - 1) + Math.abs(Math.abs(ch.eta[0] / ch.eta[1]) - 1), 1e-14));
  ch = charCoords(1, 2, 3); T = transformed(1, 2, 3, ch.xi, ch.eta);
  rows.push(P.benchRow('u_xx+4u_xy+3u_yy：特征坐标下 Ā = C̄ = 0', '代入变换公式', Math.abs(T.A) + Math.abs(T.C), 1e-13));
  ch = charCoords(1, 1, 5); T = transformed(1, 1, 5, ch.xi, ch.eta);
  rows.push(P.benchRow('u_xx+2u_xy+5u_yy：Ā = C̄, B̄ = 0', '代入变换公式', Math.abs(T.A - T.C) + Math.abs(T.B), 1e-13));
  ch = charCoords(0, 0.5, 0);
  rows.push(P.benchRow('u_xy = 0：特征线为 x = c, y = c', '方向向量 (1,0),(0,1)', Math.abs(ch.xi[0] - 1) + Math.abs(ch.xi[1]) + Math.abs(ch.eta[0]) + Math.abs(ch.eta[1] - 1), 0));
  let inv = 0; const rnd = (s => () => (s = (s * 16807) % 2147483647) / 2147483647 * 6 - 3)(7);
  for (let k = 0; k < 20; k++) {
    const A = rnd(), B = rnd(), C = rnd(), a1 = rnd(), b1 = rnd(), a2 = rnd(), b2 = rnd();
    const t = transformed(A, B, C, [a1, b1], [a2, b2]), D1 = t.B * t.B - t.A * t.C, D0 = (B * B - A * C) * t.J * t.J;
    inv = Math.max(inv, Math.abs(D1 - D0) / (Math.abs(D0) + 1));
  }
  rows.push(P.benchRow('判别式变换律 Δ̄ = J²Δ（20 组随机系数）', '独立计算两边', inv, 1e-12));
  let rr = 0;[[1, 2, 3], [4, 2, 1], [1, 1, 5], [1, 0, -1]].forEach(([A, B, C]) => { rr = Math.max(rr, residual(A, B, C, charCoords(A, B, C))); });
  rows.push(P.benchRow('样例解满足主部方程（4 个方程）', '中心差分 h=10⁻³ 的相对残差', rr, 1e-5));
  P.benchTable($('cl-bench'), 'classify', rows);
}

/* ---------- 文字内容 ---------- */
const QUIZ = [
  { type: '判断类型', q: String.raw`\(u_{xx}+4u_{xy}+3u_{yy}=0\) 属于哪一类？`, o: ['椭圆型', '抛物型', '双曲型', '无法判断'], a: 2, e: String.raw`\(2B=4\Rightarrow B=2\)，\(\Delta=4-3=1\gt0\)，双曲型。特征线 \(y-3x=c_1\)、\(y-x=c_2\)。` },
  { type: '易错：混合项', q: String.raw`方程 \(4u_{xx}+4u_{xy}+u_{yy}=0\) 中，判别式里的 \(B\) 等于？`, o: ['4', '2', '1', '8'], a: 1, e: String.raw`标准形式是 \(Au_{xx}+2Bu_{xy}+Cu_{yy}\)，混合项系数是 \(2B=4\)，所以 \(B=2\)，\(\Delta=4-4=0\)。若误用 \(B=4\) 会得到 \(\Delta=12\)，错判为双曲型。` },
  { type: '判断特征线', q: String.raw`例题 \(4u_{xx}+4u_{xy}+u_{yy}=0\) 的特征线是？`, o: ['x + 2y = c', '2x − y = c', 'x − 2y = c', 'y − 2x = c'], a: 2, e: String.raw`\(4dy^2-4dx\,dy+dx^2=(2dy-dx)^2=0\Rightarrow dy/dx=1/2\Rightarrow x-2y=c\)。` },
  { type: '变系数', q: String.raw`Tricomi 方程 \(yu_{xx}+u_{yy}=0\) 在 \(y\lt0\) 的区域是？`, o: ['椭圆型', '双曲型', '抛物型', '处处同一类型'], a: 1, e: String.raw`\(A=y,\ B=0,\ C=1\)，\(\Delta=-y\)。\(y\lt0\) 时 \(\Delta\gt0\) 为双曲型，\(y\gt0\) 为椭圆型，\(y=0\) 上退化为抛物型。类型可以逐点变化。` },
  { type: '识别推导错误', q: String.raw`某同学把一个椭圆型方程经过主轴旋转和伸缩后写成 \(u_{ss}-u_{tt}=0\)。说明什么？`, o: ['正确，椭圆型也可以这样化简', '推导有误：可逆变换不改变 Δ 的符号', '旋转角取错了，应换一个角度', '说明方程其实是抛物型'], a: 1, e: String.raw`\(\bar\Delta=J^2\Delta\)，可逆变换（\(J\neq0\)）不改变判别式的符号；矩阵 \(M\) 的特征值符号也不变（惯性定理）。椭圆型只能化为 \(u_{ss}+u_{tt}\)。` },
  { type: '物理解释', q: String.raw`热传导方程 \(u_t=a^2u_{xx}\) 的特征线是？`, o: ['x ± at = c', 't = c', 'x = c', '没有实特征线'], a: 1, e: String.raw`把 \(t\) 当作 \(y\)：\(A=a^2,\ B=C=0\)，\(\Delta=0\)，特征方程 \(a^2dt^2=0\Rightarrow t=c\)。这条"特征线"是水平线，说明扰动瞬间影响整根杆（无限传播速度）。` },
  { type: '坐标变换', q: '对抛物型方程，取特征坐标 ξ 后，η 应该怎样选？', o: ['必须取另一族特征线', '任取一个与 ξ 函数无关的光滑函数', '必须取 η = ξ²', '必须取主轴方向'], a: 1, e: String.raw`抛物型只有一族实特征线。\(\xi\) 取特征坐标保证 \(\bar A=0\)，进而 \(\bar B=0\)；\(\eta\) 只需使雅可比行列式非零。例题中取 \(\eta=x\) 或 \(\eta=y\) 都可以。` }
];
const STEPS = [
  { t: '读出系数（注意混合项）', b: String.raw`<p>对照 \(Au_{xx}+2Bu_{xy}+Cu_{yy}\)：</p><div class="formula-box">\[ 4u_{xx}+4u_{xy}+u_{yy}=0\;\Rightarrow\;A=4,\quad 2B=4\Rightarrow B=2,\quad C=1. \]</div>` },
  { t: '计算判别式，判断类型', b: String.raw`<div class="formula-box">\[ \Delta=B^2-AC=4-4=0\;\Rightarrow\;\text{抛物型（在整个平面上）}. \]</div>` },
  { t: '写出特征方程并求特征线', b: String.raw`<div class="formula-box">\[ A\,dy^2-2B\,dx\,dy+C\,dx^2=4dy^2-4dx\,dy+dx^2=(2dy-dx)^2=0\;\Rightarrow\;\frac{dy}{dx}=\frac12 . \]</div><p>积分得一族特征线 \(x-2y=c\)。</p>` },
  { t: '选取新坐标', b: String.raw`<p>取 \(\xi=x-2y\)（特征坐标），\(\eta=x\)（任选，与 \(\xi\) 无关）。雅可比 \(\partial(\xi,\eta)/\partial(x,y)=\begin{vmatrix}1&-2\\1&0\end{vmatrix}=2\neq0\)，变换可逆。</p>` },
  { t: '链式法则换元', b: String.raw`<div class="formula-box">\[ u_x=u_\xi+u_\eta,\quad u_y=-2u_\xi;\qquad u_{xx}=u_{\xi\xi}+2u_{\xi\eta}+u_{\eta\eta},\quad u_{xy}=-2u_{\xi\xi}-2u_{\xi\eta},\quad u_{yy}=4u_{\xi\xi}. \]</div>` },
  { t: '代入并化简', b: String.raw`<div class="formula-box">\[ 4(u_{\xi\xi}+2u_{\xi\eta}+u_{\eta\eta})+4(-2u_{\xi\xi}-2u_{\xi\eta})+4u_{\xi\xi}=4u_{\eta\eta}=0 . \]</div><p>\(u_{\xi\xi}\) 和 \(u_{\xi\eta}\) 的系数都消失了：这不是巧合，因为 \(\xi\) 满足 \(Q(\xi_x,\xi_y)=0\)。</p>` },
  { t: '求通解并检验', b: String.raw`<div class="formula-box">\[ u_{\eta\eta}=0\;\Rightarrow\;u=f(\xi)+\eta\,g(\xi)=f(x-2y)+x\,g(x-2y). \]</div><p>可直接代回原方程验证。下方实验中的样例解就是这种形式，页面实时计算它的有限差分残差。</p>` },
  { t: '物理意义', b: String.raw`<p>抛物型方程只有一族特征线，与热传导方程同类：沿特征方向信息"瞬时"传遍，没有有限的传播锥。</p>` }
];
function html() {
  const R = String.raw;
  return P.modHead(ID, '01', '二阶偏微分方程分类、特征线与标准型', R`对一般的二阶线性方程 \(Au_{xx}+2Bu_{xy}+Cu_{yy}\) \(+\,Du_x+Eu_y+Fu=G\)，二阶主部决定方程的类型、特征线和可以提的定解条件。本模块用判别式分类，求特征线，并用坐标变换化为标准型。`) +
  P.sec(ID, 'phys', '三类方程从哪里来', R`<div class="cols">
    <div class="card"><h3>双曲型：波动</h3><p class="small muted">弦振动 \(u_{tt}=a^2u_{xx}\)。扰动以有限速度沿两族特征线 \(x\pm at=c\) 传播。</p></div>
    <div class="card"><h3>抛物型：扩散</h3><p class="small muted">热传导 \(u_t=a^2u_{xx}\)。只有一族特征线 \(t=c\)，扰动瞬间影响全域，随时间平滑化。</p></div>
    <div class="card"><h3>椭圆型：平衡</h3><p class="small muted">稳态温度、静电势 \(u_{xx}+u_{yy}=0\)。没有实特征线，没有"时间方向"，边界上每一点都影响区域内每一点。</p></div></div>`) +
  P.sec(ID, 'model', '一般形式与主部', R`<div class="formula-box">\[ Au_{xx}+2Bu_{xy}+Cu_{yy}+Du_x+Eu_y+Fu=G,\qquad Q(\xi_1,\xi_2)=A\xi_1^2+2B\xi_1\xi_2+C\xi_2^2=\begin{pmatrix}\xi_1&\xi_2\end{pmatrix}\begin{pmatrix}A&B\\B&C\end{pmatrix}\begin{pmatrix}\xi_1\\\xi_2\end{pmatrix} \]</div>
    <div class="prose small muted"><p>只有二阶项（主部）决定类型；\(D,E,F,G\) 是低阶项，不影响分类。系数可以依赖 \(x,y\)，此时类型逐点判定。混合项写作 \(2B\) 是为了让系数矩阵对称，这是最常见的出错点。</p>
    <p>在可逆变换 \(\xi=\xi(x,y),\ \eta=\eta(x,y)\) 下，主部系数变为 \(\bar A=Q(\xi_x,\xi_y)\)，\(\bar B=(\xi_x,\xi_y)M(\eta_x,\eta_y)^T\)，\(\bar C=Q(\eta_x,\eta_y)\)，并且 \(\bar B^2-\bar A\bar C=J^2(B^2-AC)\)。所以判别式的符号（即类型）与坐标无关。</p></div>`) +
  P.sec(ID, 'cond', '类型决定什么样的定解条件合适', R`<div class="tbl-wrap"><table><thead><tr><th>类型</th><th>判别式</th><th>标准型</th><th>典型的适定问题</th></tr></thead><tbody>
    <tr><td>双曲型</td><td>\(\Delta\gt0\)</td><td>\(u_{\xi\eta}=\cdots\) 或 \(u_{\xi\xi}-u_{\eta\eta}=\cdots\)</td><td>初值问题（给 \(u\) 与 \(u_t\)）、初边值问题</td></tr>
    <tr><td>抛物型</td><td>\(\Delta=0\)</td><td>\(u_{\eta\eta}=\cdots\)</td><td>初边值问题（只给 \(u(x,0)\)）</td></tr>
    <tr><td>椭圆型</td><td>\(\Delta\lt0\)</td><td>\(u_{\xi\xi}+u_{\eta\eta}=\cdots\)</td><td>边值问题（Dirichlet、Neumann、Robin）</td></tr></tbody></table></div>
    <p class="prose small muted">对椭圆型方程给初值（Cauchy 数据）一般是不适定的（Hadamard 例子：解不连续依赖于数据）。这是分类在应用中的直接后果。</p>`) +
  P.sec(ID, 'idea', '求解思路', R`<div class="prose small muted"><p>目标是选新坐标让主部最简单。若 \(\xi(x,y)\) 满足 \(Q(\xi_x,\xi_y)=0\)，则 \(\bar A=0\)。这个一阶偏微分方程等价于沿曲线 \(\xi=c\) 有 \(A\,dy^2-2B\,dx\,dy+C\,dx^2=0\)，即特征方程，其积分曲线就是特征线。实特征线的条数由 \(\Delta\) 的符号决定：两族、一族或没有。</p></div>`) +
  P.sec(ID, 'deriv', String.raw`典型例题：\(4u_{xx}+4u_{xy}+u_{yy}=0\)`, '<div id="cl-steps"></div>') +
  P.sec(ID, 'viz', '交互实验：分类、特征线与标准化',
    P.lab('CLASSIFY', '调节主部系数，观察类型、特征线与变换后的系数',
      `<div class="ctl"><label class="ctl-label" for="cl-preset">预设方程</label><select id="cl-preset">${PRESETS.map(p => `<option value="${p.k}">${p.label}</option>`).join('')}</select></div>
       ${P.rangeHTML('cl-A', '\\(A\\)（\\(u_{xx}\\) 系数）', -5, 5, 0.5, 4)}
       ${P.rangeHTML('cl-B', '\\(B\\)（\\(u_{xy}\\) 系数为 \\(2B\\)）', -5, 5, 0.5, 2)}
       ${P.rangeHTML('cl-C', '\\(C\\)（\\(u_{yy}\\) 系数）', -5, 5, 0.5, 1)}
       <label class="small"><input type="checkbox" id="cl-sample" checked> 显示样例解的热力图</label>
       <div class="tasks"><b>试一试</b><ul><li>从例题出发把 \\(C\\) 从 1 调到 0.5 和 1.5：一族特征线分裂成两族，或消失。</li><li>选 \\(u_{xy}=0\\)：\\(A=C=0\\) 时特征线是坐标轴方向，不能用 \\(dy/dx\\) 的公式直接求。</li><li>对照"主轴旋转"一行：双曲型旋转后仍是两个二阶项，特征坐标才得到 \\(u_{\\xi\\eta}\\)。</li></ul></div>`,
      `<div class="problem" id="cl-eq"></div>
       <div class="formula-box" id="cl-disc"></div>
       <div class="stats" id="cl-stats"></div>
       <div class="formula-box" id="cl-char"></div>
       <div class="formula-box" id="cl-can"></div>
       <div class="formula-box small" id="cl-rot"></div>
       <div class="two">
         <div class="plot-block"><div class="plot-cap"><span>\\((x,y)\\) 平面：特征线与样例解</span><span class="legend" id="cl-xyLegend"></span></div><canvas id="cl-cvXY" class="square"></canvas></div>
         <div class="plot-block"><div class="plot-cap"><span>\\((\\xi,\\eta)\\) 平面：变换后的坐标线</span></div><canvas id="cl-cvXE" class="square"></canvas></div>
       </div>
       <div class="plot-block"><div class="plot-cap"><span>主部二次型 \\(Q(\\xi_1,\\xi_2)=\\pm1\\) 的水平集</span><span class="legend">${P.lg('accent', 'Q = +1')}${P.lg('warm', 'Q = −1')}${P.lg('violet', 'Q = 0：特征线的法向')}</span></div><canvas id="cl-cvQ" class="square" style="max-width:320px"></canvas>
       <p class="small muted">椭圆型：Q=1 是椭圆；双曲型：Q=±1 是双曲线，渐近方向 Q=0 给出特征线的法向量 \\((\\xi_x,\\xi_y)\\)；抛物型：Q=1（或 −1）是一对平行直线。</p></div>`)) +
  P.sec(ID, 'meaning', '特征线的物理意义', R`<div class="prose small muted"><p>特征线是"信息传播的路径"。双曲型方程中，初值数据的间断沿特征线传播而不被抹平；抛物型方程中，唯一的一族特征线是 \(t=\text{常数}\)，所有点同时受影响；椭圆型方程没有实特征线，解在区域内部总是光滑的（系数光滑时甚至解析）。</p></div>` +
    P.pitfalls([
      ['混合偏导项：', String.raw`\(u_{xy}\) 的系数是 \(2B\)，判别式用 \(B\)。例题 \(4u_{xy}\) 对应 \(B=2\)。`],
      ['特征方程的符号：', String.raw`是 \(A\,dy^2-2B\,dx\,dy+C\,dx^2=0\)，与二次型 \(Q(\xi_x,\xi_y)=0\) 等价（\(\xi_x/\xi_y=-dy/dx\)）。`],
      ['A = 0 的情形：', String.raw`不能用 \(dy/dx=(B\pm\sqrt\Delta)/A\)；应改解 \(dx/dy\)，或直接求 \(Q(\alpha,\beta)=0\) 的方向。`],
      ['特征坐标 ≠ 主轴旋转：', '旋转只能消去混合项。双曲型的特征坐标一般不正交；椭圆型根本没有实特征线，标准化用的是复特征的实部与虚部。'],
      ['变系数方程：', '类型逐点判定，同一方程可在不同区域属于不同类型（Tricomi 方程）。']
    ]) + `<div class="bench exp-only"><h3 class="sub">数值基准</h3><div id="cl-bench"></div></div>`) +
  P.sec(ID, 'quiz', '检查你的理解 <span class="score" id="cl-score"></span>', '<div id="cl-quiz"></div>');
}

let rA, rB, rC;
P.register('classify', {
  init(root) {
    root.innerHTML = html();
    P.stepper($('cl-steps'), STEPS);
    const upd = () => render();
    rA = P.bindRange('cl-A', v => { S.A = v; upd(); }, v => v.toFixed(1));
    rB = P.bindRange('cl-B', v => { S.B = v; upd(); }, v => v.toFixed(1) + '（2B=' + (2 * v).toFixed(0) + '）');
    rC = P.bindRange('cl-C', v => { S.C = v; upd(); }, v => v.toFixed(1));
    $('cl-preset').addEventListener('change', e => { const p = PRESETS.find(q => q.k === e.target.value); S.A = p.A; S.B = p.B; S.C = p.C; rA.set(p.A); rB.set(p.B); rC.set(p.C); upd(); });
    $('cl-sample').addEventListener('change', e => { S.sample = e.target.checked; upd(); });
    P.quiz($('cl-quiz'), QUIZ, $('cl-score'));
    runBench();
    $('cl-bench').addEventListener('click', e => { if (e.target.closest('[data-rerun]')) runBench(); });
  },
  render
});
})();
