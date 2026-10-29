/* app.js — 导航、学习模式、按需初始化、主题与尺寸变化时重绘；单个模块出错不影响其他模块 */
(function () {
'use strict';
const P = window.PDE;
const ORDER = [
  ['classify', '01', '方程分类与特征线', '分类'],
  ['wave', '02', '波动方程与达朗贝尔公式', '波动'],
  ['heat', '03', '热传导方程与扩散', '热传导'],
  ['separation', '04', '分离变量法与特征值问题', '分离变量'],
  ['laplace', '05', 'Laplace 方程与圆盘问题', 'Laplace'],
  ['ode', '06', '常微分方程工具箱', 'ODE']
];
const HINT = { class: '课堂模式：推导逐步展开，先自己想再看下一步；数值基准与进阶指标已收起。右上角可切换到探索模式。', explore: '探索模式：推导全部展开，显示全部参数、数值指标与基准测试。' };
const inited = {};
let current = null;

const rail = P.$('rail');
rail.innerHTML = ORDER.map(([id, n, t, sh]) => `<button type="button" data-mod="${id}" title="${t}"><span class="idx">${n}</span><span class="t">${sh}</span><span class="ok-dot" id="st-${id}"></span></button>`).join('');

function section(id) { return document.querySelector(`main .module[data-mod="${id}"]`); }
function init(id) {
  const def = P.modules[id], root = section(id);
  if (!def) { root.innerHTML = `<div class="err-card">模块脚本未能加载：${id}</div>`; return; }
  try {
    def.init(root);
    inited[id] = true;
    P.typeset(root);
    syncMode();
    const b = P.bench[id];
    const st = P.$('st-' + id), np = b ? b.filter(r => r.pass).length : 0, tb = st.parentNode;
    st.classList.toggle('on', !!b && np === b.length);
    if (b) tb.title = `${ORDER.find(o => o[0] === id)[2]} · 数值基准 ${np}/${b.length} 通过`;
    root.querySelectorAll('input[type="range"]').forEach(P.fillRange);
  } catch (e) {
    console.error(e);
    root.insertAdjacentHTML('afterbegin', `<div class="err-card">该模块初始化出错：${String(e && e.message || e)}。其他模块不受影响。</div>`);
  }
}
function render() {
  if (!current || !inited[current]) return;
  P.readColors();
  try { P.modules[current].render(); } catch (e) { console.error(e); }
}
function show(id) {
  if (!P.modules[id] && !ORDER.some(o => o[0] === id)) id = 'separation';
  current = id;
  document.querySelectorAll('main .module').forEach(s => { s.hidden = s.dataset.mod !== id; });
  rail.querySelectorAll('button').forEach(b => b.setAttribute('aria-current', b.dataset.mod === id ? 'page' : 'false'));
  try { history.replaceState(null, '', '#' + id); } catch (e) {}
  P.readColors();
  if (!inited[id]) init(id);
  render(); // 同步绘制：后台标签页或缩略图截取时 requestAnimationFrame 可能不触发
  spy();
}
// 阶段条：标出当前阅读到的环节
function spy() {
  const root = section(current); if (!root) return;
  const secs = root.querySelectorAll('.sec'), line = (parseFloat(getComputedStyle(document.documentElement).scrollPaddingTop) || 120) + 8;
  let k = 0;
  secs.forEach((s, i) => { if (s.getBoundingClientRect().top <= line) k = i; });
  if (innerHeight + scrollY >= document.documentElement.scrollHeight - 4) k = secs.length - 1; // 已到页底：最后一节可能够不到判定线
  const ch = root.querySelector('.chain'); if (ch) ch.style.setProperty('--prog', ((k + 1) / secs.length).toFixed(3));
  root.querySelectorAll('.chain button').forEach((b, i) => {
    if (i === k) b.setAttribute('aria-current', 'step'); else b.removeAttribute('aria-current');
    b.classList.toggle('past', i < k);
    if (i === k && b.parentNode.scrollWidth > b.parentNode.clientWidth) b.parentNode.scrollLeft = b.offsetLeft - b.parentNode.clientWidth / 2 + b.offsetWidth / 2;
  });
}
let sq = 0;
window.addEventListener('scroll', () => { if (!sq) sq = requestAnimationFrame(() => { sq = 0; spy(); }); }, { passive: true });
document.addEventListener('input', e => { if (e.target.type === 'range') P.fillRange(e.target); });

// 配色：跟随系统 → 浅色 → 深色
const THEMES = [['', '跟随系统'], ['light', '浅色'], ['dark', '深色']];
function setTheme(t) {
  if (t) document.documentElement.dataset.theme = t; else delete document.documentElement.dataset.theme;
  P.$('theme-btn').title = '配色：' + THEMES.find(x => x[0] === t)[1] + '（点击切换）';
}
P.$('theme-btn').addEventListener('click', () => {
  const cur = document.documentElement.dataset.theme || '', i = THEMES.findIndex(x => x[0] === cur);
  const nx = THEMES[(i + 1) % THEMES.length][0];
  setTheme(nx);
  try { localStorage.setItem('pde-lab-theme', nx); } catch (e) {}
});
try { const t = localStorage.getItem('pde-lab-theme'); setTheme(t === 'light' || t === 'dark' ? t : ''); } catch (e) { setTheme(''); }

// 公式库加载失败时给出提示（页面其余功能照常）
window.__mjFail = () => { if (!(window.MathJax && window.MathJax.typesetPromise)) P.$('mj-warn').hidden = false; };
if (window.__mjFailed) window.__mjFail();
function syncMode() {
  document.documentElement.dataset.mode = P.mode;
  document.querySelectorAll('.mode-seg button').forEach(b => b.setAttribute('aria-pressed', String(b.dataset.mode === P.mode)));
  document.querySelectorAll('.mode-hint').forEach(e => { e.textContent = HINT[P.mode]; });
  document.querySelectorAll('.stepper').forEach(s => s._apply && s._apply());
}
function setMode(m) {
  P.mode = m;
  try { localStorage.setItem('pde-lab-mode', m); } catch (e) {}
  syncMode();
  render();
}

rail.addEventListener('click', e => { const b = e.target.closest('button[data-mod]'); if (b) { show(b.dataset.mod); window.scrollTo(0, 0); } });
document.addEventListener('click', e => {
  const m = e.target.closest('[data-mode]'); if (m && m.closest('.mode-seg')) { setMode(m.dataset.mode); return; }
  if (!e.target.closest('main')) return;
  const g = e.target.closest('[data-go]');
  if (g) { const t = P.$(g.dataset.go); if (t) t.scrollIntoView({ behavior: matchMedia('(prefers-reduced-motion: reduce)').matches ? 'auto' : 'smooth' }); }
  const j = e.target.closest('[data-jump]'); if (j) { show(j.dataset.jump); window.scrollTo(0, 0); }
});

let rq = 0;
const schedule = () => { cancelAnimationFrame(rq); rq = requestAnimationFrame(render); };
new ResizeObserver(schedule).observe(P.$('main'));
matchMedia('(prefers-color-scheme: dark)').addEventListener('change', schedule);
new MutationObserver(schedule).observe(document.documentElement, { attributes: true, attributeFilter: ['data-theme'] });
if (document.fonts && document.fonts.ready) document.fonts.ready.then(schedule);

try { const m = localStorage.getItem('pde-lab-mode'); if (m === 'class' || m === 'explore') P.mode = m; } catch (e) {}
syncMode();
P.show = show; P.render = render; P.initAll = () => ORDER.forEach(([id]) => { if (!inited[id]) init(id); });
show((location.hash || '').slice(1) || 'classify');
})();
