"""Layout tool for Figure 5: place, size and turn the panel images and the GPCRdb number circles.

    python build_fig5_layout_tool.py   (writes fig4_layout_tool.html; open it from this folder)

The page draws the figure background of make_fig5.py (header and helix pills) and the panel
renders from fig4_panels/ (linked, so re-rendered panels show on reload). Drag panels and circles;
set size, angle and text in the side bar. Motif names (bold, one shared font size) can be moved,
turned, edited, added and deleted the same way; so can the double-headed distance arrows (drag the
line, its ends or its label). "Download fig4_layout.json" saves the layout; put the file
in this folder and run make_fig5.py to build fig5_structure.svg from it.
"""
import json

import make_fig5 as F

OUT = F.HERE / "fig4_layout_tool.html"
PAD = 120   # room shown around the template (dashed frame), so items placed outside stay reachable

PAGE = r"""<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1">
<title>Fig 4 Layout</title>
<style>
:root { --bg: #ffffff; --fg: #1f2328; --muted: #656d76; --card: #f6f8fa; --line: #d0d7de; --accent: #3d4a6b; --sel: #d1242f; }
@media (prefers-color-scheme: dark) { :root:not([data-theme="light"]) { --bg: #0d1117; --fg: #e6edf3; --muted: #8d96a0; --card: #161b22; --line: #30363d; --accent: #9aa8d6; --sel: #ff7b72; } }
:root[data-theme="dark"] { --bg: #0d1117; --fg: #e6edf3; --muted: #8d96a0; --card: #161b22; --line: #30363d; --accent: #9aa8d6; --sel: #ff7b72; }
* { box-sizing: border-box; }
body { margin: 0; background: var(--bg); color: var(--fg); font: 13px/1.4 Arial, sans-serif; display: flex; flex-wrap: wrap; min-height: 100vh; }
main { flex: 1 1 600px; padding: 16px; min-width: 0; }
aside { flex: 0 0 300px; max-width: 100%; padding: 16px; border-left: 1px solid var(--line); background: var(--card); overflow-y: auto; max-height: 100vh; }
svg#fig { width: 100%; height: auto; background: #fff; border: 1px solid var(--line); border-radius: 6px; touch-action: none; user-select: none; }
h1 { font-size: 16px; margin: 0 0 8px; }
h2 { font-size: 13px; margin: 16px 0 6px; color: var(--muted); text-transform: uppercase; letter-spacing: .04em; }
.bar { display: flex; flex-wrap: wrap; gap: 6px; margin-bottom: 10px; }
button, input, select { font: inherit; color: var(--fg); background: var(--bg); border: 1px solid var(--line); border-radius: 6px; padding: 4px 8px; }
button { cursor: pointer; }
button.primary { background: var(--accent); color: var(--bg); border-color: var(--accent); }
label { display: grid; grid-template-columns: 70px 1fr 64px; gap: 6px; align-items: center; margin: 4px 0; }
label.check { grid-template-columns: 70px 1fr; }
input[type=range] { width: 100%; padding: 0; }
input[type=number] { width: 64px; }
.list { display: flex; flex-wrap: wrap; gap: 4px; }
.list button { padding: 2px 6px; font-size: 12px; }
.list button.on { border-color: var(--sel); color: var(--sel); }
.list button.off { opacity: .5; text-decoration: line-through; }
.hint { color: var(--muted); font-size: 12px; }
.panel image { cursor: move; }
.badge { cursor: move; }
.outline { fill: none; stroke: var(--sel); stroke-dasharray: 4 3; stroke-width: 1; pointer-events: none; }
</style></head><body>
<main>
<h1>Figure 4 layout</h1>
<div class="bar">
  <button class="primary" id="download">Download fig4_layout.json</button>
  <button id="load">Load layout file</button><input type="file" id="file" accept=".json" hidden>
  <button id="reset">Reset to defaults</button>
  <label class="check" style="display:flex;gap:6px;margin:0"><input type="checkbox" id="outlines"> show image frames</label>
</div>
<svg id="fig" viewBox="-__PAD__ -__PAD__ __PW__ __PH__" xmlns="http://www.w3.org/2000/svg">
  <rect x="-__PAD__" y="-__PAD__" width="__PW__" height="__PH__" fill="white"/>
  <rect x="0" y="0" width="__W__" height="__H__" fill="none" stroke="#c0c0c0" stroke-dasharray="6 4"/>
  <g id="bg">__BG__</g>
  <g id="panels"></g>
  <g id="arrows"></g>
  <g id="badges"></g>
  <g id="dists"></g>
  <g id="labels"></g>
  <g id="header">__HEADER__</g>
</svg>
<p class="hint">The dashed frame is the template; the saved figure grows to include anything placed outside it. Drag panels and circles to move them. Click one to select it, then use the side bar, or the keys:
arrows move (shift: 10), + / - resize a panel, [ / ] turn a panel or name by 1 degree (shift: 10), Delete removes the selected panel, circle or name (Reset goes back to the saved fig4_layout.json).
Names share one font size (side bar).
The layout is kept in this browser; download the file into new_pipeline/step12_structure_superposition/ and run make_fig4.py.</p>
</main>
<aside>
<h2>Selected</h2>
<div id="props"><p class="hint">Nothing selected.</p></div>
<h2>Panels</h2>
<div class="list" id="plist"></div>
<h2>Number circles</h2>
<div class="list" id="blist"></div>
<div class="bar" style="margin-top:8px"><button id="addbadge">Add circle</button></div>
<label>radius<input type="range" id="br" min="8" max="30" step="0.5"><input type="number" id="brn" step="0.5"></label>
<label>font<input type="range" id="bf" min="6" max="20" step="0.5"><input type="number" id="bfn" step="0.5"></label>
<h2>Distances</h2>
<div class="list" id="dlist"></div>
<label>font<input type="range" id="df" min="6" max="20" step="0.5"><input type="number" id="dfn" step="0.5"></label>
<h2>Distance arrows</h2>
<div class="list" id="alist"></div>
<div class="bar" style="margin-top:8px"><button id="addarrow">Add arrow</button></div>
<h2>Motif names</h2>
<div class="list" id="llist"></div>
<div class="bar" style="margin-top:8px"><button id="addlabel">Add name</button></div>
<label>font<input type="range" id="lf" min="8" max="30" step="0.5"><input type="number" id="lfn" step="0.5"></label>
</aside>
<script>
const DEFAULT = __DEFAULT__;
const INITIAL = __INITIAL__;
const KEY = "fig4_layout";
const SUBRATIO = __SUBRATIO__;
const DIST = __DIST__, DIST_PT = __DIST_PT__, DIST_FILL = "__DIST_FILL__", DIST_BG_OPACITY = __DIST_BG__, DIST_PAD = __DIST_PAD__;   // distance label under each circle (make_fig4.distances)
const ARROW_W = __ARROW_W__, ARROW_HEAD = __ARROW_HEAD__, ARROW_GAP = __ARROW_GAP__;   // as make_fig4.arrow_shape, arrow_label_at
const NS = "http://www.w3.org/2000/svg";
const svg = document.getElementById("fig");
const clone = o => JSON.parse(JSON.stringify(o));
let L;
try { L = JSON.parse(localStorage.getItem(KEY)); } catch (e) {}
L = merge(L || clone(INITIAL));
let sel = null;   // {kind: "panel", key} or {kind: "badge", i}
const stamp = Date.now();

function merge(l) {   // panels added since the layout was saved get their defaults
  for (const [k, v] of Object.entries(DEFAULT.panels)) if (!l.panels[k]) l.panels[k] = clone(v);
  for (const k of Object.keys(l.panels)) if (!DEFAULT.panels[k]) delete l.panels[k];
  if (l.badge_r == null) l.badge_r = DEFAULT.badge_r;
  if (l.badge_pt == null) l.badge_pt = DEFAULT.badge_pt;
  if (!l.labels) l.labels = clone(DEFAULT.labels);
  else for (const t of DEFAULT.labels) if ((t.since || 1) > (l.version || 1)) l.labels.push(clone(t));   // labels added since
  l.version = DEFAULT.version;
  if (l.label_pt == null) l.label_pt = DEFAULT.label_pt;
  if (!l.dists)   // saved before the distance labels: put them under their circles
    l.dists = l.badges.map(b => ({ key: b.text + "|" + (b.sub || ""), x: b.x, y: b.y + l.badge_r + DIST_PT * 0.7 })).filter(d => DIST[d.key]);
  if (l.dist_pt == null) l.dist_pt = DEFAULT.dist_pt;
  if (!l.arrows) l.arrows = clone(DEFAULT.arrows);
  for (const a of l.arrows) a.y2 = a.y1;   // arrows are horizontal
  return l;
}
function save() { try { localStorage.setItem(KEY, JSON.stringify(L)); } catch (e) {} }
function el(tag, attrs) { const e = document.createElementNS(NS, tag); for (const k in attrs) e.setAttribute(k, attrs[k]); return e; }
function pt(evt) { const p = svg.createSVGPoint(); p.x = evt.clientX; p.y = evt.clientY; return p.matrixTransform(svg.getScreenCTM().inverse()); }

// drag handles for an arrow: one end (1 or 2) or the whole arrow (its label goes with it)
function handle(a, which) {
  if (which === "all") return {
    get x() { return a.x1; }, set x(v) { const d = v - a.x1; a.x1 += d; a.x2 += d; },
    get y() { return a.y1; }, set y(v) { a.y1 = a.y2 = v; } };
  // an end: the arrow stays horizontal, so an end only changes the length
  const kx = "x" + which;
  return { get x() { return a[kx]; }, set x(v) { a[kx] = v; }, get y() { return a.y1; }, set y(v) {} };
}
function drawArrows() {
  const ag = document.getElementById("arrows"); ag.innerHTML = "";
  L.arrows.forEach((a, i) => {
    const on = sel && sel.kind === "arrow" && sel.i === i, col = on ? "#d1242f" : DIST_FILL;
    const n = Math.max(Math.hypot(a.x2 - a.x1, a.y2 - a.y1), 1e-6), ux = (a.x2 - a.x1) / n, uy = (a.y2 - a.y1) / n, hw = ARROW_HEAD * 0.45;
    const g = el("g", { class: "badge" }), base = [];
    for (const [x, y, s] of [[a.x1, a.y1, 1], [a.x2, a.y2, -1]]) {
      const bx = x + s * ux * ARROW_HEAD, by = y + s * uy * ARROW_HEAD;
      base.push([bx, by]);
      g.appendChild(el("polygon", { points: `${x},${y} ${bx - uy * hw},${by + ux * hw} ${bx + uy * hw},${by - ux * hw}`, fill: col }));
    }
    g.appendChild(el("line", { x1: base[0][0], y1: base[0][1], x2: base[1][0], y2: base[1][1], stroke: col, "stroke-width": ARROW_W }));
    g.appendChild(el("line", { x1: a.x1, y1: a.y1, x2: a.x2, y2: a.y2, stroke: "transparent", "stroke-width": 10 }));   // easier to grab
    drag(g, () => handle(a, "all"), () => select({ kind: "arrow", i }));
    ag.appendChild(g);
    if (on) for (const w of [1, 2]) {
      const c = el("circle", { cx: a["x" + w], cy: a["y" + w], r: 5, fill: "none", stroke: "#d1242f", "stroke-width": 1.5, class: "badge", "pointer-events": "all" });
      drag(c, () => handle(a, w), () => {});
      ag.appendChild(c);
    }
    if (a.text) {
      const lg = el("g", { class: "badge" }), box = el("rect", { fill: "white", "fill-opacity": DIST_BG_OPACITY });
      const lx = (a.x1 + a.x2) / 2, ly = a.y1 - ARROW_W / 2 - ARROW_GAP - (0.8 * L.dist_pt + 2 * DIST_PAD) / 2;
      const t = el("text", { x: lx, y: ly + 0.35 * L.dist_pt, "text-anchor": "middle", "font-family": "Arial",
        "font-size": L.dist_pt, "font-weight": "bold", fill: col });
      t.textContent = a.text;
      lg.appendChild(box); lg.appendChild(t);
      drag(lg, () => handle(a, "all"), () => select({ kind: "arrow", i }));
      ag.appendChild(lg);
      const w = t.getBBox().width + 2 * DIST_PAD, h = 0.8 * L.dist_pt + 2 * DIST_PAD;
      Object.entries({ x: lx - w / 2, y: ly - h / 2, width: w, height: h }).forEach(([k, v]) => box.setAttribute(k, v));
    }
  });
}

function draw() {
  const pg = document.getElementById("panels"), bg = document.getElementById("badges"), lg = document.getElementById("labels"), dg = document.getElementById("dists");
  pg.innerHTML = ""; bg.innerHTML = ""; lg.innerHTML = ""; dg.innerHTML = "";
  const frames = document.getElementById("outlines").checked;
  for (const [k, v] of Object.entries(L.panels)) {
    if (v.hidden || v.deleted) continue;
    const g = el("g", { class: "panel", transform: `rotate(${v.angle} ${v.x} ${v.y})` });
    const img = el("image", { href: `fig4_panels/${encodeURIComponent(k)}.png?v=${stamp}`, x: v.x - v.size / 2, y: v.y - v.size / 2, width: v.size, height: v.size });
    g.appendChild(img);
    if (frames || (sel && sel.kind === "panel" && sel.key === k))
      g.appendChild(el("rect", { class: "outline", x: v.x - v.size / 2, y: v.y - v.size / 2, width: v.size, height: v.size }));
    drag(img, () => v, () => select({ kind: "panel", key: k }));
    pg.appendChild(g);
  }
  L.badges.forEach((b, i) => {
    const g = el("g", { class: "badge" });
    const on = sel && sel.kind === "badge" && sel.i === i;
    g.appendChild(el("circle", { cx: b.x, cy: b.y, r: L.badge_r, fill: "white", stroke: on ? "#d1242f" : "#222", "stroke-width": on ? 2 : 1.2 }));
    const t = el("text", { x: b.x, y: b.y + L.badge_pt * 0.35, "text-anchor": "middle", "font-family": "Arial", "font-size": L.badge_pt });
    t.textContent = b.text;
    if (b.sub) { const s = el("tspan", { dy: L.badge_pt * 0.3, "font-size": L.badge_pt * SUBRATIO }); s.textContent = b.sub; t.appendChild(s); }
    g.appendChild(t);
    drag(g, () => b, () => select({ kind: "badge", i }));
    bg.appendChild(g);
  });
  L.dists.forEach((d, i) => {
    if (!DIST[d.key]) return;
    const on = sel && sel.kind === "dist" && sel.i === i;
    const e = el("text", { x: d.x, y: d.y + 0.35 * L.dist_pt, "text-anchor": "middle", "font-family": "Arial",
      "font-size": L.dist_pt, "font-weight": "bold", fill: on ? "#d1242f" : DIST_FILL, class: "badge" });
    e.textContent = DIST[d.key];
    const g = el("g", { class: "badge" });
    const box = el("rect", { fill: "white", "fill-opacity": DIST_BG_OPACITY });
    g.appendChild(box); g.appendChild(e);
    drag(g, () => d, () => select({ kind: "dist", i }));
    dg.appendChild(g);
    const w = e.getBBox().width + 2 * DIST_PAD, h = 0.8 * L.dist_pt + 2 * DIST_PAD;   // as make_fig4.dists
    Object.entries({ x: d.x - w / 2, y: d.y - h / 2, width: w, height: h }).forEach(([k, v]) => box.setAttribute(k, v));
  });
  L.labels.forEach((t, i) => {
    const on = sel && sel.kind === "label" && sel.i === i;
    const e = el("text", { x: t.x, y: t.y + 0.35 * L.label_pt, "text-anchor": "middle", "font-family": "Arial",
      "font-size": L.label_pt, "font-weight": "bold", fill: on ? "#d1242f" : "#000", class: "badge",
      transform: `rotate(${t.angle} ${t.x} ${t.y})` });
    e.textContent = t.text;
    drag(e, () => t, () => select({ kind: "label", i }));
    lg.appendChild(e);
  });
  drawArrows();
  lists();
}

function drag(node, obj, onSelect) {
  node.addEventListener("pointerdown", e => {
    e.preventDefault(); e.stopPropagation(); onSelect();
    const o = obj(), p0 = pt(e), x0 = o.x, y0 = o.y;
    const move = ev => { const p = pt(ev); o.x = Math.round((x0 + p.x - p0.x) * 2) / 2; o.y = Math.round((y0 + p.y - p0.y) * 2) / 2; draw(); props(); };
    const up = () => { window.removeEventListener("pointermove", move); window.removeEventListener("pointerup", up); save(); };
    window.addEventListener("pointermove", move); window.addEventListener("pointerup", up);
  });
}

function select(s) { sel = s; draw(); props(); }
svg.addEventListener("pointerdown", () => { sel = null; draw(); props(); });

function field(parent, name, obj, key, min, max, step, onchange) {
  const lab = document.createElement("label");
  lab.innerHTML = `${name}<input type="range" min="${min}" max="${max}" step="${step}"><input type="number" step="${step}">`;
  const [r, n] = lab.querySelectorAll("input");
  r.value = n.value = obj[key];
  const set = v => { obj[key] = +v; r.value = n.value = obj[key]; draw(); save(); if (onchange) onchange(); };
  r.oninput = () => set(r.value); n.onchange = () => set(n.value);
  parent.appendChild(lab);
}
function textField(parent, name, obj, key) {
  const lab = document.createElement("label");
  lab.className = "check";
  lab.innerHTML = `${name}<input type="text">`;
  const t = lab.querySelector("input"); t.value = obj[key] || "";
  t.oninput = () => { obj[key] = t.value; draw(); save(); };
  parent.appendChild(lab);
}

function props() {
  const box = document.getElementById("props");
  box.innerHTML = "";
  if (!sel) { box.innerHTML = '<p class="hint">Nothing selected.</p>'; return; }
  if (sel.kind === "panel") {
    const v = L.panels[sel.key];
    box.insertAdjacentHTML("beforeend", `<p><b>${sel.key}</b></p>`);
    field(box, "x", v, "x", -__PAD__, __W__ + __PAD__, 0.5);
    field(box, "y", v, "y", -__PAD__, __H__ + __PAD__, 0.5);
    field(box, "size", v, "size", 40, 500, 1);
    field(box, "angle", v, "angle", -180, 180, 1);
    const lab = document.createElement("label"); lab.className = "check";
    lab.innerHTML = `hidden<input type="checkbox" ${v.hidden ? "checked" : ""}>`;
    lab.querySelector("input").onchange = e => { v.hidden = e.target.checked; if (v.hidden) sel = null; draw(); props(); save(); };
    box.appendChild(lab);
    const del = document.createElement("button"); del.textContent = "Delete panel";
    del.onclick = () => { v.deleted = true; sel = null; draw(); props(); save(); };
    box.appendChild(del);
  } else if (sel.kind === "label") {
    const t = L.labels[sel.i];
    textField(box, "text", t, "text");
    field(box, "x", t, "x", -__PAD__, __W__ + __PAD__, 0.5);
    field(box, "y", t, "y", -__PAD__, __H__ + __PAD__, 0.5);
    field(box, "angle", t, "angle", -180, 180, 1);
    const del = document.createElement("button"); del.textContent = "Delete name";
    del.onclick = () => { L.labels.splice(sel.i, 1); sel = null; draw(); props(); save(); };
    box.appendChild(del);
  } else if (sel.kind === "dist") {
    const d = L.dists[sel.i];
    box.insertAdjacentHTML("beforeend", `<p><b>${DIST[d.key]}</b> (circle ${d.key.replace("|", " ")})</p>`);
    field(box, "x", d, "x", -__PAD__, __W__ + __PAD__, 0.5);
    field(box, "y", d, "y", -__PAD__, __H__ + __PAD__, 0.5);
    const del = document.createElement("button"); del.textContent = "Delete distance";
    del.onclick = () => { L.dists.splice(sel.i, 1); sel = null; draw(); props(); save(); };
    box.appendChild(del);
  } else if (sel.kind === "arrow") {
    const a = L.arrows[sel.i];
    box.insertAdjacentHTML("beforeend", '<p class="hint">Horizontal, label centred above it. Drag the line or label to move it, a red ring to change the length; or set length below.</p>');
    textField(box, "label", a, "text");
    a.y2 = a.y1;   // horizontal
    const geo = { get x() { return a.x1; }, set x(v) { a.x2 += v - a.x1; a.x1 = v; },   // left end; keeps the length
      get y() { return a.y1; }, set y(v) { a.y1 = a.y2 = v; },
      get length() { return Math.round((a.x2 - a.x1) * 2) / 2; }, set length(v) { a.x2 = a.x1 + v; } };
    field(box, "x (left)", geo, "x", -__PAD__, __W__ + __PAD__, 0.5);
    field(box, "y", geo, "y", -__PAD__, __H__ + __PAD__, 0.5);
    field(box, "length", geo, "length", 10, __W__, 0.5);
    const del = document.createElement("button"); del.textContent = "Delete arrow";
    del.onclick = () => { L.arrows.splice(sel.i, 1); sel = null; draw(); props(); save(); };
    box.appendChild(del);
  } else {
    const b = L.badges[sel.i];
    textField(box, "number", b, "text");
    textField(box, "subscript", b, "sub");
    field(box, "x", b, "x", -__PAD__, __W__ + __PAD__, 0.5);
    field(box, "y", b, "y", -__PAD__, __H__ + __PAD__, 0.5);
    const del = document.createElement("button"); del.textContent = "Delete circle";
    del.onclick = () => { L.badges.splice(sel.i, 1); sel = null; draw(); props(); save(); };
    box.appendChild(del);
  }
}

function lists() {
  const pl = document.getElementById("plist"), bl = document.getElementById("blist"), ll = document.getElementById("llist");
  pl.innerHTML = ""; bl.innerHTML = ""; ll.innerHTML = "";
  L.labels.forEach((x, i) => {
    const b = document.createElement("button"); b.textContent = x.text || "(empty)";
    if (sel && sel.kind === "label" && sel.i === i) b.className = "on";
    b.onclick = () => select({ kind: "label", i });
    ll.appendChild(b);
  });
  for (const [k, v] of Object.entries(L.panels)) {
    if (v.deleted) continue;
    const b = document.createElement("button"); b.textContent = k;
    b.className = (sel && sel.kind === "panel" && sel.key === k ? "on " : "") + (v.hidden ? "off" : "");
    b.onclick = () => { if (v.hidden) { v.hidden = false; save(); } select({ kind: "panel", key: k }); };
    pl.appendChild(b);
  }
  L.badges.forEach((x, i) => {
    const b = document.createElement("button"); b.textContent = x.text + (x.sub ? "_" + x.sub : "");
    if (sel && sel.kind === "badge" && sel.i === i) b.className = "on";
    b.onclick = () => select({ kind: "badge", i });
    bl.appendChild(b);
  });
  const al = document.getElementById("alist"); al.innerHTML = "";
  L.arrows.forEach((x, i) => {
    const b = document.createElement("button"); b.textContent = x.text || "(arrow)";
    if (sel && sel.kind === "arrow" && sel.i === i) b.className = "on";
    b.onclick = () => select({ kind: "arrow", i });
    al.appendChild(b);
  });
  const dl = document.getElementById("dlist"); dl.innerHTML = "";
  L.dists.forEach((x, i) => {
    const b = document.createElement("button"); b.textContent = x.key.replace("|", "_").replace(/_$/, "");
    if (sel && sel.kind === "dist" && sel.i === i) b.className = "on";
    b.onclick = () => select({ kind: "dist", i });
    dl.appendChild(b);
  });
}

function globals() {
  for (const [id, key] of [["br", "badge_r"], ["bf", "badge_pt"], ["lf", "label_pt"], ["df", "dist_pt"]]) {
    const r = document.getElementById(id), n = document.getElementById(id + "n");
    r.value = n.value = L[key];
    const set = v => { L[key] = +v; r.value = n.value = L[key]; draw(); save(); };
    r.oninput = () => set(r.value); n.onchange = () => set(n.value);
  }
}

document.addEventListener("keydown", e => {
  if (!sel || ["INPUT", "TEXTAREA"].includes(document.activeElement.tagName)) return;
  const lists_ = { label: L.labels, badge: L.badges, dist: L.dists, arrow: L.arrows };
  const o = sel.kind === "panel" ? L.panels[sel.key] : sel.kind === "arrow" ? handle(L.arrows[sel.i], "all") : lists_[sel.kind][sel.i];
  const step = e.shiftKey ? 10 : 1;
  const moves = { ArrowLeft: [-step, 0], ArrowRight: [step, 0], ArrowUp: [0, -step], ArrowDown: [0, step] };
  if (moves[e.key]) { o.x += moves[e.key][0]; o.y += moves[e.key][1]; }
  else if (sel.kind === "panel" && (e.key === "+" || e.key === "=")) o.size += step;
  else if (sel.kind === "panel" && e.key === "-") o.size = Math.max(20, o.size - step);
  else if ("angle" in o && (e.key === "]" || e.key === "}")) o.angle += step;
  else if ("angle" in o && (e.key === "[" || e.key === "{")) o.angle -= step;
  else if (e.key === "Delete" || e.key === "Backspace") {
    if (sel.kind === "panel") o.deleted = true; else lists_[sel.kind].splice(sel.i, 1);
    sel = null;
  }
  else return;
  e.preventDefault(); draw(); props(); save();
});

document.getElementById("addbadge").onclick = () => { L.badges.push({ text: "1x50", sub: "", x: __W__ / 2, y: __H__ / 2 }); select({ kind: "badge", i: L.badges.length - 1 }); save(); };
document.getElementById("addlabel").onclick = () => { L.labels.push({ text: "Motif", x: __W__ / 2, y: __H__ / 2, angle: 0 }); select({ kind: "label", i: L.labels.length - 1 }); save(); };
document.getElementById("addarrow").onclick = () => { L.arrows.push({ text: "4.0 Å", x1: __W__ / 2 - 60, y1: __H__ / 2, x2: __W__ / 2 + 60, y2: __H__ / 2 }); select({ kind: "arrow", i: L.arrows.length - 1 }); save(); };
document.getElementById("outlines").onchange = draw;
document.getElementById("download").onclick = () => {
  const a = document.createElement("a");
  a.href = URL.createObjectURL(new Blob([JSON.stringify(L, null, 1)], { type: "application/json" }));
  a.download = "fig4_layout.json"; a.click();
};
document.getElementById("load").onclick = () => document.getElementById("file").click();
document.getElementById("file").onchange = e => {
  const f = e.target.files[0]; if (!f) return;
  f.text().then(t => { L = merge(JSON.parse(t)); sel = null; save(); globals(); draw(); props(); });
};
document.getElementById("reset").onclick = () => {
  if (!confirm("Discard this layout and go back to the saved defaults (fig4_layout.json)?")) return;
  L = clone(DEFAULT); sel = null; save(); globals(); draw(); props();
};
globals(); draw();
</script></body></html>
"""


def main():
    default = F.layout()      # the saved fig4_layout.json is the starting point (and what Reset returns to)
    page = (PAGE.replace("__BG__", "".join(F.pills()))
            .replace("__HEADER__", "".join(F.header()))
            .replace("__DEFAULT__", json.dumps(default))
            .replace("__INITIAL__", json.dumps(F.layout()))
            .replace("__SUBRATIO__", str(F.BADGE_SUB_PT / F.BADGE_PT))
            .replace("__DIST__", json.dumps(F.distances(), ensure_ascii=False))
            .replace("__DIST_PT__", str(F.DIST_PT)).replace("__DIST_FILL__", F.DIST_FILL)
            .replace("__DIST_BG__", str(F.DIST_BG_OPACITY)).replace("__DIST_PAD__", str(F.DIST_PAD))
            .replace("__ARROW_W__", str(F.ARROW_W)).replace("__ARROW_HEAD__", str(F.ARROW_HEAD))
            .replace("__ARROW_GAP__", str(F.ARROW_GAP))
            .replace("__PW__", str(F.W + 2 * PAD)).replace("__PH__", str(F.H + 2 * PAD)).replace("__PAD__", str(PAD))
            .replace("__W__", str(F.W)).replace("__H__", str(F.H)))
    OUT.write_text(page, encoding="utf-8")
    print("wrote", OUT, len(default["panels"]), "panels,", len(default["badges"]), "circles")


if __name__ == "__main__":
    main()
