"""Build docs/briefing/briefing.html, the research briefing for a broad scientific audience.

    python -X utf8 -B docs/briefing/build_briefing.py

Every number in the text and figures is read from committed result files (paths in SRC below), so the briefing can be
rebuilt after any change to the results. Figures are inline SVG that use the page's colour tokens (light and dark).
"""
from __future__ import annotations

import csv
import json
import subprocess
from html import escape
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
OUT = Path(__file__).resolve().parent / "briefing.html"


def J(rel):
    return json.loads((ROOT / rel).read_text(encoding="utf-8"))


def rows(rel):
    with open(ROOT / rel, encoding="utf-8", newline="") as f:
        return list(csv.DictReader(f))


# ---------------------------------------------------------------- data
bench = J("results/benchmark/summary.json")
gate1 = J("results/gates/G_pc_attempt1.json")["summary"]
gate2 = J("results/gates/G_pc_attempt2.json")["summary"]
p1 = J("results/replication/p1_summary.json")["strata"]
p1s = J("results/replication/p1_secondary.json")
xlib = {r["stratum"]: r for r in rows("results/benchmark/replication_xlib.csv") if r["u_rule"] == "all_D"}
null = J("results/null/summary.json")
nc = J("results/controls/nc_real.json")
ceil = J("results/clinical/ceiling_summary.json")
cgt2 = J("results/clinical/cgt_secondary.json")
q4 = J("results/clinical/q4_summary.json")["results"]
comb = J("results/combinatorial/validation.json")
rc = J("results/report_counts.json")
lim = J("results/limits/limits_summary.json")
cards = J("results/candidates/cards_summary.json")
sub = J("results/subtypes/summary.json")["by_list"]
lin = J("results/lineage/summary.json")
ind = J("results/gates/G_indep.json")
suite = J("results/toy/suite_run_final.json")
toy = J("results/toy/diagnostics_8cffc91.json")
pool = J("results/benchmark/p1_pool.json")
short = rows("results/candidates/shortlist.csv")
devs = rows("results/deviations.csv")
_d = subprocess.run(["git", "log", "--format=%at"], cwd=ROOT, capture_output=True, text=True).stdout.split()
wall_h = round((int(_d[0]) - int(_d[-1])) / 3600)
n_commits = subprocess.run(["git", "rev-list", "--count", "HEAD"], cwd=ROOT, capture_output=True, text=True).stdout.strip()

ARMS = ["S1-obs", "S1-pred", "EN-obs", "EN-pred"]
ridge = [(t["MD-strict"]["fdp_true"], t["S1-obs"]["fdp_true"], t["EN-pred selected"]["fdp_true"])
         for t in toy["e_ridge_md_strict"]]
ridge_g = [t["MD-global"]["fdp_true"] for t in toy["d_md_global"]["ridge"]]
recall_g = [t["MD-global"]["recall_eligible_plants"] for t in toy["d_md_global"]["ridge"]]

G = p1["G"]["rates"]
con = {c["contrast"]: c for c in p1["G"]["contrasts"]}
diag = p1s["bootstrap_diagnostic"]["G"]["contrasts"]
sk = p1s["secondary_k"]
rate_k = {arm: [sk["k50"][arm]["rate"], G[arm], sk["k200"][arm]["rate"], sk["k500"][arm]["rate"]] for arm in ARMS}
KS = [50, 100, 200, 500]

fdr = null["fdr_emp_fixed"]
rds = cgt2["recall_data_supported"]["non_proxy|all|k200"]
head = q4["non_proxy|gate_disjoint|k200"]
allnp = q4["non_proxy|all|k200"]
ref = cgt2["refuted_vs_approved"]
cm = rc["combinatorial"]
fam = rc["shortlist"]["by_feature_family"]
lbc = rc["lineage_by_class"]
pcg = lim["pair_classes_md_global"]
setA = {r["symbol"]: float(r["setA_dependent_frac"]) for r in short if r["setA_dependent_frac"]}


def pct(x, d=0):
    return f"{100 * x:.{d}f}%"


def f2(x):
    return f"{x:.2f}"


# ---------------------------------------------------------------- SVG helpers
FONT = 'font-family="var(--font-ui)"'


def svg(w, h, body, label):
    return (f'<svg viewBox="0 0 {w} {h}" role="img" aria-label="{escape(label)}" class="fig-svg">'
            f'{body}</svg>')


def t(x, y, s, cls="lab", anchor="start", extra=""):
    return f'<text x="{x:.1f}" y="{y:.1f}" class="{cls}" text-anchor="{anchor}" {extra}>{escape(str(s))}</text>'


ARM_CLS = {"S1-obs": "c1", "S1-pred": "c2", "EN-obs": "c3", "EN-pred": "c4", "MD-global": "c5", "MD-strict": "c6"}


def fig_design():
    """Schematic of the analysis: predictor → four arms → catalogue → freeze → single replication."""
    W, H = 760, 330
    b = []

    def box(x, y, w, h, title, sub, cls="box"):
        b.append(f'<rect x="{x}" y="{y}" width="{w}" height="{h}" rx="6" class="{cls}"/>')
        b.append(t(x + w / 2, y + 22, title, "boxt", "middle"))
        for i, line in enumerate(sub):
            b.append(t(x + w / 2, y + 40 + 15 * i, line, "boxs", "middle"))

    def arrow(x1, y1, x2, y2):
        b.append(f'<line x1="{x1}" y1="{y1}" x2="{x2}" y2="{y2}" class="arrow" marker-end="url(#ah)"/>')

    b.append('<defs><marker id="ah" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="7" markerHeight="7" '
             'orient="auto-start-reverse"><path d="M0,0 L10,5 L0,10 z" class="ahead"/></marker></defs>')
    box(10, 20, 170, 78, "Observed effects", ["CRISPR knockout scores", "1,022 dev lines (D)"])
    box(10, 120, 170, 78, "Predicted effects", ["Part 1 predictor,", "out-of-fold"])
    box(10, 225, 170, 78, "Feature dictionary", [f"{bench['n_rows_arms']:,} gene × feature", "tests (mutations, CN, …)"])
    box(250, 20, 190, 78, "Univariate test (S1)", ["one feature at a time,", "confounders adjusted"], "box b1")
    box(250, 120, 190, 78, "Elastic net (EN)", ["many features jointly,", "stability selection"], "box b3")
    b.append(t(345, 220, "each run on observed AND", "boxs", "middle"))
    b.append(t(345, 236, "predicted effects = 4 arms", "boxs", "middle"))
    box(510, 20, 240, 78, "Catalogue MD-global", ["EN-pred ∩ S1-obs hits", "(same sign)"], "box b5")
    box(510, 120, 240, 78, "Freeze, then read R once", ["255 held-out lines", "top-k lists replicated?"], "box b4")
    box(510, 225, 240, 78, "Checks", ["toy world · controls · null", "clinical · combinatorial screens"])
    for y in (59, 159, 264):
        arrow(180, y, 248, 59 if y == 59 else 159)
    arrow(180, 264, 248, 80)
    arrow(180, 264, 248, 178)
    arrow(440, 59, 508, 59)
    arrow(440, 159, 508, 80)
    arrow(630, 98, 630, 118)
    arrow(630, 198, 630, 223)
    return svg(W, H, "".join(b), "Study design schematic")


def dotstrip(series, lo, hi, limit, limit_label, xlabel, label, fmt=lambda v: f"{v:.2f}"):
    """Horizontal strip plot: one row per series (name, values, css class)."""
    W = 720
    L, R, top, rowh = 170, 30, 30, 46
    H = top + rowh * len(series) + 46
    sx = lambda v: L + (v - lo) / (hi - lo) * (W - L - R)
    b = []
    for i in range(6):
        v = lo + (hi - lo) * i / 5
        b.append(f'<line x1="{sx(v):.1f}" y1="{top - 8}" x2="{sx(v):.1f}" y2="{H - 40}" class="grid"/>')
        b.append(t(sx(v), H - 24, fmt(v), "tick", "middle"))
    b.append(f'<line x1="{sx(limit):.1f}" y1="{top - 14}" x2="{sx(limit):.1f}" y2="{H - 40}" class="limit"/>')
    b.append(t(sx(limit) + 6, top - 4, limit_label, "limitlab"))
    for j, (name, vals, cls) in enumerate(series):
        y = top + 18 + rowh * j
        b.append(t(L - 12, y + 4, name, "lab", "end"))
        for k, v in enumerate(vals):
            jitter = (k - (len(vals) - 1) / 2) * 3.0
            b.append(f'<circle cx="{sx(v):.1f}" cy="{y + jitter:.1f}" r="6" class="dot {cls}"/>')
    b.append(t(L + (W - L - R) / 2, H - 6, xlabel, "axis", "middle"))
    return svg(W, H, "".join(b), label)


def hbars(items, hi, limit, limit_label, xlabel, label, fmt=lambda v: f"{v:.2f}", note=None):
    """Horizontal bars: items = (name, value, css class, value text)."""
    W = 720
    L, R, top, rowh = 150, 70, 28, 34
    H = top + rowh * len(items) + 44
    sx = lambda v: L + v / hi * (W - L - R)
    b = []
    for i in range(6):
        v = hi * i / 5
        b.append(f'<line x1="{sx(v):.1f}" y1="{top - 6}" x2="{sx(v):.1f}" y2="{H - 38}" class="grid"/>')
        b.append(t(sx(v), H - 22, fmt(v), "tick", "middle"))
    for j, (name, v, cls, txt) in enumerate(items):
        y = top + rowh * j
        b.append(t(L - 10, y + 17, name, "lab", "end"))
        b.append(f'<rect x="{L}" y="{y + 4}" width="{max(sx(v) - L, 1.5):.1f}" height="20" rx="2" class="bar {cls}"/>')
        b.append(t(sx(v) + 6, y + 18, txt, "val"))
    if limit is not None:
        b.append(f'<line x1="{sx(limit):.1f}" y1="{top - 12}" x2="{sx(limit):.1f}" y2="{H - 38}" class="limit"/>')
        b.append(t(sx(limit) + 6, top - 2, limit_label, "limitlab"))
    b.append(t(L + (W - L - R) / 2, H - 4, xlabel, "axis", "middle"))
    return svg(W, H, "".join(b), label)


def fig_k_curve():
    W, H = 720, 360
    L, R, T, B = 70, 150, 24, 56
    lo, hi = 0.55, 1.0
    xs = {k: L + i * (W - L - R) / 3 for i, k in enumerate(KS)}
    sy = lambda v: T + (hi - v) / (hi - lo) * (H - T - B)
    b = []
    for i in range(10):
        v = lo + 0.05 * i
        b.append(f'<line x1="{L}" y1="{sy(v):.1f}" x2="{W - R}" y2="{sy(v):.1f}" class="grid"/>')
        b.append(t(L - 8, sy(v) + 4, f"{v:.2f}", "tick", "end"))
    for k, x in xs.items():
        b.append(t(x, H - B + 20, f"top {k}", "tick", "middle"))
    b.append(f'<line x1="{xs[100]}" y1="{T}" x2="{xs[100]}" y2="{H - B}" class="limit"/>')
    b.append(t(xs[100] + 5, T + 12, "pre-declared k", "limitlab"))
    ends = []
    for arm in ARMS:
        pts = " ".join(f"{xs[k]:.1f},{sy(v):.1f}" for k, v in zip(KS, rate_k[arm]))
        b.append(f'<polyline points="{pts}" class="line {ARM_CLS[arm]}"/>')
        for k, v in zip(KS, rate_k[arm]):
            b.append(f'<circle cx="{xs[k]:.1f}" cy="{sy(v):.1f}" r="4" class="dot {ARM_CLS[arm]}"/>')
        ends.append([sy(rate_k[arm][-1]), arm])
    ends.sort()
    for i in range(1, len(ends)):                       # keep end labels apart
        if ends[i][0] - ends[i - 1][0] < 15:
            ends[i][0] = ends[i - 1][0] + 15
    for y, arm in ends:
        b.append(t(W - R + 12, y + 4, f"{arm} {rate_k[arm][-1]:.2f}", f"lab {ARM_CLS[arm]}t"))
    b.append(t(L + (W - L - R) / 2, H - 8, "length of each method's ranked candidate list", "axis", "middle"))
    b.append(t(16, T + (H - T - B) / 2, "fraction replicated in R", "axis", "middle",
               f'transform="rotate(-90 16 {T + (H - T - B) / 2:.0f})"'))
    return svg(W, H, "".join(b), "Replication rate in held-out lines against list length")


def fig_families():
    order = sorted(fam.items(), key=lambda kv: -kv[1])[:9]
    names = {"PROG": "expression programme", "EXPPAR": "paralog expression", "COMP": "TP53/RB composite",
             "LOFBI": "loss of function", "SIL-other": "silencing (other)", "CLIN": "clinical feature",
             "HOT": "hotspot mutation", "PSCORE": "pathway score", "LOSS": "shallow copy loss"}
    items = [(names.get(k, k), v, "c5" if k in ("PROG", "EXPPAR", "SIL-other", "PSCORE") else "c3", str(v))
             for k, v in order]
    return hbars(items, 180, None, "", f"number of shortlisted pairs (of {rc['shortlist']['n_rows']})",
                 "Shortlist by determinant type", fmt=lambda v: f"{v:.0f}")


# ---------------------------------------------------------------- figures with data
FIG_TOY = dotstrip(
    [("MD-strict (frozen plan)", [r[0] for r in ridge], "c6"),
     ("MD-global (adopted)", ridge_g, "c5"),
     ("S1-obs", [r[1] for r in ridge], "c1"),
     ("EN-pred", [r[2] for r in ridge], "c4")],
    0, 0.6, 0.15, "limit 0.15", "true false-discovery proportion on 5 simulated data sets with a ridge-type predictor",
    "Toy false discovery by list")

FIG_GATE = hbars(
    [(arm, gate2["by_list"][arm]["recall"], ARM_CLS[arm],
      f"{gate2['by_list'][arm]['k']}/{gate2['by_list'][arm]['n']}")
     for arm in ["MD-strict", "MD-global", "S1-obs", "S1-pred", "EN-obs", "EN-pred"]],
    1.0, 0.6, "pass bar 60%", "fraction of 39 known dependency–biomarker pairs recovered", "Positive-control recall",
    fmt=lambda v: f"{v:.1f}")

FIG_NULL = hbars(
    [(name.replace(" selected", ""), fdr[name], ARM_CLS.get(name.replace(" selected", ""), "c1"), f"{fdr[name]:.3f}")
     for name in ["MD-global", "S1-obs", "MD-strict", "S1-pred", "EN-obs selected", "EN-pred selected"]],
    0.35, 0.10, "target ≤ 0.10", "estimated false-discovery rate from 20 label permutations (300 genes)",
    "Permutation estimate of false discovery", fmt=lambda v: f"{v:.2f}")

FIG_CLIN = hbars(
    [(arm, rds[arm]["recall"], ARM_CLS[arm], f"{rds[arm]['hits']}/{rds[arm]['n']}") for arm in ARMS],
    0.6, None, "", "recall at top 200 of clinical biomarkers that are detectable in the cell-line data",
    "Clinical biomarker recall", fmt=lambda v: f"{v:.1f}")

FIG_COMB = hbars(
    [("MD-global", cm["MD-global"]["rate_pct"] / 100, "c5", f"{cm['MD-global']['rate_pct']:.1f}%"),
     ("S1-obs", cm["S1-obs"]["rate_pct"] / 100, "c1", f"{cm['S1-obs']['rate_pct']:.1f}%"),
     ("matched expectation", cm["MD-global"]["matched_expected_rate_pct"] / 100, "c0",
      f"{cm['MD-global']['matched_expected_rate_pct']:.1f}%"),
     ("all screened pairs", cm["MD-global"]["background_rate_pct"] / 100, "c0",
      f"{cm['MD-global']['background_rate_pct']:.1f}%")],
    0.25, None, "", "fraction of tested paralog candidates that are hits in combinatorial CRISPR screens",
    "Combinatorial screen support", fmt=lambda v: f"{100 * v:.0f}%")

# ---------------------------------------------------------------- text values
V = dict(
    n_dev=1022, n_r=255,
    n_rows=f"{bench['n_rows_arms']:,}", n_genes=f"{bench['n_genes']:,}",
    s1o=f"{bench['selected']['S1-obs']:,}", eno=f"{bench['selected']['EN-obs']:,}",
    enp=f"{bench['selected']['EN-pred']:,}", s1p=f"{bench['selected']['S1-pred']:,}",
    mdg=f"{bench['n_md_global']:,}", mds=f"{bench['n_md_strict']:,}",
    gate_k=gate2["overall"]["k"], gate_n=gate2["overall"]["n"], gate_r=pct(gate2["overall"]["recall"]),
    gate1_k=gate1["overall"]["k"],
    gate_el_k=gate2["by_universe"]["eligible"]["k"], gate_el_n=gate2["by_universe"]["eligible"]["n"],
    rs_min=f2(min(r[0] for r in ridge)), rs_max=f2(max(r[0] for r in ridge)),
    rg_max=f"{max(ridge_g):.3f}", rec_min=f2(min(recall_g)), rec_max=f2(max(recall_g)),
    suite_s=suite["summary"],
    G_s1o=f2(G["S1-obs"]), G_s1p=f2(G["S1-pred"]), G_eno=f2(G["EN-obs"]), G_enp=f2(G["EN-pred"]),
    c_pred=f2(con["predictor_adds_value"]["estimate"]), c_pred_lo=f2(con["predictor_adds_value"]["ci_lo"]),
    c_pred_hi=f2(con["predictor_adds_value"]["ci_hi"]),
    c_mv=f2(con["multivariate_adds_value"]["estimate"]), c_mv_lo=f2(con["multivariate_adds_value"]["ci_lo"]),
    c_mv_hi=f2(con["multivariate_adds_value"]["ci_hi"]),
    basic_lo=f2(diag["EN-pred - EN-obs"]["basic_ci95"][0]), basic_hi=f2(diag["EN-pred - EN-obs"]["basic_ci95"][1]),
    k200_enp=f2(sk["k200"]["EN-pred"]["rate"]), k200_eno=f"{sk['k200']['EN-obs']['rate']:.3f}",
    k500_enp=f"{sk['k500']['EN-pred']['rate']:.3f}", k500_eno=f2(sk["k500"]["EN-obs"]["rate"]),
    pool=pool["pool_size"],
    x_s1=f2(float(xlib["G"]["rate_S1-obs"])), x_en=f2(float(xlib["G"]["rate_EN-obs"])),
    f_mdg=f"{fdr['MD-global']:.4f}", f_s1o=f"{fdr['S1-obs']:.3f}", f_eno=f2(fdr["EN-obs selected"]),
    f_enp=f2(fdr["EN-pred selected"]),
    nc2a=f"{nc['NC2a']['value']:.3f}", nc2a_b=nc["NC2a"]["bound"], nc3=f2(nc["NC3"]["value"]),
    nc3_b=nc["NC3"]["bound"], nc6=f"{nc['NC6']['value']:.1f}",
    ce_t=ceil["n_ceiling_testable"], ce_s=ceil["labels"]["data-supported"], ce_n=ceil["labels"]["data-null"],
    ce_all=ceil["n_tests"],
    rds_s1=f"{rds['S1-obs']['hits']}/{rds['S1-obs']['n']}", rds_enp=f"{rds['EN-pred']['hits']}/{rds['EN-pred']['n']}",
    head_n=head["n_testable"],
    np_n=allnp["n_testable"], np_s1=allnp["per_arm"]["S1-obs"]["hits"], np_enp=allnp["per_arm"]["EN-pred"]["hits"],
    mc_p=f"{allnp['mcnemar']['p']:.3f}",
    ref_n=ref["n_refuted"],
    cm_g=f"{cm['MD-global']['rate_pct']:.1f}", cm_s=f"{cm['S1-obs']['rate_pct']:.1f}",
    cm_e=f"{cm['MD-global']['matched_expected_rate_pct']:.1f}", cm_r=f"{cm['MD-global']['rate_over_matched_expected']:.1f}",
    n_short=rc["shortlist"]["n_rows"], prog=fam["PROG"], exppar=fam["EXPPAR"], hot=fam["HOT"], lofbi=fam["LOFBI"],
    c9=pcg["9"], c5=pcg["5"], c6=pcg["6"], c7=pcg["7"],
    pre9=f"{rc['pre_power_fix']['pair_classes_md_global']['9']:,}", pre_short=rc["pre_power_fix"]["n_shortlist"],
    lin_x=lin["by_source"]["MD-global"]["labels"]["cross-lineage"],
    lin_r=lin["by_source"]["MD-global"]["labels"]["lineage-restricted"],
    lin_u=lin["by_source"]["MD-global"]["labels"]["lineage-unresolvable"],
    lin_x9=pct(lbc["class_9_share_pct"]["cross-lineage"] / 100), lin_u9=pct(lbc["class_9_share_pct"]["lineage-unresolvable"] / 100),
    mods=sub["MD-global"]["n_modules"], mods_multi=sub["MD-global"]["n_multi_gene_modules"],
    sa_rab=f2(setA.get("RAB6A", float("nan"))), sa_mdm2=f2(setA.get("MDM2", float("nan"))),
    sa_fam=f2(setA.get("FAM50A", float("nan"))),
    suspects=lim["cn_cutting_suspect"]["n_suspect_genes"], susp_short=lim["cn_cutting_suspect"]["shortlist_rule_pairs_on_suspect_genes"],
    gc3=f"{lim['gene_classes']['3']:,}", gc3b=f"{lim['gene_classes']['3b']:,}",
    locus_pct=f"{rc['locus']['md_global_locus_level_pct']:.0f}%", hours=wall_h,
    fam_lin=next(r["n_powered_lineages"] for r in short if r["symbol"] == "FAM50A"),
    n_dev_rows=len(devs), n_commits=n_commits,
)

# ---------------------------------------------------------------- page
CSS = r"""
/* Layout: one reading column (~68ch) with full-width figures; numbered steps follow the real order of the work. */
:root{
  --bg:#f7f8f9; --paper:#ffffff; --ink:#1b2430; --muted:#5b6673; --rule:#d9dee4; --accent:#1f6f78; --warn:#a4561a;
  --c1:#1f6f78; --c2:#7fb4ba; --c3:#6a4c93; --c4:#b39cd0; --c5:#b5651d; --c6:#8a8f97; --c0:#c4c9cf;
  --font-display:"Newsreader", "Iowan Old Style", Georgia, serif;
  --font-body:"Source Serif 4", "Iowan Old Style", Georgia, serif;
  --font-ui:"IBM Plex Sans", "Segoe UI", system-ui, sans-serif;
}
@media (prefers-color-scheme: dark){ :root:not([data-theme="light"]){
  --bg:#12171d; --paper:#182029; --ink:#e3e8ee; --muted:#9aa6b2; --rule:#2c3742; --accent:#6cc3cc; --warn:#e3a066;
  --c1:#5fb8c1; --c2:#2f6d73; --c3:#b49be0; --c4:#6e5a91; --c5:#e3a066; --c6:#7d8894; --c0:#3d4955; color-scheme:dark } }
:root[data-theme="dark"]{
  --bg:#12171d; --paper:#182029; --ink:#e3e8ee; --muted:#9aa6b2; --rule:#2c3742; --accent:#6cc3cc; --warn:#e3a066;
  --c1:#5fb8c1; --c2:#2f6d73; --c3:#b49be0; --c4:#6e5a91; --c5:#e3a066; --c6:#7d8894; --c0:#3d4955; color-scheme:dark }
body{background:var(--bg);color:var(--ink);font-family:var(--font-body);font-size:17px;line-height:1.6;
  padding-left:max(16px, env(safe-area-inset-left, 0px));padding-right:max(16px, env(safe-area-inset-right, 0px))}
main{max-width:46rem;margin:0 auto;padding-block:2.5rem 4rem}
h1,h2,h3{font-family:var(--font-display);font-weight:600;line-height:1.2;text-wrap:balance}
h1{font-size:2.3rem;margin:.2rem 0 .6rem}
h2{font-size:1.55rem;margin:2.6rem 0 .6rem;padding-top:1.2rem;border-top:1px solid var(--rule)}
h3{font-size:1.15rem;margin:1.6rem 0 .3rem}
p{margin:.7rem 0}
.kicker{font-family:var(--font-ui);font-size:.78rem;letter-spacing:.08em;text-transform:uppercase;color:var(--accent);font-weight:600}
.byline{font-family:var(--font-ui);font-size:.88rem;color:var(--muted)}
.lede{font-size:1.15rem}
.note{font-family:var(--font-ui);font-size:.86rem;color:var(--muted);border-left:3px solid var(--rule);padding:.2rem 0 .2rem .9rem;margin:1rem 0}
.disclaimer{background:var(--paper);border:1px solid var(--rule);border-radius:6px;padding:.9rem 1.1rem;font-family:var(--font-ui);font-size:.88rem;color:var(--muted)}
.take{background:var(--paper);border:1px solid var(--rule);border-radius:6px;padding:1rem 1.3rem;margin:1.4rem 0}
.take h3{margin-top:.2rem;font-family:var(--font-ui);font-size:.82rem;letter-spacing:.07em;text-transform:uppercase;color:var(--accent)}
.take ul{margin:.3rem 0;padding-left:1.2rem}
.take li{margin:.35rem 0}
dl.gloss{display:grid;grid-template-columns:max-content 1fr;gap:.35rem 1rem;font-size:.95rem}
dl.gloss dt{font-family:var(--font-ui);font-weight:600;font-size:.88rem}
dl.gloss dd{margin:0;min-width:0}
figure{margin:1.6rem 0;background:var(--paper);border:1px solid var(--rule);border-radius:6px;padding:1rem 1rem .8rem}
figcaption{font-family:var(--font-ui);font-size:.84rem;color:var(--muted);line-height:1.5;margin-top:.5rem}
figcaption b{color:var(--ink)}
.fig-svg{width:100%;height:auto;display:block}
.fig-svg text{font-family:var(--font-ui)}
.lab{fill:var(--ink);font-size:13px}.tick{fill:var(--muted);font-size:11.5px}.val{fill:var(--ink);font-size:12px;font-variant-numeric:tabular-nums}
.axis{fill:var(--muted);font-size:12px}.limitlab{fill:var(--warn);font-size:12px;font-weight:600}
.grid{stroke:var(--rule);stroke-width:1}.limit{stroke:var(--warn);stroke-width:1.6;stroke-dasharray:5 4}
.box{fill:var(--bg);stroke:var(--rule);stroke-width:1.2}.b1{stroke:var(--c1)}.b3{stroke:var(--c3)}.b4{stroke:var(--warn)}.b5{stroke:var(--c5)}
.boxt{fill:var(--ink);font-size:13.5px;font-weight:600}.boxs{fill:var(--muted);font-size:11.5px}
.arrow{stroke:var(--muted);stroke-width:1.3}.ahead{fill:var(--muted)}
.line{fill:none;stroke-width:2.4}
.c0{fill:var(--c0)}.c1{fill:var(--c1)}.c2{fill:var(--c2)}.c3{fill:var(--c3)}.c4{fill:var(--c4)}.c5{fill:var(--c5)}.c6{fill:var(--c6)}
.line.c1{stroke:var(--c1);fill:none}.line.c2{stroke:var(--c2);fill:none}.line.c3{stroke:var(--c3);fill:none}.line.c4{stroke:var(--c4);fill:none}
.c1t{fill:var(--c1)}.c2t{fill:var(--c2)}.c3t{fill:var(--c3)}.c4t{fill:var(--c4)}
.dot{stroke:var(--paper);stroke-width:1.5}
table{border-collapse:collapse;width:100%;font-family:var(--font-ui);font-size:.88rem;font-variant-numeric:tabular-nums}
th,td{border-bottom:1px solid var(--rule);padding:.4rem .5rem;text-align:left;vertical-align:top}
th{color:var(--muted);font-weight:600}
.tw{overflow-x:auto;margin:1rem 0}
ol.refs{font-size:.9rem;padding-left:1.3rem}ol.refs li{margin:.25rem 0}
a{color:var(--accent)}
code{font-family:"IBM Plex Mono", ui-monospace, monospace;font-size:.86em}
@media (max-width:520px){ body{font-size:16px} h1{font-size:1.8rem} dl.gloss{grid-template-columns:1fr} }
"""

HTML = f"""<title>Reading Biomarkers From a Dependency Predictor</title>
<meta name="claude:viewport-fit-sides" content="cover">
<link rel="preconnect" href="https://fonts.googleapis.com"><link rel="preconnect" href="https://fonts.gstatic.com" crossorigin>
<link rel="stylesheet" href="https://fonts.googleapis.com/css2?family=IBM+Plex+Sans:wght@400;600&family=Newsreader:opsz,wght@6..72,500;6..72,600&family=Source+Serif+4:ital,opsz,wght@0,8..60,400;0,8..60,600;1,8..60,400&display=swap">
<style>{CSS}</style>
<main>
<p class="kicker">Research briefing · agentic bioinformatics</p>
<h1>Can an accurate essentiality predictor tell us why cancer cells depend on a gene?</h1>
<p class="byline">Written by the AI principal investigator (Claude Opus 5.5) of an AI-run project · 2026-10-30 · source repository: QuantitativeBiology/ai-hackathon-cancer_biomarkers</p>

<div class="disclaimer">This briefing, the analysis and all results were produced by AI agents with minimal human intervention, as part of a hackathon on agentic bioinformatics (Beltrao Lab, ETH Zurich, with INESC-ID and IST, Universidade de Lisboa). No human has yet verified the results. Every number below is read by a script from the project's committed result files.</div>

<p class="lede">A previous project built a model that predicts how much each cancer cell line depends on each gene. If a model predicts dependencies accurately, one might expect it to have learnt <em>why</em> a line depends on a gene, for example because the line carries a particular mutation. This project tested whether those reasons can be read out of the model as biomarkers that group cell lines into dependency subtypes, and whether such biomarkers hold up against simpler methods, held-out data, clinical knowledge and independent experiments.</p>

<div class="take">
<h3>Main findings</h3>
<ul>
<li><b>The predictor did not add detectable value.</b> In 255 cell lines held back until the end, the top 100 candidate biomarkers read from the predictor replicated at {V['G_enp']}, against {V['G_eno']} for the same method applied to measured data and {V['G_s1o']} for a simple one-feature-at-a-time test.</li>
<li><b>That verdict is fragile.</b> The pre-declared comparison shows no detectable difference. An alternative confidence interval says the predictor did worse, and with longer candidate lists it did slightly better.</li>
<li><b>The simplest method was the most reliable.</b> The confounder-adjusted single-feature test replicated best, was best calibrated and recovered the most clinical biomarkers.</li>
<li><b>The final catalogue is mostly a paralog story.</b> The {V['n_short']} best-supported pairs are dominated by dependencies on one member of a gene pair when its partner is weakly expressed. Several are also dependencies of non-cancerous cells.</li>
<li><b>Agents can run a pre-registered analysis end to end, but the independent checks mattered.</b> Several consequential errors were caught only by agents that had not produced the result.</li>
</ul>
</div>

<h2>1. The question and its vocabulary</h2>
<p>In a CRISPR knockout screen each gene is disrupted in turn, and the resulting loss of growth is scored. A strongly negative score means the line <em>depends</em> on that gene. The DepMap project has screened over a thousand cancer cell lines this way; we used its 26Q3 release.</p>
<p>The Part 1 predictor estimates these scores from a line's molecular profile. It averages three models: a kernel ridge regression on gene expression, a model built on paralog genes, and a gradient-boosted tree model on the residual. Its predictions are "out-of-fold": each line's prediction comes from a model that never saw that line.</p>
<p>The premise in the project brief was strong: <em>every</em> accurate prediction implies a molecular determinant that could serve as a biomarker. We treated that as a hypothesis to test, not a fact.</p>
<dl class="gloss">
<dt>Dependency</dt><dd>A gene whose knockout slows a cell line's growth.</dd>
<dt>Determinant / biomarker</dt><dd>A measurable feature of the line (a mutation, a copy-number change, methylation, expression of a gene or a gene programme) that marks the lines that depend on a gene.</dd>
<dt>Subtype</dt><dd>The group of lines sharing a dependency for the same molecular reason.</dd>
<dt>Lineage</dt><dd>The tissue of origin of a cell line (lung, bowel, blood, …). Lineage is the largest confounder: a feature and a dependency can both be common in one tissue without being linked.</dd>
<dt>Paralog</dt><dd>A related gene from an ancient duplication. Paralogs can back each other up, so losing one makes a cell depend on the other.</dd>
<dt>Replication</dt><dd>A finding from the 1,022 development lines that reappears, with the same sign and one-sided p &lt; 0.05, in 255 lines that were set aside.</dd>
</dl>

<h2>2. How the project was run</h2>
<p>The project was executed by an AI "principal investigator" that delegated work to other AI agents. Separate agents wrote the analysis code, built a simulated test world, reviewed the design, reimplemented the statistics independently and wrote this text. One human set the task, supplied data that could not be downloaded by script, and made one decision: to extend the 24-hour budget.</p>
<p>Three practices shaped the work:</p>
<ul>
<li><b>Pre-registration in stages.</b> The design was frozen after three rounds of review. Each later stage was then fixed with a git tag before the data it depended on were analysed: the clinical test table, the control bounds, and finally the catalogue, which was frozen before the held-out lines were read once. Every later change was logged as a numbered deviation; there are {V['n_dev_rows']} in total.</li>
<li><b>Pre-committed responses to failure.</b> Before any real result existed, a table fixed what would happen if each control failed. When controls did fail, there was no room to rescue them after the fact.</li>
<li><b>Authors never verify their own work.</b> The simulated test world and its answer key, the reimplementation of the statistics, and two critical reviews were all done by agents that had not written the code under test.</li>
</ul>
<p>The project ran for about {V['hours']} hours of wall time across {V['n_commits']} commits, including an overnight interruption and a laptop with 15.7 GB of memory. Several decisions below were shaped by those limits.</p>

<figure>{fig_design()}
<figcaption><b>Figure 1. Study design.</b> The same two methods (a univariate test and an elastic net) were each run on observed and on predicted dependency scores, giving four "arms". The catalogue keeps pairs found by the predictor-based elastic net that the observed univariate test also confirms. Everything was checked on simulated data with a known answer, against controls, against clinical and combinatorial-screen data, and finally in held-out cell lines read once.</figcaption></figure>

<h2>3. Step one: build the candidate biomarkers and protect a test set</h2>
<p><b>Held-out lines.</b> Part 1 had already split DepMap into 1,022 development lines and 255 test lines. We kept that split. The 255 lines were locked behind a software guard that refuses to load their scores unless the catalogue tag exists and an explicit unlock is set. Every read is logged, and the log shows one read.</p>
<p><b>Feature dictionary.</b> Each candidate biomarker is one feature: a hotspot mutation, a loss-of-function call, an amplification or deletion, a fusion, promoter methylation, a paralog's expression, an expression programme from matrix factorisation, or a clinical composite such as "TP53 mutant". Thresholds were set on the development lines only. Testing every eligible gene against its features gave {V['n_rows']} gene–feature pairs across {V['n_genes']} genes.</p>
<p><b>Eligible genes.</b> We kept genes whose dependency varies between lines and is predicted with useful accuracy beyond what lineage and screening library explain. Genes that fail are not dropped from the story. They are counted later as one kind of limit.</p>

<h2>4. Step two: four ways to find biomarkers</h2>
<p>The central comparison is a two-by-two design. One axis is the method; the other is whether it reads measured or predicted dependency scores.</p>
<ul>
<li><b>Univariate test (S1).</b> One feature at a time, a regression of the dependency score on the feature, adjusted for lineage, screening library, culture conditions and growth rate, and where appropriate for mutation burden, microsatellite instability, TP53 status and the gene's own copy number. This follows the per-gene association design used in large cell-line screens (Garnett et al. 2012; Iorio et al. 2016; Behan et al. 2019; Pacini et al. 2024), with the covariate structure modelled on Gonçalves et al. 2020. Standard errors are clustered by patient, because some lines come from the same person. They use a small-sample correction ("CR2") with an adjusted number of degrees of freedom (Bell and McCaffrey 2002; Pustejovsky and Tipton 2018). Multiple testing across the millions of tests is controlled with the Benjamini–Hochberg procedure (1995).</li>
<li><b>Elastic net (EN).</b> All features of a gene are fitted jointly, in two stages: genomic features first, then expression features. It uses a penalised regression that keeps a few features and shrinks the rest to zero (Zou and Hastie 2005). A feature counts as selected only if it is chosen in at least 60% of 50 half-sample refits, with consistent sign. This is complementary-pairs stability selection (Meinshausen and Bühlmann 2010; Shah and Samworth 2013).</li>
</ul>
<p>Run on measured scores these give <b>S1-obs</b> and <b>EN-obs</b>; run on the predictor's scores, <b>S1-pred</b> and <b>EN-pred</b>. EN-pred is the purest reading of "what the predictor learnt": it never sees a measured dependency score.</p>
<p>On the development lines the arms selected very different numbers of pairs: S1-obs {V['s1o']}, EN-obs {V['eno']}, EN-pred {V['enp']} and S1-pred {V['s1p']}. Predicted scores are smoother than measured ones, so methods reading them find many more "associations". Section 5 explains why many of those are not trustworthy.</p>

<h2>5. Step three: test the machinery on a world with a known answer</h2>
<p>Before touching real associations, a separate agent built a simulated data set with planted biomarkers of known kinds:</p>
<ul>
<li>pan-cancer and lineage-restricted effects;</li>
<li>passenger mutations near a true driver;</li>
<li>culture artefacts and artefacts of CRISPR cutting;</li>
<li>a fold-specific bias of the kind that out-of-fold predictions can carry.</li>
</ul>
<p>Our pipeline then had to recover the plants and avoid false hits. On the final engine, {V['suite_s']} of these checks passed. The three failures are reported in the full report, not hidden.</p>
<p>The simulation found the most important problem in the design.</p>
<ul>
<li><b>The planned catalogue failed its error bound.</b> The frozen plan, "MD-strict", took the predictor-based elastic-net picks and re-checked them against measured scores within that selection. When the predictor was ridge-like, as Part 1's is, the true false-discovery proportion was {V['rs_min']}–{V['rs_max']} against a limit of 0.15 (Figure 2).</li>
<li><b>The cause is a double dip.</b> A ridge predictor trained on the development lines absorbs their noise. Re-checking its picks against those same lines then confirms the noise.</li>
<li><b>The fix was decided on the simulation alone.</b> We switched the catalogue to "MD-global", which keeps a pick only if it is also significant in the observed univariate test corrected across all tests. Its false-discovery proportion was at most {V['rg_max']}, with the same recall ({V['rec_min']}–{V['rec_max']}). The switch was made before any real association result existed; the only real-data output seen by then was the count of pairs selected in a 50-gene timing pilot.</li>
</ul>
<p>The simulation also showed something this design could not solve. On a ridge-like predictor, lists built from predicted scores have a large hidden false-discovery rate. A permutation null cannot see it, because the noise lives inside the predictions.</p>

<figure>{FIG_TOY}
<figcaption><b>Figure 2. Why the catalogue definition was changed.</b> True false-discovery proportion on five simulated data sets built with a ridge-type predictor, for four lists. The pre-registered catalogue (MD-strict, grey) exceeds the 0.15 limit every time. The adopted catalogue (MD-global, amber) stays close to zero. The predictor-based elastic-net list (EN-pred) is around 0.5, which is why prediction-based lists are labelled "not calibrated" throughout. A caveat: MD-global was chosen on these same simulations, so its good performance here is not an independent test.</figcaption></figure>

<h2>6. Step four: controls on the real data</h2>
<h3>Known biomarkers: the positive-control gate failed</h3>
<p>Before any biological claim, the pipeline had to recover a pre-specified list of known dependency–biomarker pairs. The list covered:</p>
<ul>
<li>KRAS-mutant lines depending on KRAS;</li>
<li>loss of one paralog making a line depend on the other, as with SMARCA4/SMARCA2 and MTAP loss with PRMT5;</li>
<li>MDM2 dependence in TP53 wild-type lines.</li>
</ul>
<p>The pre-declared bar was 60% recovery. The pipeline recovered {V['gate_k']} of {V['gate_n']} ({V['gate_r']}), so the gate failed.</p>
<p>A bug search, the only remedy the pre-committed rules allowed, found one real bug. Control genes outside the eligible set had never been analysed. After the fix, recovery was still {V['gate_k']}/{V['gate_n']}: the formerly missing controls were genuine misses. Restricted to eligible genes, recovery was {V['gate_el_k']}/{V['gate_el_n']}, but that is not the pre-declared denominator and did not change the verdict.</p>
<p>The misses are informative:</p>
<ul>
<li>EGFR-mutant dependence has too few carriers for the conservative test.</li>
<li>Paralog pairs that need a homozygous deletion are not captured by a "shallow loss" feature.</li>
<li>MYCN amplification is almost confined to neuroblastoma, so the lineage term absorbs it.</li>
</ul>

<figure>{FIG_GATE}
<figcaption><b>Figure 3. Recovery of known dependency–biomarker pairs</b> (gate attempt 2). No list reaches the pre-declared 60% bar. The predicted-score univariate list (S1-pred) scores highest only because it is enormous: it selects {V['s1p']} pairs, so it contains many controls by size alone.</figcaption></figure>

<h3>False discoveries: a permutation check</h3>
<p>To estimate how many hits arise by chance, we shuffled the feature data between cell lines within lineage-and-library groups, 20 times, and re-ran the methods at their fixed thresholds on 300 genes.</p>
<ul>
<li><b>Within target:</b> the catalogue (MD-global, {V['f_mdg']}) and the observed univariate test (S1-obs, {V['f_s1o']}) were well inside the 0.10 target.</li>
<li><b>Failed:</b> the elastic net on measured scores was not ({V['f_eno']}), and is labelled "false discovery not controlled".</li>
<li><b>Not calibrated:</b> the prediction-based lists cannot be judged by this check at all, as Figure 2 explained.</li>
</ul>

<figure>{FIG_NULL}
<figcaption><b>Figure 4. Estimated false-discovery rate from permutations.</b> Only the catalogue and the observed univariate list are both within the 0.10 target and interpretable. The values for the prediction-based lists (S1-pred, EN-pred) and for MD-strict are not reliable estimates, because a ridge predictor's hidden false discoveries are invisible to this check.</figcaption></figure>

<h3>Negative controls</h3>
<p>Three further pre-declared negative controls also failed:</p>
<ul>
<li><b>Low-tier genes:</b> the predictor's fidelity was {V['nc2a']}, against a bound of {V['nc2a_b']} set from the simulation. The bug search found no bug; the simulated bound was too strict for real data.</li>
<li><b>Self-null features:</b> a gene's own deletion should not make lines more dependent on it, yet {pct(nc['NC3']['value'])} of such pairs showed dependent carriers, against a bound of {V['nc3_b']}. This was driven by deletion calls. The {V['suspects']} affected genes are flagged as copy-number artefact suspects.</li>
<li><b>Lineage transcription factors:</b> these should be classified as "lineage-driven", but their fraction in that class was {V['nc6']}. The classification rule is therefore flagged as unreliable.</li>
</ul>
<p>We report all three as failures rather than adjusting their bounds.</p>

<h2>7. Step five: the test that counts, replication in held-out lines</h2>
<p>With the catalogue and the candidate lists frozen, the 255 held-out lines were read once. For each arm we took its top 100 genomic candidates (k = 100, fixed beforehand from a pool of {V['pool']:,} powered candidates). We then asked how many show the same effect in the new lines.</p>
<div class="tw"><table>
<tr><th>Arm</th><th>Replicated (top 100)</th><th>What it reads</th></tr>
<tr><td>S1-obs</td><td>{V['G_s1o']}</td><td>measured scores, one feature at a time</td></tr>
<tr><td>S1-pred</td><td>{V['G_s1p']}</td><td>predicted scores, one feature at a time</td></tr>
<tr><td>EN-obs</td><td>{V['G_eno']}</td><td>measured scores, joint model</td></tr>
<tr><td>EN-pred</td><td>{V['G_enp']}</td><td>predicted scores, joint model</td></tr>
</table></div>
<p>The two pre-declared contrasts:</p>
<ul>
<li><b>Does the predictor add value?</b> EN-pred − EN-obs = {V['c_pred']} (95% CI {V['c_pred_lo']} to {V['c_pred_hi']}): no detectable difference.</li>
<li><b>Does the joint model add value over single features?</b> EN-obs − S1-obs = {V['c_mv']} ({V['c_mv_lo']} to {V['c_mv_hi']}): the single-feature test is better.</li>
</ul>
<p>Replication across screening libraries within the development lines points the same way: S1-obs {V['x_s1']}, EN-obs {V['x_en']}.</p>
<p>Two features of the data make the first verdict fragile.</p>
<ul>
<li><b>The bootstrap is biased.</b> Its replicate rates sit below the observed rates, so the standard percentile interval is off. A "basic" bootstrap interval for the same contrast is {V['basic_lo']} to {V['basic_hi']}, which would say the predictor does worse.</li>
<li><b>The sign flips at deeper lists.</b> With 200 and 500 candidates, EN-pred replicates better than EN-obs ({V['k200_enp']} vs {V['k200_eno']}, and {V['k500_enp']} vs {V['k500_eno']}; Figure 5).</li>
</ul>
<p>The fair summary: the predictor neither clearly helps nor clearly hurts, and its value depends on how deep one reads the list. The single-feature test is never beaten for lists up to 200 candidates.</p>
<p>Absolute rates are high, partly for reasons unrelated to biology:</p>
<ul>
<li>The lists favour the strongest effects.</li>
<li>The held-out lines share the same processing pipeline as the development lines, so replication tests generalisation across cell lines, not across laboratories.</li>
<li>Part 1 had already scored its predictor on these lines once. They were unused by this project, but not untouched.</li>
</ul>

<figure>{fig_k_curve()}
<figcaption><b>Figure 5. Replication in 255 held-out lines, by list length.</b> The pre-declared comparison is at the top 100 (dashed line); 50, 200 and 500 are pre-declared secondary list lengths computed from the same single read. The observed univariate test (S1-obs) leads at short lists. The two elastic-net arms swap order between 100 and 200 candidates, which is why the "no added value" verdict is called fragile.</figcaption></figure>

<h2>8. Step six: clinical biomarkers and combinatorial screens</h2>
<h3>Clinical ground truth</h3>
<p>Drugs that reached trials with a companion biomarker are the strongest available ground truth: BRAF V600 for BRAF/MEK inhibitors, MTAP loss for PRMT5 inhibitors, and so on. Before joining any dependency data we built a table of {V['ce_all']} biomarker–target tests from CIViC, the Molecular Oncology Almanac (Reardon et al. 2021) and the FDA list of companion diagnostics. For each test we recorded the expected direction of the knockout effect.</p>
<p>Three results:</p>
<ul>
<li><b>Most clinical biomarkers are visible in cell lines.</b> Of {V['ce_t']} tests that can be checked in the data, {V['ce_s']} show the expected dependency directly ("data-supported") and {V['ce_n']} do not. EGFR-activating mutations in lung cancer and ALK fusions are among the misses, partly because of very few carriers.</li>
<li><b>The predictor-derived lists recover almost none of them.</b> Among the detectable, non-proxy tests, S1-obs recovers {V['rds_s1']} in its top 200 and EN-pred {V['rds_enp']} (Figure 6). Across the full non-proxy set S1-obs recovers {V['np_s1']} of {V['np_n']} and EN-pred {V['np_enp']}; the exact McNemar test for EN-pred vs S1-obs gives p = {V['mc_p']}.</li>
<li><b>The primary clinical comparison is underpowered.</b> It was restricted to tests that do not overlap the positive-control gate, which left only {V['head_n']} tests, below the pre-declared minimum of 15. It is reported as descriptive.</li>
</ul>
<p>No method recovered any of the {V['ref_n']} biomarkers that trials have refuted, which is the right answer. With 3 tests that says little.</p>

<figure>{FIG_CLIN}
<figcaption><b>Figure 6. Clinical biomarkers recovered.</b> Recall in each arm's top 200 of the clinical biomarker–target pairs that are detectable in the cell-line data (non-proxy tests). The numbers are small, so the intervals are wide, but no prediction-based arm recovers any.</figcaption></figure>

<h3>Combinatorial CRISPR screens</h3>
<p>Many catalogue pairs are paralog pairs. Those can be tested independently, because several studies knocked out both genes of a pair together (Dede et al. 2020; Thompson et al. 2021; Parrish et al. 2021; Esmaeili Anvar et al. 2024).</p>
<ul>
<li>Of the catalogue's tested paralog candidates, {V['cm_g']}% were hits in those screens.</li>
<li>A matched set of paralog pairs, with similar sequence identity and family size, would give {V['cm_e']}%: about {V['cm_r']} times enrichment.</li>
<li>The observed univariate list did equally well ({V['cm_s']}%).</li>
</ul>
<p>The candidates are real paralog dependencies, but the predictor does not add to them. The two lists are nested, and the screens are not fully independent of DepMap, so this is supporting evidence rather than proof.</p>

<figure>{FIG_COMB}
<figcaption><b>Figure 7. Support from combinatorial CRISPR screens.</b> Fraction of tested paralog candidates that are hits in published double-knockout screens, compared with paralog pairs matched for sequence identity and family size, and with all screened pairs. The per-cell-line hit rule (effect size and minimum effect) was an implementation choice made before the catalogue existed.</figcaption></figure>

<h2>9. Step seven: what the subtypes look like</h2>
<p><b>How the catalogue sorts out.</b> The catalogue contains {V['mdg']} gene–biomarker pairs, grouped into {V['mods']} modules, {V['mods_multi']} of them shared by more than one gene. After replication and a power check, {V['c9']} pairs are classed as "valid candidates": replicated in held-out lines with at least 50% power. The rest:</p>
<ul>
<li>{V['c5']} are too rare or under-powered to judge;</li>
<li>{V['c7']} were powered but did not replicate;</li>
<li>{V['c6']} are confounded.</li>
</ul>
<p>An earlier count of {V['pre9']} valid pairs was wrong. Power had not been computed for continuous features, so pairs that had not replicated stayed classed as valid. An independent reviewer caught this.</p>
<p><b>Within or across cancer types.</b> Of the catalogue pairs, {V['lin_x']} behave the same way across several lineages ("cross-lineage"), {V['lin_r']} are confined to one, and {V['lin_u']} cannot be resolved because too few lineages have enough lines. Cross-lineage pairs end up valid far more often ({V['lin_x9']} vs {V['lin_u9']} for unresolvable ones). But the cross-lineage label itself requires several well-powered lineages, so part of that gap is simply statistical power.</p>
<p><b>What kind of biology.</b> The {V['n_short']} valid pairs are dominated by expression-level determinants: {V['prog']} expression programmes and {V['exppar']} paralog-expression features. There are only {V['hot']} hotspot-mutation and {V['lofbi']} loss-of-function features (Figure 8). Typical examples:</p>
<ul>
<li><b>FAM50A</b> dependence where its paralog FAM50B is weakly expressed, seen in {V['fam_lin']} lineages;</li>
<li><b>DDX3X</b> dependence where DDX3Y is absent, which tracks sex and Y-chromosome loss, and whose gene is also flagged as a copy-number artefact suspect;</li>
<li><b>MDM2</b> dependence in TP53 wild-type lines;</li>
<li><b>KRAS</b> dependence in KRAS-mutant lines, the cleanest genomic case;</li>
<li><b>PRMT5</b> dependence with MTAP loss, the one example with a matching clinical programme. It entered the catalogue through the univariate test, not through anything the predictor revealed.</li>
</ul>
<p>Several top hits are also dependencies in non-cancerous cell lines: RAB6A in {pct(float(V['sa_rab']))} of them, MDM2 {pct(float(V['sa_mdm2']))}, FAM50A {pct(float(V['sa_fam']))}. This warns that they may not be tumour-selective. It is not a therapeutic-window measurement.</p>

<figure>{fig_families()}
<figcaption><b>Figure 8. What the {V['n_short']} best-supported pairs are made of.</b> Determinant types of the shortlisted pairs (classed valid and replicated in held-out lines). Amber marks expression-level determinants; violet marks genomic and composite determinants. Most of the signal is in expression and paralog buffering, not in the kinds of mutations a diagnostic test usually reads.</figcaption></figure>

<h2>10. Which dependencies cannot serve as biomarkers, and why</h2>
<p>Every gene and pair was assigned a "failure class". The large classes:</p>
<ul>
<li><b>No usable contrast or not predictable.</b> Many dependencies barely vary between lines, or the predictor does not capture them, so there is nothing to explain.</li>
<li><b>Explained only by tissue of origin.</b> {V['gc3']} genes are predicted accurately only through lineage or screening library. The predictor knows they are blood- or brain-specific, but not why.</li>
<li><b>Predictable but with no nameable determinant.</b> {V['gc3b']} genes are predicted well, yet no catalogue feature explains them. This alone contradicts the premise that every accurate prediction implies a readable biomarker.</li>
<li><b>Determinants that are not genes.</b> {V['locus_pct']} of copy-number determinants resolve only to a chromosomal region containing several genes, not to one gene. Expression programmes are not assayable biomarkers in the clinical sense.</li>
<li><b>Too rare to judge.</b> Biomarkers carried by a handful of lines, such as EGFR-activating mutations, cannot be confirmed with the conservative statistics needed to avoid false positives in small groups.</li>
</ul>

<h2>11. What this says about agentic science</h2>
<p><b>What worked.</b> The agents wrote a reviewed design, froze it, built a simulated world with a known answer, and ran a multi-day analysis to a single, logged read of held-out data. Failures were reported as failures: the gate, four negative controls, a calibration test and one tolerance of the independent reimplementation. Pre-committing the response to each failure removed the temptation, and the opportunity, to rescue results after seeing them.</p>
<p><b>What the independent checks caught.</b> Several consequential errors were found only by agents that had not produced the result:</p>
<ul>
<li>the catalogue's double-dip problem, found by the simulation agent;</li>
<li>an anti-conservative statistical test, which triggered a pre-declared fallback, found by the process-check agent;</li>
<li>missing statistical power for continuous features, which had inflated the valid set by about a third, found by the post-replication reviewer;</li>
<li>a miscount in the permutation summary;</li>
<li>a claim by the PI agent that seven controls were "missed by every method", which was false.</li>
</ul>
<p>Without these second readers, the report would have overstated its results.</p>
<p><b>What went wrong.</b></p>
<ul>
<li>Memory limits on a laptop forced three rewrites of the orchestration code.</li>
<li>Long jobs without checkpoints were lost to an overnight interruption.</li>
<li>One shortcut was logged in error and withdrawn: skipping the dry run before the held-out read.</li>
<li>The budget was extended by the human.</li>
</ul>
<p>These lessons are recorded as process rules in the repository.</p>

<div class="take">
<h3>Take-away messages</h3>
<ul>
<li>Reading biomarkers out of an accurate dependency predictor did not beat simply testing each feature against measured dependencies, with careful confounder control. On this evidence, prediction accuracy did not translate into better biomarkers.</li>
<li>Prediction-based association lists can carry false discoveries that standard permutation checks cannot see. Any such list needs a fresh-data test.</li>
<li>The robust signal in this catalogue is paralog buffering and expression state, supported by combinatorial screens. Most of it is not the mutation-type biomarker that drives clinical practice, and some of it is shared with normal cells.</li>
<li>The clinical comparison was too small to be decisive. That limitation is itself a finding about how few clinical biomarkers can be tested in cell-line data.</li>
<li>For agentic science, staged pre-registration and independent verification were what made the result trustworthy enough to report. They should be standard.</li>
</ul>
</div>

<p class="note">Not done or cut under the pre-declared rules: patient-level prevalence from TCGA, an attribution panel for the predictor's components, a selection-adaptive permutation re-run, and some secondary replication analyses. One held-out subset analysis would have required a second read of the held-out data and was not done. All departures are listed in the repository's deviations log.</p>

<h2>Credits: ideas taken from the literature</h2>
<p>The final method that held up, a confounder-adjusted per-feature test intersected with stability-selected joint models, assembles established ideas rather than inventing new ones:</p>
<ul>
<li><b>Per-gene association testing with context adjustment</b> follows GDSC and DepMap biomarker work.</li>
<li><b>Variable selection with error control</b> comes from stability selection.</li>
<li><b>Small-sample robust inference</b> comes from the econometrics literature on cluster-robust errors.</li>
<li><b>The confounder list</b> comes from studies of CRISPR screen artefacts: copy-number cutting toxicity, culture conditions, passenger and co-deleted genes.</li>
<li><b>The validation sources</b> are curated clinical knowledge bases and published combinatorial screens.</li>
</ul>
<ol class="refs">
<li>Tsherniak A. et al. 2017, Cell: Defining a cancer dependency map. · Dempster J.M. et al. 2021, Genome Biology: Chronos.</li>
<li>Garnett M.J. et al. 2012, Nature; Iorio F. et al. 2016, Cell: pharmacogenomic association testing with tissue adjustment.</li>
<li>Behan F.M. et al. 2019, Nature; Pacini C. et al. 2024, Cancer Cell: CRISPR-based target and biomarker prioritisation.</li>
<li>Gonçalves E. et al. 2020, Molecular Systems Biology: covariate-adjusted drug–gene association models.</li>
<li>Zou H. and Hastie T. 2005, J R Stat Soc B: the elastic net.</li>
<li>Meinshausen N. and Bühlmann P. 2010, J R Stat Soc B; Shah R.D. and Samworth R.J. 2013, J R Stat Soc B: stability selection.</li>
<li>Bell R.M. and McCaffrey D.F. 2002, Survey Methodology; Pustejovsky J.E. and Tipton E. 2018, J Bus Econ Stat: small-sample cluster-robust inference.</li>
<li>Benjamini Y. and Hochberg Y. 1995, J R Stat Soc B: false discovery rate.</li>
<li>Aguirre A.J. et al. 2016 and Munoz D.M. et al. 2016, Cancer Discovery: copy-number cutting toxicity in CRISPR screens.</li>
<li>Marjon K. et al. 2016, Cell Reports; Mavrakis K.J. et al. 2016, Science: MTAP loss and PRMT5/MAT2A dependence.</li>
<li>Chan E.M. et al. 2019, Nature: WRN dependence in microsatellite-unstable cancers.</li>
<li>De Kegel B. and Ryan C.J. 2019, PLoS Genetics; De Kegel B. et al. 2021, Cell Systems: paralog buffering of dependencies.</li>
<li>Dede M. et al. 2020, Genome Biology; Thompson N.A. et al. 2021, Nature Communications; Parrish P.C.R. et al. 2021, Cell Reports; Esmaeili Anvar N. et al. 2024, Nature Communications: combinatorial paralog screens.</li>
<li>Griffith M. et al. 2017, Nature Genetics (CIViC); Reardon B. et al. 2021, Nature Cancer (Molecular Oncology Almanac); US FDA list of companion diagnostics.</li>
</ol>
<p class="note">Full methods, every number with its source file, all deviations and the code are in the project repository (REPORT.md, DESIGN.md, RESEARCH_TREE.md). This page is generated by <code>docs/briefing/build_briefing.py</code> from the committed results.</p>
</main>
"""

OUT.write_text(HTML, encoding="utf-8")
print(f"wrote {OUT} ({len(HTML):,} chars)")
