"""M5 analyses on the frozen arms outputs (DESIGN v3.2 §8; amendment A7: the catalogue is MD-global, MD-strict is
reported alongside): per-gene subtypes (§8.1), biomarker modules (§8.2), within vs across lineage (§8.3, through
``audit.cross_lineage_stats`` / ``cross_lineage_labels``), co-location and cutting (§8.4-§8.5, ``audit.locus_resolution``
/ ``audit.amplicon_test``), culture (§8.6) and hypermutation / TP53 (§8.7).

Shared infrastructure for every M5-M8 script:
- ``Workspace``: where frozen inputs are read and outputs written (the project's results/ cache/ configs/ external/,
  or a toy root with the same layout);
- ``RealSource`` / ``ToySource``: the data-level loaders (``pipeline.Data``, D effects and expression, controls,
  gene info, CGT, symbols, paralogs). Real effects are read only through ``guard.load_effects(scope='dev')`` (via
  ``audit._real_loaders`` and ``pipeline.Data.from_real``); nothing here reads R.
- ``prepare_toy``: runs the unchanged engine on the toy (``toysim.make_toy`` + ``pipeline.run_pipeline``) and writes
  the M4 outputs in the real layout, so the M5-M8 scripts run end to end on it.

Numbers come from design_params.json (through ``params``) and configs/m5_m8_params.json (``cfg()``).
"""
from __future__ import annotations

import functools
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats as sps

from . import arms as A
from . import audit as AU
from . import params as PR
from . import paths as P
from . import stats as ST


@functools.lru_cache(maxsize=1)
def cfg() -> dict:
    """configs/m5_m8_params.json."""
    return P.load_config("m5_m8_params")


def lists() -> list:
    return list(cfg()["catalogue"]["lists"])


PRIMARY = "MD-global"
SECONDARY = "MD-strict"


# ------------------------------------------------------------------------------------------------
# Workspace (file layout) and sources (data loaders)
# ------------------------------------------------------------------------------------------------


class Workspace:
    """Frozen inputs and outputs in the project layout (``root=None``) or under a toy root with the same layout."""

    def __init__(self, root=None):
        self.root = Path(root) if root is not None else P.PROJECT
        self.is_toy = root is not None
        self.results = self.root / "results"
        self.cache = self.root / "cache"
        self.configs = self.root / "configs"
        self.external = self.root / "external" if self.is_toy else P.EXTERNAL

    # ---- paths ---------------------------------------------------------------------------------
    @property
    def timing_dir(self):
        return self.results / "timing" if self.is_toy else None

    @property
    def registry(self):
        return self.results / "registry.csv"

    def rel(self, path) -> str:
        return Path(path).resolve().relative_to(self.root.resolve()).as_posix()

    # ---- frozen M4 inputs ------------------------------------------------------------------------
    def arms(self, arms=None, columns=None, genes=None, selected_only: bool = False,
             universe: str = "eligible") -> pd.DataFrame:
        """Rows of the arms table, never the whole table (real: ~26 M rows): only ``columns`` and the rows matching
        the filters are materialised. Real workspace: ``endpoint_io.read_arms(columns=, filters=)`` when the blocked
        dataset exists (cache/benchmark/arms/part-*.parquet; the toy workspace writes the same layout under its own
        cache/), else a pyarrow read of a legacy arms.parquet with the same column / row pushdown.
        ``universe='ctrl'``: the U_ctrl dataset cache/benchmark/arms_ctrl/ (DESIGN §2; empty when it does not exist);
        its rows never enter a top-k list, the null, P1 or P2."""
        from . import endpoint_io as IO
        bench = self.cache / "benchmark"
        if universe != "eligible":
            d = bench / f"arms_{universe}"
            if not any(d.glob("part-*.parquet")):
                return pd.DataFrame(columns=list(columns or ["gene", "feature"]))
            a = IO.read_arms(columns=columns, arms=arms, genes=None if genes is None else [str(g) for g in genes],
                             selected_only=selected_only, dataset_dir=d)
        elif hasattr(IO, "read_arms") and any((bench / "arms").glob("part-*.parquet")):
            a = IO.read_arms(columns=columns, arms=arms, genes=None if genes is None else [str(g) for g in genes],
                             selected_only=selected_only, dataset_dir=bench / "arms")
        else:                                   # legacy single file / reader absent: same pushdown via pyarrow
            filters = []
            if arms is not None:
                filters.append(("arm", "in", list(arms)))
            if genes is not None:
                filters.append(("gene", "in", sorted({str(g) for g in genes})))
            if selected_only:
                filters.append(("selected", "==", True))
            a = read_arms_dataset(bench, columns, filters)
        if "gene" in a:
            a["gene"] = a["gene"].astype(str)
        return a.reset_index(drop=True)

    def in_universe(self, pairs: pd.DataFrame, universe: str = "eligible") -> set:
        """(gene, feature) keys of ``pairs`` that are rows of U (S1-obs rows of the arms table; ``universe='ctrl'``:
        of U_ctrl), read for the pairs' genes only."""
        if not len(pairs):
            return set()
        a = self.arms(arms=["S1-obs"], columns=["gene", "feature"], genes=pairs["gene"].unique(), universe=universe)
        return set(zip(a["gene"].astype(str), a["feature"])) & set(zip(pairs["gene"].astype(str), pairs["feature"]))

    def ctrl_genes(self) -> set:
        """U_ctrl genes (configs/u_ctrl_genes.txt of this workspace; empty if absent)."""
        f = self.configs / "u_ctrl_genes.txt"
        return {x.strip() for x in f.read_text(encoding="utf-8").splitlines() if x.strip()} if f.exists() else set()

    def has_ctrl(self) -> bool:
        """True when the U_ctrl run's outputs exist (cache/benchmark/arms_ctrl parts and results/md_ctrl lists)."""
        return (bool(self.ctrl_genes()) and any((self.cache / "benchmark" / "arms_ctrl").glob("part-*.parquet"))
                and all((self.results / "md_ctrl" / f).exists() for f in ("md_global.csv", "md_strict.csv")))

    def md(self, name: str, universe: str = "eligible") -> pd.DataFrame:
        d = self.results / ("md" if universe == "eligible" else f"md_{universe}")
        f = d / ("md_global.csv" if name == PRIMARY else "md_strict.csv")
        return pd.read_csv(f, dtype={"gene": str})

    def lam(self) -> pd.DataFrame:
        return self.read_csv("results/benchmark/lambda.csv", gene=str)

    def read_csv(self, rel: str, required: bool = True, **dtypes) -> pd.DataFrame | None:
        f = self.root / rel
        if not f.exists():
            if required:
                raise FileNotFoundError(f"{rel} missing (run the script that writes it first)")
            return None
        d = pd.read_csv(f, dtype=dtypes or {"gene": str}, keep_default_na=True)
        if "gene" in d:
            d["gene"] = d["gene"].astype(str)
        return d

    def read_json(self, rel: str) -> dict | None:
        f = self.root / rel
        return json.loads(f.read_text(encoding="utf-8")) if f.exists() else None

    # ---- outputs ---------------------------------------------------------------------------------
    def write_csv(self, df: pd.DataFrame, rel: str) -> str:
        f = self.root / rel
        f.parent.mkdir(parents=True, exist_ok=True)
        df.to_csv(f, index=False)
        return rel

    def write_json(self, obj, rel: str) -> str:
        from . import endpoint_io as IO
        IO.write_json(self.root / rel, obj)
        return rel

    def timing(self, stem: str, rec: dict, start, workers: int) -> None:
        from . import analyses as AN
        AN.write_timing(stem, rec, self.timing_dir, start_utc=start, workers=workers)

    def register(self, script: str, outputs) -> None:
        from . import analyses as AN
        AN.append_registry(script, "toy" if self.is_toy else "ok", outputs, registry=self.registry)


def read_arms_dataset(bench: Path, columns, filters) -> pd.DataFrame:
    """pyarrow dataset read of cache/benchmark/arms/ (blocked dataset) or arms.parquet with column projection and
    row-filter pushdown (filters: list of (column, op, value), op in '==', 'in')."""
    import pyarrow.dataset as pds
    src = bench / "arms" if (bench / "arms").is_dir() else bench / "arms.parquet"
    d = pds.dataset(str(src), format="parquet")
    expr = None
    for col, op, val in filters or []:
        f = pds.field(col)
        e = f.isin(list(val)) if op == "in" else (f == val)
        expr = e if expr is None else (expr & e)
    cols = None if columns is None else list(dict.fromkeys(columns))
    return d.to_table(columns=cols, filter=expr).to_pandas()


def require_tag(check=None) -> None:
    """Every M5-M8 script refuses to run unless HEAD descends from tag freeze-bounds (pattern of run_gate.py)."""
    from . import endpoints as EP
    EP.require_descends("freeze-bounds", check=check)


def source_for(ws: Workspace):
    return ToySource(ws) if ws.is_toy else RealSource()


class RealSource:
    """Real D data through the project loaders (guard.load_effects scope 'dev'; p1 loaders). Never R."""

    is_toy = False

    def data(self, genes):
        from . import pipeline as PL
        return PL.Data.from_real(genes=sorted({str(g) for g in genes}))

    def loaders(self, data):
        return AU._real_loaders(data)

    @functools.cached_property
    def _controls_raw(self) -> pd.DataFrame:
        return pd.read_csv(P.CONFIGS / "positive_controls.csv", dtype=str)

    def controls(self) -> pd.DataFrame:
        return control_pairs(self._controls_raw)

    @functools.cached_property
    def gene_map(self) -> pd.Series:
        from . import p1
        gm = p1.gene_map()["symbol"]
        gm.index = gm.index.astype(str)
        return gm

    def symbols(self) -> pd.Series:
        return self.gene_map

    def symbol_to_entrez(self, syms) -> dict:
        from . import genesets as GS
        s = GS.symbol_to_entrez(list(dict.fromkeys(syms))).dropna().astype(str)
        return s.to_dict()

    def gene_info(self) -> pd.DataFrame:
        """Per gene: symbol, tier, variable, selective, r_within, eligible, sd (D effect SD of Part 1 gene_classes,
        same Chronos D effects), common_essential (DepMap common-essential flag = the pan-lethal label)."""
        from . import p1
        el = pd.read_csv(P.CONFIGS / "eligible_genes.csv", dtype={"entrez": str}).set_index("entrez")
        gc = p1.gene_classes()
        gc.index = gc.index.astype(str)
        out = pd.DataFrame(index=el.index.union(gc.index))
        out["symbol"] = el["symbol"].reindex(out.index).fillna(gc["symbol"].reindex(out.index))
        out["tier"] = el["tier"].reindex(out.index)
        for c in ("variable", "selective", "eligible"):
            out[c] = _as_bool(el[c].reindex(out.index))
        out["r_within"] = el["r_beyond"].reindex(out.index)
        out["sd"] = gc["sd"].reindex(out.index)
        out["common_essential"] = _as_bool(gc["common_essential"].reindex(out.index))
        return out

    def cgt(self) -> pd.DataFrame:
        from . import endpoints as EP
        c = pd.read_csv(P.CONFIGS / "clinical_ground_truth.csv")
        ez = self.symbol_to_entrez(sorted({x for v in c["target_set"] for x in EP._targets(v)}))
        c["targets"] = [[str(ez.get(x, x)) for x in EP._targets(v)] for v in c["target_set"]]
        return c

    def tmb(self, data) -> np.ndarray:
        from . import endpoint_io as IO
        f = cfg()["hypermutation"]["tmb_feature"]
        return IO.line_values(list(data.lines), [f])[f].to_numpy(dtype=np.float64)

    def paralogs(self) -> pd.DataFrame:
        from . import genesets as GS
        return GS.paralog_top()

    def spec_partners(self) -> set:
        from . import genesets as GS
        cp = GS.control_partner_genes()
        ez = self.symbol_to_entrez(list(cp["dependency_symbol"]) + list(cp["partner_symbol"]))
        return {(ez.get(a), ez.get(b)) for a, b in zip(cp["dependency_symbol"], cp["partner_symbol"])
                if ez.get(a) and ez.get(b)}

    def curated_pairs(self) -> set:
        """(dependency, biomarker gene or feature id) keys of positive controls, CGT, P1 biomarker_pairs and
        oncogene_pairs (DESIGN §12.3 not_in_curated_lists)."""
        keys = set()
        for _, r in self.controls().iterrows():
            keys |= {(r["gene"], r["feature"])} | {(r["gene"], h) for h in feature_genes_of_id(r["feature"])}
        for _, t in self.cgt().iterrows():
            f = t["feature_D"]
            if isinstance(f, str):
                for g in t["targets"]:
                    keys |= {(g, f)} | {(g, h) for h in feature_genes_of_id(f)}
        op = pd.read_csv(P.P1_ONCOGENE_PAIRS, dtype=str)
        keys |= set(zip(op["dependency_entrez"], op["biomarker_entrez"]))
        bp = pd.read_csv(P.P1_BIOMARKER_PAIRS, dtype=str)
        first = bp["biomarker"].fillna("").str.split().str[0].fillna("")
        ez = self.symbol_to_entrez(list(bp["dependency"]) + list(first))
        keys |= {(ez.get(d), ez.get(b)) for d, b in zip(bp["dependency"], first) if ez.get(d) and ez.get(b)}
        return keys

    def meth_profiled(self) -> int | None:
        """D lines with RRBS methylation (METH_absent == 0; DESIGN §12.2 coverage report)."""
        from . import p1
        from . import pipeline as PL
        c = PL.real_covariates(p1.dev_ids())
        return int((c["METH_absent"].astype(float) == 0).sum()) if "METH_absent" in c else None

    def fmeta(self) -> pd.DataFrame:
        from . import limits as LM
        return LM.fmeta_from_dictionary()


class ToySource:
    """The toy of ``prepare_toy`` (same seed / predictor), through the same interface as ``RealSource``."""

    is_toy = True

    def __init__(self, ws: Workspace):
        from . import toysim as TS
        meta = ws.read_json("cache/toy/meta.json")
        if meta is None:
            raise FileNotFoundError(f"{ws.root}: no toy prepared (m5.prepare_toy)")
        self.meta = meta
        self.toy = TS.make_toy(seed=int(meta["seed"]), predictor=meta["predictor"])

    def data(self, genes=None):
        """The toy's eligible genes; requested genes outside them (U_ctrl controls) are added in toy gene order."""
        from . import pipeline as PL
        if genes is None:
            return PL.Data.from_toy(self.toy)
        el = PL.toy_eligibility(self.toy)["eligible"]
        extra = set(map(str, genes)) - set(el.index[el].astype(str))
        if not extra:
            return PL.Data.from_toy(self.toy)
        keep = [g for g in self.toy.genes.index if bool(el.get(g, False)) or str(g) in extra]
        return PL.Data.from_toy(self.toy, genes=keep)

    def loaders(self, data):
        t = self.toy
        return (lambda gs: t.effects.reindex(index=data.lines, columns=list(gs)),
                lambda gs: t.expression.reindex(index=data.lines, columns=list(gs)))

    def controls(self) -> pd.DataFrame:
        return control_pairs(self.toy.controls)

    def symbols(self) -> pd.Series:
        return self.toy.genes["symbol"]

    def symbol_to_entrez(self, syms) -> dict:
        inv = {s: g for g, s in self.toy.genes["symbol"].items()}
        return {s: inv[s] for s in syms if s in inv}

    def gene_info(self) -> pd.DataFrame:
        from . import pipeline as PL
        gi = self.toy.genes
        el = PL.toy_eligibility(self.toy)
        out = pd.DataFrame(index=gi.index)
        out["symbol"] = gi["symbol"]
        out["tier"] = gi["tier"]
        out["variable"] = gi["variable"].astype(bool)
        out["selective"] = gi["selective"].astype(bool)
        out["eligible"] = el["eligible"].reindex(gi.index).fillna(False).astype(bool)
        out["r_within"] = el["r_within"].reindex(gi.index)
        out["sd"] = self.toy.effects.std()
        out["common_essential"] = False
        return out

    def cgt(self) -> pd.DataFrame:
        from . import endpoints as EP
        c = self.toy.cgt.copy()
        c["targets"] = [EP._targets(v) for v in c["target_set"]]
        return c

    def tmb(self, data) -> np.ndarray:
        col = cfg()["hypermutation"]["toy_tmb_covariate"]
        return self.toy.covariates[col].reindex(data.lines).to_numpy(dtype=np.float64)

    def paralogs(self) -> pd.DataFrame:
        """paralog_top-like table from the toy's paralog lists (no sequence identity in the toy: scores NaN)."""
        rows = {}
        k = int(PR.dp()["dictionary_s3"]["n_paralogs"])
        for g, p in self.toy.genes["paralogs"].items():
            par = [x for x in str(p).split("|") if x]
            rec = {f"par{r + 1}": (par[r] if r < len(par) else None) for r in range(k)}
            rec.update({f"par{r + 1}_score": np.nan for r in range(k)})
            rec["family_size"] = len(par) + 1
            rows[g] = rec
        return pd.DataFrame.from_dict(rows, orient="index")

    def spec_partners(self) -> set:
        return set()

    def curated_pairs(self) -> set:
        keys = set()
        for _, r in self.controls().iterrows():
            keys |= {(r["gene"], r["feature"])} | {(r["gene"], h) for h in feature_genes_of_id(r["feature"])}
        for _, t in self.cgt().iterrows():
            for g in t["targets"]:
                keys |= {(g, t["feature_D"])} | {(g, h) for h in feature_genes_of_id(t["feature_D"])}
        return keys

    def meth_profiled(self) -> int | None:
        return None

    def fmeta(self) -> pd.DataFrame:
        return self.data().fmeta


def _as_bool(s: pd.Series) -> pd.Series:
    return s.astype(str).str.lower().isin(["true", "1"])


def feature_genes_of_id(f) -> list:
    """Gene ids encoded in a dictionary feature id ``FAMILY:<entrez>`` (named features such as SIG:MSI: none)."""
    if not isinstance(f, str) or ":" not in f:
        return []
    tail = f.split(":", 1)[1]
    return [tail] if tail.isdigit() else []


def control_pairs(ctl: pd.DataFrame) -> pd.DataFrame:
    """Long (control_id, gene, feature, expected_sign, scope, family, gate_eligible): recoverable feature ids per
    control (positive_controls.csv: endpoints.control_features; toy: primary + secondary ids); feature None for a
    control without one (class_check rows), so every control is listed."""
    from . import endpoints as EP
    rows = []
    for _, c in ctl.iterrows():
        if "feature_primary_ids" in ctl.columns:
            feats = EP.control_features(c)
        else:
            feats = [x for x in (c.get("feature_primary"), c.get("feature_secondary")) if isinstance(x, str) and x]
        for f in feats or [None]:               # a control without a 𝒟 feature (e.g. LINEAGE-CHECK) keeps one row
            rows.append({"control_id": c["control_id"], "gene": str(c["dependency_entrez"]), "feature": f,
                         "expected_sign": c.get("expected_sign", ""), "scope": c.get("scope", ""),
                         "family": c.get("family", ""), "gate_eligible": str(c.get("gate_eligible", "")).lower()
                         in ("true", "1")})
    return pd.DataFrame(rows, columns=["control_id", "gene", "feature", "expected_sign", "scope", "family",
                                       "gate_eligible"])


# ------------------------------------------------------------------------------------------------
# Catalogue, pair sets, signs, shrinkage
# ------------------------------------------------------------------------------------------------

CAT_COLS = ["gene", "feature", "beta", "se", "p", "n_obs", "n_carriers", "collinear_with_design", "dropped_cov",
            "stage", "selection_freq"]


def catalogue(ws: Workspace, universe: str = "eligible") -> pd.DataFrame:
    """MD-global ∪ MD-strict pairs with membership flags; values from MD-global where present (β = S1-obs β, the
    D-stage test of both lists). ``universe='ctrl'``: the U_ctrl lists (results/md_ctrl/), used only for controls."""
    parts = []
    for name in lists():
        m = ws.md(name, universe)
        m = m[[c for c in CAT_COLS if c in m]].copy()
        m["list"] = name
        parts.append(m)
    c = pd.concat(parts, ignore_index=True)
    flags = c.groupby(["gene", "feature"])["list"].agg(lambda s: set(s))
    c = c.drop_duplicates(["gene", "feature"], keep="first").drop(columns="list").reset_index(drop=True)
    key = list(zip(c["gene"], c["feature"]))
    c["in_md_global"] = [PRIMARY in flags[k] for k in key]
    c["in_md_strict"] = [SECONDARY in flags[k] for k in key]
    c["lists"] = ["|".join(n for n in lists() if n in flags[k]) for k in key]
    return c


def sign_map(cat: pd.DataFrame) -> dict:
    return {(g, f): float(np.sign(b)) for g, f, b in zip(cat["gene"], cat["feature"], cat["beta"])}


def in_universe(pairs: pd.DataFrame, U: set) -> pd.DataFrame:
    return pairs[[k in U for k in zip(pairs["gene"], pairs["feature"])]]


def pair_set(cat: pd.DataFrame, s1_hits: pd.DataFrame, ctl: pd.DataFrame, U: set, top_s1: int | None = None) -> pd.DataFrame:
    """MD-global ∪ MD-strict ∪ S1-obs BH hits (optionally the top ``top_s1`` by rank) ∪ controls in U, with a
    '|'-joined ``source``."""
    rows = [(g, f, n) for g, f, ls in zip(cat["gene"], cat["feature"], cat["lists"]) for n in ls.split("|") if n]
    s = s1_hits.sort_values("rank") if top_s1 is not None else s1_hits
    if top_s1 is not None:
        s = s.head(int(top_s1))
    rows += [(g, f, "S1-obs BH") for g, f in zip(s["gene"], s["feature"])]
    rows += [(g, f, "control") for g, f in zip(ctl["gene"], ctl["feature"])]
    p = pd.DataFrame(rows, columns=["gene", "feature", "source"])
    p = in_universe(p, U)
    return p.groupby(["gene", "feature"], as_index=False, sort=False)["source"].agg(
        lambda v: "|".join(sorted(set(v))))


def shrunk_beta_map(s1: pd.DataFrame, keep: set | None = None) -> dict:
    """§6.1 EB shrinkage from an in-memory S1-obs table (toy-size inputs; ``stats.shrunk_beta`` unchanged)."""
    ok = s1["beta"].notna() & s1["se"].notna() & s1["sd_dtilde"].notna()
    s = s1[ok]
    b = ST.shrunk_beta(s["beta"].to_numpy(), s["se"].to_numpy(), s["sd_dtilde"].to_numpy())
    out = dict(zip(zip(s["gene"], s["feature"]), b))
    return out if keep is None else {k: v for k, v in out.items() if k in keep}


def eb_fit_streamed(b: np.ndarray, s: np.ndarray, chunk: int) -> tuple:
    """The EM of ``stats.eb_scale_mixture`` (same grid, start, tolerance and iteration cap) with the n x K
    likelihood recomputed in row chunks, so memory stays O(chunk x K) on the ~6.5 M S1-obs rows of real U.
    Returns (sigma grid, weights)."""
    rp = PR.replication_params()
    tiny = PR.stats_params()["tiny"]
    ok = np.isfinite(b) & np.isfinite(s) & (s > 0)
    bo, so = b[ok], s[ok]
    lo = rp["eb_sigma_min_frac"] * so.min()
    hi = rp["eb_sigma_max_mult"] * np.sqrt(max(np.max(bo ** 2 - so ** 2), (2 * lo) ** 2))
    sig = np.geomspace(lo, max(hi, lo * 2), rp["eb_n_sigma"])
    pi = np.full(sig.size, 1.0 / sig.size)
    n = bo.size
    for _ in range(rp["eb_max_iter"]):
        acc = np.zeros(sig.size)
        for i in range(0, n, chunk):
            w = _eb_lik(bo[i:i + chunk], so[i:i + chunk], sig) * pi
            w /= np.maximum(w.sum(1, keepdims=True), tiny)
            acc += w.sum(0)
        new = acc / n
        done = np.max(np.abs(new - pi)) < rp["eb_tol"]
        pi = new
        if done:
            break
    return sig, pi


def _eb_lik(b, s, sig):
    var = sig[None, :] ** 2 + s[:, None] ** 2
    return np.exp(-b[:, None] ** 2 / var / 2) / np.sqrt(var)


def eb_posterior(b: np.ndarray, s: np.ndarray, sig: np.ndarray, pi: np.ndarray) -> np.ndarray:
    tiny = PR.stats_params()["tiny"]
    out = np.full(b.shape, np.nan)
    ok = np.isfinite(b) & np.isfinite(s) & (s > 0)
    if ok.any():
        var = sig[None, :] ** 2 + s[ok][:, None] ** 2
        w = _eb_lik(b[ok], s[ok], sig) * pi
        w /= np.maximum(w.sum(1, keepdims=True), tiny)
        out[ok] = (w * (sig[None, :] ** 2 / var)).sum(1) * b[ok]
    return out


def shrunk_beta_ws(ws: Workspace, keep: set) -> dict:
    """§6.1 shrunk β for the ``keep`` pairs, with the EB prior fitted once on every S1-obs U row (incl. F2): pass 1
    reads only (beta, se, sd_dtilde) of S1-obs rows, pass 2 the rows of the keep pairs' genes."""
    a = ws.arms(arms=["S1-obs"], columns=["beta", "se", "sd_dtilde"])
    sd = a["sd_dtilde"].to_numpy(dtype=np.float64)
    sig, pi = eb_fit_streamed(a["beta"].to_numpy(dtype=np.float64) * sd, a["se"].to_numpy(dtype=np.float64) * sd,
                              int(cfg()["memory"]["eb_chunk_rows"]))
    del a
    genes = sorted({g for g, _ in keep})
    cols = ["gene", "feature", "beta", "se", "sd_dtilde"]
    ctrl = ws.ctrl_genes() & set(genes)                 # U_ctrl controls: their rows come from the U_ctrl dataset
    ge = [g for g in genes if g not in ctrl]
    parts = [ws.arms(arms=["S1-obs"], columns=cols, genes=ge)] if ge or not ctrl else []
    if ctrl:
        parts.append(ws.arms(arms=["S1-obs"], columns=cols, genes=sorted(ctrl), universe="ctrl"))
    parts = [x for x in parts if len(x)] or parts[:1]
    r = parts[0] if len(parts) == 1 else pd.concat(parts, ignore_index=True)
    r = r[[k in keep for k in zip(r["gene"].astype(str), r["feature"])]]
    sd = r["sd_dtilde"].to_numpy(dtype=np.float64)
    with np.errstate(divide="ignore", invalid="ignore"):
        post = eb_posterior(r["beta"].to_numpy(dtype=np.float64) * sd, r["se"].to_numpy(dtype=np.float64) * sd,
                            sig, pi) / sd
    return dict(zip(zip(r["gene"], r["feature"]), post))


def analysis_pairs(ws: Workspace, src, top_s1: int | None = None, gate_only: bool = False,
                   u_ctrl_controls: bool = False) -> tuple:
    """(catalogue, S1-obs BH hits, pair set): MD-global ∪ MD-strict ∪ S1-obs BH (top ``top_s1``) ∪ positive
    controls in U, read without loading U (only S1-obs selected rows and the control genes' rows).
    ``u_ctrl_controls``: also the controls whose gene is in U_ctrl (DESIGN §2), taken where they are rows of U_ctrl
    (§8.3 uses every control); the pair set then has a ``universe`` column ('eligible' / 'ctrl')."""
    cat = catalogue(ws)
    s1h = ws.arms(arms=["S1-obs"], columns=["gene", "feature", "beta", "rank"], selected_only=True)
    ctl = src.controls()
    if gate_only:
        ctl = ctl[ctl["scope"] == "gate"]
    cg = ws.ctrl_genes() if u_ctrl_controls else set()
    ctl_e = ctl[~ctl["gene"].astype(str).isin(cg)]
    U = set(zip(cat["gene"], cat["feature"])) | set(zip(s1h["gene"], s1h["feature"])) | ws.in_universe(ctl_e)
    pairs = pair_set(cat, s1h, ctl_e, U, top_s1=top_s1)
    if not u_ctrl_controls:
        return cat, s1h, pairs
    pairs = pairs.assign(universe="eligible")
    if cg:
        ctl_c = ctl[ctl["gene"].astype(str).isin(cg)]
        Uc = ws.in_universe(ctl_c, "ctrl")
        if Uc:
            pc = pair_set(cat.iloc[:0], s1h.iloc[:0], ctl_c, Uc).assign(universe="ctrl")
            pairs = pd.concat([pairs, pc], ignore_index=True)
    return cat, s1h, pairs


class RankLookup:
    """``endpoints.RankIndex.rank_of`` semantics (§10.5 rank_of_feature_D) without the whole arms table: per arm the
    sorted ranks of the de-duplicated rows (two float columns) plus the rows of the queried genes."""

    def __init__(self, ws: Workspace, arms_list, genes, universe: str = "eligible"):
        """``universe='ctrl'``: ranks within U_ctrl (its own de-duplicated rows), for U_ctrl controls / CGT targets."""
        self.dr, self.by = {}, {}
        for arm in arms_list:
            r = ws.arms(arms=[arm], columns=["rank", "rank_dedup"], universe=universe)
            self.dr[arm] = np.sort(r.loc[r["rank_dedup"].notna(), "rank"].to_numpy(dtype=float))
        t = ws.arms(arms=list(arms_list), columns=["arm", "gene", "feature", "rank", "rank_dedup", "variant_group",
                                                   "beta", "access_class"], genes=genes, universe=universe)
        t["gene"] = t["gene"].astype(str)
        for arm, s in t.groupby("arm"):
            self.by[arm] = {g: x for g, x in s.groupby("gene")}
        s1 = t[t["arm"] == "S1-obs"]
        self.access = dict(zip(zip(s1["gene"], s1["feature"]), s1["access_class"]))

    def rank_of(self, arm: str, gene: str, feature: str) -> tuple:
        G = self.by.get(arm, {}).get(str(gene))
        if G is None:
            return np.nan, np.nan
        me = G[G["feature"] == feature]
        if me.empty:
            return np.nan, np.nan
        me = me.iloc[0]
        r = me["rank"]
        if pd.isna(r):
            return np.nan, float(me["beta"])
        better = int(np.searchsorted(self.dr[arm], r, side="left"))
        d = G[G["rank_dedup"].notna()]
        vg = me["variant_group"]
        same = (d["feature"] == feature) | ((d["variant_group"] == vg) if (isinstance(vg, str) and vg) else False)
        better -= int((same & (d["rank"] < r)).sum())
        return float(1 + better), float(me["beta"])


# ------------------------------------------------------------------------------------------------
# Parallel map over pair chunks (payload staged through workerstore, as in run_confounder_audit)
# ------------------------------------------------------------------------------------------------

_POOL: dict = {}


def _task_cross_lineage(pl, chunk):
    return AU.cross_lineage_stats(pl["data"], chunk, pl["bsh"], pl["sign_of"])


def _task_locus(pl, chunk):
    return AU.locus_resolution(pl["data"], chunk, pl["sign_of"])


TASKS = {"cross_lineage": _task_cross_lineage, "locus": _task_locus}


def _pool_init(path):
    from . import workerstore as WS
    _POOL.clear()
    _POOL.update(WS.load(path))
    try:
        from threadpoolctl import threadpool_limits
        threadpool_limits(int(PR.null_params()["blas_threads"]))
    except Exception:  # noqa: BLE001 - optional
        pass


def _pool_task(args):
    name, chunk = args
    return TASKS[name](_POOL, chunk)


def chunks_of(pairs: pd.DataFrame, workers: int) -> list:
    n = max(1, int(PR.null_params()["chunk_genes"]))
    if workers > 1:
        n = max(1, min(n, -(-len(pairs) // workers)))
    return [pairs.iloc[i:i + n] for i in range(0, len(pairs), n)] or [pairs]


def map_pairs(task: str, data, pairs: pd.DataFrame, workers: int, **extra) -> list:
    if workers > 1:
        try:                                    # pre-spawn memory guard, when the runtime config provides one
            from . import analyses as AN
            workers = int(AN.guard_workers(workers, task)[0])
        except Exception:  # noqa: BLE001 - guard optional
            pass
    chunks = chunks_of(pairs, workers)
    if workers <= 1 or len(chunks) <= 1:
        pl = {"data": data, **extra}
        return [TASKS[task](pl, c) for c in chunks]
    from concurrent.futures import ProcessPoolExecutor

    from . import workerstore as WS
    with WS.staged({"data": data.for_workers(include_obs=True), **extra}, task) as path, \
            ProcessPoolExecutor(max_workers=workers, initializer=_pool_init, initargs=(path,)) as ex:
        return list(ex.map(_pool_task, [(task, c) for c in chunks]))


# ------------------------------------------------------------------------------------------------
# §8.1 per-gene subtypes
# ------------------------------------------------------------------------------------------------


def frame(data, M: np.ndarray) -> pd.DataFrame:
    return pd.DataFrame(M, index=data.lines, columns=data.genes)


def dependency_status(data) -> tuple:
    """(dependent, accurate-dependent) D lines x genes (NaN where the effect is NaN), DESIGN §2."""
    from . import eligibility as ELG
    Y = frame(data, data.target("obs")).where(frame(data, data.rows_ok).astype(bool))
    dep = ELG.dependent(Y, ELG.dependent_rule(Y))
    acc = ELG.accurate_dependent_D(dep, frame(data, data.pred))
    return dep, acc


def line_contributions(data, g: str, feats: list, beta: np.ndarray) -> tuple:
    """β_gd'·(x̃_ℓd' − x̄_d') on g's rows for the given features: x̃ = 𝒟 standardised on D, residualised on the EN
    unpenalised block (design blocks + NA indicators) estimated on g's rows (the full-D EN-pred fit, §4.1).
    Returns (row mask, n_rows x k contribution matrix)."""
    gj = data._gpos[g]
    rows = data.rows_ok[:, gj].astype(bool)
    pos = data.fmeta.index.get_indexer(list(feats))
    U = np.hstack([data.cov.Z(rows, data.cov.design_names), data.na_flags[rows]])
    Xt = ST.resid(ST.orth(U), data.Xstd[rows][:, pos])
    return rows, (Xt - Xt.mean(0)) * np.asarray(beta, dtype=np.float64)[None, :]


def subtypes(data, en_pred: pd.DataFrame, cat: pd.DataFrame, dep: pd.DataFrame, acc: pd.DataFrame,
             symbols: pd.Series) -> dict:
    """§8.1. S(g, d) = accurate-dependent lines of g whose line share s_ℓd (β = full-D EN-pred coefficient over the
    EN-pred selected features of g, both stages; lines with a negative sum only) is >= line_share_min; binary d also
    needs carrier status. Returns {'summary', 'lines', 'genes'}."""
    s8 = PR.dp()["subtypes_s8"]
    smin = float(s8["line_share_min"])
    nl = int(cfg()["subtypes"]["max_lineages_listed"])
    lin = data.cov_df.get("lineage_raw", data.cov_df["lineage"]).astype(str).to_numpy()
    lines = np.asarray(data.lines)
    en = en_pred[en_pred["selected"].astype(bool)]
    by_gene = {g: s for g, s in en.groupby("gene")}
    summ, long, member_sets = [], [], {}
    for g, cg in cat.groupby("gene", sort=False):
        sel = by_gene.get(g)
        gj = data._gpos[g]
        dcol, acol = dep[g].to_numpy(), acc[g].to_numpy()
        base = {"gene": g, "symbol": symbols.get(g, ""), "n_dependent": int(np.nansum(dcol == 1)),
                "n_accurate_dependent": int(np.nansum(acol == 1))}
        rows, C = (None, None)
        if sel is not None and len(sel):
            rows, C = line_contributions(data, g, list(sel["feature"]), sel["beta_full"].fillna(0).to_numpy())
            total = C.sum(1)
        for _, r in cg.iterrows():
            f = r["feature"]
            rec = dict(base, feature=f, family=data.fmeta.at[f, "family"], d_type=data.fmeta.at[f, "d_type"],
                       lists=r["lists"], beta_s1=r["beta"], beta_full_en=np.nan, n_members=0, note="")
            if sel is None or f not in set(sel["feature"]):
                rec["note"] = "not EN-pred selected"
                summ.append(rec)
                continue
            j = list(sel["feature"]).index(f)
            rec["beta_full_en"] = float(sel["beta_full"].iloc[j])
            if not np.isfinite(rec["beta_full_en"]) or rec["beta_full_en"] == 0:
                rec["note"] = "full-D EN-pred coefficient is 0 (selected by stability only)"
            with np.errstate(divide="ignore", invalid="ignore"):
                share = np.where(total < 0, C[:, j] / total, np.nan)
            ok = (share >= smin) & (acol[rows] == 1)
            x = data.X[rows, data.fmeta.index.get_loc(f)]
            if rec["d_type"] == "binary":
                ok &= x == 1
                rec["n_carriers"] = int(np.nansum(data.X[data.rows_ok[:, gj], data.fmeta.index.get_loc(f)] == 1))
            idx = np.flatnonzero(rows)[ok]
            rec["n_members"] = int(idx.size)
            rec["member_frac_of_accurate"] = idx.size / base["n_accurate_dependent"] if base["n_accurate_dependent"] else np.nan
            vc = pd.Series(lin[idx]).value_counts()
            rec["n_lineages"] = int(vc.size)
            rec["lineages"] = ";".join(f"{k}:{v}" for k, v in vc.head(nl).items())
            summ.append(rec)
            member_sets.setdefault(g, []).append(set(idx))
            long += [{"gene": g, "feature": f, "ModelID": lines[i], "share": float(share[k]), "lineage": lin[i]}
                     for k, i in zip(np.flatnonzero(ok), idx)]
    genes = []
    for g in cat["gene"].unique():
        dcol, acol = dep[g].to_numpy(), acc[g].to_numpy()
        sets = member_sets.get(g, [])
        assigned = set().union(*sets) if sets else set()
        cnt = pd.Series([i for s in sets for i in s]).value_counts() if sets else pd.Series(dtype=int)
        genes.append({"gene": g, "symbol": symbols.get(g, ""), "n_dependent": int(np.nansum(dcol == 1)),
                      "n_accurate_dependent": int(np.nansum(acol == 1)),
                      "n_unexplained": int(np.nansum((dcol == 1) & (acol != 1))),
                      "n_accurate_unassigned": int(np.nansum(acol == 1)) - len(assigned),
                      "n_multi_label": int((cnt > 1).sum()), "n_catalogue_pairs": int((cat["gene"] == g).sum())})
    return {"summary": pd.DataFrame(summ), "lines": pd.DataFrame(long, columns=["gene", "feature", "ModelID", "share",
                                                                                "lineage"]),
            "genes": pd.DataFrame(genes)}


# ------------------------------------------------------------------------------------------------
# §8.2 biomarker modules
# ------------------------------------------------------------------------------------------------


def carriers(data, f: str, sign: float) -> np.ndarray:
    """Binary d: x == 1. Continuous d: the dependency side of the D median (configs/m5_m8_params subtypes)."""
    x = data.X[:, data.fmeta.index.get_loc(f)]
    if data.fmeta.at[f, "d_type"] == "binary":
        return x == 1
    with np.errstate(invalid="ignore"):
        return np.isfinite(x) & (sign * (x - np.nanmedian(x)) < 0)


def effective_clusters(Y: pd.DataFrame) -> int:
    """Number of gene clusters: average linkage on 1 − |r| of D effects, cut at 1 − gene_cluster_abs_r."""
    from scipy.cluster.hierarchy import fcluster, linkage
    from scipy.spatial.distance import squareform
    if Y.shape[1] < 2:
        return int(Y.shape[1])
    s8 = PR.dp()["subtypes_s8"]
    D = (1 - Y.corr().abs()).fillna(1).to_numpy()
    np.fill_diagonal(D, 0)
    Z = linkage(squareform(np.clip(D, 0, None), checks=False), method=s8["gene_cluster_linkage"])
    return int(np.unique(fcluster(Z, t=1 - float(s8["gene_cluster_abs_r"]), criterion="distance")).size)


def modules(data, cat: pd.DataFrame, acc: pd.DataFrame, symbols: pd.Series) -> dict:
    """§8.2 per list (MD-global; MD-strict alongside): M(d, sign) = genes with (g, d, same sign) in the list;
    member lines = carriers accurate-dependent on >= module_line_frac_min of the module genes."""
    s8 = PR.dp()["subtypes_s8"]
    fmin = float(s8["module_line_frac_min"])
    nl = int(cfg()["subtypes"]["max_lineages_listed"])
    lin = data.cov_df.get("lineage_raw", data.cov_df["lineage"]).astype(str).to_numpy()
    Y = frame(data, data.target("obs")).where(frame(data, data.rows_ok).astype(bool))
    rows, long = [], []
    for name in lists():
        L = cat[cat["lists"].str.split("|").map(lambda v: name in v)]
        for (f, s), sub in L.assign(sign=np.sign(L["beta"])).groupby(["feature", "sign"], sort=True):
            G = list(dict.fromkeys(sub["gene"]))
            car = carriers(data, f, s)
            a = acc[G].to_numpy(dtype=np.float64)
            with np.errstate(invalid="ignore"):
                frac = np.nanmean(np.where(np.isnan(a), np.nan, a == 1), axis=1) if a.size else np.full(len(car), np.nan)
            mem = car & (np.nan_to_num(frac, nan=-1) >= fmin)
            idx = np.flatnonzero(mem)
            vc = pd.Series(lin[idx]).value_counts()
            rows.append({"list": name, "feature": f, "sign": int(s), "family": data.fmeta.at[f, "family"],
                         "n_genes": len(G), "genes": "|".join(G), "symbols": "|".join(str(symbols.get(g, "")) for g in G),
                         "n_effective_clusters": effective_clusters(Y[G]), "n_carriers": int(car.sum()),
                         "n_member_lines": int(idx.size), "n_lineages": int(vc.size),
                         "lineages": ";".join(f"{k}:{v}" for k, v in vc.head(nl).items())})
            long += [{"list": name, "feature": f, "sign": int(s), "ModelID": data.lines[i], "lineage": lin[i]}
                     for i in idx]
    return {"modules": pd.DataFrame(rows), "lines": pd.DataFrame(long, columns=["list", "feature", "sign", "ModelID",
                                                                               "lineage"])}


# ------------------------------------------------------------------------------------------------
# §8.6 culture
# ------------------------------------------------------------------------------------------------


def cramers_v(member: np.ndarray, values: np.ndarray) -> float:
    """Cramér's V of a binary membership vector with a categorical vector (2 x K table, levels with lines only)."""
    t = pd.crosstab(pd.Series(member, dtype=bool), pd.Series(values, dtype=str))
    t = t.loc[:, t.sum(0) > 0]
    if t.shape[0] < 2 or t.shape[1] < 2:
        return np.nan
    from scipy.stats.contingency import association
    return float(association(t.to_numpy(), method="cramer"))


def culture_values(data, var: dict) -> np.ndarray:
    v = data.cov_df[var["col"]].astype(str).to_numpy()
    return (v == var["level"]).astype(str) if "level" in var else v


def culture_audit(data, summary: pd.DataFrame, lines: pd.DataFrame, sign_of: dict) -> pd.DataFrame:
    """§8.6: Cramér's V of subtype membership (members vs the other lines with a g effect) with each culture /
    library variable; V >= culture_flag_cramers_v flags the subtype, unless the confounder model within the stratum
    of the flagging variable that holds most members keeps the catalogue sign at p < rescue_p_max."""
    s8 = PR.dp()["subtypes_s8"]
    vthr = float(s8["culture_flag_cramers_v"])
    cc = cfg()["culture"]
    pmax = float(cc["rescue_p_max"])
    vars_ = cc["variables"]
    vals = {v["name"]: culture_values(data, v) for v in vars_}
    pos = {m: i for i, m in enumerate(data.lines)}
    mem = lines.groupby(["gene", "feature"])["ModelID"].apply(list).to_dict() if len(lines) else {}
    out = []
    for _, r in summary.iterrows():
        g, f = r["gene"], r["feature"]
        rows = data.rows_ok[:, data._gpos[g]].astype(bool)
        m = np.zeros(len(data.lines), dtype=bool)
        m[[pos[x] for x in mem.get((g, f), [])]] = True
        rec = {"gene": g, "feature": f, "lists": r.get("lists", ""), "n_members": int(m.sum())}
        best, bv = None, -np.inf
        for v in vars_:
            cv = cramers_v(m[rows], vals[v["name"]][rows]) if m.any() else np.nan
            rec[f"V_{v['name']}"] = cv
            if np.isfinite(cv) and cv > bv:
                best, bv = v["name"], cv
        rec["V_max"] = bv if best else np.nan
        rec["V_max_variable"] = best or ""
        rec["culture_flag"] = bool(best and bv >= vthr)
        rec["rescue_stratum"], rec["rescue_beta"], rec["rescue_p"], rec["rescued"] = "", np.nan, np.nan, False
        if rec["culture_flag"]:
            vv = vals[best]
            lev = pd.Series(vv[m]).value_counts().idxmax()
            fit = AU.pair_fit(data, g, f, extra_rows=vv == lev)
            s = sign_of.get((g, f), np.sign(r.get("beta_s1", np.nan)))
            rec.update(rescue_stratum=f"{best}={lev}", rescue_beta=fit["beta"], rescue_p=fit["p"],
                       rescued=bool(np.sign(fit["beta"]) == s and fit["p"] < pmax))
        rec["culture_confounded"] = bool(rec["culture_flag"] and not rec["rescued"])
        out.append(rec)
    return pd.DataFrame(out)


def culture_overlap(mods: pd.DataFrame, universe_genes, sym2ez: dict) -> pd.DataFrame:
    """§8.6: module genes vs the culture gene sets (design_params subtypes_s8.culture_genesets), one-sided
    hypergeometric over the arms' gene universe, BH over every (module, set) test."""
    s8 = PR.dp()["subtypes_s8"]
    Ug = set(map(str, universe_genes))
    N = len(Ug)
    out = []
    for name, syms in s8["culture_genesets"].items():
        S = {sym2ez[x] for x in syms if x in sym2ez} & Ug
        for _, m in mods.iterrows():
            G = set(str(m["genes"]).split("|")) & Ug
            k = len(G & S)
            p = float(sps.hypergeom.sf(k - 1, N, len(S), len(G))) if N and len(S) and len(G) else np.nan
            out.append({"list": m["list"], "feature": m["feature"], "sign": m["sign"], "gene_set": name,
                        "n_module_genes": len(G), "n_set_in_universe": len(S), "overlap": k, "p": p})
    o = pd.DataFrame(out, columns=["list", "feature", "sign", "gene_set", "n_module_genes", "n_set_in_universe",
                                   "overlap", "p"])
    o["q"] = ST.bh(o["p"].to_numpy()) if len(o) else []
    o["flag"] = o["q"] < float(s8["culture_overlap_bh_q"])
    return o


# ------------------------------------------------------------------------------------------------
# §8.7 hypermutation and TP53
# ------------------------------------------------------------------------------------------------


def hypermutated(data, tmb: np.ndarray) -> np.ndarray:
    """MSI ∪ TMB in the top dictionary_s3.hypermutated_top_frac of D (NaN TMB: not hypermutated by TMB)."""
    hc = cfg()["hypermutation"]
    top = float(PR.dp()["dictionary_s3"]["hypermutated_top_frac"])
    t = np.asarray(tmb, dtype=np.float64)
    cut = np.nanquantile(t, 1 - top) if np.isfinite(t).any() else np.inf
    with np.errstate(invalid="ignore"):
        return (data.cov_df[hc["msi_column"]].astype(str).to_numpy() == hc["msi_level"]) | (t >= cut)


def hypermut_tp53(data, pairs: pd.DataFrame, tmb: np.ndarray) -> pd.DataFrame:
    """§8.7: HOT/ALLELE/LOFBI determinants keep sign with |β| >= persist_beta_frac·|β_full| in non-hypermutated
    lines; SIG determinants keep sign within TP53-mutant lines. ``passed`` NaN = not testable in the subset."""
    hc = cfg()["hypermutation"]
    frac = float(hc["persist_beta_frac"])
    hyp = hypermutated(data, tmb)
    tp53m = data.cov_df[hc["tp53_column"]].astype(str).to_numpy() == hc["tp53_mut_level"]
    out = []
    for g, f in zip(pairs["gene"], pairs["feature"]):
        fam = data.fmeta.at[f, "family"]
        if fam in hc["persist_families"]:
            test, sub, rule = "hypermutation", ~hyp, "non-hypermutated lines"
        elif fam in hc["tp53_families"]:
            test, sub, rule = "TP53-mutant", tp53m, "TP53-mutant lines"
        else:
            continue
        full = AU.pair_fit(data, g, f)
        r = AU.pair_fit(data, g, f, extra_rows=sub)
        ok_fit = np.isfinite(r["beta"]) and np.isfinite(full["beta"])
        same = bool(np.sign(r["beta"]) == np.sign(full["beta"])) if ok_fit else None
        if test == "hypermutation" and ok_fit:
            passed = bool(same and abs(r["beta"]) >= frac * abs(full["beta"]))
        else:
            passed = same
        out.append({"gene": g, "feature": f, "family": fam, "test": test, "subset": rule,
                    "n_subset": int(r.get("n", 0)), "n_carriers_subset": int(r.get("n_carriers", 0)),
                    "beta_full": full["beta"], "beta_subset": r["beta"], "p_subset": r["p"],
                    "passed": passed, "failed": passed is False})
    return pd.DataFrame(out, columns=["gene", "feature", "family", "test", "subset", "n_subset", "n_carriers_subset",
                                      "beta_full", "beta_subset", "p_subset", "passed", "failed"])


# ------------------------------------------------------------------------------------------------
# §8.3 cross-lineage (parallel stats, one BH family)
# ------------------------------------------------------------------------------------------------


def cross_lineage(data, pairs: pd.DataFrame, bsh: dict, sign_of: dict, workers: int = 1) -> tuple:
    """audit.cross_lineage_stats over pair chunks, then audit.cross_lineage_labels on the whole set (F4 / F5 BH over
    every pair), plus the A7-iii caveat column. Returns (labels, per-lineage table)."""
    keep = set(zip(pairs["gene"], pairs["feature"]))
    parts = map_pairs("cross_lineage", data, pairs[["gene", "feature"]], workers,
                      bsh={k: v for k, v in bsh.items() if k in keep}, sign_of=sign_of)
    R = pd.concat([p[0] for p in parts], ignore_index=True)
    P5 = pd.concat([p[1] for p in parts if len(p[1])], ignore_index=True) if any(len(p[1]) for p in parts) \
        else pd.DataFrame(columns=["gene", "feature", "lineage", "beta", "p", "power"])
    lab = AU.cross_lineage_labels(R, P5)
    pmax = float(cfg()["cross_lineage"]["single_lineage_p_max"])
    if len(P5):
        P5 = P5.copy()
        P5["q"] = ST.bh(P5["p"].to_numpy())
    info = {(g, f): (set(str(pw).split("|")) if isinstance(pw, str) else set(), np.sign(b))
            for g, f, pw, b in zip(lab["gene"], lab["feature"], lab.get("lineages_powered", pd.Series(dtype=object)),
                                   lab.get("beta_pooled", pd.Series(np.nan, index=lab.index)))}
    n_sig = {}
    groups = P5.groupby(["gene", "feature"]) if len(P5) else []
    for (g, f), s in groups:
        pl, sg = info.get((g, f), (set(), 0))
        s = s[s["lineage"].isin(pl)]
        n_sig[(g, f)] = int(((s["p"] < pmax) & (np.sign(s["beta"]) == sg)).sum())
    lab["n_powered_nominal_same_sign"] = [n_sig.get((g, f), 0) for g, f in zip(lab["gene"], lab["feature"])]
    lab["single_lineage_driven"] = (lab["label"] == "cross-lineage") & (lab["n_powered_nominal_same_sign"] <= 1)
    lab = lab.merge(pairs, on=["gene", "feature"], how="left")
    return lab, P5


# ------------------------------------------------------------------------------------------------
# Toy preparation: the unchanged engine on the toy, M4 outputs in the real layout
# ------------------------------------------------------------------------------------------------


def prepare_toy(root, seed: int | None = None, predictor: str | None = None, workers: int = 1) -> Workspace:
    """Run ``pipeline.run_pipeline`` (four arms, MD-strict / MD-global, no null, no built-in audits) on
    ``toysim.make_toy(seed, predictor)`` and write what run_arms.py and run_confounder_audit.py write for real D:
    cache/benchmark/arms/part-00000.parquet, results/md/{md_global,md_global_all,md_strict,md_strict_all}.csv,
    results/benchmark/lambda.csv, results/confounders/delta.csv (+ delta_flags.csv). Idempotent."""
    from . import pipeline as PL
    from . import toysim as TS
    ws = Workspace(root)
    if ws.read_json("cache/toy/meta.json") is not None:
        return ws
    toy = TS.make_toy(seed=seed, predictor=predictor)
    data = PL.Data.from_toy(toy)
    out = PL.run_pipeline(data, arms=A.ARMS, null_B=0, workers=workers, audits=False)
    (ws.cache / "benchmark" / "arms").mkdir(parents=True, exist_ok=True)
    a = out["arms"].copy()
    a["gene"] = a["gene"].astype(str)
    a.to_parquet(ws.cache / "benchmark" / "arms" / "part-00000.parquet")       # layout of endpoint_io.read_arms
    for k, fn in (("md_global", "md_global.csv"), ("md_global_all", "md_global_all.csv"),
                  ("md_strict", "md_strict.csv"), ("md_strict_all", "md_strict_all.csv")):
        ws.write_csv(out[k], f"results/md/{fn}")
    ws.write_csv(out["lam"], "results/benchmark/lambda.csv")
    cat = catalogue(ws)
    s1h = a[(a["arm"] == "S1-obs") & a["selected"].astype(bool)]
    U = set(zip(a.loc[a["arm"] == "S1-obs", "gene"], a.loc[a["arm"] == "S1-obs", "feature"]))
    ctl = control_pairs(toy.controls)
    ctl = ctl[ctl["scope"] == "gate"]
    pairs = pair_set(cat, s1h, ctl, U)
    delta = AU.confounding_delta(data, pairs).merge(pairs, on=["gene", "feature"], how="left")
    ws.write_csv(delta, "results/confounders/delta.csv")
    ws.write_json({"seed": int(toy.seed), "predictor": toy.predictor, "scenario": toy.scenario,
                   "n_genes": len(data.genes), "n_md_global": int(len(out["md_global"])),
                   "n_md_strict": int(len(out["md_strict"]))}, "cache/toy/meta.json")
    return ws


# ------------------------------------------------------------------------------------------------
# Script plumbing
# ------------------------------------------------------------------------------------------------


def parse_args(argv, doc: str, workers: bool = False):
    """Common CLI: ``--toy DIR`` runs on a toy workspace (prepared on first use, engine unchanged) instead of the
    project's results/; ``--workers N`` for the pair-parallel steps (default compute_s14 foreground workers)."""
    import argparse
    ap = argparse.ArgumentParser(description=doc.splitlines()[0] if doc else None)
    ap.add_argument("--toy", default=None, help="toy workspace directory (never the project's results/)")
    if workers:
        ap.add_argument("--workers", type=int, default=None)
    a = ap.parse_args(argv)
    if workers:
        a.workers = int(a.workers or PR.null_params()["fg_workers"])
    return a


def open_workspace(a) -> tuple:
    """(Workspace, source) after the freeze-bounds check; a toy workspace is prepared on first use."""
    require_tag()
    if a.toy:
        prepare_toy(a.toy, workers=int(getattr(a, "workers", 1) or 1))
    ws = Workspace(a.toy)
    return ws, source_for(ws)
