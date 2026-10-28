"""Part 2 replication guard. Grant paths use a SYNTHETIC effects parquet, so no real replication or
non-cancerous effect value is ever read by the test suite."""
import numpy as np
import pandas as pd
import pytest

from biomark import guard, p1


@pytest.fixture
def tmp_log(tmp_path, monkeypatch):
    log = tmp_path / "replication_access_log.csv"
    monkeypatch.setattr(guard, "ACCESS_LOG", log)
    monkeypatch.setattr(guard, "LOCK_FILE", tmp_path / "log.lock")
    return log


@pytest.fixture
def fake_effects(tmp_path, monkeypatch):
    ids = list(p1.test_ids()[:4]) + list(p1.noncancerous_ids()[:2]) + list(p1.dev_ids()[:2])
    df = pd.DataFrame(np.arange(len(ids) * 3, dtype=np.float32).reshape(len(ids), 3),
                      index=pd.Index(ids, name="ModelID"), columns=["1", "2", "3"])
    path = tmp_path / "fake_crispr.parquet"
    df.to_parquet(path)
    monkeypatch.setattr(guard, "CRISPR_PARQUET", path)
    return df


@pytest.fixture
def unblinded(monkeypatch):
    monkeypatch.setenv(guard.UNBLIND_ENV, "1")
    monkeypatch.setattr(guard, "tag_exists", lambda: True)
    monkeypatch.setattr(guard, "dirty_entries", lambda: [])


def test_dev_allowed_real_data_no_log(tmp_log):
    ids = list(p1.dev_ids()[:5])
    y = guard.load_effects(ids, genes=["7157", "4609", "3845"])
    assert y.shape == (5, 3) and list(y.index) == ids and y.dtypes.eq(np.float32).all()
    assert not tmp_log.exists()


def test_dev_default_all_dev(tmp_log):
    y = guard.load_effects(genes=["7157"])
    assert y.shape == (1022, 1)
    assert not tmp_log.exists()


def test_test_refused_without_tag(tmp_log, fake_effects, monkeypatch):
    monkeypatch.setenv(guard.UNBLIND_ENV, "1")
    monkeypatch.setattr(guard, "tag_exists", lambda: False)
    monkeypatch.setattr(guard, "dirty_entries", lambda: [])
    with pytest.raises(guard.ReplicationGuardError, match="freeze-catalogue"):
        guard.load_effects(list(p1.test_ids()[:2]), scope="replication", purpose="x")
    assert not tmp_log.exists()


@pytest.mark.skipif(guard.tag_exists(), reason="freeze-catalogue already exists")
def test_test_refused_now_real_state(tmp_log, monkeypatch):
    monkeypatch.setenv(guard.UNBLIND_ENV, "1")
    with pytest.raises(guard.ReplicationGuardError):
        guard.load_effects(list(p1.test_ids()[:1]), scope="replication", purpose="probe")
    with pytest.raises(guard.ReplicationGuardError):
        guard.p1_test_metrics("", purpose="probe")
    assert not tmp_log.exists()


@pytest.mark.parametrize("missing", ["env", "clean", "purpose"])
def test_test_refused_each_condition(tmp_log, fake_effects, unblinded, monkeypatch, missing):
    purpose = "replication of catalogue v1"
    if missing == "env":
        monkeypatch.setenv(guard.UNBLIND_ENV, "0")
    elif missing == "clean":
        monkeypatch.setattr(guard, "dirty_entries", lambda: [" M src/biomark/guard.py"])
    else:
        purpose = "  "
    with pytest.raises(guard.ReplicationGuardError):
        guard.load_effects(list(p1.test_ids()[:2]), scope="replication", purpose=purpose)
    assert not tmp_log.exists()


def test_test_granted_is_logged(tmp_log, fake_effects, unblinded):
    ids = list(p1.test_ids()[:3]) + list(p1.dev_ids()[:1])
    y = guard.load_effects(ids, scope="replication", purpose="unit test", genes=["2", "3"])
    assert list(y.index) == ids and list(y.columns) == ["2", "3"]
    pd.testing.assert_frame_equal(y, fake_effects.loc[ids, ["2", "3"]])
    log = pd.read_csv(tmp_log)
    assert len(log) == 1 and list(log.columns) == guard.LOG_COLS
    r = log.iloc[0]
    assert r.purpose == "unit test" and r.n_ids == 4 and r.scope == "replication" and int(r.n_genes) == 2
    assert "test:3" in r.roles and r.utc_time.endswith("+00:00") and len(str(r.commit)) == 40
    guard.load_effects(ids[:1], scope="replication", purpose="second")
    assert len(pd.read_csv(tmp_log)) == 2


def test_scope_must_cover_roles(tmp_log, fake_effects, unblinded):
    with pytest.raises(guard.ReplicationGuardError):
        guard.load_effects(list(p1.test_ids()[:1]), scope="dev", purpose="x")
    with pytest.raises(guard.ReplicationGuardError):
        guard.load_effects(list(p1.noncancerous_ids()[:1]), scope="replication", purpose="limits: x")
    assert not tmp_log.exists()


def test_noncancerous_rules(tmp_log, fake_effects, monkeypatch):
    nc = list(p1.noncancerous_ids()[:2])
    monkeypatch.setattr(guard, "tag_exists", lambda: False)
    with pytest.raises(guard.ReplicationGuardError):
        guard.load_effects(nc, scope="limits", purpose="limits: normal-tissue toxicity")
    monkeypatch.setattr(guard, "tag_exists", lambda: True)
    with pytest.raises(guard.ReplicationGuardError):
        guard.load_effects(nc, scope="limits", purpose="toxicity")
    assert not tmp_log.exists()
    y = guard.load_effects(nc, scope="limits", purpose="limits: normal-tissue toxicity")
    assert y.shape == (2, 3)
    assert len(pd.read_csv(tmp_log)) == 1


def test_unassigned_refused(tmp_log):
    sc_ids = guard.screened_ids()
    in_pop = set(p1.splits().index)
    unassigned = [i for i in sc_ids if i not in in_pop]
    assert len(unassigned) == 69
    with pytest.raises(guard.ReplicationGuardError):
        guard.load_effects(unassigned[:1], scope="dev")
    assert not tmp_log.exists()


def test_p1_loaders_refuse_effects_and_test_metrics():
    with pytest.raises(PermissionError):
        p1.raw_matrix("crispr")
    with pytest.raises(PermissionError):
        p1.read_parquet(guard.CRISPR_PARQUET)
    with pytest.raises(PermissionError):
        p1.read_parquet(p1.P.P1_TEST_RESULTS / "anything.parquet")


def test_p1_test_metrics_refused_and_path_confined(tmp_log, monkeypatch):
    monkeypatch.setattr(guard, "tag_exists", lambda: False)
    with pytest.raises(guard.ReplicationGuardError):
        guard.p1_test_metrics("", purpose="x")
    with pytest.raises(ValueError):
        guard.p1_test_metrics("../runs", purpose="x")
    assert not tmp_log.exists()


def test_dirty_entries_rules(monkeypatch):
    lines = [" M src/biomark/p1.py", "?? results/foo.csv", "?? scripts/new.py",
             " M results/replication_access_log.csv"]
    monkeypatch.setattr(guard.gitutil, "status_porcelain", lambda repo, untracked=True: lines)
    assert guard.dirty_entries() == [" M src/biomark/p1.py", "?? scripts/new.py"]
