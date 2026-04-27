from pathlib import Path
from tempfile import TemporaryDirectory
from unittest import SkipTest

import numpy as np
import MDAnalysis as mda


TOP = "tests/data/trex1.tpr"
GMX_DIR = Path("data_gmx")
MDA_DIR = Path("data_mda")


def _window_ids(base):
    return sorted(
        int(p.name) for p in base.iterdir()
        if p.is_dir() and p.name.isdigit())


def _require_split_outputs():
    missing = [
        str(base) for base in (GMX_DIR, MDA_DIR)
        if not base.exists()
    ]
    if missing:
        raise SkipTest(
            "Split trajectory comparison outputs are missing: {}".format(
                ", ".join(missing)))

    gmx_windows = _window_ids(GMX_DIR)
    mda_windows = _window_ids(MDA_DIR)
    if not gmx_windows or not mda_windows:
        raise SkipTest("Split trajectory comparison outputs are empty.")

    return gmx_windows, mda_windows


def _frame_times(topology, trajectory):
    u = mda.Universe(topology, str(trajectory))
    return np.array([float(ts.time) for ts in u.trajectory], dtype=float)


def test_window_ids_returns_sorted_numeric_directories_only():
    with TemporaryDirectory() as tmpdir:
        base = Path(tmpdir)
        for name in ("10", "notes", "2", "1"):
            (base / name).mkdir()
        (base / "3.txt").write_text("")

        assert _window_ids(base) == [1, 2, 10]


def test_splittraj_gmx_mda_all_frame_times_match():
    gmx_windows, mda_windows = _require_split_outputs()

    assert gmx_windows == mda_windows, (
        f"Window folders differ: GMX={gmx_windows}, MDA={mda_windows}"
    )

    for i in gmx_windows:
        gmx_xtc = GMX_DIR / str(i) / "aa.xtc"
        mda_xtc = MDA_DIR / str(i) / "aa.xtc"

        assert gmx_xtc.exists(), f"Missing GMX output: {gmx_xtc}"
        assert mda_xtc.exists(), f"Missing MDA output: {mda_xtc}"

        gmx_times = _frame_times(TOP, gmx_xtc)
        mda_times = _frame_times(TOP, mda_xtc)

        assert len(gmx_times) == len(mda_times), (
            f"Frame count mismatch in window {i}: "
            f"GMX={len(gmx_times)}, MDA={len(mda_times)}"
        )

        np.testing.assert_allclose(
            gmx_times,
            mda_times,
            rtol=0,
            atol=1e-6,
            err_msg=f"Frame-time mismatch in window {i}",
        )
