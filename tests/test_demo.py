import json
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]


@pytest.mark.integration
def test_complete_demo_generates_valid_outputs(tmp_path):
    output_dir = tmp_path / "outputs"
    completed = subprocess.run(
        [
            sys.executable,
            str(REPOSITORY_ROOT / "examples" / "run_core_demo.py"),
            "--output",
            str(output_dir),
        ],
        cwd=REPOSITORY_ROOT,
        check=True,
        capture_output=True,
        text=True,
        timeout=60,
    )

    summary = json.loads((output_dir / "summary.json").read_text())
    assert summary["cells"] == 48
    assert summary["frames"] == 3000
    assert summary["active_samples"] == 7779
    assert summary["network_nodes"] == 47
    assert summary["network_edges"] == 163
    assert "Saved numerical results" in completed.stdout

    with np.load(output_dir / "demo_results.npz") as results:
        assert results["binary_activity"].shape == (48, 3000)
        assert results["pca_embedding"].shape == (300, 48)
        assert np.isfinite(results["dff"]).all()
        assert np.isfinite(results["correlation_matrix"]).all()

    for name in ("calcium_binarization", "population_geometry"):
        png_path = output_dir / f"{name}.png"
        svg_path = output_dir / f"{name}.svg"
        assert png_path.stat().st_size > 100_000
        assert svg_path.stat().st_size > 100_000
        assert all(line == line.rstrip() for line in svg_path.read_text().splitlines())
