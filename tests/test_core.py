import hashlib
from pathlib import Path

import networkx as nx
import numpy as np

from bmitools import Calcium
from bmitools.utils import compute_dff0
from bmitools.utils.network.network import (
    generate_graph_from_connected_nodes,
    get_degree_distribution,
    load_correlation_matrix,
)


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
PLANE0 = (
    REPOSITORY_ROOT
    / "test_data"
    / "DON-014371"
    / "20230411Rec"
    / "calibration"
    / "suite2p"
    / "plane0"
)


def test_bundled_data_integrity():
    expected = {
        "F.npy": "850a828ff6059178f39da0270b9f709ac29232871f438b2565651bb502234bf6",
        "Fneu.npy": "b97b2fff07448c3f91828ffea3ba9e2541a7f87371bcabf7226bf1d6d914e55f",
        "iscell.npy": "c307344615887ef5270b238a03fb73df83d25ae800d6c67ed8593c32f9ddcc70",
        "ops.npy": "a671846075ed8437b68e7c5cb727f11334d9b8d458e33cfe5b100285f8611892",
        "source_cell_ids.npy": "3d099d9c6bf16cb77a7f2b6c43d5a3e665b07404337521535b64515a6e27e105",
        "spks.npy": "2569d1019ae1668a550a98460221fa0460d374d3f10b8dc9860ba69c2062d58e",
        "stat.npy": "0b35ad34c64fd95d742ffa606e505ef4e7bfd0a381cc1312fe97211faf799451",
    }
    for name, checksum in expected.items():
        digest = hashlib.sha256((PLANE0 / name).read_bytes()).hexdigest()
        assert digest == checksum


def test_public_calcium_loader_reads_bundled_suite2p():
    calcium = Calcium()
    calcium.data_dir = str(PLANE0)
    calcium.load_suite2p()

    assert calcium.F.shape == (48, 3000)
    assert calcium.Fneu.shape == calcium.F.shape
    assert calcium.spks.shape == calcium.F.shape
    assert calcium.stat.shape == (48,)
    assert np.isfinite(calcium.F).all()
    assert np.all(calcium.iscell[:, 0] == 1)


def test_compute_dff0_uses_per_cell_median_baseline():
    time_by_cell = np.array(
        [
            [10.0, 20.0],
            [12.0, 18.0],
            [14.0, 22.0],
        ]
    )
    baseline, dff = compute_dff0(time_by_cell)

    np.testing.assert_allclose(baseline, [12.0, 20.0])
    np.testing.assert_allclose(dff[1], [0.0, -0.1])


def test_calcium_binarize_saves_binary_activity(tmp_path):
    calcium = Calcium()
    calcium.data_dir = str(tmp_path)
    traces = np.array(
        [
            [0.0, 0.0, 4.0, 0.0, 0.0],
            [0.0, 3.0, 0.0, 0.0, 0.0],
        ],
        dtype=np.float32,
    )

    binary, _ = calcium.binarize(traces, thresh=1.5)

    assert set(np.unique(binary)) == {0.0, 1.0}
    assert binary[0, 2] == 1
    assert binary[1, 1] == 1
    np.testing.assert_array_equal(np.load(tmp_path / "binarized.npy"), binary)


def test_network_thresholding_and_degree_distribution(tmp_path):
    pairs = np.array(
        [
            [0, 1, 0.45, 0.01],
            [1, 0, 0.45, 0.01],
            [1, 2, 0.20, 0.01],
            [2, 1, 0.20, 0.01],
            [0, 2, 0.80, 0.50],
            [2, 0, 0.80, 0.50],
        ],
        dtype=float,
    )
    pair_path = tmp_path / "pairs.npy"
    np.save(pair_path, pairs)

    adjacency = load_correlation_matrix(str(pair_path), corr_thresh=0.30)
    graph = nx.Graph(generate_graph_from_connected_nodes(adjacency))
    graph.remove_nodes_from(list(nx.isolates(graph)))
    degree_bins, degree_counts = get_degree_distribution(graph)

    assert set(graph.edges()) == {(0, 1)}
    assert degree_bins[0] == 0
    assert degree_counts[0] == 2
