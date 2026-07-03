# Bundled Real Suite2p Test Data

## Provenance

- Source animal: `DON-014371`
- Source session: `20230411Rec`
- Recording: `calibration`
- Source directory: `/media/cat/8TB/donato/bmi/DON-014371/20230411Rec/calibration/suite2p/plane0`
- Source classified cells: 123
- Included cells: 48 highest-variance classified cells, restored to source order
- Included frames: 0–2,999
- Sampling rate: 30.947 Hz
- Duration: 96.9 seconds

The subset preserves the real Suite2p `F`, `Fneu`, `spks`, `iscell`, `stat`, and mean-image data required by the core demonstration. `ops.npy` contains only the fields needed by the demo: `Ly`, `Lx`, `fs`, `nframes`, and `meanImg`.

Original cell indices are stored in `source_cell_ids.npy` and repeated in `metadata.json`. No synthetic fluorescence values were introduced.

## Integrity

| File | SHA-256 |
|---|---|
| `F.npy` | `850a828ff6059178f39da0270b9f709ac29232871f438b2565651bb502234bf6` |
| `Fneu.npy` | `b97b2fff07448c3f91828ffea3ba9e2541a7f87371bcabf7226bf1d6d914e55f` |
| `iscell.npy` | `c307344615887ef5270b238a03fb73df83d25ae800d6c67ed8593c32f9ddcc70` |
| `ops.npy` | `a671846075ed8437b68e7c5cb727f11334d9b8d458e33cfe5b100285f8611892` |
| `source_cell_ids.npy` | `3d099d9c6bf16cb77a7f2b6c43d5a3e665b07404337521535b64515a6e27e105` |
| `spks.npy` | `2569d1019ae1668a550a98460221fa0460d374d3f10b8dc9860ba69c2062d58e` |
| `stat.npy` | `0b35ad34c64fd95d742ffa606e505ef4e7bfd0a381cc1312fe97211faf799451` |

## Run

From the repository root:

```bash
python examples/run_core_demo.py
```
