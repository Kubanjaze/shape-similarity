# shape-similarity — Phase 33

3D shape similarity between compounds using the Ultrafast Shape Recognition (USR) algorithm.
Identifies shape-similar compounds across different scaffold families (scaffold hops).

## Usage

```bash
PYTHONUTF8=1 python main.py --input data/compounds.csv
```

## Outputs

| File | Description |
|---|---|
| `output/usr_similarity_matrix.csv` | N×N pairwise USR scores |
| `output/usr_heatmap.png` | Clustered heatmap sorted by family |
| `output/shape_hop_candidates.csv` | Pairs with USR≥0.7 and different family |

## Algorithm

USR computes 12 descriptors (3 moments × 4 reference points) from 3D atom coordinates.
Similarity = 1 / (1 + mean absolute difference of descriptors).
Score: 0 (different shape) → 1.0 (identical shape).
