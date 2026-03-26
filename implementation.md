# Phase 33 — Shape Similarity (USR Algorithm)

**Version:** 1.1 | **Tier:** Standard | **Date:** 2026-03-26

## Goal
Compute 3D shape similarity between all compound pairs using the Ultrafast Shape Recognition (USR) algorithm.
Identify shape-similar compounds across different scaffold families — potential scaffold hops.

CLI: `python main.py --input data/compounds.csv --query benz_001_H`

Outputs: usr_similarity_matrix.csv, usr_heatmap.png, shape_hop_candidates.csv

## Algorithm — USR
- Generate 3D conformer: ETKDGv3 + MMFF optimization
- Compute 4 reference points: centroid, closest-to-centroid, farthest-from-centroid, farthest-from-farthest
- For each reference point: compute 3 moments (μ1, μ2, μ3) of distance distribution → 12 descriptors total
- USR similarity = 1 / (1 + (1/12) * Σ|m_i - m_j|)  for i,j in range(12)
- Score range: 0 (different shape) to 1.0 (identical shape)

## Logic
- Embed all compounds with ETKDGv3, seed=42
- MMFF minimize (maxIters=200)
- Compute USR descriptors via RDKit `AllChem.GetUSRScore` / `AllChem.GetUSR`
- Pairwise USR scores → N×N matrix
- Scaffold hop candidates: USR > 0.7 AND different scaffold family

## Outputs
- usr_similarity_matrix.csv: N×N pairwise USR scores
- usr_heatmap.png: clustered heatmap sorted by family, colored by USR score
- shape_hop_candidates.csv: pairs with USR>0.7 + different family (query_compound, hit_compound, usr_score, query_family, hit_family)

## Actual Results (v1.1)

| Family | Mean intra-family USR |
|---|---|
| benz | 0.806 |
| naph | 0.860 |
| ind | 0.732 |
| quin | 0.725 |
| pyr | 0.741 |
| bzim | 0.741 |

**228 cross-family shape hop candidates** at USR ≥ 0.7
**Top pair:** ind_001_H / bzim_001_H — USR=0.979 (nearly identical 3D shape despite different 2D scaffolds)
**Key insight:** ind and quin families are highly shape-similar to each other and to bzim — all bicyclic systems with similar dimensions. naph/benz show highest intra-family shape consistency (rigid planar scaffolds).
**45/45 valid** — all conformers generated successfully
