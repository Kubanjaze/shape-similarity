import sys
if sys.platform == "win32":
    import io
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")
    sys.stderr = io.TextIOWrapper(sys.stderr.buffer, encoding="utf-8", errors="replace")

import argparse, os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem
RDLogger.DisableLog("rdApp.*")

FAMILY_COLORS = {"benz": "#4C72B0", "naph": "#DD8452", "ind": "#55A868",
                 "quin": "#C44E52", "pyr": "#8172B2", "bzim": "#937860", "other": "#808080"}

def load_compounds(path):
    df = pd.read_csv(path)
    records, n_bad = [], 0
    for _, row in df.iterrows():
        mol = Chem.MolFromSmiles(str(row["smiles"]))
        if mol is None:
            n_bad += 1
            continue
        fam = str(row["compound_name"]).split("_")[0]
        records.append({
            "compound_name": str(row["compound_name"]),
            "family": fam if fam in FAMILY_COLORS else "other",
            "pic50": float(row.get("pic50", float("nan"))),
            "mol": mol,
        })
    print(f"  {len(records)} valid ({n_bad} skipped)")
    return pd.DataFrame(records)

def generate_conformers(df):
    records = []
    params = AllChem.ETKDGv3()
    params.randomSeed = 42
    n_fail = 0
    for _, row in df.iterrows():
        mol = Chem.AddHs(row["mol"])
        cid = AllChem.EmbedMolecule(mol, params)
        if cid < 0:
            n_fail += 1
            mol = Chem.RemoveHs(mol)
            records.append({**row.to_dict(), "mol_3d": None})
            continue
        AllChem.MMFFOptimizeMolecule(mol, maxIters=200)
        mol = Chem.RemoveHs(mol)
        records.append({**row.to_dict(), "mol_3d": mol})
    if n_fail:
        print(f"  WARNING: {n_fail} conformer generation failures (excluded from USR)")
    return pd.DataFrame(records)

def compute_usr_matrix(df):
    valid = df[df["mol_3d"].notna()].reset_index(drop=True)
    n = len(valid)
    usr_desc = []
    for _, row in valid.iterrows():
        desc = AllChem.GetUSR(row["mol_3d"])
        usr_desc.append(desc)

    scores = np.zeros((n, n))
    for i in range(n):
        for j in range(i, n):
            s = AllChem.GetUSRScore(usr_desc[i], usr_desc[j])
            scores[i, j] = s
            scores[j, i] = s

    names = valid["compound_name"].tolist()
    return pd.DataFrame(scores, index=names, columns=names), valid

def plot_heatmap(sim_matrix, valid_df, output_path):
    order = valid_df.sort_values(["family", "compound_name"])["compound_name"].tolist()
    mat = sim_matrix.loc[order, order]
    families = valid_df.set_index("compound_name")["family"]

    fig = plt.figure(figsize=(12, max(10, len(order) * 0.25)))
    ax_h = fig.add_axes([0.16, 0.05, 0.75, 0.88])
    ax_c = fig.add_axes([0.04, 0.05, 0.06, 0.88])

    sns.heatmap(mat, ax=ax_h, cmap="YlOrRd", vmin=0, vmax=1,
                xticklabels=False, yticklabels=order,
                linewidths=0, cbar_kws={"shrink": 0.5})
    ax_h.tick_params(axis="y", labelsize=6)
    ax_h.set_title("USR Shape Similarity Matrix", fontsize=13, fontweight="bold")

    for i, name in enumerate(order):
        fam = families.get(name, "other")
        ax_c.barh(i, 1, color=FAMILY_COLORS.get(fam, "#808080"), edgecolor="none")
    ax_c.set_xlim(0, 1); ax_c.set_ylim(-0.5, len(order) - 0.5)
    ax_c.invert_yaxis(); ax_c.axis("off")

    plt.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close()

def find_shape_hops(sim_matrix, valid_df, threshold=0.7):
    families = valid_df.set_index("compound_name")["family"]
    rows = []
    names = sim_matrix.index.tolist()
    for i, q in enumerate(names):
        for j, h in enumerate(names):
            if i >= j:
                continue
            score = sim_matrix.loc[q, h]
            if score >= threshold and families.get(q) != families.get(h):
                rows.append({
                    "query": q, "hit": h,
                    "usr_score": round(float(score), 4),
                    "query_family": families.get(q, "other"),
                    "hit_family": families.get(h, "other"),
                })
    df = pd.DataFrame(rows)
    if len(df):
        df = df.sort_values("usr_score", ascending=False).reset_index(drop=True)
    return df

def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--input", required=True)
    parser.add_argument("--threshold", type=float, default=0.7,
                        help="USR score threshold for shape hop candidates")
    parser.add_argument("--output-dir", default="output")
    args = parser.parse_args()
    os.makedirs(args.output_dir, exist_ok=True)

    print(f"\nLoading: {args.input}")
    df = load_compounds(args.input)

    print("Generating 3D conformers (ETKDGv3 + MMFF)...")
    df3d = generate_conformers(df)

    print("Computing USR similarity matrix...")
    sim_matrix, valid_df = compute_usr_matrix(df3d)

    csv_path = os.path.join(args.output_dir, "usr_similarity_matrix.csv")
    sim_matrix.to_csv(csv_path)
    print(f"Saved: {csv_path}")

    plot_heatmap(sim_matrix, valid_df, os.path.join(args.output_dir, "usr_heatmap.png"))
    print(f"Saved: {args.output_dir}/usr_heatmap.png")

    hops = find_shape_hops(sim_matrix, valid_df, threshold=args.threshold)
    hops_path = os.path.join(args.output_dir, "shape_hop_candidates.csv")
    hops.to_csv(hops_path, index=False)
    print(f"Saved: {hops_path}")

    print(f"\n--- Shape hop candidates (USR >= {args.threshold}, different family) ---")
    print(f"  Found: {len(hops)} pairs")
    if len(hops):
        print(hops.head(10).to_string(index=False))

    print(f"\n--- Mean self-similarity per family ---")
    name_to_fam = valid_df.set_index("compound_name")["family"].to_dict()
    self_scores = []
    for name in sim_matrix.index:
        fam = name_to_fam.get(name, "other")
        same = [n for n in sim_matrix.index if name_to_fam.get(n) == fam and n != name]
        if same:
            mean_s = sim_matrix.loc[name, same].mean()
            self_scores.append({"compound_name": name, "family": fam, "mean_intra_usr": mean_s})
    if self_scores:
        ss_df = pd.DataFrame(self_scores)
        print(ss_df.groupby("family")["mean_intra_usr"].mean().round(3).to_string())
    print("\nDone.")

if __name__ == "__main__":
    main()
