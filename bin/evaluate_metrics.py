#!/usr/bin/env python3

"""
Stage 7 Validation: evaluate_metrics.py

Calculates Precision, Recall, and F1-score per confidence tier by comparing
the pipeline's tier_resolution.py output against the closed-genome ground truth.

Logic:
Since the ground-truth consists of closed genomes (chromosome + resolved plasmids):
- If a predicted ARG is listed in the PLSDB ground-truth for that genome -> True = Plasmid
- If a predicted ARG is NOT listed in the PLSDB ground-truth -> True = Chromosome

Metrics are computed independently for each tier to ensure no tier hides systemic issues.
"""

import argparse
import pandas as pd
import numpy as np
import sys
import re

__version__ = "0.1.0"

def clean_gene_name(gene: str) -> str:
    """Normalize gene names for safer matching (e.g. 'blaZ' vs 'blaZ-1')."""
    if pd.isna(gene):
        return ""
    # Lowercase and remove common suffixes that might mismatch
    return str(gene).strip().lower()

def main():
    parser = argparse.ArgumentParser(
        prog="evaluate_metrics",
        description="Calculate Precision/Recall/F1 per tier against ground truth.",
    )
    parser.add_argument("--predictions", required=True, help="Output TSV from tier_resolution.py")
    parser.add_argument("--ground-truth", required=True, help="Curated ground-truth CSV (e.g., table_summary.csv)")
    parser.add_argument("--genome-id", required=False, help="Genome Accession (if processing a single genome prediction file)")
    parser.add_argument("--output", required=True, help="Output metrics TSV path")
    args = parser.parse_args()

    # 1. Load Ground Truth
    try:
        gt_df = pd.read_csv(args.ground_truth, dtype=str)
    except Exception as e:
        sys.exit(f"Error reading ground truth: {e}")

    # Build a dictionary of true plasmid ARGs per genome
    # Format: { 'GCF_...': {'blaz', 'meca', ...} }
    true_plasmid_args = {}
    for _, row in gt_df.iterrows():
        genome = row.get("Assembly Accession", "")
        if pd.isna(genome) or not genome:
            continue
        
        if genome not in true_plasmid_args:
            true_plasmid_args[genome] = set()
            
        genes_str = row.get("Documented AMR Content", "")
        if pd.notna(genes_str) and genes_str.strip() != "None":
            # Split by comma and clean
            genes = [clean_gene_name(g) for g in genes_str.split(",")]
            true_plasmid_args[genome].update(genes)

    # 2. Load Predictions
    try:
        pred_df = pd.read_csv(args.predictions, sep="\t", dtype=str)
    except Exception as e:
        sys.exit(f"Error reading predictions: {e}")

    if "tier" not in pred_df.columns or "gene_symbol" not in pred_df.columns:
        sys.exit("Error: Predictions file must contain 'tier' and 'gene_symbol' columns.")

    # 3. Evaluate Predictions
    results = []
    for _, row in pred_df.iterrows():
        gene = clean_gene_name(row["gene_symbol"])
        tier = row["tier"]
        
        # Determine genome_id for this row
        genome_id = args.genome_id
        if not genome_id:
            if "genome_id" in row:
                genome_id = row["genome_id"]
            elif "Assembly Accession" in row:
                genome_id = row["Assembly Accession"]
            else:
                # Fallback: try to extract GCF_... from input_sequence_id
                seq_id = str(row.get("input_sequence_id", ""))
                match = re.search(r'(GCF_\d+\.\d+)', seq_id)
                if match:
                    genome_id = match.group(1)
                else:
                    sys.exit("Error: No --genome-id provided and could not find genome ID in columns.")
        
        # True location logic
        genome_plasmid_genes = true_plasmid_args.get(genome_id, set())
        
        # Substring/fuzzy match since PLSDB might say 'blaZ' and tool might say 'blaZ-1'
        is_plasmid = any(gene in pt or pt in gene for pt in genome_plasmid_genes if pt)
        true_location = "Plasmid" if is_plasmid else "Chromosome"
        
        results.append({
            "genome_id": genome_id,
            "gene_symbol": row["gene_symbol"],
            "predicted_tier": tier,
            "true_location": true_location
        })

    eval_df = pd.DataFrame(results)
    
    if eval_df.empty:
        sys.exit("No predictions to evaluate.")

    # 4. Calculate Metrics per Tier
    tiers = [
        "High-confidence plasmid",
        "Moderate-confidence plasmid",
        "Chromosomal",
        "Ambiguous"
    ]
    
    metrics = []
    total_true_plasmids = len(eval_df[eval_df["true_location"] == "Plasmid"])
    total_true_chromosomes = len(eval_df[eval_df["true_location"] == "Chromosome"])
    
    for t in tiers:
        tier_preds = eval_df[eval_df["predicted_tier"] == t]
        
        if t in ["High-confidence plasmid", "Moderate-confidence plasmid"]:
            target_true = "Plasmid"
            total_target = total_true_plasmids
        elif t == "Chromosomal":
            target_true = "Chromosome"
            total_target = total_true_chromosomes
        else: # Ambiguous
            # Ambiguous is a rejection class; it has no 'true' equivalent.
            # We just report how many were routed here.
            metrics.append({
                "Tier": t,
                "TP": np.nan, "FP": np.nan, "FN": np.nan,
                "Precision": np.nan, "Recall": np.nan, "F1_Score": np.nan,
                "Total_Calls": len(tier_preds)
            })
            continue
            
        tp = len(tier_preds[tier_preds["true_location"] == target_true])
        fp = len(tier_preds[tier_preds["true_location"] != target_true])
        
        # FN for a specific tier = True genes of that target type that were NOT predicted as this tier
        fn = total_target - tp
        
        precision = tp / (tp + fp) if (tp + fp) > 0 else 0.0
        recall = tp / (tp + fn) if (tp + fn) > 0 else 0.0
        f1 = 2 * (precision * recall) / (precision + recall) if (precision + recall) > 0 else 0.0
        
        metrics.append({
            "Tier": t,
            "TP": tp,
            "FP": fp,
            "FN": fn,
            "Precision": round(precision, 4),
            "Recall": round(recall, 4),
            "F1_Score": round(f1, 4),
            "Total_Calls": len(tier_preds)
        })
        
    metrics_df = pd.DataFrame(metrics)
    
    # Save results
    metrics_df.to_csv(args.output, sep="\t", index=False)
    print(f"Metrics successfully computed and saved to {args.output}")

if __name__ == "__main__":
    main()