#!/usr/bin/env python3
# Script to simulate MetaPhlAn4 and HUMAnN3 outputs for demonstration/testing
# Author: Dawid X (dawidx1233)

import argparse
import os
import random

def simulate_profile(sample_id, output_path):
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    
    # Simple simulated taxonomic profile (MetaPhlAn4 style)
    taxa = [
        ("k__Bacteria|p__Bacteroidetes|c__Bacteroidia|o__Bacteroidales|f__Bacteroidaceae|g__Bacteroides", 40.0 + random.uniform(-5, 5)),
        ("k__Bacteria|p__Firmicutes|c__Clostridia|o__Clostridiales|f__Lachnospiraceae|g__Coprococcus", 15.0 + random.uniform(-3, 3)),
        ("k__Bacteria|p__Firmicutes|c__Clostridia|o__Clostridiales|f__Ruminococcaceae|g__Faecalibacterium", 20.0 + random.uniform(-4, 4)),
        ("k__Bacteria|p__Proteobacteria|c__Gammaproteobacteria|o__Enterobacterales|f__Enterobacteriaceae|g__Escherichia", 5.0 + random.uniform(-2, 2)),
        ("k__Bacteria|p__Actinobacteria|c__Actinomycetia|o__Bifidobacteriales|f__Bifidobacteriaceae|g__Bifidobacterium", 10.0 + random.uniform(-3, 3)),
        ("k__Bacteria|p__Verrucomicrobia|c__Verrucomicrobiae|o__Verrucomicrobiales|f__Akkermansiaceae|g__Akkermansia", 10.0 + random.uniform(-2, 2)),
    ]
    
    # Normalize abundances to sum to 100
    total = sum(val for _, val in taxa)
    taxa = [(t, (val / total) * 100.0) for t, val in taxa]
    
    with open(output_path, 'w') as f:
        f.write("#clade_name\trelative_abundance\n")
        for taxon, abundance in taxa:
            f.write(f"{taxon}\t{abundance:.4f}\n")
    print(f"Simulated MetaPhlAn4 profile for {sample_id} written to {output_path}")

def simulate_pathways(sample_id, output_path):
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    
    pathways = [
        ("PWY-3781: aerobic respiration I (cytochrome c)", 1200.0 + random.uniform(-100, 100)),
        ("GLUCOSE1PMETAB-PWY: glucose-1-phosphate degradation", 800.0 + random.uniform(-80, 80)),
        ("PWY-5097: L-lysine biosynthesis II", 450.0 + random.uniform(-40, 40)),
        ("PWY-6163: chorismate biosynthesis from 3-dehydroquinate", 300.0 + random.uniform(-30, 30)),
        ("FA-SYNTHESIS-PWY: fatty acid biosynthesis", 1500.0 + random.uniform(-150, 150)),
    ]
    
    with open(output_path, 'w') as f:
        f.write("# Pathway\tAbundance\n")
        for pwy, val in pathways:
            f.write(f"{pwy}\t{val:.4f}\n")
    print(f"Simulated HUMAnN3 pathways for {sample_id} written to {output_path}")

def merge_profiles(input_profiles, output_path):
    import pandas as pd
    profiles = input_profiles.split()
    merged_data = {}
    
    for prof in profiles:
        sample_name = os.path.basename(prof).replace("_profile.txt", "")
        df = pd.read_csv(prof, sep="\t", comment="#", names=["clade_name", "relative_abundance"])
        for _, row in df.iterrows():
            clade = row["clade_name"]
            if clade not in merged_data:
                merged_data[clade] = {}
            merged_data[clade][sample_name] = row["relative_abundance"]
            
    df_merged = pd.DataFrame(merged_data).T.fillna(0.0)
    df_merged.index.name = "clade_name"
    df_merged.to_csv(output_path, sep="\t")
    print(f"Merged MetaPhlAn4 profiles written to {output_path}")

def merge_pathways(input_pathways, output_path):
    import pandas as pd
    pathways = input_pathways.split()
    merged_data = {}
    
    for pwy in pathways:
        sample_name = os.path.basename(pwy).replace("_pathabundance.tsv", "")
        df = pd.read_csv(pwy, sep="\t", comment="#", names=["Pathway", "Abundance"])
        for _, row in df.iterrows():
            p = row["Pathway"]
            if p not in merged_data:
                merged_data[p] = {}
            merged_data[p][sample_name] = row["Abundance"]
            
    df_merged = pd.DataFrame(merged_data).T.fillna(0.0)
    df_merged.index.name = "Pathway"
    df_merged.to_csv(output_path, sep="\t")
    print(f"Merged HUMAnN3 pathways written to {output_path}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Simulate MetaPhlAn4/HUMAnN3 outputs")
    parser.add_argument("--sample", type=str, help="Sample ID")
    parser.add_argument("--output_profile", type=str, help="Output taxonomic profile path")
    parser.add_argument("--output_pathway", type=str, help="Output pathway abundance path")
    parser.add_argument("--merge_profiles", action="store_true", help="Merge taxonomic profiles")
    parser.add_argument("--merge_pathways", action="store_true", help="Merge pathway abundances")
    parser.add_argument("--input_profiles", type=str, help="Space-separated list of taxonomic profiles")
    parser.add_argument("--input_pathways", type=str, help="Space-separated list of pathway abundances")
    parser.add_argument("--output_merged", type=str, help="Output path for merged table")
    
    args = parser.parse_args()
    
    if args.output_profile:
        simulate_profile(args.sample, args.output_profile)
    elif args.output_pathway:
        simulate_pathways(args.sample, args.output_pathway)
    elif args.merge_profiles:
        merge_profiles(args.input_profiles, args.output_merged)
    elif args.merge_pathways:
        merge_pathways(args.input_pathways, args.output_merged)
