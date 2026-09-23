#!/bin/bash --login
# ---------------------------------------------------------------------------
# Get MLST (ST, clonal complex) + provenance for all Peru isolates from PubMLST,
# including private records (uses the BIGSdb_downloader token in ~/.bigsdb_tokens),
# then merge with the campypump results into one analysis-ready CSV.
# Submit with:   sbatch metadata_job.slurm
# Resumable: if it stops, submit again and it continues where it left off.
# ---------------------------------------------------------------------------
#SBATCH --job-name=pubmlst_meta
#SBATCH --account=YOUR_GROUP                  # <-- EDIT (run: va)
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=06:00:00
#SBATCH --output=metadata_%j.log

PROJECT=/xdisk/YOUR_PI_NETID/YOUR_NETID/campypump_project   # <-- EDIT

cd "$PROJECT" || { echo "Project folder not found: $PROJECT"; exit 1; }
module load micromamba
micromamba activate campypump
mkdir -p results

echo "Step 1: MLST + provenance from PubMLST  $(date)"
campypump metadata --ids peru_ccoli_ids.txt peru_cjejuni_ids.txt \
    --token-dir ~/.bigsdb_tokens -o results/pubmlst_metadata.tsv

echo "Step 2: merging with campypump results  $(date)"
campypump annotate results/campypump_all_peru.tsv \
    --metadata results/pubmlst_metadata.tsv \
    -o results/campypump_all_peru_annotated.csv

echo "Finished  $(date)"
