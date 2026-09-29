#!/bin/bash --login
# ---------------------------------------------------------------------------
# Retrieve PubMLST metadata (ST, clonal complex, provenance) and merge it into the results.
# Only relevant when the genomes come from PubMLST and you have their isolate ids.
# ---------------------------------------------------------------------------
#SBATCH --job-name=pubmlst_meta
#SBATCH --account=YOUR_ACCOUNT          # <-- EDIT
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=06:00:00
#SBATCH --output=metadata_%j.log

PROJECT=/path/to/your/project           # <-- EDIT
IDS="$PROJECT/ids.txt"                  # <-- EDIT
TOKENS="$HOME/.bigsdb_tokens"           # only needed for records that are not publicly released

cd "$PROJECT" || exit 1
module load micromamba 2>/dev/null || true
micromamba activate campypump 2>/dev/null || conda activate campypump 2>/dev/null || true

campypump metadata --ids "$IDS" ${TOKENS:+--token-dir "$TOKENS"} \
    -o results/pubmlst_metadata.tsv

# species_ani.tsv (skani) and exclude_ids.txt are used only if present
campypump annotate results/campypump_results.tsv \
    --metadata results/pubmlst_metadata.tsv \
    $([ -s species_ani.tsv ] && echo "--ani species_ani.tsv") \
    $([ -s exclude_ids.txt ] && echo "--exclude exclude_ids.txt") \
    -o results/campypump_annotated.csv
