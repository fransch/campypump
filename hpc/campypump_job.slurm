#!/bin/bash --login
# ---------------------------------------------------------------------------
# campypump on a Slurm cluster: type a folder of genome assemblies.
# Edit the lines marked <-- EDIT, then:   sbatch campypump_job.slurm
# ---------------------------------------------------------------------------
#SBATCH --job-name=campypump
#SBATCH --account=YOUR_ACCOUNT          # <-- EDIT: the account/allocation to charge
#SBATCH --partition=standard            # <-- EDIT if your cluster names it differently
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --time=04:00:00
#SBATCH --output=campypump_%j.log

PROJECT=/path/to/your/project           # <-- EDIT: holds assemblies/ and will hold results/
ASSEMBLIES="$PROJECT/assemblies"        # <-- EDIT: one FASTA per isolate
SPECIES=""                              # optional: "C. jejuni" or "C. coli", labels the output

cd "$PROJECT" || { echo "Project folder not found: $PROJECT"; exit 1; }
module load micromamba 2>/dev/null || true
micromamba activate campypump 2>/dev/null || conda activate campypump 2>/dev/null || true
command -v campypump >/dev/null || { echo "campypump not on PATH: activate its environment"; exit 1; }
mkdir -p results results/figures

echo "campypump $(campypump --version | awk '{print $2}')  $(date)"

campypump run "$ASSEMBLIES" -t "${SLURM_NTASKS:-8}" ${SPECIES:+--species "$SPECIES"} \
    -o results/campypump_results.tsv \
    --alleles results/cmeB_alleles.fasta \
    --segment-plots results/segment_plots

campypump figure results/campypump_results.tsv -o results/figures
campypump diversity results/cmeB_alleles.fasta \
    -o results/figures/cmeB_allele_diversity.svg --table results/figures/alleles.tsv

echo "finished $(date)"
