#!/bin/bash --login
# ---------------------------------------------------------------------------
# campypump batch job for UArizona HPC (Puma): Peru C. coli + C. jejuni
# Edit the two lines marked <-- EDIT before submitting.
# Submit with:   sbatch campypump_job.slurm
# ---------------------------------------------------------------------------
#SBATCH --job-name=campypump
#SBATCH --account=YOUR_GROUP                  # <-- EDIT (run: va)
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --time=10:00:00
#SBATCH --output=campypump_%j.log

PROJECT=/xdisk/YOUR_PI_NETID/YOUR_NETID/campypump_project   # <-- EDIT

cd "$PROJECT" || { echo "Project folder not found: $PROJECT"; exit 1; }
module load micromamba
micromamba activate campypump
mkdir -p genomes/ccoli genomes/cjejuni results

# If C. coli genomes were downloaded earlier straight into genomes/, move them into genomes/ccoli/
if ls genomes/id-*.fasta >/dev/null 2>&1; then
    echo "Moving previously downloaded genomes into genomes/ccoli/"
    mv genomes/id-*.fasta genomes/ccoli/
fi

echo "Step 1: downloading genomes from PubMLST  $(date)"
campypump fetch --ids peru_ccoli_ids.txt   --outdir genomes/ccoli/
campypump fetch --ids peru_cjejuni_ids.txt --outdir genomes/cjejuni/

# paste-ready lists of genomes that could not be downloaded (private on PubMLST)
cut -f1 genomes/ccoli/failed_ids.txt   > results/pending_ccoli_ids.txt
cut -f1 genomes/cjejuni/failed_ids.txt > results/pending_cjejuni_ids.txt
echo "Pending (not public): $(wc -l < results/pending_ccoli_ids.txt) C. coli, $(wc -l < results/pending_cjejuni_ids.txt) C. jejuni"

echo "Step 2: typing cmeB  $(date)"
campypump run genomes/ccoli/ -t "$SLURM_NTASKS" --species "C. coli" \
    -o results/campypump_ccoli.tsv --alleles results/cmeB_alleles_ccoli.fasta \
    --segment-plots results/segment_plots_ccoli
campypump run genomes/cjejuni/ -t "$SLURM_NTASKS" --species "C. jejuni" \
    -o results/campypump_cjejuni.tsv --alleles results/cmeB_alleles_cjejuni.fasta \
    --segment-plots results/segment_plots_cjejuni

echo "Step 3: combining both species into one table  $(date)"
campypump combine results/campypump_ccoli.tsv results/campypump_cjejuni.tsv \
    -o results/campypump_all_peru.tsv

# if PubMLST metadata has been retrieved (metadata_job.slurm), refresh the annotated table
if [ -s results/pubmlst_metadata.tsv ]; then
    echo "Step 4: merging with PubMLST metadata  $(date)"
    campypump annotate results/campypump_all_peru.tsv \
        --metadata results/pubmlst_metadata.tsv \
        -o results/campypump_all_peru_annotated.csv
fi

echo "Finished  $(date)"
