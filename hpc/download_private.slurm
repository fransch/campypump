#!/bin/bash --login
# ---------------------------------------------------------------------------
# Download the PRIVATE (logged-in only) Peru genomes from PubMLST with
# BIGSdb_downloader, using the access token created once with --setup.
# Reads the pending lists written by the campypump job:
#   results/pending_ccoli_ids.txt   -> genomes/ccoli/
#   results/pending_cjejuni_ids.txt -> genomes/cjejuni/
# Submit with:   sbatch download_private.slurm
# Afterwards:    sbatch campypump_job.slurm   (types everything)
# ---------------------------------------------------------------------------
#SBATCH --job-name=pubmlst_private
#SBATCH --account=YOUR_GROUP                  # <-- EDIT (run: va)
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=08:00:00
#SBATCH --output=download_private_%j.log

PROJECT=/xdisk/YOUR_PI_NETID/YOUR_NETID/campypump_project   # <-- EDIT
TOKENS=$HOME/.bigsdb_tokens
URL=https://rest.pubmlst.org/db/pubmlst_campylobacter_isolates/isolates

cd "$PROJECT" || { echo "Project folder not found: $PROJECT"; exit 1; }
module load micromamba
micromamba activate campypump

download_list () {
    local list=$1 outdir=$2 ok=0 skip=0 fail=0 n=0 total
    mkdir -p "$outdir"
    total=$(grep -c '[0-9]' "$list")
    : > "$outdir/private_failed_ids.txt"
    while IFS= read -r id; do
        id=$(echo "$id" | tr -cd '0-9')          # keep digits only (also strips Windows line endings)
        [ -z "$id" ] && continue
        n=$((n + 1))
        if [ -s "$outdir/id-$id.fasta" ]; then   # already downloaded
            skip=$((skip + 1)); continue
        fi
        tmp="$outdir/id-$id.fasta.part"
        err=$(bigsdb-downloader --cron --key_name PubMLST --site PubMLST --token_dir "$TOKENS" \
              --url "$URL/$id/contigs_fasta" --output_file "$tmp" 2>&1)
        if [ -s "$tmp" ] && [ "$(head -c 1 "$tmp")" = ">" ]; then
            mv "$tmp" "$outdir/id-$id.fasta"
            ok=$((ok + 1))
        else
            fail=$((fail + 1))
            echo "$id" >> "$outdir/private_failed_ids.txt"
            echo "[$n/$total] id $id: FAILED  $(echo "$err" | tail -n 1) $(head -c 150 "$tmp" 2>/dev/null)"
            rm -f "$tmp"
        fi
        if [ $((n % 100)) -eq 0 ]; then
            echo "[$n/$total] $outdir: downloaded $ok, already present $skip, failed $fail  $(date +%H:%M)"
        fi
        sleep 1
    done < "$list"
    echo "DONE $outdir: downloaded $ok, already present $skip, failed $fail"
}

echo "Start $(date)"
download_list results/pending_ccoli_ids.txt   genomes/ccoli
download_list results/pending_cjejuni_ids.txt genomes/cjejuni
echo "Finished $(date)"
