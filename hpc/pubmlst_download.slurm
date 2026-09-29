#!/bin/bash --login
# ---------------------------------------------------------------------------
# Download assemblies from PubMLST by isolate id (one id per line in ids.txt).
# Public records need no account; records that are not publicly released need BIGSdb_downloader
# and a personal API key — see docs/hpc_quickstart.md.
# ---------------------------------------------------------------------------
#SBATCH --job-name=pubmlst_fetch
#SBATCH --account=YOUR_ACCOUNT          # <-- EDIT
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=06:00:00
#SBATCH --output=pubmlst_fetch_%j.log

PROJECT=/path/to/your/project           # <-- EDIT
IDS="$PROJECT/ids.txt"                  # <-- EDIT: one PubMLST isolate id per line

cd "$PROJECT" || exit 1
module load micromamba 2>/dev/null || true
micromamba activate campypump 2>/dev/null || conda activate campypump 2>/dev/null || true

# resumable: genomes already present are skipped, failures are listed in assemblies/failed_ids.txt
campypump fetch --ids "$IDS" --outdir assemblies/
