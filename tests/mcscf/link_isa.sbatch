#!/usr/bin/env bash
#SBATCH --job-name=link-isa
#SBATCH --partition=CPU-192C768GB
#SBATCH --qos=qos_cpu_192c768gb
#SBATCH --account=xiatao
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --time=00:02:00
#SBATCH --output=link-isa-%j.out
#SBATCH --error=link-isa-%j.err
set -eo pipefail
[[ -n ${SLURM_JOB_ID:-} ]] || exit 64
. /etc/profile
module purge
module load gcc/14.1.0
cd "$SLURM_SUBMIT_DIR"
binary=build/src/xmvb-cpp.exe
hostname
sha256sum "$binary"
# Static inspection only: no new build, workload or compiler-flag change.
# Inspect the linked functions from the actual benchmark executable.
stem="link-isa-$SLURM_JOB_ID"
nm --defined-only --format=posix "$binary" | awk '/build_spin_connections/ {print $1}' > "$stem.symbols"
while IFS= read -r symbol; do
  objdump -d --disassemble="$symbol" "$binary"
done < "$stem.symbols" > "$stem.asm"
grep -q -E '^[[:space:]]*[[:xdigit:]]+:' "$stem.asm" || exit 65
if ! grep -n -E -B3 -A3 'popcnt|popcount|tzcnt|bsf|pdep|pext|cmovp|cmovnp|setp|setnp' "$stem.asm"; then
  printf 'No matching instructions; inspect the retained symbol/assembly files.\n'
fi
