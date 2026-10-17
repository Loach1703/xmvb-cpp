#!/usr/bin/env bash
#SBATCH --job-name=C12-spin-freeze
#SBATCH --partition=CPU-192C768GB
#SBATCH --qos=qos_cpu_192c768gb
#SBATCH --account=xiatao
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --output=C12-spin-freeze-%j.out
#SBATCH --error=C12-spin-freeze-%j.err
set -eo pipefail
[[ -n ${SLURM_JOB_ID:-} && $# == 3 ]] || exit 64
. /etc/profile
module purge
module load gcc/14.1.0
set -u
source_dir="$SLURM_SUBMIT_DIR"
original=$(realpath -e "$1")
qualified=$(realpath -e "$2")
bundle="$source_dir/c12-spaces-bundle-$SLURM_JOB_ID"
mkdir -p "$bundle"
grep -q '100% tests passed, 0 tests failed out of 11' "$qualified/ctest.out"
sha256sum --check "$qualified/sources.sha256" > "$bundle/source-check.txt"
(cd "$original" && sha256sum --check hashes.sha256) > "$bundle/input-bundle-check.txt"
# Reuse the exact five already measured inputs. Only the qualified current
# executable changes; no old executable or old benchmark is executed here.
cp -a "$original/cases" "$bundle/"
cp "$original/"{cases.json,original.xmi,cc-pvdz.gbs,c12_spaces.py,c12_spaces.sbatch,summarize_scf_step.py} "$bundle/"
cp "$qualified/"{CMakeCache.txt,sources.sha256,ctest.out} "$bundle/"
cp "$source_dir/build/src/xmvb-cpp.exe" "$bundle/"
printf '%s\n' "$qualified" > "$bundle/qualification.txt"
printf '%s\n' "$3" > "$bundle/revision.txt"
(cd "$bundle" && sha256sum xmvb-cpp.exe original.xmi cc-pvdz.gbs \
    c12_spaces.py c12_spaces.sbatch summarize_scf_step.py cases.json > hashes.sha256)
