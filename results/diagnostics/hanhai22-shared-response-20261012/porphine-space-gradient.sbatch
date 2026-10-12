#!/usr/bin/env bash
#SBATCH --job-name=porphine-spin-gradient
#SBATCH --partition=CPU-192C768GB
#SBATCH --qos=qos_cpu_192c768gb
#SBATCH --account=xiatao
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=192G
#SBATCH --time=01:00:00
#SBATCH --output=porphine-space-gradient-%j.out
#SBATCH --error=porphine-space-gradient-%j.err
set -eo pipefail
[[ -n ${SLURM_JOB_ID:-} ]] || exit 64
. /etc/profile
module purge
module load gcc/14.1.0 cmake/3.30.5
root="$SLURM_SUBMIT_DIR"
source_dir=/home/guqqgroup/xiatao/work/xmvb-shared-response-1eb1b50a
artifact="$root/porphine-space-gradient-$SLURM_JOB_ID"
mkdir -p "$artifact"
export LC_ALL=C OMP_NUM_THREADS=32 OMP_PLACES=cores OMP_PROC_BIND=close OMP_MAX_ACTIVE_LEVELS=1
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1
export LD_LIBRARY_PATH="/home/guqqgroup/xiatao/software/xmvb-deps-20261012/lib:${LD_LIBRARY_PATH:-}"
cmake --build "$source_dir/build" --target check_structure_spin audit_orbital_gradient -j8 > "$artifact/build.log" 2>&1
srun --ntasks=1 --cpus-per-task=32 --cpu-bind=cores "$source_dir/build/src/check_structure_spin" \
 "$root/inputs/singlet/input.xmi" "$root/inputs/triplet/input.xmi" > "$artifact/spin.csv" 2> "$artifact/spin.err"
cd "$artifact"
cp "$root/inputs/triplet/cc-pvdz.gbs" .
python3 "$root/molden_endpoint.py" "$root/inputs/triplet/input.xmi" \
 "$root/porphine-ciah-1247745/triplet/input.molden" endpoint.xmi > conversion.txt
sha256sum "$source_dir/build/src/audit_orbital_gradient" endpoint.xmi > hashes.txt
hostname > node.txt
# Two symmetric displacements at each of the existing diagnostic probe sizes.
# This is a directional derivative check, not another fixed-point energy
# comparison against PySCF and not a continuation optimization.
srun --ntasks=1 --cpus-per-task=32 --cpu-bind=cores /usr/bin/time -v -o resource.time \
 "$source_dir/build/src/audit_orbital_gradient" endpoint.xmi davidson 1e-7 1e-3 gradient \
 > gradient.out 2> gradient.err
