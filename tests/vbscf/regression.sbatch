#!/usr/bin/env bash
#SBATCH --job-name=vb-adapter
#SBATCH --partition=CPU-192C768GB
#SBATCH --qos=qos_cpu_192c768gb
#SBATCH --account=xiatao
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --threads-per-core=1
#SBATCH --mem=32G
#SBATCH --time=00:20:00
#SBATCH --output=vb-adapter-%j.out
#SBATCH --error=vb-adapter-%j.err
set -eo pipefail
[[ -n ${SLURM_JOB_ID:-} && -x ${1:-} ]] || exit 64
case ${2:-} in
  240) deck=240_vb.xmi ;;
  241) deck=241_VBSCF.xmi ;;
  FeCl2) deck=FeCl2.xmi ;;
  *) exit 64 ;;
esac
. /etc/profile
module purge
module load gcc/14.1.0 cmake/3.30.5
source_dir="$SLURM_SUBMIT_DIR"
artifact="$source_dir/vb-adapter-$SLURM_JOB_ID"
mkdir -p "$artifact"
cp "$1" "$artifact/xmvb-cpp.exe"
cp "$source_dir/testdata/vbscf/$deck" "$artifact/input.xmi"
ln -s "$source_dir/basis" "$artifact/basis"
export LC_ALL=C OMP_NUM_THREADS=8 OMP_PLACES=cores OMP_PROC_BIND=close
export OMP_MAX_ACTIVE_LEVELS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1
export LD_LIBRARY_PATH="/home/guqqgroup/xiatao/software/xmvb-deps-20261012/lib:${LD_LIBRARY_PATH:-}"
unset XMVB_PHASE_TIMING XMVB_CIAH_TRACE
hostname > "$artifact/node.txt"
printf '%s\n' "$2" > "$artifact/case.txt"
cd "$artifact"
sha256sum input.xmi xmvb-cpp.exe > hashes.sha256
srun --exact --ntasks=1 --cpus-per-task=8 --cpu-bind=cores \
  /usr/bin/time -v -o resource.time \
  ./xmvb-cpp.exe input.xmi --optimizer-backend ciah --eigensolver davidson \
  --orbital-preconditioner hessian-diagonal --gradient-tolerance 1e-3 \
  --energy-tolerance 1e-7 --ciah-trace ciah.tsv > output.out 2> run.err
grep -q '^Termination reason *: ciah_projected_gradient_tolerance' output.out
