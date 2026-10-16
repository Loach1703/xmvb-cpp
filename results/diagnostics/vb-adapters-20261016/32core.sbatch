#!/usr/bin/env bash
#SBATCH --job-name=vb-adapter32
#SBATCH --partition=CPU-192C768GB
#SBATCH --qos=qos_cpu_192c768gb
#SBATCH --account=xiatao
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --threads-per-core=1
#SBATCH --mem=32G
#SBATCH --time=00:20:00
#SBATCH --array=0-1
#SBATCH --output=vb-adapter32-%A_%a.out
#SBATCH --error=vb-adapter32-%A_%a.err
set -eo pipefail
[[ -n ${SLURM_JOB_ID:-} && ${SLURM_CPUS_PER_TASK:-0} == 32 ]] || exit 64
. /etc/profile
module purge
module load gcc/14.1.0 cmake/3.30.5
set -u
source_dir=/home/guqqgroup/xiatao/work/xmvb-source-tiles-20261015
case "$SLURM_ARRAY_TASK_ID" in
  0) name=240; deck=240_vb.xmi; input_sha=61d1035fd059c52a7361d6ed8f49f00a125286c644c155127bf6f8ee865eba48 ;;
  1) name=241; deck=241_VBSCF.xmi; input_sha=5dbd9c72c4ab5b8dec9b1bda5ef8c8d5d64741cc181f25cc8604e2464f9f07c8 ;;
  *) exit 64 ;;
esac
artifact="$SLURM_SUBMIT_DIR/32core-$SLURM_ARRAY_JOB_ID/$name"
mkdir -p "$artifact"
cp "$source_dir/build/src/xmvb-cpp.exe" "$artifact/xmvb-cpp.exe"
cp "$source_dir/testdata/vbscf/$deck" "$artifact/input.xmi"
ln -s "$source_dir/basis" "$artifact/basis"
cd "$artifact"
printf '%s  xmvb-cpp.exe\n%s  input.xmi\n' \
  27966b5ea621c92392457a56ecd6e25d87167e3f31c56039333b64d963d70b57 \
  "$input_sha" | sha256sum --check > hash-check.txt
sha256sum input.xmi xmvb-cpp.exe basis/cc-pvdz.gbs > hashes.sha256
printf '%s\n' 48c145baea76357abb5105e38004f48dfdc8645b > revision.txt
printf '%s\n' "$SLURM_JOB_ID" > job.txt
hostname > node.txt
lscpu > cpu.txt
module list > modules.txt 2>&1
export LC_ALL=C OMP_NUM_THREADS=32 OMP_PLACES=cores OMP_PROC_BIND=close
export OMP_MAX_ACTIVE_LEVELS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1
export LD_LIBRARY_PATH="/home/guqqgroup/xiatao/software/xmvb-deps-20261012/lib:${LD_LIBRARY_PATH:-}"
unset XMVB_PHASE_TIMING XMVB_CIAH_TRACE
env | sort | sed -n '/^OMP_/p; /^OPENBLAS_/p; /^MKL_/p; /^BLIS_/p; /^SLURM_/p' > environment.txt
set +e
srun --exact --ntasks=1 --cpus-per-task=32 --cpu-bind=cores bash -c '
  grep Cpus_allowed_list /proc/self/status > affinity.txt
  exec /usr/bin/time -v -o resource.time \
    ./xmvb-cpp.exe input.xmi --optimizer-backend ciah --eigensolver davidson \
    --orbital-preconditioner hessian-diagonal --gradient-tolerance 1e-3 \
    --energy-tolerance 1e-7 --ciah-trace ciah.tsv
' > output.out 2> run.err
status=$?
printf '%s\n' "$status" > run.exit
set -e
[[ "$status" == 0 ]]
grep -q '^Termination reason *: ciah_projected_gradient_tolerance' output.out
