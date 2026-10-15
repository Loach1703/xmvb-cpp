#!/usr/bin/env bash
#SBATCH --job-name=240-spaces
#SBATCH --partition=CPU-192C768GB
#SBATCH --qos=qos_cpu_192c768gb
#SBATCH --account=xiatao
#SBATCH --array=0-4
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time=15-00:00:00
#SBATCH --output=240-spaces-%A_%a.out
#SBATCH --error=240-spaces-%A_%a.err

# Run the nested CAS(8,8) structure spaces to convergence, not a fixed prefix.
# Argument: frozen bundle containing binary, original input and basis file.
set -eo pipefail
[[ -n ${SLURM_JOB_ID:-} && -d ${1:-} ]] || exit 64
. /etc/profile
module purge
module load gcc/14.1.0
bundle=$(realpath "$1")
classes=('COV' 'ION(0:1)' 'ION(0:2)' 'ION(0:3)' 'FULL')
labels=(cov ion01 ion02 ion03 full)
index=$SLURM_ARRAY_TASK_ID
artifact="$SLURM_SUBMIT_DIR/240-spaces-$SLURM_ARRAY_JOB_ID/${labels[index]}"
mkdir -p "$artifact"
cd "$artifact"
# Only the structure class and optimizer change; the original GUS is intact.
sed -E -e "s/^str=.*/str=${classes[index]}/I" -e 's/^iscf=.*/iscf=9/I' \
 "$bundle/240_vb.xmi" > input.xmi
cp "$bundle/cc-pvdz.gbs" .
export LD_LIBRARY_PATH="/home/guqqgroup/xiatao/software/xmvb-deps-20261012/lib:${LD_LIBRARY_PATH:-}"
export LC_ALL=C OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OMP_PLACES=cores OMP_PROC_BIND=close OMP_MAX_ACTIVE_LEVELS=1
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1
hostname > node.txt
lscpu > cpu.txt
module list > modules.txt 2>&1
env | sort | grep -E '^(OMP_|OPENBLAS_|MKL_|BLIS_|SLURM_|XMVB_)' > environment.txt
sha256sum "$bundle/xmvb-cpp.exe" input.xmi cc-pvdz.gbs > hashes.txt
command=("$bundle/xmvb-cpp.exe" input.xmi --optimizer-backend neo
 --eigensolver davidson --orbital-preconditioner hessian-diagonal
 --gradient-tolerance 1e-3 --energy-tolerance 1e-7 --neo-trace neo.tsv)
printf '%q ' "${command[@]}" > command.txt
printf '\n' >> command.txt
srun --exact --ntasks=1 --cpus-per-task="$SLURM_CPUS_PER_TASK" --cpu-bind=cores \
 bash -c 'grep Cpus_allowed_list /proc/self/status' > affinity.txt
status=0
srun --exact --ntasks=1 --cpus-per-task="$SLURM_CPUS_PER_TASK" --cpu-bind=cores \
 /usr/bin/time -v -o resource.time "${command[@]}" > output.out 2> run.err || status=$?
printf '%s\n' "$status" > run.exit
exit "$status"
