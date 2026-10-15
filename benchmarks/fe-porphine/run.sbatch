#!/usr/bin/env bash
#SBATCH --job-name=fe-porphine-prepare
#SBATCH --partition=CPU-256C768GB
#SBATCH --qos=qos_cpu-256c768gb
#SBATCH --account=taoxia
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=192G
#SBATCH --time=04:00:00
#SBATCH --output=slurm-%j.out
#SBATCH --error=slurm-%j.err

set -eo pipefail
[[ -n ${SLURM_JOB_ID:-} ]] || exit 64
. /etc/profile
module purge
module load anaconda3/25.06
set -u
source_dir="$1"
dependencies=/home/guqqgroup/taoxia/software/conda-envs/xmvb-cpp-89c74a6
for command in srun sha256sum timeout mkdir cp date; do
    command -v "$command" >/dev/null || { printf 'Missing required command: %s\n' "$command" >&2; exit 69; }
done
for executable in /usr/bin/time /home/guqqgroup/taoxia/conda-envs/pyscf/bin/python \
    "$dependencies/bin/cc" "$dependencies/bin/c++" "$dependencies/bin/gfortran" \
    "$dependencies/bin/cmake" "$dependencies/bin/ninja"; do
    [[ -x "$executable" ]] || { printf 'Missing required executable: %s\n' "$executable" >&2; exit 69; }
done
run_dir="$SLURM_SUBMIT_DIR/job-$SLURM_JOB_ID"
mkdir -p "$run_dir"
cp "$SLURM_SUBMIT_DIR/prepare.py" "$SLURM_SUBMIT_DIR/geometry.xyz" \
    "$SLURM_SUBMIT_DIR/source-revision.txt" "$SLURM_SUBMIT_DIR/benchmark.sha256" "$run_dir/"
cd "$run_dir"
export LC_ALL=C OMP_NUM_THREADS="$SLURM_CPUS_PER_TASK"
export OMP_PLACES=cores OMP_PROC_BIND=close OMP_MAX_ACTIVE_LEVELS=1
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1 GOTO_NUM_THREADS=1
export CCACHE_DIR=/home/guqqgroup/taoxia/.cache/ccache-xmvb-cpp
printf '%s\n' "$SLURMD_NODENAME" > node.txt
if command -v lscpu >/dev/null; then
    lscpu > cpu.txt
else
    printf 'lscpu unavailable; node recorded in node.txt\n' > cpu.txt
fi
module list > modules.txt 2>&1
for variable in OMP_NUM_THREADS OMP_PLACES OMP_PROC_BIND OMP_MAX_ACTIVE_LEVELS \
    OPENBLAS_NUM_THREADS MKL_NUM_THREADS BLIS_NUM_THREADS GOTO_NUM_THREADS \
    SLURM_JOB_ID SLURM_CPUS_PER_TASK SLURM_JOB_NODELIST; do
    printf '%s=%s\n' "$variable" "${!variable-}"
done > allocation.txt
printf '%s BUILD START\n' "$(date -Is)" >> progress.txt
(
    export PATH="$dependencies/bin:$PATH" CONDA_PREFIX="$dependencies"
    export LD_LIBRARY_PATH="$dependencies/lib:${LD_LIBRARY_PATH:-}"
    cmake -S "$source_dir" -B "$source_dir/build" -G Ninja \
        -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTING=OFF \
        -DCMAKE_PREFIX_PATH="$dependencies" \
        -DCMAKE_C_COMPILER="$dependencies/bin/cc" \
        -DCMAKE_CXX_COMPILER="$dependencies/bin/c++" \
        -DCMAKE_Fortran_COMPILER="$dependencies/bin/gfortran" \
        -DXMVB_ENABLE_EIGEN_BLAS_BACKEND=ON
    cmake --build "$source_dir/build" --target xmvb -j16
) > build.log 2>&1
printf '%s PREPARE START\n' "$(date -Is)" >> progress.txt
srun --cpu-bind=cores /usr/bin/time -v -o prepare.time \
    /home/guqqgroup/taoxia/conda-envs/pyscf/bin/python -u prepare.py "$source_dir" \
    > prepare.out 2> prepare.err
printf '%s PREPARE END\n' "$(date -Is)" >> progress.txt
sha256sum "$source_dir/build/src/xmvb-cpp.exe" > binary.sha256
date -Is > finished.txt
