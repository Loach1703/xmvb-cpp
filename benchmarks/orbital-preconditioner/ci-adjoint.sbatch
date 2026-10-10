#!/usr/bin/env bash
#SBATCH --job-name=ci-adjoint
#SBATCH --partition=CPU-256C768GB
#SBATCH --qos=qos_cpu-256c768gb
#SBATCH --account=taoxia
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=16G
#SBATCH --time=15-00:00:00
#SBATCH --output=ci-adjoint-%j.out
#SBATCH --error=ci-adjoint-%j.err
set -eo pipefail
[[ -n ${SLURM_JOB_ID:-} ]] || exit 64
. /etc/profile
module purge
module load anaconda3/25.06
set -u
source_dir="$1"
reference_sigma="$2"
dependencies=/home/guqqgroup/taoxia/software/conda-envs/xmvb-cpp-89c74a6
export PATH="$dependencies/bin:$PATH" CONDA_PREFIX="$dependencies"
export LD_LIBRARY_PATH="$dependencies/lib:${LD_LIBRARY_PATH:-}"
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1 GOTO_NUM_THREADS=1
export OMP_PLACES=cores OMP_PROC_BIND=close OMP_MAX_ACTIVE_LEVELS=1 LC_ALL=C
unset OMP_WAIT_POLICY GOMP_SPINCOUNT
root="$SLURM_SUBMIT_DIR/ci-adjoint-$SLURM_JOB_ID"
mkdir -p "$root"
hostname > "$root/node.txt"
lscpu > "$root/cpu.txt"
flags=(-O3 -DNDEBUG -g1 -fno-omit-frame-pointer -std=gnu++17 -fopenmp
       -DEIGEN_USE_BLAS=1 -DEIGEN_USE_LAPACKE=1
       -I"$source_dir/src" -I"$dependencies/include/eigen3")
libraries=("$source_dir/build/src/libxmvb-vbscf.a" -L"$dependencies/lib"
           -lcint -lopenblas -lm -Wl,--allow-shlib-undefined)
"$dependencies/bin/c++" "${flags[@]}" \
    "$source_dir/benchmarks/orbital-preconditioner/ci-adjoint.cpp" \
    "${libraries[@]}" -o "$root/current.exe" > "$root/build.log" 2>&1
# Rename only the reference class; do not introduce a production legacy path.
"$dependencies/bin/c++" "${flags[@]}" -DDirectCiSigmaAction=BaselineDirectCiSigmaAction \
    "$reference_sigma" "$source_dir/benchmarks/orbital-preconditioner/ci-adjoint.cpp" \
    "${libraries[@]}" -o "$root/baseline.exe" >> "$root/build.log" 2>&1
sha256sum "$reference_sigma" "$root/current.exe" "$root/baseline.exe" > "$root/hashes.txt"
for threads in 8 16 32; do
    for variant in current baseline; do
        mkdir -p "$root/$variant-$threads"
        cd "$root/$variant-$threads"
        srun --exact --cpus-per-task=32 --cpu-bind=cores \
            env OMP_NUM_THREADS="$threads" \
            python3 "$source_dir/benchmarks/orbital-preconditioner/affinity_probe.py" \
            --threads "$threads" --compact --record affinity.json -- \
            /usr/bin/time -v -o resource.time "$root/$variant.exe" > output.out 2> run.err
    done
done
