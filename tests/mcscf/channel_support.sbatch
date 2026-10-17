#!/usr/bin/env bash
#SBATCH --job-name=channel-support
#SBATCH --partition=CPU-192C768GB
#SBATCH --qos=qos_cpu_192c768gb
#SBATCH --account=xiatao
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --output=channel-support-%j.out
#SBATCH --error=channel-support-%j.err
set -eo pipefail
[[ -n ${SLURM_JOB_ID:-} ]]
. /etc/profile
module purge
module load gcc/14.1.0 cmake/3.30.5
set -u
cd "$SLURM_SUBMIT_DIR"
artifact="$SLURM_SUBMIT_DIR/channel-support-$SLURM_JOB_ID"
mkdir -p "$artifact"
export LD_LIBRARY_PATH="/home/guqqgroup/xiatao/software/xmvb-deps-20261012/lib:${LD_LIBRARY_PATH:-}"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1
cmake -S . -B build > "$artifact/configure.log" 2>&1
cmake --build build --target audit_channel_support -j2 > "$artifact/build.log" 2>&1
hostname > "$artifact/node.txt"
srun --exact --ntasks=1 --cpus-per-task=2 --cpu-bind=cores \
 ctest --test-dir build --output-on-failure \
 -R '^(structure_tiles|representation_overlap|selected_metric|ci_transport|selected_ci_transport|coupled_structure|keyframe|orthogonal_direct_ci|ciah|integrals|channel_support)$' \
 > "$artifact/ctest.out" 2>&1
sha256sum build/src/audit_channel_support tests/mcscf/audit_channel_support.cpp \
 tests/mcscf/targets.cmake src/mcscf/ci/product_support.hpp \
 src/mcscf/ci/product_paths.hpp src/mcscf/ci/covalent_support.hpp > "$artifact/hashes.txt"
# Serial combinatorial census, not a production CI timing or a scaling test.
srun --exact --ntasks=1 --cpus-per-task=1 --cpu-bind=cores \
 /usr/bin/time -v -o "$artifact/resource.time" build/src/audit_channel_support \
 > "$artifact/census.txt" 2> "$artifact/run.err"
