#!/usr/bin/env bash
#SBATCH --job-name=ao-fock-entry
#SBATCH --partition=CPU-192C768GB
#SBATCH --qos=qos_cpu_192c768gb
#SBATCH --account=xiatao
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=16G
#SBATCH --time=00:30:00
#SBATCH --output=ao-fock-entry-%j.out
#SBATCH --error=ao-fock-entry-%j.err
set -eo pipefail
. /etc/profile
module purge
module load gcc/14.1.0 cmake/3.30.5
root="$SLURM_SUBMIT_DIR"
source_dir=/home/guqqgroup/xiatao/work/xmvb-shared-response-1eb1b50a
run="$root/ao-fock-entry-$SLURM_JOB_ID"
mkdir -p "$run"
cd "$run"
cmake -S "$source_dir" -B "$source_dir/build" > configure.log 2>&1
cmake --build "$source_dir/build" --target audit_ao_fock_entry -j16 > build.log 2>&1
cp "$root/porphine-core-fd-1253527/endpoint.xmi" .
cp "$root/porphine-core-fd-1253527/cc-pvdz.gbs" .
cp "$root/porphine-core-fd-1253527/audit_orbital_gradient_inactive_density_f64.bin" density.bin
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export LD_LIBRARY_PATH="/home/guqqgroup/xiatao/software/xmvb-deps-20261012/lib:${LD_LIBRARY_PATH:-}"
binary="$source_dir/build/src/audit_ao_fock_entry"
sha256sum "$binary" endpoint.xmi density.bin > hashes.txt
hostname > node.txt
srun --exact --ntasks=1 --cpus-per-task=1 --cpu-bind=cores \
  /usr/bin/time -v -o resource.time "$binary" endpoint.xmi density.bin 185 454 > audit.out 2> audit.err
