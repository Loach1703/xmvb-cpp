#!/usr/bin/env bash
#SBATCH --job-name=shell-screening-fix
#SBATCH --partition=CPU-192C768GB
#SBATCH --qos=qos_cpu_192c768gb
#SBATCH --account=xiatao
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=192G
#SBATCH --time=01:00:00
#SBATCH --output=shell-screening-fix-%j.out
#SBATCH --error=shell-screening-fix-%j.err
set -eo pipefail
. /etc/profile
module purge
module load gcc/14.1.0 cmake/3.30.5
root="$SLURM_SUBMIT_DIR"
source_dir=/home/guqqgroup/xiatao/work/xmvb-shared-response-1eb1b50a
run="$root/shell-screening-fix-$SLURM_JOB_ID"
mkdir -p "$run/entry" "$run/gradient"
cd "$run"
hostname > node.txt
cmake -S "$source_dir" -B "$source_dir/build" > configure.log 2>&1
cmake --build "$source_dir/build" --target test_libcint_shell_screening test_libcint_ri \
  audit_ao_fock_entry audit_orbital_gradient xmvb -j16 > build.log 2>&1
export LD_LIBRARY_PATH="/home/guqqgroup/xiatao/software/xmvb-deps-20261012/lib:${LD_LIBRARY_PATH:-}"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1
srun --exact --ntasks=1 --cpus-per-task=1 --cpu-bind=cores \
  ctest --test-dir "$source_dir/build" -R '^libcint_(shell_screening|ri)$' --output-on-failure > tests.out 2> tests.err
export OMP_NUM_THREADS=4
srun --exact --ntasks=1 --cpus-per-task=4 --cpu-bind=cores \
  "$source_dir/build/src/test_libcint_shell_screening" > tests-4core.out 2> tests-4core.err
for stage in entry gradient; do
  cp "$root/porphine-core-fd-1253527/endpoint.xmi" "$run/$stage/"
  cp "$root/porphine-core-fd-1253527/cc-pvdz.gbs" "$run/$stage/"
done
cd "$run/entry"
export OMP_NUM_THREADS=1
srun --exact --ntasks=1 --cpus-per-task=1 --cpu-bind=cores \
  "$source_dir/build/src/audit_ao_fock_entry" endpoint.xmi \
  "$root/porphine-core-fd-1253527/audit_orbital_gradient_inactive_density_f64.bin" \
  185 454 > audit.out 2> audit.err
cd "$run/gradient"
export OMP_NUM_THREADS=32 OMP_PLACES=cores OMP_PROC_BIND=close OMP_MAX_ACTIVE_LEVELS=1
sha256sum "$source_dir/build/src/xmvb-cpp.exe" "$source_dir/build/src/audit_orbital_gradient" \
  endpoint.xmi > hashes.txt
srun --exact --ntasks=1 --cpus-per-task=32 --cpu-bind=cores \
  /usr/bin/time -v -o resource.time "$source_dir/build/src/audit_orbital_gradient" \
  endpoint.xmi davidson 1e-7 1e-3 gradient coupled-keyframe > audit.out 2> audit.err
