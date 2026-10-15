#!/usr/bin/env bash
#SBATCH --job-name=porphine-molcas
#SBATCH --partition=CPU-192C768GB
#SBATCH --qos=qos_cpu_192c768gb
#SBATCH --account=xiatao
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=192G
#SBATCH --time=15-00:00:00
#SBATCH --output=porphine-molcas-%j.out
#SBATCH --error=porphine-molcas-%j.err
set -eo pipefail
. /etc/profile
module purge
module load gcc/14.1.0 mkl/2024.1 anaconda/3/22.05
root="$SLURM_SUBMIT_DIR"
artifact="$root/porphine-molcas-$SLURM_JOB_ID"
mkdir -p "$artifact"
cd "$artifact"
export PYTHONPATH="$root/porphine-reference-audit-1247977/python" PYTHONNOUSERSITE=1
export MOLCAS=/home/guqqgroup/xiatao/software/build/OpenMolcas-v26.06
export LD_LIBRARY_PATH="/gpfs/opt/gcc/14.1.0/lib64:/home/guqqgroup/xiatao/software/libxc-7/lib:/opt/Anaconda3/2022.05/lib:${LD_LIBRARY_PATH:-}"
export MOLCAS_NPROCS=1 MOLCAS_MEM=120000 MOLCAS_THREADS=32
export MOLCAS_PROJECT=setup
export MOLCAS_WORKDIR
MOLCAS_WORKDIR=$(mktemp -d /tmp/porphine-molcas.XXXXXX)
printf '%s\n' "$MOLCAS_WORKDIR" > scratch-path.txt
export OMP_NUM_THREADS=32 OMP_PROC_BIND=false OMP_MAX_ACTIVE_LEVELS=1
unset OMP_PLACES GOMP_CPU_AFFINITY
export MKL_NUM_THREADS=32 MKL_DYNAMIC=false OPENBLAS_NUM_THREADS=1
sha256sum "$MOLCAS/bin/rasscf.exe" "$root/porphine-molcas.py" "$root/porphine-cross-restart.py" > hashes.txt
srun --ntasks=1 --cpus-per-task=32 --cpu-bind=cores \
  python3 "$root/porphine-molcas.py" "$root" prepare > prepare.out 2> prepare.err
srun --ntasks=1 --cpus-per-task=32 --cpu-bind=cores \
  /usr/bin/time -v -o setup.time python3 "$MOLCAS/pymolcas" setup.input > setup.out 2> setup.err
srun --ntasks=1 --cpus-per-task=32 --cpu-bind=cores \
  python3 "$root/porphine-cross-restart.py" "$root" export-frame > export.out 2> export.err
srun --ntasks=1 --cpus-per-task=32 --cpu-bind=cores \
  python3 "$root/porphine-molcas.py" "$root" convert > conversion.out 2> conversion.err
sha256sum initial.RasOrb setup.input optimize.input >> hashes.txt
status=0
srun --ntasks=1 --cpus-per-task=32 --cpu-bind=cores \
  /usr/bin/time -v -o optimize.time python3 "$MOLCAS/pymolcas" optimize.input > optimize.out 2> optimize.err || status=$?
printf '%s\n' "$status" > run.exit
exit "$status"
