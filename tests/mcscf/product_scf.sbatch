#!/usr/bin/env bash
#SBATCH --job-name=product-scf
#SBATCH --partition=CPU-192C768GB
#SBATCH --qos=qos_cpu_192c768gb
#SBATCH --account=xiatao
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --sockets-per-node=1
#SBATCH --cores-per-socket=32
#SBATCH --threads-per-core=1
#SBATCH --mem=48G
#SBATCH --time=00:20:00
#SBATCH --output=product-scf-%j.out
#SBATCH --error=product-scf-%j.err
set -eo pipefail
[[ -n ${SLURM_JOB_ID:-} && -d ${1:-} ]] || exit 64
. /etc/profile
module purge
module load gcc/14.1.0 cmake/3.30.5
export LC_ALL=C OMP_NUM_THREADS="$SLURM_CPUS_PER_TASK" OMP_PLACES=cores OMP_PROC_BIND=close
export OMP_MAX_ACTIVE_LEVELS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1
export LD_LIBRARY_PATH="/home/guqqgroup/xiatao/software/xmvb-deps-20261012/lib:${LD_LIBRARY_PATH:-}"
source_dir="$SLURM_SUBMIT_DIR"
artifact="$source_dir/selected-transport-$SLURM_JOB_ID"
mkdir -p "$artifact"
# Reject an unsuitable allocation before measuring, not a slow numerical result.
srun --exact --ntasks=1 --cpus-per-task="$SLURM_CPUS_PER_TASK" --cpu-bind=cores \
 python3 -c 'import json, os, pathlib
cpus=sorted(os.sched_getaffinity(0))
sockets=sorted({int(pathlib.Path(f"/sys/devices/system/cpu/cpu{c}/topology/physical_package_id").read_text()) for c in cpus})
print(json.dumps(dict(cpus=cpus,sockets=sockets)))
raise SystemExit(0 if len(cpus)==int(os.environ["SLURM_CPUS_PER_TASK"]) and len(sockets)==1 else 77)' \
 > "$artifact/topology.json"
cmake --build "$source_dir/build" --target xmvb -j8 > "$artifact/build.log" 2>&1
cp "$source_dir/build/CMakeCache.txt" "$artifact/"
sha256sum "$source_dir/src/vbscf/structures/assembly/"{action.cpp,action.hpp,product_jet.cpp,column_tiles.hpp} \
 "$source_dir/src/mcscf/ci/"{sigma.cpp,sigma.hpp,product_index.hpp,product_support.hpp} > "$artifact/jet-sources.sha256"
bash "$source_dir/tests/mcscf/selected_transport.sbatch" --240 "$1"
# Native results above are separate from this bounded phase attribution run.
mkdir -p "$artifact/phases/full"
cp "$artifact/native/full/"{input.xmi,cc-pvdz.gbs} "$artifact/phases/full/"
cd "$artifact/phases/full"
srun --exact --ntasks=1 --cpus-per-task="$SLURM_CPUS_PER_TASK" --cpu-bind=cores \
 env XMVB_PHASE_TIMING=1 /usr/bin/time -v -o resource.time "$artifact/xmvb-cpp.exe" input.xmi \
 --optimizer-backend neo --eigensolver davidson --orbital-preconditioner hessian-diagonal \
 --gradient-tolerance 1e-3 --energy-tolerance 1e-7 --neo-trace neo.tsv > output.out 2> run.err
