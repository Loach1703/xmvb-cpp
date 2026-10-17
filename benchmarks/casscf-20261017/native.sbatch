#!/usr/bin/env bash
#SBATCH --job-name=casscf-xmi
#SBATCH --partition=CPU-256C768GB
#SBATCH --qos=qos_cpu-256c768gb
#SBATCH --account=taoxia
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --threads-per-core=1
#SBATCH --mem=96G
#SBATCH --time=01:00:00
#SBATCH --output=native-%j.out
#SBATCH --error=native-%j.err
set -eo pipefail
[[ -n ${SLURM_JOB_ID:-} && $# == 2 ]] || exit 64
. /etc/profile
module purge
module load gcc/14.3.0 anaconda3/25.06
set -u
source_dir=$(realpath -e "$1")
case_name=$2
case "$case_name" in
  F2) deck=F2_OEO_PYSCF_LOW.xmi ;;
  FeCl2) deck=FeCl2.xmi ;;
  *) exit 64 ;;
esac
bundle=/home/guqqgroup/taoxia/xmvb-runs/porphine-current-20261017-c491230b/bundle
export LD_LIBRARY_PATH="$bundle/lib:${LD_LIBRARY_PATH:-}"
export OMP_NUM_THREADS=32 OMP_PLACES=cores OMP_PROC_BIND=close
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1
run_dir="$SLURM_SUBMIT_DIR/native-${case_name}-${SLURM_JOB_ID}"
mkdir -p "$run_dir"
cd "$run_dir"
hostname > node.txt
lscpu > cpu.txt
cp "$source_dir/source.sha256" .
python3 - "$source_dir/testdata/vbscf/$deck" input.xmi <<'PY'
from pathlib import Path
import re
import sys
source=Path(sys.argv[1]).read_text()
deck,count=re.subn(r'(?im)^\s*vbscf\s*$', 'casscf', source)
if count != 1:
    raise ValueError('expected exactly one VBSCF method flag')
Path(sys.argv[2]).write_text(deck)
PY
sha256sum "$source_dir/build/src/xmvb-cpp.exe" input.xmi > hashes.txt
srun --exact --ntasks=1 --cpus-per-task=32 --cpu-bind=cores bash -c '
  grep Cpus_allowed_list /proc/self/status > affinity.txt
  exec /usr/bin/time -v -o resource.time "$1" input.xmi
' bash "$source_dir/build/src/xmvb-cpp.exe" > output.out 2> run.err
grep -q '^Converged: yes$' output.out
[[ -s input.orb ]]
! grep -Eq '[0-9][eE][+-][0-9]' input.orb
