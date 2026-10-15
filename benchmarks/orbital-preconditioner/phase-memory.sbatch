#!/usr/bin/env bash
#SBATCH --job-name=orbital-memory-report
#SBATCH --partition=CPU-256C768GB
#SBATCH --qos=qos_cpu-256c768gb
#SBATCH --account=taoxia
#SBATCH --array=0-1%45
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=00:10:00
#SBATCH --output=phase-memory-%A_%a.out
#SBATCH --error=phase-memory-%A_%a.err
set -eo pipefail
[[ -n ${SLURM_JOB_ID:-} ]] || exit 64
set -u
captures_root="$1"
cases=(241_vb-iscf8 FeCl2-iscf9)
case_name="${cases[SLURM_ARRAY_TASK_ID]}"
analyzer=/home/guqqgroup/taoxia/software/bin/heaptrack_print
export LC_ALL=C
for variant in baseline updated; do
  cd "$captures_root/$case_name/$variant/heap"
  [[ -f memory.zst ]] || exit 66
  srun --exact --cpu-bind=cores "$analyzer" memory.zst \
    --print-allocators=0 --print-peaks=0 --print-temporary=0 \
    --print-massif memory.massif --massif-detailed-freq=0 \
    --print-flamegraph memory-peak.stacks --flamegraph-cost-type peak \
    > memory-detail-report.txt 2> detail-analyzer.err
done
