#!/usr/bin/env bash
#SBATCH --job-name=molcas-build-22
#SBATCH --partition=CPU-192C768GB
#SBATCH --qos=qos_cpu_192c768gb
#SBATCH --account=xiatao
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=96G
#SBATCH --time=12:00:00
#SBATCH --output=molcas-build-22-%j.out
#SBATCH --error=molcas-build-22-%j.err
set -eo pipefail
. /etc/profile
module purge
module load gcc/14.1.0 cmake/3.30.5 mkl/2024.1 anaconda/3/22.05
software=/home/guqqgroup/xiatao/software
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export LD_LIBRARY_PATH="/gpfs/opt/gcc/14.1.0/lib64:$software/libxc-7/lib:/opt/Anaconda3/2022.05/lib:${LD_LIBRARY_PATH:-}"
cmake -S "$software/src/OpenMolcas-v26.06" -B "$software/build/OpenMolcas-v26.06" \
  -DCMAKE_BUILD_TYPE=Release -DCMAKE_C_COMPILER=gcc -DCMAKE_Fortran_COMPILER=gfortran \
  -DMPI=OFF -DGA=OFF -DOPENMP=ON -DHDF5=ON -DHDF5_ROOT=/opt/Anaconda3/2022.05 \
  -DLINALG=MKL -DMKLROOT="$MKLROOT" -DEXTERNAL_LIBXC="$software/libxc-7" \
  -DTOOLS=ON -DNECI=OFF -DEFPLIB=OFF -DMSYM=OFF -DWFA=OFF -DBUILD_TESTING=OFF
cmake --build "$software/build/OpenMolcas-v26.06" --parallel 32
test -x "$software/build/OpenMolcas-v26.06/bin/rasscf.exe"
test -x "$software/build/OpenMolcas-v26.06/bin/seward.exe"
test -x "$software/build/OpenMolcas-v26.06/bin/guessorb.exe"
ldd "$software/build/OpenMolcas-v26.06/bin/rasscf.exe"
