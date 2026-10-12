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
# Match the successful Hanhai25 bundled dependency, but build it offline.
# Its Fortran interface uses baseline flags (not Molcas's global I8 flags).
cmake -S "$software/src/libwignernj-0.6.0" -B "$software/build/libwignernj-0.6.0" \
  -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX="$software/libwignernj-0.6.0" \
  -DCMAKE_INSTALL_LIBDIR=lib -DCMAKE_C_COMPILER=gcc -DCMAKE_Fortran_COMPILER=gfortran \
  -DCMAKE_POSITION_INDEPENDENT_CODE=ON -DBUILD_SHARED_LIBS=OFF -DBUILD_FORTRAN=ON \
  -DBUILD_TESTS=OFF -DBUILD_CXX_TESTS=OFF -DBUILD_EXAMPLES=OFF -DBUILD_PYTHON=OFF -DBUILD_LTO=OFF
cmake --build "$software/build/libwignernj-0.6.0" --parallel 32
cmake --install "$software/build/libwignernj-0.6.0"
cmake -S "$software/src/OpenMolcas-v26.06" -B "$software/build/OpenMolcas-v26.06" \
  -DCMAKE_BUILD_TYPE=Release -DCMAKE_C_COMPILER=gcc -DCMAKE_Fortran_COMPILER=gfortran \
  -DMPI=OFF -DGA=OFF -DOPENMP=ON -DHDF5=ON -DHDF5_ROOT=/opt/Anaconda3/2022.05 \
  -DLINALG=MKL -DMKLROOT="$MKLROOT" -DEXTERNAL_LIBXC="$software/libxc-7" \
  -DCMAKE_PREFIX_PATH="$software/libwignernj-0.6.0" \
  -DTOOLS=ON -DNECI=OFF -DEFPLIB=OFF -DMSYM=OFF -DWFA=OFF -DBUILD_TESTING=OFF
cmake --build "$software/build/OpenMolcas-v26.06" --parallel 32
test -x "$software/build/OpenMolcas-v26.06/bin/rasscf.exe"
test -x "$software/build/OpenMolcas-v26.06/bin/seward.exe"
test -x "$software/build/OpenMolcas-v26.06/bin/guessorb.exe"
ldd "$software/build/OpenMolcas-v26.06/bin/rasscf.exe"
