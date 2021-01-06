#!/usr/bin/env bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=10000
#SBATCH --partition=short
#SBATCH --job-name=DRUID.valgrind
#SBATCH --output=valgrind.out.%j
#SBATCH --mail-user=yh362@cornell.edu
#SBATCH --mail-type=ALL

host=$(hostname)
echo "On host $host"
date

if [ ! -d /fs/cbsubscb09/storage/yilei/simulate/chrom ]; then
  # need to mount cbsubscb09 storage
  /programs/bin/labutils/mount_server cbsubscb09 /storage
fi

# command to run test on SAMFAS dataset
#prefix="/fs/cbsubscb09/storage/yilei/simulate/SAMAFS"
#valgrind --track-origins=yes ./DRUID -i $prefix/safs.seg --bim $prefix/safs.bim --Ne $prefix/safs.ibdne-ped2.ne -o safs --max 10 

# command to run on my simulated dataset using ukb
prefix="/fs/cbsubscb09/storage/yilei/simulate/DRUID_cpp"
valgrind --track-origins=yes --leak-check=yes ./DRUID -i $prefix/ped.seg --bim $prefix/ped.bim -o test --max 10 
