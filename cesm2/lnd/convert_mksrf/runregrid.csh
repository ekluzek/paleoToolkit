#!/bin/tcsh
#
# This is setup to run the regrid script on Caspar the Data Visualization system at NCAR
#
#SBATCH -J regridding
#SBATCH -n 1
#SBATCH -t 02:30:00
#SBATCH -A P93300606
#SBATCH -p dav
#SBATCH --mem=8G
#SBATCH -e %x.err.%J
#SBATCH -o %x.out.%J
#
setenv DIN_LOC_ROOT /gpfs/fs1/p/cesmdata/cseg/inputdata
set echo
module load ncl
module load nco
ncl regrid_GLCMEC_n_PFT.ncl
ncl regrid_GIA_glc.ncl
ncl lake_dataset_21ka.ncl
