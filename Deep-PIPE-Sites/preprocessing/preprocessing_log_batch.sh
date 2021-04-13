#!/bin/bash
#SBATCH --account=def-jrgreen
#SBATCH --cpus-per-task=32
#SBATCH --time=00:40:00
#SBATCH --mem-per-cpu=4G 

RUNPATH=/home/williamm/scratch/final-wrap-up/PIPE/Deep-PIPE-Sites/preprocessing
cd $RUNPATH
source ~/SYSC4907-pt/bin/activate
python preprocessing_log.py
