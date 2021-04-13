#!/bin/bash
# perfoms swissprot psiblast
infile="$1"
outfile="$(basename "${infile%.*}").pssm"
database="/home/williamm/scratch/final-wrap-up/PIPE/Deep-PIPE-Sites/preprocessing/pssm/swissprot/swissprot"
module load StdEnv/2020
module load gcc/9.3.0
module load blast+/2.11.0 

psiblast \
 -query "${infile}" \
 -db "${database}" \
 -inclusion_ethresh 0.001 \
 -num_iterations 2 \
 -out_ascii_pssm "../pssms/${outfile}" \
 -save_pssm_after_last_round
