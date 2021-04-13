#!/bin/bash
#SBATCH --account=def-jrgreen
#SBATCH --cpus-per-task=32
#SBATCH --time=2:00:00
#SBATCH --mem-per-cpu=2G    # Total memory for all tasks

parallel --joblog parallel.log < ./blastall.cmds
