#!/bin/bash
# Sample batchscript to run a parallel python job on HPC using 10 CPU cores
 
#SBATCH --partition=bch-compute                                # queue to be used
#SBATCH --time=9-20:05:00                         # Running time (in hours-minutes-seconds)
#SBATCH --job-name=perturbationsim                        # Job name
#SBATCH --mail-type=BEGIN,END,FAIL              # send and email when the job begins, ends or fails
#SBATCH --mail-user=lakindu.pathirakankanamge@childrens.harvard.edu          # Email address to send the job status
#SBATCH --output=output_%j.txt                  # Name of the output file
#SBATCH --nodes=1                            # Number of compute nodes
#SBATCH --ntasks=10                         # Number of cpu cores on one node
#SBATCH --mem=900G

echo "Hello BCH user $(whoami)!"
echo "This is a message from E2 node $(hostname) sent at $(date)"



cd /temp_work/ch250798/Scenic

python NFYA1c.py 