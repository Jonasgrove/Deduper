#!/bin/bash
#SBATCH --partition=bgmp          ## Partition (like a queue in PBS)
#SBATCH --job-name=deduper         ## Job Name
#SBATCH --output=deduper.er         ## File in which to store job output
#SBATCH --error=deduper.out          ## File in which to store job error messages
#SBATCH --time=10-24:00:00         ## Wall clock time limit in Days-HH:MM:SS
#SBATCH --nodes=1                 ## Number of nodes needed for the job
#SBATCH --ntasks-per-node=1       ## Number of tasks to be launched per Node
#SBATCH --account=bgmp            ## Account used for job
#SBATCH --cpus-per-task=40
#SBATCH --mem=179G

data_base=/projects/bgmp/jonasg/bioinfo/bi624/Deduper/package/
input=/projects/bgmp/jonasg/bioinfo/bi624/Deduper/in/Dataset3.sam
#input=/projects/bgmp/jonasg/bioinfo/bi624/Deduper/in/test_14000.sam
output=/projects/bgmp/jonasg/bioinfo/bi624/Deduper/out/Dataset3_deduped.sam

python ./deduper_v2.py -db $data_base -i $input -u /projects/bgmp/jonasg/bioinfo/bi624/Deduper/in/STL96.txt -o $output -p True -t 40 -s 4999546

# 1 757912
# 2 874147
# 3 4999546