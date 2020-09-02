#!/bin/bash
#SBATCH --mail-user=labellal@vanderbilt.edu  # email address
#SBATCH --mail-type=END  # Alerts sent when job begins, ends, or aborts
#SBATCH --nodes=2   # Number of nodes required
#SBATCH --ntasks=2   # Number of nodes required
#SBATCH --mem=16G  # Total memory_RAM_required, per node
#SBATCH --array=1-9
#SBATCH --time=08:00:00  # Wall Clock time_dd-hh:mm:ss_[max of 14 days]
#SBATCH --output=NAMEHERE_perm_test_%A_%a.out  # output and error messages go to this file

module load GCC/6.4.0-2.28
module load OpenMPI/2.1.1
module load R/3.4.3

#wi_perm.pl NAMEHERE_output_wi_file.txt 10000



for i in NAMEHERE*wi_file.$SLURM_ARRAY_TASK_ID*.txt
do
n=${i/.txt/}
n=${n/*wi_file./}
tAI_calc_perm.pl $i NAMEHERE.fas
spec=${i/_output_wi_file.*.txt/}
s_test_combine_perm.pl $spec.tAI.$n.txt
Rscript s_test_calc.R $spec.stest.$n.txt >> NAMEHERE.stest.all.txt
rm $spec.tAI.$n.txt
rm $spec.stest.$n.txt
rm $i
done
