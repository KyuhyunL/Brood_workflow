#!/bin/bash
#SBATCH -N 1
#SBATCH -n 16
#SBATCH -J Brood_test
#SBATCH -p dev,dev48,dev96,cadd
#SBATCH -o Brood_test_%j.log
#SBATCH -e Brood_test_%j.err

MYCODE='/db2/users/kyuhyunlee/git_repos/Brood_workflow/src/brood_workflow.py'
MYCORE='COc1cc(Nc2cc(ncn2)-c2cc(cc3ccoc23)C#N)ccc1N1CCOCC1'
MYSCAFFOLD='N(c1ccccc1)c1ccncn1'

python $MYCODE $MYCORE $MYSCAFFOLD --nproc 16
