# load modules
module load boost/1.81.0-m2umk6c
module load cmake/3.26.2-gbzqcn2
module load autoconf/2.69-cd4mshe
module load automake/1.16.5-omfl6zl

# add rb to the path
export PATH=$PATH:/work/LAS/phylo-lab/petrucci/revbayes/projects/cmake

# go to the correct directory
cd /work/LAS/phylo-lab/petrucci/ssetests_chap1/

# get the rep value
rep=$SLURM_ARRAY_TASK_ID

# create a file to hold the rep
touch aux/aux_${1}_$rep.Rev

# echo the definitions on it
printf "rep <- " >> aux/aux_${1}_$rep.Rev
echo $rep >> aux/aux_${1}_$rep.Rev

# source it, the parameter combination, and the actual script
timeout 121h rb aux/aux_${1}_$rep.Rev analysis/refs/refs_${1}.Rev analysis/master.Rev

# timeout and requeue if it takes more than 24h, MCMC is probably stuck
if [[ $? == 124 ]]; then
  scontrol requeue ${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}
fi

# remove file
rm aux/aux_${1}_$rep.Rev
