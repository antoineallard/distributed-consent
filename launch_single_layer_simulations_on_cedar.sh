#!/bin/bash

# Antoine Allard
# antoineallard.info
# Octobre 2021

# Ce code soumet une tache sur cedar.

# =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
# Ecriture du script de soumission et des commandes a executer.
cat <<END_OF_SCRIPT > model_script.pbs
#!/bin/bash
# ---------------------------------------------------------------------
# SLURM script for job submission on a Compute Canada cluster.
# ---------------------------------------------------------------------
#SBATCH --job-name=NAME1
#SBATCH --account=def-aallard
#SBATCH --time=TIME
#SBATCH --output=log_files/%x-%j.txt
# ---------------------------------------------------------------------
echo ""
echo "Current working directory: \`pwd\`"
echo "Starting run at: \`date\`"
echo "Submitted by launch_single_layer_simulations_on_cedar.sh"
# ---------------------------------------------------------------------
echo ""
echo "Job Name: \$SLURM_JOB_NAME"
echo "Job ID: \$SLURM_JOB_ID"
echo ""
# ---------------------------------------------------------------------

# adoption_within_private_profiles_values=\$(seq START STEP STOP)
for adoption_within_private_profiles in \${adoption_within_private_profiles_values[*]}; do
  echo \$adoption_within_private_profiles
  ./bin/single_layer NAME1.txt OBSERVATION_DEPTH_L APP_COVERAGE FRACTION_OF_PRIVATE_PROFILES \$adoption_within_private_profiles NB_SIMULATIONS NAME1 > \$HOME/scratch/\$SLURM_JOB_ID.txt
done

cat \$HOME/scratch/\$SLURM_JOB_ID.txt >> OUTPUT_FILENAME
rm \$HOME/scratch/\$SLURM_JOB_ID.txt

# ---------------------------------------------------------------------
echo "Job finished with exit code \$? at: \`date\`"
echo ""
# ---------------------------------------------------------------------

exit

END_OF_SCRIPT

network_name=$1
observation_depth_L=$2
app_coverage=$3
fraction_of_private_profiles=$4
nb_simulations=$5
output_filename=$6
start=$7
step=$8
stop=$9

time="11:00:00"


#  Copies the model script (see above) file and changes the variables.
sed -i 's,NAME1,'"${network_name}"',g'                                        model_script.pbs
sed -i 's,OBSERVATION_DEPTH_L,'"${observation_depth_L}"',g'                   model_script.pbs
sed -i 's,APP_COVERAGE,'"${app_coverage}"',g'                                 model_script.pbs
sed -i 's,FRACTION_OF_PRIVATE_PROFILES,'"${fraction_of_private_profiles}"',g' model_script.pbs
sed -i 's,NB_SIMULATIONS,'"${nb_simulations}"',g'                             model_script.pbs
sed -i 's,OUTPUT_FILENAME,'"${output_filename}"',g'                           model_script.pbs
sed -i 's,START,'"${start}"',g'                                               model_script.pbs
sed -i 's,STEP,'"${step}"',g'                                                 model_script.pbs
sed -i 's,STOP,'"${stop}"',g'                                                 model_script.pbs
sed -i 's,TIME,'"${time}"',g'                                                 model_script.pbs

#  Submits the task.
sbatch model_script.pbs

# Deletes the script.
# rm model_script.pbs
