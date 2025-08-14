#!/bin/bash

#SBATCH --cpus-per-task=1

#SBATCH --time=14-00:00:00

#SBATCH --job-name=photon_AuAu_pTHat3_20_out_alldecayoff_prompt_1M_final_state_hadrons

#SBATCH --output=output_AuAu_hydro_direct_photon3_20_1M.log

#SBATCH --error=error_AuAu_hydro_direct_photon3_20_1M.log

#SBATCH --cpus-per-task=1

#SBATCH --ntasks=1

#SBATCH --mem=64G

#SBATCH --account=rke_group

#ml GCC/10.2.0

#SEED=$SLURM_ARRAY_TASK_ID
#export LD_LIBRARY_PATH=/data/rke_group/Software/LHAPDF6/local/lib:$LD_LIBRARY_PATH

#export LHAPATH=/data/rke_group/Software/LHAPDF6/local/share/LHAPDF/

#job_id=$1
#EVENTS_PER_JOB=10
#seed=$((job_id * EVENTS_PER_JOB))
#output_file="photon_pp_pTHat35_100_out_alldecayoff_noprompt_10_test_${job_id}_final_state_hadrons.dat"

#export PATH=/data/rke_group/Software/JEWEL/jewel-2.4.0/:$PATH
#sed -i "s|<outputFilename>.*</outputFilename>|<outputFilename>${output_file}</outputFilename>|" ../config/jetscape_pp_frag_photons.xml

BASE_DIR="/data/rke_group/sleethmr"

# Array of pTHat ranges to process
pTHat_ranges=(
    "3_4" "4_5" "5_6" "6_7" "7_8" "8_9" "9_10" "10_11" "11_12" "12_13" "13_14" "14_15" "15_16"
)

FILE_PATHS=()

for RANGE in "${pTHat_ranges[@]}"; do
    FILE_NAME="photon_AuAu_pTHat${RANGE}_out_alldecayoff_prompt_nohydro_1M_final_state_hadrons.dat"
    
    for i in $(seq -w 000 005); do	
	# Construct the full path to the .dat file
	DAT_FILE="${BASE_DIR}/job_direct${RANGE}_${i}/${FILE_NAME}"
    #DAT_FILE="${BASE_DIR}/${FILE_NAME}"
    
	# Check if the file exists
	if [ -f "$DAT_FILE" ]; then
	    FILE_PATHS+=("$DAT_FILE")

	else
            echo "File $DAT_FILE does not exist."
	fi
    done
done

# Check if we found any files to process
if [ ${#FILE_PATHS[@]} -eq 0 ]; then
    echo "No files found to process."
    exit 1
fi

echo "Processing $DAT_FILE"
export APPTAINER_BINDPATH="/panfs,/data,/nobackup,/data/rke_group/sleethmr"
 
# Run the singularity command
apptainer exec  /home/sleethmr/jetbox/jetbox ./direct_photon_cross_sec_slurm "${FILE_PATHS[@]}"
        

#singularity exec --bind /data/rke_group/sleethmr:/data/rke_group/sleethmr /home/sleethmr/jetbox/jet.sif ./photonfrag_func_deltaR_slurm /data/rke_group/sleethmr/job_0*/photon_AuAu_pTHat15_100_out_alldecayoff_prompt_hydro_10M_final_state_hadrons.dat

#singularity exec --bind /data/rke_group/sleethmr:/mnt /home/sleethmr/jetbox/jet.sif /home/sleethmr/jetbox/myX-SCAPE/X-SCAPE/build/Analysis/photonfrag_func_deltaR_slurm_norm_integral /mnt/photon_pp_pTHat35_100_out_alldecayoff_noprompt_10M_1_final_state_hadrons.dat /mnt/photon_pp_pTHat35_100_out_alldecayoff_noprompt_10M_2_final_state_hadrons.dat /mnt/photon_pp_pTHat35_100_out_alldecayoff_noprompt_10M_3_final_state_hadrons.dat /mnt/photon_pp_pTHat35_100_out_alldecayoff_noprompt_10M_4_final_state_hadrons.dat /mnt/photon_pp_pTHat35_100_out_alldecayoff_noprompt_10M_5_final_state_hadrons.dat /mnt/photon_pp_pTHat35_100_out_alldecayoff_noprompt_10M_6_final_state_hadrons.dat /mnt/photon_pp_pTHat35_100_out_alldecayoff_noprompt_10M_7_final_state_hadrons.dat /mnt/photon_pp_pTHat35_100_out_alldecayoff_noprompt_10M_8_final_state_hadrons.dat /mnt/photon_pp_pTHat35_100_out_alldecayoff_noprompt_10M_9_final_state_hadrons.dat /mnt/photon_pp_pTHat35_100_out_alldecayoff_noprompt_10M_10_final_state_hadrons.dat

#    --output photon_pp_pTHat35_100_out_alldecayoff_noprompt_10M_final_state_hadrons_$SLURM_ARRY_TASK_ID.dat \
#    --seed $SEED

#jewel-2.4.0-simple /home/shengj2/WORK/jewel/simple/RHIC/params_PbPb_ppyq.dat


echo "Job Finished"
