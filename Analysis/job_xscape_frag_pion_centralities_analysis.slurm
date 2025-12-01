#!/bin/bash

#SBATCH --cpus-per-task=1

#SBATCH --time=14-00:00:00

#SBATCH --job-name=pion_AuAu_pTHat15_100_hard_sector_hydro_z_deltaR_100k_events_centr0_10_analysis

#SBATCH --output=output_AuAu_pTHat15_100_hard_sector_hydro_z_deltaR_centr0_10_100k_events.log

#SBATCH --error=error_AuAu_pTHat15_100_hard_sector_hydro_z_deltaR_centr0_10_100k_events.log

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
#output_file="pion_pp_pTHat35_100_out_alldecayoff_noprompt_10_test_${job_id}_final_state_hadrons.dat"

#export PATH=/data/rke_group/Software/JEWEL/jewel-2.4.0/:$PATH
#sed -i "s|<outputFilename>.*</outputFilename>|<outputFilename>${output_file}</outputFilename>|" ../config/jetscape_pp_frag_pions.xml

BASE_DIR="/data/rke_group/AuAuOnlyHydroFiles/Simulations"

# Array of pTHat ranges to process
centrality_ranges=(
    "0" "1" "2" "3" "4" "5" "6" "7" "8" "9" 
)

pTHat_ranges=(
    "15_100"
)

FILE_PATHS=()
for pTHat in "${pTHat_ranges[@]}"; do
    #FILE_NAME="jetscape_AuAu_pTHat5_10_hydro_1k_centr${RANGE}/test_out_pTHat5_10_hydro_AuAu_200_centr_${RANGE}_1k_events_final_state_hadrons.dat"
    #FILE_NAME="job_AuAu_only_soft_sector_hydro_1_event_centr${RANGE}_test_000/event_${RANGE}_final_state_hadrons.dat"
    #FILE_NAME="event_${RANGE}_final_state_hadrons.dat"
    for i in $(seq -w 0 9); do
		for centrality in "${centrality_ranges[@]}"; do
	    	FILE_NAME="jetscape_AuAu_pTHat${pTHat}_frag_hydro_file_10k_centr${centrality}_${i}_final_state_hadrons.dat"
    #photon_pp_pTHat5_10_out_alldecayoff_prompt_nohydro_10k_final_state_hadrons.dat
    #photon_pp_pTHat1_3_out_alldecayoff_1M_final_state_hadrons.dat
    #photon_pp_pTHat19_20_out_alldecayoff_1M_fixedPID_final_state_hadrons.dat
		
	# Construct the full path to the .dat file
	#DAT_FILE="${BASE_DIR}/job_direct${RANGE}_${i}/${FILE_NAME}"
	    	DAT_FILE="${BASE_DIR}/jetscape_AuAu_pTHat${pTHat}_frag_hydro_file_10k_centr${centrality}_00${i}/${FILE_NAME}"
    
	# Check if the file exists
		if [ -f "$DAT_FILE" ]; then
	    	FILE_PATHS+=("$DAT_FILE")
	    
		else
	    	echo "File $DAT_FILE does not exist."
		fi
		done   
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
apptainer exec  /home/sleethmr/jetbox/jetbox ./test_jet_etaphi "${FILE_PATHS[@]}"

#singularity exec --bind /data/rke_group/sleethmr:/data/rke_group/sleethmr /home/sleethmr/jetbox/jet.sif ./pionfrag_func_deltaR_slurm /data/rke_group/sleethmr/job_0*/pion_pp_pTHat15_100_out_alldecayoff_prompt_hydro_10M_final_state_hadrons.dat

#singularity exec --bind /data/rke_group/sleethmr:/mnt /home/sleethmr/jetbox/jet.sif /home/sleethmr/jetbox/myX-SCAPE/X-SCAPE/build/Analysis/pionfrag_func_deltaR_slurm_norm_integral /mnt/pion_pp_pTHat35_100_out_alldecayoff_noprompt_10M_1_final_state_hadrons.dat /mnt/pion_pp_pTHat35_100_out_alldecayoff_noprompt_10M_2_final_state_hadrons.dat /mnt/pion_pp_pTHat35_100_out_alldecayoff_noprompt_10M_3_final_state_hadrons.dat /mnt/pion_pp_pTHat35_100_out_alldecayoff_noprompt_10M_4_final_state_hadrons.dat /mnt/pion_pp_pTHat35_100_out_alldecayoff_noprompt_10M_5_final_state_hadrons.dat /mnt/pion_pp_pTHat35_100_out_alldecayoff_noprompt_10M_6_final_state_hadrons.dat /mnt/pion_pp_pTHat35_100_out_alldecayoff_noprompt_10M_7_final_state_hadrons.dat /mnt/pion_pp_pTHat35_100_out_alldecayoff_noprompt_10M_8_final_state_hadrons.dat /mnt/pion_pp_pTHat35_100_out_alldecayoff_noprompt_10M_9_final_state_hadrons.dat /mnt/pion_pp_pTHat35_100_out_alldecayoff_noprompt_10M_10_final_state_hadrons.dat

#    --output pion_pp_pTHat35_100_out_alldecayoff_noprompt_10M_final_state_hadrons_$SLURM_ARRY_TASK_ID.dat \
#    --seed $SEED

#jewel-2.4.0-simple /home/shengj2/WORK/jewel/simple/RHIC/params_PbPb_ppyq.dat


echo "Job Finished"
