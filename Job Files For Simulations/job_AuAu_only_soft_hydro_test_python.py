import os
import random
import subprocess
import shutil
import time
import threading
import glob

def create_jetscape_jobs(num_jobs=30,
                         main_config='/home/sleethmr/jetbox/myJETSCAPE/JETSCAPE/config/jetscape_main.xml',
                         user_config='/home/sleethmr/jetbox/myJETSCAPE/JETSCAPE/config/jetscape_user_AA_2011.01430.xml'):
    """
    Create job directories, generate Slurm jobs, and submit them

    Args:
        num_jobs (int): Number of job directories to create
        main_config (str): Path to main XML configuration file
        user_config (str): Path to user XML configuration file
    """
    # Use current directory (build directory) as base path
    build_path = os.getcwd()
    data_rkegroup_path = '/data/rke_group/sleethmr'
    os.makedirs(data_rkegroup_path, exist_ok=True)

    EOS_path = '/home/sleethmr/jetbox/myJETSCAPE/JETSCAPE/build/EOS'
    iSS_tables_path = '/home/sleethmr/jetbox/myJETSCAPE/JETSCAPE/build/iSS_tables'
    iSS_parameters_path = '/home/sleethmr/jetbox/myJETSCAPE/JETSCAPE/build/iSS_parameters.dat'
    music_path = '/home/sleethmr/jetbox/myJETSCAPE/JETSCAPE/build/music_input'
    #LBT_tables_path = '/home/sleethmr/jetbox/myJETSCAPE/JETSCAPE/build/LBT-tables'
    freestream_path = '/home/sleethmr/jetbox/myJETSCAPE/JETSCAPE/build/freestream_input'
    smash_config_file = '/home/sleethmr/jetbox/myJETSCAPE/JETSCAPE/external_packages/smash/smash_config.yaml'
    smash_particles_file = '/home/sleethmr/jetbox/myJETSCAPE/JETSCAPE/external_packages/smash/smash_code/input/particles.txt'
    smash_decaymodes_file = '/home/sleethmr/jetbox/myJETSCAPE/JETSCAPE/external_packages/smash/smash_code/input/decaymodes.txt'
    
    #hydro_files_path = '/data/rke_group/AuAuOnlyHydroFiles/Results'
    
    # Get absolute paths to config files
    abs_main_config = os.path.abspath(main_config)
    abs_user_config = os.path.abspath(user_config)

    # Initialize the submitted_jobs list
    submitted_jobs = []

    # Create a thread-local random number generator
    thread_local_random = threading.local()

    for job_id in range(0,1):
    #For job_id == 9:
    # Create unique directory for each job
        dir_path = os.path.join(data_rkegroup_path, f"jetscape_AuAu_only_soft_hydro_centr0_10_10events_gridsmall_{job_id:03d}")
        os.makedirs(dir_path, exist_ok=True)

        # Initialize thread-local random generator if not exists
        if not hasattr(thread_local_random, 'generator'):
            thread_local_random.generator = random.Random()

        # Generate random seed using thread-local generator
        random_seed = thread_local_random.generator.randint(1, 1000000)

        job_copy_dir_EOS = os.path.join(dir_path, os.path.basename(EOS_path))
        shutil.copytree(EOS_path, job_copy_dir_EOS)

        job_copy_dir_iSS_tables = os.path.join(dir_path, os.path.basename(iSS_tables_path))
        shutil.copytree(iSS_tables_path, job_copy_dir_iSS_tables)

        job_copy_iSS_parameters_file = os.path.join(dir_path, os.path.basename(iSS_parameters_path))
        shutil.copy2(iSS_parameters_path, job_copy_iSS_parameters_file)

        #job_copy_dir_LBT_tables = os.path.join(dir_path, os.path.basename(LBT_tables_path))
        #shutil.copytree(LBT_tables_path, job_copy_dir_LBT_tables)
        
        job_copy_music_file = os.path.join(dir_path, os.path.basename(music_path))
        shutil.copy2(music_path, job_copy_music_file)

        job_copy_freestream_file = os.path.join(dir_path, os.path.basename(freestream_path))
        shutil.copy2(freestream_path, job_copy_freestream_file)

        job_copy_smash_config_file = os.path.join(dir_path, os.path.basename(smash_config_file))
        shutil.copy2(smash_config_file, job_copy_smash_config_file)

        job_copy_smash_particles_file = os.path.join(dir_path, os.path.basename(smash_particles_file))
        shutil.copy2(smash_particles_file, job_copy_smash_particles_file)

        job_copy_smash_decaymodes_file = os.path.join(dir_path, os.path.basename(smash_decaymodes_file))
        shutil.copy2(smash_decaymodes_file, job_copy_smash_decaymodes_file)
        
        # Add error handling and retry mechanism
        max_attempts = 3
        for attempt in range(max_attempts):
            try:
                # Create Slurm job file
                with open(os.path.join(dir_path, "jetscape_AuAu_only_soft_hydro_centr0_10_10events_gridsmall.slurm"), "w") as job_file:
                    job_file.write(f"""#!/bin/bash
#SBATCH --job-name=jetscape_AuAu_only_soft_hydro_centr0_10_10events_gridsmall_{job_id}
#SBATCH --output={dir_path}/slurm_output.log
#SBATCH --error={dir_path}/slurm_error.log
#SBATCH --time=14-00:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1
#SBATCH --mem=64G  # Add this line to request more memory


# Export additional diagnostic environment variables
export APPTAINER_BINDPATH="/panfs,/data,/nobackup,/data/rke_group/sleethmr"

cd {dir_path}

#apptainer overlay create --size 1024 /tmp/my_overlay.img

# Run with additional error checking
apptainer exec \
  /home/sleethmr/jetbox/jetbox \
  /home/sleethmr/jetbox/myJETSCAPE/JETSCAPE/build/runJetscape \
  {abs_user_config} \
  {abs_main_config} \
  || echo "Jetscape simulation failed" >&2
""")

                # Make job file executable
                os.chmod(os.path.join(dir_path, "jetscape_AuAu_only_soft_hydro_centr0_10_10events_gridsmall.slurm"), 0o755)

                # Submit the job to Slurm
                result = subprocess.run(['sbatch', os.path.join(dir_path, 'jetscape_AuAu_only_soft_hydro_centr0_10_10events_gridsmall.slurm')],
                                        capture_output=True,
                                        text=True,
                                        timeout=300,  # 5-minute timeout
                                        check=True)
                
                # Store the job ID
                submitted_jobs.append({
                    'dir': dir_path,
                    'job_id': result.stdout.strip().split()[-1]
                })
                print(f"Submitted job in {dir_path}")
                
                # If successful, break retry loop
                break
            
            except subprocess.TimeoutExpired:
                print(f"Job submission timed out on attempt {attempt + 1}")
                if attempt == max_attempts - 1:
                    print(f"Failed to submit job in {dir_path} after {max_attempts} attempts")
            
            except subprocess.CalledProcessError as e:
                print(f"Job submission failed on attempt {attempt + 1}: {e}")
                if attempt == max_attempts - 1:
                    print(f"Failed to submit job in {dir_path} after {max_attempts} attempts")
                
                # Add a longer delay between retry attempts
                #time.sleep(random.uniform(1, 3))

        # Add a small random delay between job submissions
        #time.sleep(random.uniform(0, 0.5))

    # Print summary of submitted jobs
    print("\nJob Submission Summary:")
    print(f"Total Jobs Submitted: {len(submitted_jobs)}")
    #for job in submitted_jobs:
    #    print(f"Fragory: {job_frag119_1150['dir']}, Slurm Job ID: {job_frag['job_id']}")

# Run the job creation and submission script
create_jetscape_jobs()
