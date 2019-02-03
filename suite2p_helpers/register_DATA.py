import os

def mkdir_p(dir):
    '''make a directory (dir) if it doesn't exist'''
    if not os.path.exists(dir):
        os.mkdir(dir)
    

job_directory = "%s/register_jobs" %os.environ['SCRATCH']
scratch = os.environ['SCRATCH']
#data_dir = os.path.join(scratch, '/registration/')

# Make top level directories
mkdir_p(job_directory)
#mkdir_p(data_dir)


mice=["AA_190110_022", "AA_190111_026","AA_190111_027","AA_190111_029"]
root_dir=os.path.join(os.environ['OAK'],'attialex','DATA')
for mouse in mice:
    log_dir =os.path.join(job_directory,mouse)
    mkdir_p(log_dir)
    job_file = os.path.join(job_directory,"%s.job" %mouse)
    #mouse_data = os.path.join(data_dir, mouse)

    # Create mouse directories
    #mkdir_p(mouse_data)

    with open(job_file,'w') as fh:
        fh.writelines("#!/bin/bash\n")
        fh.writelines("#SBATCH -p giocomo\n")
        fh.writelines("#SBATCH --ntasks-per-node=1\n")
        fh.writelines("#SBATCH --cpus-per-task=4\n")
        fh.writelines("#SBATCH --job-name=%s.job\n" % mouse)
        fh.writelines("#SBATCH -o "+log_dir+"/slurm.%N.%j.out\n")
        fh.writelines("#SBATCH -e " + log_dir+"/slurm.%N.%j.err\n")
        fh.writelines("#SBATCH --time=12:00:00\n")
        fh.writelines("#SBATCH --mem-per-cpu=40G\n")
        fh.writelines("#SBATCH --mail-type=BEGIN,END,FAIL\n")
        fh.writelines("#SBATCH --mail-user=attialex@stanford.edu\n")
        fh.writelines("ml python/3.6.1\n")
        fh.writelines("python3 $HOME/AlexA_Library/suite2p_helpers/run_single_file.py " + root_dir +' ' + mouse)

    os.system("sbatch %s" %job_file)