import os
import glob
def mkdir_p(dir):
    '''make a directory (dir) if it doesn't exist'''
    if not os.path.exists(dir):
        os.mkdir(dir)
    

job_directory = "%s/attialex/speed_jobs_scripts" %os.environ['OAK']
scratch = os.environ['SCRATCH']
#data_dir = os.path.join(scratch, '/registration/')

# Make top level directories
mkdir_p(job_directory)
#mkdir_p(data_dir)

data_dir = '/oak/stanford/groups/giocomo/attialex'

files = glob.glob('/oak/stanford/groups/giocomo/attialex/NP_DATA/np*_gain_*.mat')
session_names = files

for sn in session_names:
    sn_path = os.path.basename(sn)
    sn_path = os.path.splitext(sn_path)[0]
    log_dir =os.path.join(job_directory,sn_path)
    mkdir_p(log_dir)
    job_file = os.path.join(job_directory,"%s.job" %sn_path)
    #mouse_data = os.path.join(data_dir, mouse)

    # Create mouse directories
    #mkdir_p(mouse_data)

    with open(job_file,'w') as fh:
        fh.writelines("#!/bin/bash\n")
        fh.writelines("#SBATCH -p owners\n")
        fh.writelines("#SBATCH --ntasks-per-node=1\n")
        fh.writelines("#SBATCH --cpus-per-task=1\n")
        fh.writelines("#SBATCH --job-name=%s.job\n" % sn)
        fh.writelines("#SBATCH -o "+log_dir+"/slurm.%N.%j.out\n")
        fh.writelines("#SBATCH -e " + log_dir+"/slurm.%N.%j.err\n")
        fh.writelines("#SBATCH --time=0:30:00\n")
        fh.writelines("#SBATCH --mem-per-cpu=2G\n")
        fh.writelines("#SBATCH --mail-type=FAIL\n")
        fh.writelines("#SBATCH --mail-user=attialex@stanford.edu\n")
        fh.writelines("ml matlab\n")
        fh.writelines("cd $HOME/AlexA_Library/NeuroPixel/speed_sorting\n")
        fh.writelines("matlab -nosplash -nodisplay -nojvm -r \'fit_speed(\""  + sn +"\");exit\'")

    os.system("sbatch %s" %job_file)
    print('submitted job for: ' + sn_path)