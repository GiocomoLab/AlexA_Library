import os

def mkdir_p(dir):
    '''make a directory (dir) if it doesn't exist'''
    if not os.path.exists(dir):
        os.mkdir(dir)
    

job_directory = "%s/.job" %os.environ['SCRATCH']
scratch = os.environ['SCRATCH']
data_dir = os.path.join(scratch, '/registration/')

# Make top level directories
mkdir_p(job_directory)
mkdir_p(data_dir)

mice=["AA_190110_022"]
root_dir=os.path.join(os.environ['OAK'],'attialex','DATA')
for mouse in mice:

    job_file = os.path.join(job_directory,"%s.job" %mouse)
    #mouse_data = os.path.join(data_dir, mouse)

    # Create mouse directories
    #mkdir_p(mouse_data)

    with open(job_file) as fh:
        fh.writelines("#!/bin/bash\n")
        fh.writelines("#SBATCH --job-name=%s.job\n" % mouse)
        fh.writelines("#SBATCH --output=.out/%s.out\n" % mouse)
        fh.writelines("#SBATCH --error=.out/%s.err\n" % mouse)
        fh.writelines("#SBATCH --time=4:00\n")
        fh.writelines("#SBATCH --mem=1\n")
        fh.writelines("#SBATCH --qos=normal\n")
        fh.writelines("#SBATCH --mail-type=ALL\n")
        fh.writelines("#SBATCH --mail-user=$USER@stanford.edu\n")
        fh.writelines("ml python/3.6.1")
        fh.writelines("python3 $HOME/AlexA_Library/suite2p_helpers/run_single_file.py " + root_dir +' ' + mouse)

    #os.system("sbatch %s" %job_file)