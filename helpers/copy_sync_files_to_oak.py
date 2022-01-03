import glob
import os
from distutils.dir_util import copy_tree, remove_tree
import errno
from distutils.file_util import copy_file
from distutils import log
log.set_verbosity(log.INFO)
log.set_threshold(log.INFO)

files = glob.glob(r'H:\catGT2\catgt_AA_210114_5*\*\*SY_384_6_500.txt')

for iF in range(0,len(files)):
    (root,tmp) = os.path.split(files[iF])
    #json_files = glob.glob(os.path.join(root,'*.json'))
    sync_file = glob.glob(os.path.join(root,'*SY_384_6_500.txt'))
    #cluster_file = [os.path.join(os.path.join(files[iF],'cluster_group.tsv'))]
    sync_nidaq = glob.glob(os.path.join(os.path.split(root)[0],'*XA_0_500.txt'))
    file_list = sync_file+sync_nidaq
    print(sync_nidaq)

    parts = files[iF].split('\\')

    session = parts[2][6::]
    print('session: '+session)
    animal = session[0:11]
    print('animal: '+animal)
    oak = r'Z:\giocomo\export\data\Projects\AlexA_NP'

    oak_animal = os.path.join(oak,animal)
    oak_animal_ks_data = os.path.join(oak_animal,'ks_data')

    target_dir = os.path.join(oak_animal_ks_data)
    target_dir = os.path.join(target_dir,session)
    #print('copy {} to {}'.format(sync_file[0],target_dir))
    for jfi in file_list:
        (a,b)=copy_file(jfi,target_dir,verbose=1,dry_run=0)