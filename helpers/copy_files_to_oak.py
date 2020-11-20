import glob
import os
from distutils.dir_util import copy_tree
import errno
from distutils.file_util import copy_file
from distutils import log
log.set_verbosity(log.INFO)
log.set_threshold(log.INFO)

files = glob.glob(r'F:\CatGT\cat*\*\imec0_ks2')

for iF in range(2,len(files)):
    (root,tmp) = os.path.split(files[iF])
    json_files = glob.glob(os.path.join(root,'*.json'))
    #print(json_files)

    parts = files[iF].split('\\')

    session = parts[2][6::]
    print('session: '+session)
    animal = session[0:11]
    print('animal: '+animal)
    oak = r'Z:\giocomo\export\data\Projects\AlexA_NP'

    oak_animal = os.path.join(oak,animal)
    if not os.path.isdir(oak_animal):
        print(oak_animal)
        print('does not exist, aborting')
        exit()

    oak_animal_ks_data = os.path.join(oak_animal,'ks_data')

    if not os.path.isdir(oak_animal_ks_data):
        os.mkdir(oak_animal_ks_data)

    target_dir = os.path.join(oak_animal_ks_data)
    target_dir = os.path.join(target_dir,session)
    if not os.path.isdir(target_dir):
        os.mkdir(target_dir)


    copy_tree(files[iF],target_dir,verbose=1)
    for jfi in json_files:
        (a,b)=copy_file(jfi,target_dir,verbose=1,dry_run=0)
        #print(a)
        #print(b)
    #print('from: '+files[iF])
    #print('to: ' +target_dir)
    #print('#######################')


