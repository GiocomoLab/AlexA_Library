import numpy as np
import glob
import os
def single_reg_file(opspath, regpath):
    ops=np.load(opspath)
    ops=ops.item()
    ops['reg_file']=regpath
    ops.pop('reg_file_raw',None) #kick out the raw file path otherwise suite2p will crash looking for it
    np.save(opspath,ops)

def do_for_animal(animaldir):
    for filename in glob.iglob(animaldir + r'\**\plane*\ops.npy', recursive=True):
        opspath = filename
        exp = filename.split('\\')[8]
        plane= filename.split('\\')[10]
        regpath = '//oak-smb-giocomo.stanford.edu/groups/giocomo/attialex/REGISTRATION/'+ exp +'/suite2p/' + plane +'/data.bin'
        #print(regpath)
        #print(os.path.isfile(regpath))
        single_reg_file(opspath,regpath)


if __name__=='__main__':
    #opspath = 'Z:/giocomo/attialex/DATA/AA_190111_029/AA_190111_029_000_001/suite2p/plane1/ops.npy'
    #regpath = 'Z:/giocomo/attialex/REGISTRATION/AA_190111_029_000_001/suite2p/plane1/data.bin'
    #single_reg_file(opspath,regpath)
    do_for_animal(r'\\oak-smb-giocomo.stanford.edu\groups\giocomo\attialex\DATA\AA_190111_027')

