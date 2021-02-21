import numpy as np
import glob

fi = glob.glob(r"F:\CatGT2\*\*\imec0_ks2\cluster_KSLabel.tsv")

for fp in fi:
    data=np.recfromcsv(fp,delimiter='\t')
    print(fp)
    cluster = np.array([clu[0] for clu in data])
    label = np.array([clu[1] for clu in data])

    sid=np.lexsort((cluster,label))

    ranking = sid.argsort()
    out = np.vstack((cluster,ranking))
    out_name = fp.replace('KSLabel','order')
    np.savetxt(out_name,out.T,delimiter='\t',header = 'cluster_id\torder',fmt = '%i',comments='')
