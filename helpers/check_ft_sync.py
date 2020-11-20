import glob
import numpy as np
import matplotlib.pyplot as plt

ft_files = glob.glob(r'Z:\giocomo\export\data\Projects\AlexA_NP\AA_2009*\*.ft')
time_diff = np.empty((len(ft_files),))
time_diff.fill(np.nan)
for idx,fi in enumerate(ft_files):
    try:
        ft_data = np.genfromtxt(fi)
        vr_name = fi.replace('.ft','.log')
        vr_data = np.genfromtxt(vr_name,skip_header=1)

        vid_duration = ft_data[-1,1]-ft_data[0,1]
        vr_duration = vr_data[-1,0] - vr_data[0,0]

        time_diff[idx]=vid_duration-vr_duration
        if abs(time_diff[idx])>2:
            print(fi)
            print(time_diff[idx])
    except Exception as e:
        print(e)


print(np.nanmean(time_diff))
print(np.nanmedian(time_diff))
plt.hist(time_diff,range=(-.2,.2),bins=30)
plt.xlabel('video - vr')
plt.show()