'''Concatenate multiple files into a single virtual dataset
    TODO: discard last frame if odd number of total frames
'''
import h5py
import numpy as np

file_names_to_concatenate = ['AA_190111_027_002_027.h5', 'AA_190111_027_002_027.h5']
entry_key = 'data' # where the data is inside of the source files.
total_frames=0
for filename in file_names_to_concatenate:
    sh = h5py.File(filename, 'r')[entry_key].shape # get the first ones shape.
    n_frames = sh[0]
    total_frames = total_frames + n_frames




layout = h5py.VirtualLayout(shape=(total_frames,sh[1],sh[2]),
                            dtype=np.int16)
total_frames=0
for i, filename in enumerate(file_names_to_concatenate):
    sh = h5py.File(filename, 'r')[entry_key].shape
    vsource = h5py.VirtualSource(filename, entry_key, shape=sh)
    layout[total_frames:total_frames+sh[0], :, :] = vsource
    total_frames=total_frames+sh[0]

with h5py.File("VDS2.h5", 'w', libver='latest') as f:
    f.create_virtual_dataset(entry_key, layout, fillvalue=0)