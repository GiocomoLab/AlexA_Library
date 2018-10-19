# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

"""Adding a column in the cluster view.
This file should be saved in `~/.phy/plugins/custom_stats.py`.
In addition, you need to edit `~/.phy/phy_config.py` and add the following
line:
```
c.KwikController.plugins = ['MyKwikPlugin']
```
"""

import numpy as np
from phy import IPlugin


import pdb

class AmplitudePlugin(IPlugin):
    def attach_to_controller(self, controller):
        """This method is called when a controller is created.
        The controller emits a few events:
        * `init()`: when the controller is created
        * `create_gui(gui)`: when the controller creates a GUI
        * `add_view(gui, view)`: when a view is added to a GUI
        You can register callback functions to these events.
        """

        # The controller defines several objects for the GUI.

        # The ManualClustering instance is responsible for the manual
        # clustering logic and the cluster views.
        mc = controller.manual_clustering
        # The context provides `cache()` and `memcache()` methods to cache
        # functions on disk or in memory, respectively.
        ctx = controller.context

        # We add a column in the cluster view and set it as the default.
        @mc.add_column(default=True)
        # We memcache it.
        @ctx.memcache
        #print('plugin')
        
        
        def Amp20(cluster_id):
            # This function takes a cluster id as input and returns a scalar.
    
            # We retrieve the spike_ids and waveforms for that cluster.
            # waveforms is a (n_spikes, n_samples, n_channels) array.
            # data = controller.get_waveforms(cluster_id)[0]
            #pdb.set_trace()
            bc = controller.get_best_channel(cluster_id)
            #waveforms_b = controller._select_data(cluster_id,controller.all_waveforms,500)
            #wf=waveforms_b.data
            spike_ids = controller._select_spikes(cluster_id, 500, batch_size=None)
            wf=controller.all_waveforms[spike_ids]
            wf_m=wf.mean(axis=0)
            #plt.figure()
            #plt.subplot(2,1,1)
            #plt.imshow(wf_m.transpose())


            trace=wf_m[:,bc]
            amp=trace.max()-trace.min()
            #amp2=controller.get_waveforms_amplitude(cluster_id)
            #plt.subplot(2,1,2)
            #plt.plot(trace)
            #plt.show()
            #data = controller.get_amplitudes(cluster_id)
            #noise=controller.get_channel_noise()
            #wfa=controller.get_waveforms_amplitude(cluster_id)
            #pdb.set_trace()
            return amp
        
