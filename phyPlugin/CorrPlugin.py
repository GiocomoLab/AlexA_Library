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


class CorrPlugin(IPlugin):
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
        mc = controller.supervisor
        # The context provides `cache()` and `memcache()` methods to cache
        # functions on disk or in memory, respectively.
        ctx = controller.context

        # We add a column in the cluster view and set it as the default.
        @mc.add_column(default=True)
        # We memcache it.
        @ctx.memcache
        
        def VioRate132(clusterID):
            #n_bins = 100
            #middle_bin = 50
            #violation = 2
            thresh=0.002
            #to = middle_bin+violation+1
            clu_idx=controller.spike_clusters==clusterID
            spike_times=controller.spike_times[clu_idx]
            diffST=np.diff(spike_times)
            n_violated=np.count_nonzero(diffST<thresh);
            total_spikes=len(spike_times)*1.0
            #data = controller.get_correlograms((clusterID,),1,n_bins);
            #total_spikes=data[0,0,middle_bin:].sum();
            #violation_spikes=data[0,0,middle_bin:to].sum();
            #pdb.set_trace()
            violation_rate = n_violated/total_spikes
            return violation_rate