# -*- coding: utf-8 -*-
"""
Created on Mon Oct 22 15:20:41 2018

@author: giocomolab
"""

from phy import IPlugin
import numpy as np
import matplotlib.pyplot as plt


class CustomView(IPlugin):
    def attach_to_controller(self, c):

        # Create the figure when initializing the GUI.
        f, ax = plt.subplots()

        @c.connect
        def on_create_gui(gui):
            # Called when the GUI is created.

            # We add the matplotlib figure to the GUI.
            gui.add_view(f, name='ISI')

            # We connect this function to the "select" event triggered
            # by the controller at every cluster selection change.
            @c.connect
            def on_select(clusters, **kwargs):
                # We clear the figure.
                ax.clear()

                # We compute the ISI.
                spikes = c.spikes_per_cluster(clusters[0])
                ax.hist(np.diff(spikes), bins=50)

                # We update the figure.
                f.canvas.draw()