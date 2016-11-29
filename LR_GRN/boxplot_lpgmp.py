#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 29 15:18:01 2016

@author: Jiajia
"""
import sys
import matplotlib.pyplot as plt
#import numpy as np

#np.random.seed(10)
#collectn_1 = np.random.normal(100, 10, 200)
#collectn_2 = np.random.normal(80, 30, 200)
#collectn_3 = np.random.normal(90, 20, 200)
#collectn_4 = np.random.normal(70, 25, 200)

#data_to_plot = [collectn_1, collectn_2, collectn_3, collectn_4]

def boxplot_lpgmp(data_to_plot,xlabels,width,height,file_path):
    try:
        # Create a figure instance
        fig = plt.figure(1, figsize=(width , height))

        # Create an axes instance
        ax = fig.add_subplot(111)

        # Create the boxplot
        bp = ax.boxplot(data_to_plot)

        ## add patch_artist=True option to ax.boxplot() 
        ## to get fill color
        bp = ax.boxplot(data_to_plot, patch_artist=True)
        
        ## change outline color, fill color and linewidth of the boxes
        for box in bp['boxes']:
            # change outline color
            box.set( color='#252525', linewidth=2)
            # change fill color
            box.set( facecolor = '#bdbdbd' )
    
        ## change color and linewidth of the whiskers
        for whisker in bp['whiskers']:
            whisker.set(color='#636363', linewidth=2)

        ## change color and linewidth of the caps
        for cap in bp['caps']:
            cap.set(color='#636363', linewidth=2)

        ## change color and linewidth of the medians
        for median in bp['medians']:
            median.set(color='#252525', linewidth=2)
 
        ## change the style of fliers and their fill
        for flier in bp['fliers']:
            flier.set(marker='o', color='#636363', alpha=0.5)


        ## Custom x-axis labels
        ax.set_xticklabels(xlabels)   

        ## Remove top axes and right axes ticks
        ax.get_xaxis().tick_bottom()
        ax.get_yaxis().tick_left()


        # Save the figure
        fig.savefig(file_path, bbox_inches='tight')
    except AttributeError:
        sys.exit()

