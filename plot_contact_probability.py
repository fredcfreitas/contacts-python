#!/usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np


def getspacedelements(input_vector, num=4):
    """Function to return #num elements from an 1D array, evenly distributed if\
    possible.
    Input: 1-D list, #num of elements.
    Ouput: 1-D list with #num elements.
    """
    out = input_vector[np.round(np.linspace(0, len(input_vector)-1, \
                                                             num)).astype(int)]
    return out

def customize_axis_ticks(old_label, num=4):
    """
    Function to return the new labels and their respective positions.
    Input:
     - list of available ticks.
     - number of desired ticks.
    Output:
     - ticks indexes.
     - chosen ticks.
     """
    # What is going on?
    chosen = getspacedelements(old_label, num=num)
    idx = np.isin(old_label, chosen)
    location = np.arange(np.shape(old_label)[0])[idx]
    return location, chosen


def plot_contact_probability(xlabel, ylabel, plot_title, datapoints, xvalues, \
                             num_xticks, yvalues, num_yticks, colorbar_label, \
                             colormap='rainbow', saveit=True):
    """
    Function to create a 2-D plot, similar to a heatmap.
    Input: One 2-D numpy array.
    """
    fig, ax = plt.subplots()
    img = ax.imshow(datapoints, cmap=colormap, aspect='auto', origin='lower')
    fig.colorbar(img, ax=ax, orientation='vertical', label=colorbar_label)
    ax.set_title(str(plot_title))
    ax.set_ylabel(str(ylabel))
    ax.set_xlabel(str(xlabel))
    x_ax_idx, x_ax_ticks = customize_axis_ticks(xvalues, num_xticks)
    plt.xticks(x_ax_idx, x_ax_ticks.astype(int))
    y_ax_idx, y_ax_ticks = customize_axis_ticks(yvalues, num_yticks)
    plt.yticks(y_ax_idx, y_ax_ticks.astype(int))
    if saveit:
        plt.savefig("plot_" + str(plot_title) + ".png", format="png", dpi=300)
    plt.show()
    return
