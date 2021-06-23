#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import InsetPosition


def test_normalized(vector):
    """
    Function to test if a vector is normalized. If this is true, a formatted\
    output with 3 decimal places will be returned.
    """
    if np.equal(np.max(vector), 1):
        output = ["{:.3f}".format(i) for i in vector]
    else:
        output = np.asarray(vector).astype(int)
    return output


def getspacedelements(input_vector, num=4):
    """Function to return #num elements from an 1D array, evenly distributed if\
    possible.
    Input: 1-D list, #num of elements.
    Output: 1-D list with #num elements.
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


def format_output_name(input_string):
    """Function to get a string and format it if it is necessary"""
    split_string = str(input_string).split()
    if np.less_equal(np.size(split_string), 1):
        output = split_string
    else:
        two_firsts = split_string[0] + "-" + split_string[1]
        output = two_firsts
    return output


def plot_contact_probability(name, xlabel, ylabel, plot_title, datapoints, \
                             xvalues, num_xticks, yvalues, num_yticks, \
                             colorbar_label, colormap='rainbow', saveit=True):
    """
    Function to create a 2-D plot, similar to a heatmap.
    Input: One 2-D numpy array.
    """
    fig, ax = plt.subplots()
    img = ax.imshow(datapoints, cmap=colormap, aspect='auto', origin='lower', \
                    interpolation='hanning')
    fig.colorbar(img, ax=ax, orientation='vertical', label=colorbar_label)
    ax.set_title(str(plot_title))
    ax.set_ylabel(str(ylabel))
    ax.set_xlabel(str(xlabel))
    x_ax_idx, x_ax_ticks = customize_axis_ticks(xvalues, num_xticks)
    x_ax_ticks_f = test_normalized(x_ax_ticks)
    plt.xticks(x_ax_idx, x_ax_ticks_f)
    y_ax_idx, y_ax_ticks = customize_axis_ticks(yvalues, num_yticks)
    y_ax_ticks_f = test_normalized(y_ax_ticks)
    plt.yticks(y_ax_idx, y_ax_ticks_f)
    if saveit:
        plt.savefig("plot_" + str(name) + ".png", format="png", dpi=300)
    plt.show()
    pass




def plot_three_aside_one_colorbar(datapoints_1, xvalues_1, yvalues_1, datapoints_2, xvalues_2, yvalues_2, datapoints_3, xvalues_3, yvalues_3, x_label, y_label, colorbar_label, fig_name,  num_xticks = 5, num_yticks = 8, colormap="rainbow"):
    """
    Function to plot three graphs with the same y and x axis and same colorbar \
    in the right side. A lot of features can (and maybe should) be changed.
    In this version all x values will be normalized. This should be corrected for a generalization.
    """
    colormap = colormap
    colorbar_label = colorbar_label
    x_label = x_label
    num_xticks = 5
    y_label = y_label
    num_yticks = 8

    fig, (ax, ax1, ax2, cax) = plt.subplots(ncols=4, figsize=(10, 4),  gridspec_kw={"width_ratios":[1, 1, 1, 0.05]})
    # fig.subplots_adjust(wspace=0.4)

    fig.subplots_adjust(top=0.9, bottom=0.146, left=0.12, right=0.98, hspace=0.13, wspace=0.35)


    im  = ax.imshow(datapoints_1, cmap=colormap, aspect='auto', origin='lower', interpolation='hanning')
    im1 = ax1.imshow(datapoints_2, cmap=colormap, aspect='auto', origin='lower', interpolation='hanning')
    im2 = ax2.imshow(datapoints_3, cmap=colormap, aspect='auto', origin='lower', interpolation='hanning')
    #ax.set_ylabel("atom")
    # ax.set_xlabel("Q")
    # ax1.set_xlabel("Q")
    # ax2.set_xlabel("Q")
    #fig.suptitle(title, fontsize=16)

    fig.text(0.5, 0.04, x_label, ha='center', va='center', fontsize=16)
    fig.text(0.06, 0.5, y_label, ha='center', va='center', rotation='vertical', fontsize=16)

    plt.sca(ax)
    contacts = xvalues_1
    atoms = yvalues_1

    xvalues = np.divide(contacts, contacts.max())
    x_ax_idx, x_ax_ticks = customize_axis_ticks(xvalues, num_xticks)
    x_ax_ticks_f = test_normalized(x_ax_ticks)
    plt.xticks(x_ax_idx, x_ax_ticks_f)

    yvalues = atoms
    y_ax_idx, y_ax_ticks = customize_axis_ticks(yvalues, num_yticks)
    y_ax_ticks_f = test_normalized(y_ax_ticks)
    plt.yticks(y_ax_idx, y_ax_ticks_f)

    plt.sca(ax1)
    contacts = xvalues_2
    atoms = yvalues_2

    xvalues = np.divide(contacts, contacts.max())
    x_ax_idx, x_ax_ticks = customize_axis_ticks(xvalues, num_xticks)
    x_ax_ticks_f = test_normalized(x_ax_ticks)
    plt.xticks(x_ax_idx, x_ax_ticks_f)

    yvalues = atoms
    y_ax_idx, y_ax_ticks = customize_axis_ticks(yvalues, num_yticks)
    y_ax_ticks_f = test_normalized(y_ax_ticks)
    plt.yticks(y_ax_idx, y_ax_ticks_f)


    plt.sca(ax2)
    contacts = xvalues_3
    atoms = yvalues_3

    xvalues = np.divide(contacts, contacts.max())
    x_ax_idx, x_ax_ticks = customize_axis_ticks(xvalues, num_xticks)
    x_ax_ticks_f = test_normalized(x_ax_ticks)
    plt.xticks(x_ax_idx, x_ax_ticks_f)

    yvalues = atoms
    y_ax_idx, y_ax_ticks = customize_axis_ticks(yvalues, num_yticks)
    y_ax_ticks_f = test_normalized(y_ax_ticks)
    plt.yticks(y_ax_idx, y_ax_ticks_f)

    ip = InsetPosition(ax2, [1.05, 0, 0.05, 1])
    cax.set_axes_locator(ip)
    fig.colorbar(im, cax=cax, ax=[ax, ax1, ax2], orientation='vertical', label=colorbar_label)

    # fig.tight_layout(pad=0.02, h_pad = 0.03, w_pad = 0.03, rect = [0.08, 0.0146, 0.0985, 9])

    plt.savefig(fig_name + ".png", format="png", dpi=500)

    plt.show()
    pass



def generating_color_axis(axis_ticks, value):
    """Function to generate a list with the colors before and after the value"""
    used_shape = np.shape(axis_ticks)
    b = np.full(used_shape, 'b')[np.less(axis_ticks, value)]
    r = np.full(used_shape, 'r')[np.greater_equal(axis_ticks, value)]
    colors = np.concatenate((b, r))
    return colors


def plot_three_aside_one_colorbar_mod(datapoints_1, xvalues_1, yvalues_1, datapoints_2, xvalues_2, yvalues_2, datapoints_3, xvalues_3, yvalues_3, x_label, y_label, colorbar_label, fig_name,  num_xticks = 5, num_yticks = 8, colormap="rainbow"):
    """
    Function to plot three graphs with the same y and x axis and same colorbar \
    in the right side. A lot of features can (and maybe should) be changed.
    """
    colormap = colormap
    colorbar_label = colorbar_label
    x_label = x_label
    num_xticks = num_xticks
    y_label = y_label
    num_yticks = num_yticks

    fig, (ax, ax1, ax2, cax) = plt.subplots(ncols=4, figsize=(10, 4), \
                                            gridspec_kw={"width_ratios":[1, 1, \
                                                                      1, 0.05]})
    # fig.subplots_adjust(wspace=0.4)
    fig.subplots_adjust(top=0.9, bottom=0.146, left=0.12, right=0.985, \
                        hspace=0.13, wspace=0.35)
    im  = ax.imshow(datapoints_1, cmap=colormap, aspect='auto', origin='lower',\
                    interpolation='hanning')
    im1 = ax1.imshow(datapoints_2, cmap=colormap, aspect='auto', \
                     origin='lower', interpolation='hanning')
    im2 = ax2.imshow(datapoints_3, cmap=colormap, aspect='auto', \
                     origin='lower', interpolation='hanning')
    #ax.set_ylabel("atom")
    # ax.set_xlabel("Q")
    # ax1.set_xlabel("Q")
    # ax2.set_xlabel("Q")
    #fig.suptitle(title, fontsize=16)
    fig.text(0.5, 0.04, x_label, ha='center', va='center', fontsize=16)
    fig.text(0.06, 0.5, y_label, ha='center', va='center', \
             rotation='vertical', fontsize=16)
    plt.sca(ax)
    contacts = xvalues_1
    atoms = yvalues_1
    xvalues = np.divide(contacts, 1)
    x_ax_idx, x_ax_ticks = customize_axis_ticks(xvalues, num_xticks)
    x_ax_ticks_f = test_normalized(x_ax_ticks)
    plt.xticks(x_ax_idx, x_ax_ticks_f)
    yvalues = atoms
    y_ax_idx, y_ax_ticks = customize_axis_ticks(yvalues, num_yticks)
    y_ax_ticks_f = test_normalized(y_ax_ticks)
    plt.yticks(y_ax_idx, y_ax_ticks_f)
    for ticklabel, tickcolor in zip(plt.gca().get_yticklabels(), \
                                    generating_color_axis(y_ax_ticks_f, 307)):
        ticklabel.set_color(tickcolor)
    plt.sca(ax1)
    contacts = xvalues_2
    atoms = yvalues_2
    xvalues = np.divide(contacts, 1)
    x_ax_idx, x_ax_ticks = customize_axis_ticks(xvalues, num_xticks)
    x_ax_ticks_f = test_normalized(x_ax_ticks)
    plt.xticks(x_ax_idx, x_ax_ticks_f)
    yvalues = atoms
    y_ax_idx, y_ax_ticks = customize_axis_ticks(yvalues, num_yticks)
    y_ax_ticks_f = test_normalized(y_ax_ticks)
    plt.yticks(y_ax_idx, y_ax_ticks_f)
    for ticklabel, tickcolor in zip(plt.gca().get_yticklabels(), \
                                    generating_color_axis(y_ax_ticks_f, 333)):
        ticklabel.set_color(tickcolor)
    plt.sca(ax2)
    contacts = xvalues_3
    atoms = yvalues_3
    xvalues = np.divide(contacts, 1)
    x_ax_idx, x_ax_ticks = customize_axis_ticks(xvalues, num_xticks)
    x_ax_ticks_f = test_normalized(x_ax_ticks)
    plt.xticks(x_ax_idx, x_ax_ticks_f)
    yvalues = atoms
    y_ax_idx, y_ax_ticks = customize_axis_ticks(yvalues, num_yticks)
    y_ax_ticks_f = test_normalized(y_ax_ticks)
    plt.yticks(y_ax_idx, y_ax_ticks_f)
    for ticklabel, tickcolor in zip(plt.gca().get_yticklabels(), \
                                    generating_color_axis(y_ax_ticks_f, 333)):
        ticklabel.set_color(tickcolor)
    ip = InsetPosition(ax2, [1.05,0,0.05,1])
    cax.set_axes_locator(ip)
    fig.colorbar(im, cax=cax, ax=[ax, ax1, ax2], orientation='vertical', \
                 label=colorbar_label)
    # fig.tight_layout(pad=0.02, h_pad = 0.03, w_pad = 0.03, rect = [0.08, 0.0146, 0.0985, 9])
    plt.savefig(fig_name + ".png", format="png", dpi=500)
    plt.show()
    pass




colormap = 'rainbow'
colorbar_label = '$p(Q,i)$'
num_xticks = 5
num_yticks = 8

title = "$Q_{bind}$"

fig, (ax, ax1, ax2, cax) = plt.subplots(ncols=4, figsize=(10, 4),  gridspec_kw={"width_ratios":[1, 1, 1, 0.05]})
# fig.subplots_adjust(wspace=0.4)

fig.subplots_adjust(top=0.9, bottom=0.146, left=0.12, right=0.985, hspace=0.13, wspace=0.35)


im  = ax.imshow(cov_1_raw_prob_AB.T, cmap=colormap, aspect='auto', origin='lower', interpolation='hanning')
im1 = ax1.imshow(cov_2_raw_prob_AB.T, cmap=colormap, aspect='auto', origin='lower', interpolation='hanning')
im2 = ax2.imshow(cov_2g_raw_prob_AB.T, cmap=colormap, aspect='auto', origin='lower', interpolation='hanning')
#ax.set_ylabel("atom")
# ax.set_xlabel("Q")
# ax1.set_xlabel("Q")
# ax2.set_xlabel("Q")

#fig.suptitle(title, fontsize=16)

fig.text(0.5, 0.04, "Q$_{bind}$", ha='center', va='center', fontsize=16)
fig.text(0.06, 0.5, "atom", ha='center', va='center', rotation='vertical', fontsize=16)


plt.sca(ax)
contacts = cov_1_contacts_AB
atoms = cov_1_atoms_AB

xvalues = np.divide(contacts, contacts.max())
x_ax_idx, x_ax_ticks = customize_axis_ticks(xvalues, num_xticks)
x_ax_ticks_f = test_normalized(x_ax_ticks)
plt.xticks(x_ax_idx, x_ax_ticks_f)

yvalues = atoms
y_ax_idx, y_ax_ticks = customize_axis_ticks(yvalues, num_yticks)
y_ax_ticks_f = test_normalized(y_ax_ticks)
plt.yticks(y_ax_idx, y_ax_ticks_f)

plt.sca(ax1)
contacts = cov_2_contacts_AB
atoms = cov_2_atoms_AB

xvalues = np.divide(contacts, contacts.max())
x_ax_idx, x_ax_ticks = customize_axis_ticks(xvalues, num_xticks)
x_ax_ticks_f = test_normalized(x_ax_ticks)
plt.xticks(x_ax_idx, x_ax_ticks_f)

yvalues = atoms
y_ax_idx, y_ax_ticks = customize_axis_ticks(yvalues, num_yticks)
y_ax_ticks_f = test_normalized(y_ax_ticks)
plt.yticks(y_ax_idx, y_ax_ticks_f)


plt.sca(ax2)
contacts = cov_2g_contacts_AB
atoms = cov_2g_atoms_AB

xvalues = np.divide(contacts, contacts.max())
x_ax_idx, x_ax_ticks = customize_axis_ticks(xvalues, num_xticks)
x_ax_ticks_f = test_normalized(x_ax_ticks)
plt.xticks(x_ax_idx, x_ax_ticks_f)

yvalues = atoms
y_ax_idx, y_ax_ticks = customize_axis_ticks(yvalues, num_yticks)
y_ax_ticks_f = test_normalized(y_ax_ticks)
plt.yticks(y_ax_idx, y_ax_ticks_f)

ip = InsetPosition(ax2, [1.05,0,0.05,1])
cax.set_axes_locator(ip)
fig.colorbar(im, cax=cax, ax=[ax,ax1,ax2], orientation='vertical', label=colorbar_label)

# fig.tight_layout(pad=0.02, h_pad = 0.03, w_pad = 0.03, rect = [0.08, 0.0146, 0.0985, 9])

plt.savefig("plot_Q_bind_3_v1.png", format="png", dpi=500)

plt.show()
