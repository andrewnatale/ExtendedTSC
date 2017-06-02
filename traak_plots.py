import sys, os, math
import numpy as np
import matplotlib.pyplot as plt

# def convert_dihedral_array(rad_array):
#     # convert to degrees
#     deg_array = np.rad2deg(rad_array)
#     # wrap values >180 and <-180
#     too_high = deg_array > 180
#     deg_array[too_high] -= 360
#     too_low = deg_array < -180
#     deg_array[too_low] += 360
#     return deg_array

# def timeseries_stack_plot(datalist1,datalist2,time_array,ylimits=False):
#     fig = plt.subplot(211)
#     colors = ['red', 'blue', 'green', 'yellow']
#     for idx,elem in enumerate(datalist1):
#         fig.plot(time_array, elem, c=colors[idx])
#     fig.set_xlim([0,time_array[-1]])
#     if ylimits:
#         fig.set_ylim(ylimits)
#     fig = plt.subplot(212)
#     for idx,elem in enumerate(datalist2):
#         fig.plot(time_array, elem, c=colors[idx])
#     fig.set_xlim([0,time_array[-1]])
#     if ylimits:
#         fig.set_ylim(ylimits)
#     plt.show()

# def angle_space_plot(datalist):
#     fig = plt.subplot(111)
#     colors = ['red', 'blue', 'green', 'yellow']
#     for idx,elem in enumerate(datalist):
#         fig.scatter(convert_dihedral_array(elem[0]), convert_dihedral_array(elem[1]), c=colors[idx])
#     fig.set_xlim([-180,180])
#     fig.set_ylim([-180,180])
#     plt.show()

def dihedral_AB_space(plotname,datalist,save=False):
    f, ax = plt.subplots()
    colors = ['gold', 'purple']
    for idx,elem in enumerate(datalist):
        ax.scatter(convert_dihedral_array(elem[0]), convert_dihedral_array(elem[1]), c=colors[idx])
    ax.set_xlim([-180,180])
    ax.set_ylim([-180,180])
    ax.set_xlabel('Chi1')
    ax.set_ylabel('Chi2')
    ax.set_title(plotname)
    ax.text(0.5,0.05,'W262, Subunit A',color='gold',transform=ax.transAxes,verticalalignment='bottom',horizontalalignment='center')
    ax.text(0.5,0.01,'W262, Subunit B',color='purple',transform=ax.transAxes,verticalalignment='bottom',horizontalalignment='center')
    if save:
        plt.savefig(plotname+'.dihedral_space.png', bbox_inches='tight')
    else:
        plt.show()

def timeseries_basic(plotname,data,time,ylimits=False,save=False):
    f, ax = plt.subplots()
    colors = ['red', 'blue', 'green', 'purple', 'silver', 'gold']
    for idx,elem in enumerate(data):
        ax.plot(time, elem[1], c=colors[idx])
        ax.text(0.99,0.99-0.08*idx, elem[0], color=colors[idx],
                transform=ax.transAxes,verticalalignment='top',horizontalalignment='right')
    ax.set_xlim([time[0],time[-1]])
    if ylimits:
        ax.set_ylim(ylimits)
    ax.set_title(plotname)
    ax.set_xlabel('time (ns)')
    ax.set_ylabel('Relative indole ring position (angstroms)')
    if save:
        plt.savefig(plotname+'.W262_rotamers.png', bbox_inches='tight')
    else:
        plt.show()

def timeseries_AB_stack(plotname,data,time,ylimits=False,save=False):
    # data should be a list of tuples, each with 3 elements:
    # [('name1', arr1A, arr1B), ('name2', arr2A, arr2B), ...]
    f, (ax1, ax2) = plt.subplots(2, sharex=True, sharey=True)
    colors = ['red', 'blue', 'green', 'purple', 'silver', 'gold']
    for idx,elem in enumerate(data):
        ax1.plot(time, elem[1], c=colors[idx])
        ax2.plot(time, elem[2], c=colors[idx])
        ax1.text(0.99,0.99-0.1*idx, elem[0], color=colors[idx],
                 transform=ax1.transAxes,verticalalignment='top',horizontalalignment='right')
    f.subplots_adjust(hspace=.1)
    if ylimits:
        ax1.set_ylim(ylimits)
        ax2.set_ylim(ylimits)
    ax1.set_xlim([time[0],time[-1]])
    ax1.set_title(plotname)
    # hardcoded labels
    ax2.set_xlabel('time (ns)')
    ax1.set_ylabel('Subunit A distance')
    ax2.set_ylabel('Subunit B distance')
    ax2.set_xlim([time[0],time[-1]])
    plt.setp([a.get_xticklabels() for a in f.axes[:-1]], visible=False)
    ax1.grid(True)
    ax2.grid(True)
    gridlines = ax1.get_xgridlines() + ax1.get_ygridlines() + ax2.get_xgridlines() + ax2.get_ygridlines()
    for line in gridlines:
        line.set_linestyle('-.')
    if save:
        plt.savefig(plotname+'.ts_stack.png', bbox_inches='tight')
    else:
        plt.show()

def timeseries_dual_yscale(plotname,data_y1,data_y2,time,y1_limits=False,y2_limits=False,save=False):
    # data should be a list of tuples, each with 2 elements:
    # [('name1', arr1), ('name2', arr2), ...]
    f, ax1 = plt.subplots()
    colors = ['red', 'blue', 'purple', 'gold']
    labelcount = 0
    for idx,elem in enumerate(data_y1):
        ax1.plot(time, elem[1], c=colors[labelcount])
        ax1.text(0.01,0.99-0.1*idx, elem[0], color=colors[labelcount],
                 transform=ax1.transAxes,verticalalignment='top',horizontalalignment='left')
        labelcount += 1
    ax1.set_xlim([time[0],time[-1]])
    ax2 = ax1.twinx()
    for idx,elem in enumerate(data_y2):
        ax2.plot(time, elem[1], c=colors[labelcount])
        ax1.text(0.99,0.99-0.1*idx, elem[0], color=colors[labelcount],
                 transform=ax1.transAxes,verticalalignment='top',horizontalalignment='right')
        labelcount += 1
    if y1_limits:
        ax1.set_ylim(y1_limits)
    if y2_limits:
        ax2.set_ylim(y2_limits)
    ax1.set_title(plotname)
    # hardcoded labels
    ax1.set_xlabel('time (ns)')
    ax1.set_ylabel('Distance')
    ax2.set_ylabel('Rotamer')
    ax2.set_xlim([time[0],time[-1]])
    if save:
        plt.savefig(plotname+'.ts_dualy.png', bbox_inches='tight')
    else:
        plt.show()
