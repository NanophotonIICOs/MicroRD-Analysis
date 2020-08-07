import glob as glob
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.patches as patches
import sys
from matplotlib import cm
from matplotlib import gridspec
from mpl_toolkits.axes_grid1 import make_axes_locatable
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from ipywidgets import interact, interactive, fixed, interact_manual
import ipywidgets as widgets
from tqdm.auto import tqdm, trange
from IPython.display import display
from IPython.display import HTML
from IPython import display
import matplotlib.animation as animation
from matplotlib.animation import PillowWriter
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.animation
import SPECTRA_ANALYSIS as spa
from Binning_Class import Binning, Image_Binning


#matplotlib options
plt.rcParams['font.family']         = 'Times New Roman'
#plt.rcParams['font.serif']          = 'Times'
plt.rcParams['xtick.labelsize']     = 17
plt.rcParams['ytick.labelsize']     = 17
plt.rcParams['axes.linewidth']      = 2
plt.rcParams["xtick.minor.visible"] =  False
plt.rcParams["xtick.major.size"]    =  10
plt.rcParams["xtick.minor.size"]    =  5
plt.rcParams["xtick.major.width"]   =  2
plt.rcParams["xtick.minor.width"]   =  2
plt.rcParams["xtick.direction"]     =  'in'
plt.rcParams["ytick.minor.visible"] =  False
plt.rcParams["ytick.major.size"]    =  10
plt.rcParams["ytick.minor.size"]    =  5
plt.rcParams["ytick.major.width"]   =  2
plt.rcParams["ytick.minor.width"]   =  2
plt.rcParams["ytick.direction"]     =  'in'
plt.rcParams['text.usetex']         = True
plt.rcParams['legend.frameon']      = False
import WinspecUtils

def mapcell(imagen,datab,xpix,ypix,peak,yoff1,yoff2,**labels):
        fig  = plt.figure(figsize=(15, 9))
        spec = gridspec.GridSpec(ncols=2,nrows=1, width_ratios=[2,1])
        spec.update(wspace=0.1, hspace=0.05) # set the spacing between axes. 
        ax1=fig.add_subplot(spec[0])
        im=ax1.imshow(imagen,cmap=cm.viridis)
        ax1.plot([xpix],[ypix],marker = 's',color='w',markersize='13',markerfacecolor='none')
        ax1.axhline([ypix], ls = ':',color='w', lw=2)
        ax1.axvline([xpix], ls = ':',color='w', lw=2)
        ax1.set_title("%s Binning Size: %s  at %.1f nm" %(labels['Experiment'],labels['binning'],labels['waveselect']), fontsize=25),
        ax1.tick_params(direction='out',which='both')
        ax1.minorticks_on(),
        # SPECTRA PLOT
        ax2=fig.add_subplot(spec[1])
        ax2.set_title("%s CCD " %(labels['Experiment']), fontsize=25),
        ax2.plot(datab[xpix,ypix][:,0],datab[xpix,ypix][:,1],
                 lw = 1,
                 color = 'b')
        ax2.axvline([peak], ls = '--',color='k', lw=1,label="Peak at %.1f nm"%(peak)),
        ax2.set_xlabel('Wavelength [nm]',fontsize = 25)
        ax2.set_ylabel('CCD accounts',fontsize =25)
        ax2.set_ylim([yoff1,yoff2])
        ax2.minorticks_on()
        ax2.ticklabel_format(style='sci', axis='y')
        ax2.yaxis.offsetText.set_fontsize(15)
        ax2.legend(fontsize=18,frameon=False,loc=1)
        ax2.text(0.6,0.9, 'x Pixel:%d'%xpix,color='red', fontsize=18,transform=ax2.transAxes)
        ax2.text(0.6,0.85, 'y Pixel:%d'%ypix,color='red', fontsize=18,transform=ax2.transAxes)
        ax2.yaxis.tick_right()
        ax2.yaxis.set_label_position("right")
        ax2.set_xlim([labels['li']-30,labels['lf']+30])
        axins1 = inset_axes(ax2,
                    width="5%",  # width = 50% of parent_bbox width
                    height="50%",  # height : 5%
                    loc="upper left",
                    borderpad=0)
        plt.colorbar(im, cax=axins1, orientation="vertical")
        axins1.yaxis.tick_left()
        plt.show()

## Create animation
def map_video(bindata,image,ax01,ax02,xpix,ypix,x1,y1,yoffmin,yoffmax,sel,**labels):
    # Show image PR
    ax01.clear()
    im = ax01.imshow(image,cmap=cm.viridis)
    ax01.plot([xpix],[ypix],marker = 's',color='w',markersize='13',markerfacecolor='none')
    ax01.axhline([ypix], ls = ':',color='w', lw=1)
    ax01.axvline([xpix], ls = ':',color='w', lw=1)
    ax01.plot([x1,xpix],[y1,ypix],'red',lw=3)
    ax01.tick_params(direction='out',which='both')
    ax01.minorticks_on()
    ax01.set_title("%s Binning Size: %s at %.1f nm"%(labels['Experiment'],labels['binning'],labels['waveselect']), fontsize=labels['fsTitles'])
    #ax02.text(0.6,0.97, 'x Pixel:%d'%xpix,color='red', fontsize=13,transform=ax02.transAxes)
    #ax02.text(0.6,0.97-0.04, 'y Pixel:%d'%ypix,color='red', fontsize=13,transform=ax02.transAxes)
    #plots
    ax02.clear()   
    ax02.plot(bindata[xpix,ypix][:,0],bindata[xpix,ypix][:,1],
                lw = 1,
                color = 'r',
                label = 'xpix:%d-pix:%d'%(xpix,ypix))
    ax02.set_title("%s CCD " %(labels['Experiment']), fontsize=labels['fsTitles']),
    ax02.set_xlabel('Wavelength (nm)',fontsize=labels['fsTicks'])
    ax02.set_ylabel('CCD accounts ',fontsize=labels['fsTicks'],rotation=-90, labelpad=25)
    ax02.set_ylim([yoffmin,yoffmax])
    #if sel=="div":
    #        ax02.set_ylim([np.min(image)+yoffmin*1e-4,np.max(image)+yoffmax*1e-4])
    #else:
    #        ax02.set_ylim([np.min(image)+yoffmin,np.max(image)+yoffmax])
    ax02.minorticks_on()
    ax02.yaxis.tick_right()
    ax02.yaxis.set_label_position("right")
    ax02.set_xlim([labels['li']-30,labels['lf']+30])
    axins1 = inset_axes(ax02,
                    width="5%",  # width = 50% of parent_bbox width
                    height="50%",  # height : 5%
                    loc="upper left",
                    borderpad=0)
    plt.colorbar(im, cax=axins1, orientation="vertical")
    axins1.yaxis.tick_left()