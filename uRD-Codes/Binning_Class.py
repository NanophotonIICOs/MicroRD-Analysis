import sys
import os
import numpy as np
import tables
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib import cm

import WinspecUtils
import glob
from tqdm.auto import tqdm, trange
from IPython.display import display


class Binning( ):


    def im2area (file):
        #specific dimensions of VersArray CCD
        rows  = 512
        cols  = 512
        file  = WinspecUtils.readSpe(file)
        spedata = file['data'][0]
        return spedata


       
    def read_spe(spefiles):
        #specific dimensions of VersArray CCD
        rows  = 512
        cols  = 512
        lenfiles = len(spefiles)
        spedata = np.zeros((lenfiles,rows,cols)) 
        for i in tqdm(range(lenfiles),desc = "Load SPE Files"):
            files = WinspecUtils.readSpe(spefiles[i])
            spedata[i,:,:] = files['data'][0]
        return spedata

    # image binning
    def binning(binsize,image):
        while (image.shape[1]%binsize) != 0 and (image.shape[0]%binsize) != 0 :
            try :
                print("This Binning is incorrect, check this.")
                break
            except ValueError:
                 print("This Binning is incorrect, check this.")
        
        b  = binsize
        n = int(image.shape[0]/b) 
        m = int(image.shape[1]/b) 
        imagebin = np.zeros((n,m))
        for i in range(n):
            for j in range(m):
                suma = 0
                for k in range((i+1)*b-b,(i*b)+b):
                    for l in range((j+1)*b-b,(j*b)+b):
                        suma = suma + image[k,l]
                imagebin[i,j] = suma/(b*b)
                    
        return imagebin

    


class Image_Binning():
    def __init__(self,pathfiles,spefiles,binsize,xi,xf,yi,yf,l,li,lf,step):
        self.path     = pathfiles
        self.spefiles = spefiles
        self.binsize  = binsize
        # Area Parameters
        self.xi = xi
        self.xf = xf
        self.yi = yi
        self.yf = yf
        self.l  = l    # Select image that corresponds to Wavelength
        self.li = li   # Initial Wavelength in the experiment
        self.lf = lf   # Final Wavelength in the experiment
        self.expfiles = glob.glob(self.path+self.spefiles+'*.SPE')
        self.expfiles = sorted(self.expfiles)
        
    

        while (self.l < li) or (self.l > lf):
            try :
                print("Error in Image Wavlength ")
                self.l = int(input("Image Wavlength: "))
            except ValueError:
                print("Error in Image Wavlength ")
                exit

        
        self.lvec  = np.linspace(li,lf,len(self.expfiles))



        self.contl = 0
        for i in self.lvec:
            if i==self.l:
                break
            else:
                self.contl+=1
            

        self.expdata0 = Binning.im2area(self.expfiles[self.contl])




        self.lenspefiles = len(self.expfiles)


        self.image0 = self.expdata0
        self.n = self.image0[self.yi:self.yi+self.yf, self.xi:self.xi+self.xf].shape[0]
        self.m = self.image0[self.yi:self.yi+self.yf, self.xi:self.xi+self.xf].shape[1]
        self.image0out = self.expdata0[self.yi:self.yi+self.yf, self.xi:self.xi+self.xf]

        self.nb = int(self.n/self.binsize)
        self.mb = int(self.m/self.binsize)








    def draw_area(self):
            figure, (ax1,ax2) = plt.subplots(1,2,figsize=(10,5))
            ax1.imshow(self.image0,cmap=cm.viridis)
            rect = patches.Rectangle((self.xi,self.yi),self.xf,self.yf, edgecolor='k',lw=3, facecolor="none")
            ax1.add_patch(rect)
            ax1.set_title("Original Image")

            ax2.imshow(self.image0[self.yi:self.yi+self.yf, self.xi:self.xi+self.xf],cmap=cm.viridis)
            ax2.set_title("Area")
            figure.tight_layout()
            plt.show()

            

    def Binning(self,spefiles,bins):
        self.data2bin = np.zeros((self.lenspefiles,self.nb,self.mb))

        for bim in tqdm(range(self.lenspefiles),desc= "IMAGE BINNING"):
            self.data2bin[bim,:,:] = Binning.binning(bins,spefiles[bim,self.yi:self.yi+self.yf,self.xi:self.xi+self.xf])

        return self.data2bin

    def Array2Plot(self,spefiles,bins,li,lf):
        
        self.rd = self.Binning(spefiles,bins)
        self.nb = int(self.n/bins)
        self.mb = int(self.m/bins)
        self.lendim   = self.nb*self.mb
        self.nfiles   = self.rd.shape[0]
        self.rd_vec   = np.zeros((self.nfiles,self.lendim))
        self.bind     = np.zeros((self.nfiles,self.lendim))

        # Flaten Binning array
        self.rdfl   = self.rd.flatten()

        cont1  = 0
        for i in tqdm(range(self.nfiles),desc="Rearrange array"):
            for j in range(self.lendim):
                self.rd_vec[i,j] = self.rdfl[cont1]
                cont1 += 1

        #CREATE CELL
        # prcell = {}
        # for i in range(n):
        #     for j in range(m):
        #         prcell[i, j] = np.zeros((nfiles, 2))
        #Fill cell 
        self.col = self.nb
        self.row = self.mb
        self.lvec  = np.linspace(li,lf,self.nfiles)
        cont_cell = 0
        mdumm = np.zeros((self.nfiles,2))





        self.rdcell = {}
        for i in tqdm(range(self.col),desc="CELL ARRAY"):
            for j in range(self.row):
                for k in range(self.nfiles):
                #Binnign per wavelength
                    ddum = self.rd_vec[k,cont_cell]
                    mdumm[k,:] = [self.lvec[k] , ddum]
            
                # create cell
                self.rdcell[j,i] = mdumm
        
                # Counters
                cont_cell+=1
                mdumm = np.zeros((self.nfiles,2))
        
        return {"cell" : self.rdcell,
                "array": self.rd[self.contl],
                "binning": self.rd}

    def execute(self):
        self.expdata = Binning.read_spe(self.expfiles)
        self.results = self.Array2Plot(self.expdata,self.binsize,self.li,self.lf)

        return {"cell": self.results["cell"],
                "array":self.results["array"],
                 "binning":self.results["binning"]}


    
        



