# -*- coding: utf-8 -*-
"""
Created on Wed Jan 31 19:16:06 2018

@author: rfetick
"""

from numpy import zeros, meshgrid, arange, sqrt, fft, random, where
import matplotlib.pyplot as plt

class ShackHart(object):
    """
    ShackHart class
    ###   ATTRIBUTES   ###
    Npix : Number of pixels for the Shack-Hart detector
    Nsp : Linear number of sub-pupils
    NpixSP : Number of pixels per sub-pupil
    wavefront : 
    RON : Read-out noise
    photonNoise : Enable photon noise
    MLA : List of all the micro-lenses
    """
    
    def __init__(self,Npix,Nsp,RON=0,photonNoise=False):
        """
        Constructor
        Npix : Number of pixels for your Shack-Hartman detector
        Nsp : Linear number of sub-pupils in your Shack-Hart
        RON : Eventual read-out noise
        photonNoise : Enable photon noise
        """
        
        self.Npix = int(Npix)
        self.Nsp = int(Nsp)
        self.NpixSP = int(Npix/Nsp)
        self.wavefront = zeros((Npix,Npix))
        self.RON = RON
        self.photonNoise = photonNoise
        #Create an array of microlenses
        self.MLA=[]
        count = 0
        
        if self.NpixSP*self.Nsp != self.Npix:
            print "Npix is not a multiple of Nsp: partial microlenses are ignored"
            
        for i in arange(0, Nsp): #from bottom to top
            for j in arange(0, Nsp): #from left to right
                ml = Microlens(i,j,self.NpixSP,valid=True,ID=count)
                count+=1
                #define relative position of microlenses
                if j!=0:
                    ml.bot = self.MLA[len(self.MLA)-1]
                    ml.bot.top = ml
                if i!=0:
                    ml.left = self.MLA[len(self.MLA)-Nsp]
                    ml.left.right = ml
                #append microlens to array of microlenses
                self.MLA.append(ml)
        
    
    def __repr__(self):
        return "Shack-Hartmann (size="+str(self.Npix)+") (Nsp="+str(self.Nsp)+") (pix per SP="+str(self.NpixSP)+")"
    
    def acquire(self,WF=None):
        """
        Acquire a wavefront onto the detector
        """
        if WF != None:
            self.wavefront = WF
        PSF = zeros((self.Npix,self.Npix))
        for ml in self.MLA:
            # if amplitude of wavefront is nul on microlens, it is not a valid one
            ml.valid = (ml.proj(self.wavefront)).all()
            ml.acquire(self.wavefront)
            PSF[ml.XX,ml.YY] = ml.field
            if self.photonNoise:
                PSF[ml.XX,ml.YY] = random.poisson(lam=PSF[ml.XX,ml.YY],size=(self.NpixSP,self.NpixSP))
        PSF += self.RON*random.randn(self.Npix,self.Npix)
        
        return PSF
    
    def getSlopes(self):
        """
        Find the position of center
        """
        self.acquire()
        for ml in self.MLA:
            ml.valid = (ml.proj(self.wavefront)).all()            
            if ml.valid:
                ml.getSlopes(self.wavefront)
      
    def estimateWF(self):
        """
        Il faut aussi assurer la continuité du piston!
        A faire!!!
        """
        self.getSlopes()
        y,x = meshgrid(arange(self.NpixSP)-self.NpixSP/2,arange(self.NpixSP)-self.NpixSP/2)
        WF = zeros((self.Npix,self.Npix))
        for ml in self.MLA:
            if ml.valid:
                ### Method to retrieve piston
                if ml.bot!=None:
                    ml.piston = ml.bot.piston + ml.bot.tty*self.NpixSP/2 + ml.tty*self.NpixSP/2
                ### End method
                WF[ml.XX,ml.YY] = ml.piston + ml.ttx*x+ml.tty*y
        ### Horrible method, just for visu
        for i in arange(1,self.Nsp):
            WF[i*self.NpixSP:(i+1)*self.NpixSP,:] += WF[i*self.NpixSP-1,self.Npix/2]-WF[i*self.NpixSP,self.Npix/2]
        
        ### Remove average piston
        av = WF[where(self.getValidMask())]
        WF -= av.mean()
        
        return WF
                    
    def getValidMask(self):
        """
        Mask of valid micro-lenses
        """
        mask = zeros((self.Npix,self.Npix))
        for ml in self.MLA:
            if ml.valid:
                mask[ml.XX,ml.YY] = 1.
        return mask
    
    def look(self,func=None,cmap='inferno'):
        """
        Give a quick look on your Shack-Hartmann
        'func' defines the color scaling
        """
        sh = self.acquire()
        pixSP = self.Npix/self.Nsp
        # Plot lines to delimitate microlenses
        for i in arange(0,self.Npix,pixSP):
                sh[i,:] = sh.max()/10.
                sh[:,i] = sh.max()/10.
        # Choice of color scaling
        if func==None:
            plt.pcolormesh(sh.T)
        else:
            plt.pcolormesh(func(sh.T-sh.min()+1.))
        plt.set_cmap(cmap)
        plt.title('Shack-Hartmann detector')
        plt.axis('equal')
        plt.xlabel('X pixels',fontsize=15)
        plt.ylabel('Y pixels',fontsize=15)
        
    def lookNeighbours(self,ID,cmap='inferno'):
        """
        For debugging
        """
        a = zeros((self.Npix,self.Npix))
        IDleft = ''
        IDright = ''
        IDtop = ''
        IDbot = ''
        for ml in self.MLA:
            if ml.ID == ID:
                a[ml.XX,ml.YY] = 5
                if ml.top != None:
                    a[ml.top.XX,ml.top.YY] = 1
                    IDtop = ' #'+str(ml.top.ID)
                if ml.bot != None:
                    a[ml.bot.XX,ml.bot.YY] = 2
                    IDbot = ' #'+str(ml.bot.ID)
                if ml.left != None:
                    a[ml.left.XX,ml.left.YY] = 3
                    IDleft = ' #'+str(ml.left.ID)
                if ml.right != None:
                    a[ml.right.XX,ml.right.YY] = 4
                    IDright = ' #'+str(ml.right.ID)
        cax = plt.pcolormesh(a.T)
        plt.set_cmap(cmap)
        plt.title('Microlens number '+str(ID)+' and its neighbours')
        plt.axis('equal')
        plt.xlabel('X pixels',fontsize=15)
        plt.ylabel('Y pixels',fontsize=15)
        cbar = plt.colorbar(cax,ticks=[0, 1, 2, 3, 4, 5])
        cbar.ax.set_yticklabels(['Other', 'Top'+IDtop, 'Bottom'+IDbot, 'Left'+IDleft, 'Right'+IDright, 'Microlens #'+str(ID)])
        

#############################################################################
#############################################################################
#############################################################################
   
class Microlens(object):
    """
    Microlens class
    ###   ATTRIBUTES   ###
    col : column number of this microlens in the Shack-Hart
    row : row number of this microlens in the Shack-Hart
    X : position of the lower left corner of this microlens in the Shack-Hart
    Y : position of the lower left corner of this microlens in the Shack-Hart
    XX : meshgrid of this microlens
    YY : meshgrid of this microlens
    Npix : number of pixels for the micro-lens
    valid : enable/disable micro-lens
    ID : a unique ID for the micro-lens
    cX : position of the center of this microlens in the Shack-Hart
    cY : position of the center of this microlens in the Shack-Hart
    field : the image on the detector
    mask : a mask of 1 and 0 for the microlens
    left : the micro-lens situated on the left of this one
    right : the micro-lens situated on the right of this one
    top : the micro-lens situated at the top of this one
    bot : the micro-lens situated ot the bottom of this one
    ttx : tip-tilt on X of the incoming wavefront (in unit of pixels displacement)
    tty : tip-tilt on Y of the incoming wavefront (in unit of pixels displacement)
    piston : usually a Shack-Hart doesn't detect piston
             however this attribute is useful for wavefront reconstruction
    """
    
    
    def __init__(self,row,col,Npix,valid=True,ID=-1):
        """
        row : row index
        col : column index
        Npix : size in pixels
        valid : boolean
        ID : a unique number defining this microlens
        """
        self.col = col
        self.row = row
        X = Npix*col
        Y = Npix*row
        self.X = X #position of the lower left corner in the big array
        self.Y = Y
        YY,XX = meshgrid(arange(X,X+Npix),arange(Y,Y+Npix))
        self.XX = XX #array of indices
        self.YY = YY #array of indices
        self.Npix = Npix
        self.valid = valid
        self.ID = ID #unique ID number
        self.cX = X+Npix/2.
        self.cY = Y+Npix/2.
        self.field = zeros((Npix,Npix))
        yy,xx = meshgrid(arange(Npix)-Npix/2,arange(Npix)-Npix/2)
        r = sqrt((xx**2) + (yy**2))/Npix*2.
        self.mask = r < 1.
        # Reference to neighbours
        self.left = None
        self.right = None
        self.top = None
        self.bot = None
        # Tip-tilt and piston for phase reconstruction
        self.ttx = 0
        self.tty = 0
        self.piston = 0
    
    def __repr__(self):
        return "Microlens (ID="+str(self.ID)+") (status="+str(self.valid)+") (size="+str(self.Npix)+"pix) (row="+str(self.row)+", col="+str(self.col)+")"
    
    def proj(self,phase):
        """
        Project a big array onto the microlens
        """
        s = phase.shape
        if self.X+self.Npix>s[0]:
            raise ValueError("Microlens exceeds the X dimension of array to be projected on")
        if self.Y+self.Npix>s[1]:
            raise ValueError("Microlens exceeds the Y dimension of array to be projected on")
        return phase[self.XX,self.YY]
    
    def acquire(self,wf):
        """
        Project a complex wavefront onto the lens detector surface
        """
        if len(wf)>self.Npix:
            #user passed the big array
            wf = self.proj(wf)
        if self.valid:
            #Fourier transform, normalized to conserve energy according to Parseval
            self.field = (abs(fft.fftshift(fft.fft2(self.mask*wf)))**2)/(self.Npix**2)
      
    def getSlopes(self,wf):
        """
        Get the local slope of the wavefront, in unit of pixels displacement
        The real wavefront slope may include lambda and the focal length
        """
        if len(wf)>self.Npix:
            #user passed the big array
            wf = self.proj(wf)
        self.acquire(wf)
        x=arange(self.Npix)
        self.ttx = ((x*self.field.sum(axis=1)).sum())/(self.field.sum()) - self.Npix/2.
        self.tty = ((x*self.field.sum(axis=0)).sum())/(self.field.sum()) - self.Npix/2.
        