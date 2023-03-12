#plot the pv files for all models
import numpy as np
from matplotlib import pyplot as plt
import matplotlib.gridspec as gridspec
from astropy.io import fits
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.wcs import WCS
from pvextractor import extract_pv_slice
from pvextractor.geometry import Path
DTR = np.pi/180.

workpath = "/data/users/li/extragas_grid/"
cat_filename = workpath+"catalogue.txt"
from os import walk

f = []
#find all the model names
for (dirpath, dirnames, filenames) in walk("../../model_cubes/"):
    f.extend(filenames)
    break


cat = np.loadtxt(cat_filename, dtype={'names': ('name','database','RA','DEC','VSYS','D','INC','PA','noise'),'formats':(np.chararray,np.chararray,np.chararray,np.chararray,np.float,np.float,np.float,np.float,np.float)},skiprows=1)

class pv:
    def __init__(self,filename,RA,DEC,PA,VSYS):
        self.hdu = fits.open(filename)[0]
        c = SkyCoord(RA,DEC,frame='icrs')
        if(self.hdu.header['OBJECT']=='ngc2403'): #deal with VLA header. Uff.
            self.vel_conv = 1.
            self.hdu.header['CUNIT1']='DEGREE'
            self.hdu.header['CUNIT2']='DEGREE'
            self.hdu.header['CUNIT3']='km/s'
            wmap = WCS(self.hdu.header)
            wmap_channel = wmap.dropaxis(3)
            wmap_channel = wmap_channel.dropaxis(2)
        elif(self.hdu.header['OBJECT']=='ngc0949'):
            self.vel_conv = 1e3
            wmap = WCS(self.hdu.header)
            wmap_channel = wmap.dropaxis(3)
            wmap_channel = wmap_channel.dropaxis(2)
        else:
            self.vel_conv = 1e3
            wmap = WCS(self.hdu.header)
            if wmap.naxis==4:
                wmap_channel = wmap.dropaxis(3)
                wmap_channel = wmap_channel.dropaxis(2)
            if wmap.naxis==3:
                wmap_channel = wmap.dropaxis(2)
        self.xc,self.yc = wmap_channel.wcs_world2pix(c.ra.degree,c.dec.degree,0) 
        self.PA = PA
        self.VSYS = VSYS
        self.arcmin_to_pix = 1./60./np.fabs(self.hdu.header['CDELT1'])
    def load_mask(self,filename):
        mask = fits.getdata(filename)
        ind = np.where((np.isnan(mask)))
        self.hdu.data[ind] = np.nan
    def slice_maj(self,len_arcmin,offset_arcmin,thickness_pix):
        sinPA = np.sin(self.PA*DTR)
        cosPA = np.cos(self.PA*DTR)
        len_pix = len_arcmin*self.arcmin_to_pix
        offset_pix = offset_arcmin*self.arcmin_to_pix
        x0 = self.xc + offset_pix*cosPA
        y0 = self.yc + offset_pix*sinPA 
        x1 = x0 - len_pix*sinPA
        y1 = y0 + len_pix*cosPA
        x2 = x0 + len_pix*sinPA
        y2 = y0 - len_pix*cosPA
        xy = Path([(x1,y1),(x2,y2)],width=thickness_pix)
        return extract_pv_slice(self.hdu, xy).data
    def slice_min(self,len_arcmin,offset_arcmin,thickness_pix):
        sinPA = np.sin(self.PA*DTR)
        cosPA = np.cos(self.PA*DTR)
        len_pix = len_arcmin*self.arcmin_to_pix
        offset_pix = offset_arcmin*self.arcmin_to_pix
        x0 = self.xc - offset_pix*sinPA
        y0 = self.yc + offset_pix*cosPA
        x1 = x0 - len_pix*cosPA
        y1 = y0 - len_pix*sinPA
        x2 = x0 + len_pix*cosPA
        y2 = y0 + len_pix*sinPA
        xy = Path([(x1,y1),(x2,y2)],width=thickness_pix)
        return extract_pv_slice(self.hdu, xy).data
        
        
#GALAXY PROPERTIES
galaxy = "NGC2403"
model_type='fountain + corona accretion'
mock_data = False
overlay_model = True
isys = np.where(cat['name']==galaxy)[0][0]
#galaxy global properties
RA  = cat['RA'][isys]
DEC = cat['DEC'][isys]
Distance  = cat['D'][isys]           #Mpc
rms_noise = cat['noise'][isys]
#rms_noise=0.0023
#read other properties from BB file
BB_file = workpath +galaxy+"_rotcur.txt"
BB = np.loadtxt(BB_file, dtype={'names': ('RAD_kpc','RAD_arcsec','VROT','DISP','INC','PA','Z0_pc','Z0_arcsec','SIG','XPOS','YPOS','VSYS','VRAD'),'formats':(np.float,np.float,np.float,np.float,np.float,np.float,np.float,np.float,np.float,np.float,np.float,np.float,np.float)},skiprows=1)
VSYS = np.median(BB['VSYS'])
INC  = np.median(BB['INC'])
PA   = np.median(BB['PA'])
#rms_noise = 1.55e-4 # usato per NGC 5055

#SLICE PROPERTIES
length_arcmin_maj = 17
length_arcmin_min = 0.7*length_arcmin_maj
offset_arcmin = [-4,-2,0,2,4]
width_pixels = 3

#for ngc2403
#length_arcmin_maj = 25
#length_arcmin_min = 0.7*length_arcmin_maj
#offset_arcmin = [-4,-2,0,2,4]
#width_pixels = 3

conto = (2.0*rms_noise)*(2**np.arange(20))
conto = (2.0*rms_noise)*(2**np.arange(15))
conto_data = (2.0*rms_noise)*(2**np.arange(15))
conto_neg = -conto[0]
if(mock_data):
    obs_file  = workpath+"Halogas/datacubes/observations/"+galaxy+"_mock.fits"
    mask_file = workpath+"Halogas/datacubes/observations/"+galaxy+"_mock_mask.fits.gz"

else:
    for mod_file in f:
	# pv filenames
        pv_file=("../figures/pv_"+mod_file).replace(".fits",".pdf")
        mod_file="../../model_cubes/"+mod_file
        obs_file  = workpath+"cubes/"+galaxy+".fits"
    #obs_file='datacubes/observations/.fits'
        mask_file = workpath+"cubes/"+galaxy+"_mask.fits"
        rratio=1.0
        ftt=fits.open(mod_file)
        hhh=fits.getheader(mod_file)
        hhh['OBJECT']='NGC2403'
        hhh['CUNIT3']='m/s'
        hhh['CRVAL3']=-5152.86
        ftt[0].header=hhh
        ftt.writeto(mod_file,overwrite='True')
        #mod_file='datacubes/models/mk_halo.fits'
    
        obs_mypv = pv(obs_file,RA,DEC,PA,VSYS)
        obs_mypv.load_mask(mask_file)
        mask_mypv=pv(mask_file,RA,DEC,PA,VSYS)
    
        if(overlay_model): 
            mod_mypv = pv(mod_file,RA,DEC,PA,VSYS)
            mod_mypv.load_mask(mask_file)
        vel_conv = obs_mypv.vel_conv
        myvmax = np.nanpercentile(obs_mypv.hdu.data[np.where(obs_mypv.hdu.data>3*rms_noise)],95.0)
        
        fig,axes = plt.subplots(figsize=(16,6),sharey=    True)
        grid = [gridspec.GridSpec(1,len(offset_arcmin)),gridspec.GridSpec(1,len(offset_arcmin))]
        grid[0].update(top=0.95, bottom=0.48, left=0.06, right=0.97, wspace=0.0, hspace=0.0)
        grid[1].update(top=0.43, bottom=0.08, left=0.06, right=0.97, wspace=0.0, hspace=0.0)
        
        cv3 = obs_mypv.hdu.header['CRVAL3']
        cd3 = obs_mypv.hdu.header['CDELT3']
        cp3 = obs_mypv.hdu.header['CRPIX3']
        n3  = obs_mypv.hdu.header['NAXIS3']
        vini = cv3  + (1-cp3)*cd3
        vfin = vini + (n3-1)*cd3
        extent_maj = [-length_arcmin_maj,length_arcmin_maj,vini/vel_conv-VSYS,vfin/vel_conv-VSYS]
        extent_min = [-length_arcmin_min,length_arcmin_min,vini/vel_conv-VSYS,vfin/vel_conv-VSYS]
        
        #major axes pv
        for i in range(len(offset_arcmin)):
            ax = plt.subplot(grid[0][i])
            if i==0:
                ax.tick_params(axis='both',which='both',bottom=True,top=False,left=True,right=False,labelbottom=True,labelleft=True,labelsize=20)
                ax.set_ylabel("V$_{\\rm HEL}$-V$_{\\rm SYS}$ [km/s]",fontsize=20)
            else:
                ax.tick_params(axis='both',which='both',bottom=True,top=False,left=True,right=False,labelbottom=True,labelleft=False,labelsize=20)
            if i==(len(offset_arcmin)-1)/2:
                ax.text(0.92,0.92,"P.A.="+np.str(np.int(PA+0.5))+"$^\circ$",fontsize=13,ha='right',transform=ax.transAxes)
                if(mock_data):
                    ax.text(0.5,1.04,"mock data",fontsize=14,ha='center',transform=ax.transAxes)
                else:
                    ax.text(0.5,1.04,model_type,fontsize=14,ha='center',transform=ax.transAxes)
            obs_pv = obs_mypv.slice_maj(length_arcmin_maj,offset_arcmin[i],width_pixels)
            if(overlay_model): mod_pv = mod_mypv.slice_maj(length_arcmin_maj,offset_arcmin[i],width_pixels)
            ax.imshow(obs_pv,cmap="Blues",extent=extent_maj,origin='lower',vmin=conto_neg,vmax=myvmax,aspect='auto')
            ax.contour(obs_pv,conto_data,origin='lower',extent=extent_maj,linewidths=1,colors='black')
            ax.contour(obs_pv,[conto_neg],origin='lower',extent=extent_maj,linewidths=1,colors='grey')
            mask_pv=mask_mypv.slice_maj(length_arcmin_maj,offset_arcmin[i],width_pixels)
            #ax.contour(mask_pv,[conto_neg],origin='lower',extent=extent_maj,linewidths=1,colors='green',levels=[0])
            if(overlay_model&(np.nanmax(mod_pv)>conto[0])): ax.contour(mod_pv*rratio,conto,origin='lower',extent=extent_maj,extend='neither',linewidths=1.5,colors='red')
            ax.axhline(0,ls=':',color='k',lw=0.75)
            ax.axvline(0,ls=':',color='k',lw=0.75)
            ax.text(0.07,0.92,np.str(offset_arcmin[i])+"$^\prime$",fontsize=12,ha='left',transform=ax.transAxes)
        
        for i in range(len(offset_arcmin)):
            ax = plt.subplot(grid[1][i])
            if i==0:
                ax.tick_params(axis='both',which='both',bottom=True,top=False,left=True,right=False,labelbottom=True,labelleft=True,labelsize=20)
                ax.set_ylabel("V$_{\\rm HEL}$-V$_{\\rm SYS}$ [km/s]",fontsize=20)
            else:
                ax.tick_params(axis='both',which='both',bottom=True,top=False,left=True,right=False,labelbottom=True,labelleft=False,labelsize=20)
                
            if i==(len(offset_arcmin)-1)/2:
                ax.set_xlabel("offset [$^\prime$]",fontsize=20)    
            obs_pv = obs_mypv.slice_min(length_arcmin_min,offset_arcmin[i],width_pixels)
            mask_pv=mask_mypv.slice_min(length_arcmin_maj,offset_arcmin[i],width_pixels)
            if(overlay_model): mod_pv = mod_mypv.slice_min(length_arcmin_min,offset_arcmin[i],width_pixels)
            ax.imshow(obs_pv,cmap="Blues",extent=extent_min,origin='lower',vmin=conto_neg,vmax=myvmax,aspect='auto')
            ax.contour(obs_pv,conto_data,origin='lower',extent=extent_min,linewidths=1,colors='black')
            #ax.contour(mask_pv,origin='lower',extent=extent_min,linewidth=1,colors='green',levels=[0])
            ax.contour(obs_pv,[conto_neg],origin='lower',extent=extent_min,linewidths=1,colors='grey')
            if(overlay_model): ax.contour(mod_pv*rratio,conto,origin='lower',extent=extent_min,linewidths=1.,colors='red')
            ax.axhline(0,ls=':',color='k',lw=0.75)
            ax.axvline(0,ls=':',color='k',lw=0.75)
            ax.text(0.07,0.92,np.str(offset_arcmin[i])+"$^\prime$",fontsize=12,ha='center',transform=ax.transAxes)
           
        #plt.tight_layout()
        plt.savefig(pv_file,bbox_inches='tight')
