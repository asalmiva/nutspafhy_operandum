# -*- coding: utf-8 -*-
"""
Created on Thu Oct 19 10:21:51 2017

@author: lauren
"""

import pandas as pd
import numpy as np
import matplotlib.pylab as plt
from scipy import optimize
from scipy.misc import derivative
import pickle
import xarray as xr
#import Spathy_resp as stand
#import Spathy_utils
#from Spathy_utils import read_setup
#from datetime import datetime
import datetime
from netCDF4 import Dataset, date2num, num2date
from scipy.optimize import fsolve, minimize_scalar 
import os
from scipy.interpolate import interp1d
from scipy.optimize import curve_fit
from scipy.interpolate import InterpolatedUnivariateSpline as interS
from spafhy_io import read_AsciiGrid

def get_Motti_results(fm, latS, lonS):
    print ('****** Reading predefined Motti results*****************')    
    print (fm)
    #import sys; sys.exit()
    #motti = pickle.load(open(fm, 'rb'))
    with open(fm, 'rb') as f:
        u = pickle._Unpickler(f)
        u.encoding = 'latin1'
        motti = u.load()
    minDist=1e5; mClose = None    
    for m in motti.keys():
        d = ((motti[m]['Lat'] - latS/1000.)**2 + (motti[m]['Lon'] - lonS/1000.)**2)**0.5
        #print m, np.round(d,0)        
        if  d < minDist: 
            mClose=m; minDist=d
    print ('    + Closest Motti-simulations from ', mClose) 
    print ('    + Distance ', np.round(minDist,0), ' km')
    return motti[mClose]

def get_spafhy_results(ff=None):
    """

    Parameters
    ----------
    ff : TYPE,string 
        DESCRIPTION. path to netCDF4 file containing the spafhy results 

    Returns
    -------
    res : netCDF4 file
        DESCRIPTION. Daily SpaFHy results in gridded time series.

    """    
    #************Get Spathy results for the computation************************    
    print ('****** Reading SpaFHy results *********************')
    print (ff)        
    res=Dataset(ff, mode='r')                                                   # water netCDF, open in reading mode
    return res

def get_nutbal_results(f=None):
    print ('****** Reading NutSpaFHy results *********************')
    print (f)
    nres = Dataset(f, mode='r')
    return nres
    

def get_measured_data(f,fc,start_date,end_date):
    #*********** Get measured load and concentration ****************************    
    print ('****** Reading measured load and concentration **************')
    mit = pd.read_csv(f)
    utctime=pd.to_datetime(mit['utctime'],format='%Y-%m-%d %H:%M:%S')
    mit.index=utctime    
    mit=mit[(mit.index >= start_date) & (mit.index <= end_date)]
    #n_meas = mit['NTOT_g/ha/d'][:-1] 
    #n_meas = n_meas/1000.
    #p_meas = mit['PTOT_g/ha/d'][:-1]    
    #p_meas = p_meas/1000.
    #meastimes = pd.to_datetime(mit['utctime'][:-1])
    #concentration
    mitc = pd.read_csv(fc)
    utctime=pd.to_datetime(mitc['utctime'],format='%Y-%m-%d %H:%M:%S')
    mitc.index=utctime
    mitc=mitc[(mitc.index >= start_date) & (mitc.index <= end_date)]

    #measured concentration here mitc['NTOT_ugL']
    mit = mit.resample('M', convention='end').sum()
    mitc=  mitc.resample('M', convention='end').mean()

    return mit, mitc
    #return n_meas, p_meas, meastimes, mitc    
    
def process_Motti():
    """
    Constructs a dict variable having dimensions Municipality/Site fertility (1...6)/Soiltype (mineral/peat)/species
    picks Motti simulations (saved in xls format) from given folder and locates age, height, and yield to dict named 'motti'.
    Saves lat, lon, and temp sum.
    Fits s-shaped Schumacher model to h(a), and y(a), and saves the parameter values to the dict
    Saves the dict 'motti' to picke file
    """
    #******* Input **********
    folder ='C:\\Apps\\WinPython-64bit-2.7.10.3\\SAMULI_ARI\\Motti\\'
    fi = 'Stands.xlsx'   #contains the list of motti-files and the relevant metadata
    dfStands = pd.read_excel(folder+fi)
    fPickle = folder + 'Stands_bio.p'

    #****** Create dict ***************
    motti ={}    
    
    for m in set(dfStands['Municipality']):
        motti[m]={}
        motti[m]['Lat']=0
        motti[m]['Lon']=0
        motti[m]['Tsum']=0
        for s in set(dfStands['Site']):
            motti[m][s]={}
            for minpe in set(dfStands['Min/Peat']):
                motti[m][s][minpe]={}
                for spe in set(dfStands['Species']):
                    motti[m][s][minpe][spe]={}
                    motti[m][s][minpe][spe]['age']=[]
                    motti[m][s][minpe][spe]['hg']=[]
                    motti[m][s][minpe][spe]['yi']=[]

    #************   Read motti-simulations ***********
    for m, lat, lon, ts, s, minpe, spe, fname in \
        zip(dfStands['Municipality'], dfStands['Lat'], dfStands['Lon'], dfStands['TempSum'], \
        dfStands['Site'], dfStands['Min/Peat'], dfStands['Species'], dfStands['Filename']):    
    
        motti[m]['Lat'] = lat
        motti[m]['Lon'] = lon
        motti[m]['Tsum'] =ts
        dfSite = pd.read_excel(folder + fname, header=None, skiprows=1)

        #column 2: age, 5:hg, 12:Yield 
        age = dfSite[2].values[:-1]; age = np.insert(age,0,0)
        hg = dfSite[5].values[:-1]; hg = np.insert(hg, 0,0)
        yi = dfSite[12].values[:-1]; yi = np.insert(yi,0,0)
        idx = np.argwhere(np.diff(age)>0)    
        age = age[idx].flatten(); hg= hg[idx].flatten(); yi = yi[idx].flatten()
        motti[m][s][minpe][spe]['age']=age
        motti[m][s][minpe][spe]['hg']=hg
        motti[m][s][minpe][spe]['yi']=yi
        
        #*******  Fit Schumacher model****************
        a = age
        h = hg
        y =yi
        p =[0.05, 7.0]    
        idx = len(a)/2
        SI = h[idx]; Iage=a[idx]
        fh = lambda  age,  *p: SI*((1.0-np.exp(-1.0*p[0]*age))/(1.0-np.exp(-1.0*p[0]*Iage)))**p[1]        
        hopt, hcov = optimize.curve_fit(fh, a, h, p, bounds=([0.015, 1.],[0.5, 10.]))
        print (hopt)
        p=[0.2, 5.]
        SI = y[idx]
        fy = lambda  age, *p: SI*((1.0-np.exp(-1.0*p[0]*age))/(1.0-np.exp(-1.0*p[0]*Iage)))**p[1]        
        yopt, ycov = optimize.curve_fit(fy, a, y, p, bounds=([0.015, 1.],[0.5, 10.]))
        print (yopt)
        
        f = lambda age, SI, Iage, B1, B2: SI*((1.0-np.exp(-1.0*B1*age))/(1.0-np.exp(-1.0*B1*Iage)))**B2 
        motti[m][s][minpe][spe]['hpara']= {}
        motti[m][s][minpe][spe]['hpara']['Eq']='f = lambda age, SI, Iage, B1, B2: SI*((1.0-np.exp(-1.0*B1*age))/(1.0-np.exp(-1.0*B1*Iage)))**B2'
        motti[m][s][minpe][spe]['hpara']['desrc']= 'Schumacher model for height development, parameters from p0...p3: age [yr], SI h at index age [m], B1, B2 shape params'
        motti[m][s][minpe][spe]['hpara']['para']= [h[idx], a[idx], hopt[0], hopt[1]]

        motti[m][s][minpe][spe]['ypara']= {}
        motti[m][s][minpe][spe]['ypara']['Eq']='f = lambda age, SI, Iage, B1, B2: SI*((1.0-np.exp(-1.0*B1*age))/(1.0-np.exp(-1.0*B1*Iage)))**B2'
        motti[m][s][minpe][spe]['ypara']['desrc']= 'Schumacher model for yield, parameters from p0...p3: age [yr], SI y at index age [m], B1, B2 shape params'
        motti[m][s][minpe][spe]['ypara']['para']= [y[idx], a[idx], yopt[0], yopt[1]]
        
        #******** print figures for quality control**************
        printOpt = False
        if printOpt == True:
            fig= plt.figure(num = 'Stands', facecolor=(232/255.0, 243/255.0, 245.0/255), edgecolor='k',figsize=(18.0,12.0))   #Figsize(w,h), tuple inches 
            st = m + ' ' + str(minpe) +' ' + spe       
            plt.suptitle(st, fontsize = 18)
    
            filnam = folder + '\\Figs\\' + st + '.png'
        
            aext = np.arange(0,150,5)        
            fit =  f(aext, h[idx], Iage, hopt[0], hopt[1])  
            fity = f(aext, y[idx], Iage, yopt[0], yopt[1])  
            
            dh =derivative(f, aext, args=(h[idx], Iage, hopt[0], hopt[1]))
            dy =derivative(f, aext, args=(y[idx], Iage, yopt[0], yopt[1]))
            
            textstr = 'SI = ' +str(h[idx])[0:4] + '\n' + 'Iage = ' + str(Iage)[0:2] + '\n' + 'b1 = ' + str(hopt[0])[0:6] + '\n'+ 'b2 = ' + str(hopt[1])[0:6]
            plt.subplot(2,2,1); plt.plot(a, h, 'ro'); plt.plot(aext, fit, 'b-');plt.text(max(aext)*0.7,max(fit[np.isfinite(fit)])*0.8, textstr, fontsize =14)
            textstr = 'SI = ' +str(y[idx])[0:4] + '\n' + 'Iage = ' + str(Iage)[0:2] + '\n' +'b1 = ' + str(yopt[0])[0:6] + '\n'+ 'b2 = ' + str(yopt[1])[0:6]   
            plt.subplot(2,2,2); plt.plot(a, y, 'ro'); plt.plot(aext, fity, 'b-'); plt.text(max(aext)*0.7,max(fity[np.isfinite(fity)])*0.8, textstr, fontsize =14)
            plt.subplot(2,2,3); plt.plot(aext, dh, 'ro'); plt.plot(aext, dh, 'b-')
            plt.subplot(2,2,4); plt.plot(aext, dy, 'ro'); plt.plot(aext, dy, 'b-'); 
            plt.savefig(filnam)
            plt.close('all')
    #******* save to pickle ***************
    with open(fPickle, 'wb') as fp:
        pickle.dump(motti, fp)
    fp.close()

"""
def get_measured():
    #Retrieves measured datasets from Web database. 
    
    fi = 'C:\Apps\WinPython-64bit-2.7.10.3\SAMULI_ARI\work\Ini\spathy_sve.ini'
    ff= 'C:\Apps\WinPython-64bit-2.7.10.3\SAMULI_ARI\work\Spathy_runs\Spathy_ch_3.nc'    
    res=Dataset(ff, mode='r')                                                   # open in reading mode
    times=res['time']                                                           # get time variable from the file, manipulate to get dates
    dates=num2date(times[:],units=times.units,calendar=times.calendar)
    
    pgen,_,pbu,_=read_setup(fi) 
    start_date = datetime.datetime.strptime(pgen['start_date'], '%Y-%m-%d')
    end_date = datetime.datetime.strptime(pgen['end_date'], '%Y-%m-%d')


    dat=stand.vdataQuery(3,start_date, end_date,'dload')
"""
def make_video(outvid, images, outimg=None, fps=2, size=None,
               is_color=True, format="XVID"):
    """ 
    Create a video from a list of images.
 
    @param      outvid      output video
    @param      images      list of images to use in the video
    @param      fps         frame per second
    @param      size        size of each frame
    @param      is_color    color
    @param      format      see http://www.fourcc.org/codecs.php
    @return                 see http://opencv-python-tutroals.readthedocs.org/en/latest/py_tutorials/py_gui/py_video_display/py_video_display.html
 
    The function relies on http://opencv-python-tutroals.readthedocs.org/en/latest/.
    By default, the video will have the size of the first image.
    It will resize every image to this size before adding them to the video.
    """
    from cv2 import VideoWriter, VideoWriter_fourcc, imread, resize
    fourcc = VideoWriter_fourcc(*format)
    vid = None
    for image in images:
        #if not os.path.exists(image):
        #    raise FileNotFoundError(image)
        img = imread(image)
        print (image)
        if vid is None:
            if size is None:
                size = img.shape[1], img.shape[0]
            vid = VideoWriter(outvid, fourcc, float(fps), size, is_color)
        if size[0] != img.shape[1] and size[1] != img.shape[0]:
            img = resize(img, size)
        vid.write(img)
    vid.release()
    return vid

def construct_imgs():
    
    #ff= 'C:\Apps\WinPython-64bit-2.7.10.3\SAMULI_ARI\work\Spathy_runs\NutSpathy_ch_3.nc'
    folder = ofolder #'C:\Apps\WinPython-64bit-2.7.10.3\SAMULI_ARI\work\Spathy_runs/figs/'
    outvid = ofolder + 'fresh.avi' #'C:\\Apps\\WinPython-64bit-2.7.10.3\\SAMULI_ARI\\work\\Spathy_runs\\figs\\fresh.avi'

    res=Dataset(fn, mode='r') #open in reading mode
    times=res['time'] #get time variable from the file, manipulate to get dates
    dates=num2date(times[:],units=times.units,calendar=times.calendar)
    nd =  len(res['time'])  #number of days in the simulation
    stepIm = 10 #step between images in the vid, days
    #print np.shape(res['nut']['nsto'])    
    #print res['nut']['nsto'][0]    
    images = []
    
    for n in range(0,nd,stepIm):
        fname = folder + 'time' + str(n) + '.png'
        fig= plt.figure(num = 'Nitroen map', facecolor=(232/255.0, 243/255.0, 245.0/255), edgecolor='k',figsize=(24.0,12.0))   #Figsize(w,h), tuple inches 
        print (dates[n], n)
        plt.suptitle(res.description + ', date: ' +str(dates[n])[:10], fontsize = 18)
        fig.add_subplot(2,3,1); plt.imshow(res['nut']['nrel'][n,:,:], vmin = 0.0, vmax=0.3); plt.colorbar(); plt.title(str(res['nut']['nrel'].units))       
        fig.add_subplot(2,3,2); plt.imshow(res['nut']['nup'][n,:,:], vmin = 0.0, vmax=0.3); plt.colorbar(); plt.title(str(res['nut']['nup'].units))       
        fig.add_subplot(2,3,3); plt.imshow(res['nut']['nsto'][n,:,:], vmin = 0.0, vmax=30.0); plt.colorbar(); plt.title(str(res['nut']['nsto'].units))       
        fig.add_subplot(2,3,5); plt.imshow(res['nut']['nconc'][n,:,:], vmin = 0.0, vmax=1.2); plt.colorbar(); plt.title(str(res['nut']['nconc'].units))       
        _,c,r=  np.shape(res['nut']['nsto'])
        c, r = c/2,r/2
        fig.add_subplot(2,3,4); plt.plot(dates, res['nut']['nsto'][:,c,r]); plt.title(str(res['nut']['nsto'].units))       
        #ax = plt.twinx(); ax.plot(dates, res['nut']['nconc'][:,c,r])        
        fig.add_subplot(2,3,6); plt.plot(dates, res['nut']['nup'][:,c,r], 'r-'); plt.ylabel('Uptake, release [kg ha-1 day-1]')
        fig.add_subplot(2,3,6); plt.plot(dates, res['nut']['nrel'][:,c,r], 'g-')
        
        plt.savefig(fname)
        plt.close('all')
        images.append(fname)

    optVid = True
    if optVid == True:
        vid = make_video(outvid, images, fps=3)
        for im in images: os.remove(im)

    optPrT=False
    if optPrT==True:    
        fig= plt.figure(num = 'Times', facecolor=(232/255.0, 243/255.0, 245.0/255), edgecolor='k',figsize=(24.0,12.0))   #Figsize(w,h), tuple inches 
        fi =range(1,17); co=np.linspace(0,60,16, dtype=int); ro=np.linspace(0,80,16, dtype=int)
        for f, c, r in zip(fi,co,ro):    
            fig.add_subplot(4,4,f); plt.plot(res['nut']['nsto'][:,c,r])       
    
    res.close()



def extract_measured():    
    import zipfile    
    fol ='C:\Apps\WinPython-64bit-2.7.10.3\SAMULI_ARI\measured conc\\'
    zips = os.listdir(fol)
    #zf=zipfile.ZipFile(fol + zips[0], 'r')
    for z in zips:
        zf=zipfile.ZipFile(fol + z, 'r')
        name = zf.namelist()[0]
        n = z[:-4]
        outfol = fol +'\\' +n
        if not os.path.exists(outfol):
            os.makedirs(outfol)
            zf.extract(name, outfol)    
#extract_measured()
def initialize_netCDF_nut_bck(ID, gis, forc, fpath='C:\Apps\WinPython-64bit-2.7.10.3\SAMULI_ARI\work\Spathy_runs\\', fname=None):
 
    "netCDF file for NutSpathy outputs"    
    
    #from netCDF4 import Dataset, date2num 
    #from datetime import datetime
    
    #dimensions
    dlat,dlon=np.shape(gis['cmask'])
    dtime=None 
    
    if fname: ff=fpath + fname 
    else: ff=fpath +'NutSpathy_ch_'+str(ID)+'.nc'
    ncf= Dataset(ff,'w')
    ncf.description = 'SpatHy results. Catchment : ' +str(ID)
    print (ncf.description)
    ncf.history = 'created ' + datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    ncf.source = 'SpatHy -model v.0.99'
    
    ncf.createDimension('dtime',dtime)
    ncf.createDimension('dlon',dlon)
    ncf.createDimension('dlat',dlat)
    
    #create variables into base and groups 'forc','eval','cpy','bu','top'
    
    #call as createVariable(varname,type,(dimensions))
    time=ncf.createVariable('time','f8',('dtime',)); 
    time.units="days since 0001-01-01 00:00:00.0"; time.calendar='standard'
    
    lat=ncf.createVariable('lat','f4',('dlat',)); lat.units='ETRS-TM35FIN';
    lon=ncf.createVariable('lon','f4',('dlon',)); lon.units='ETRS-TM35FIN'; 

    #tvec=[k.to_datetime() for k in forc.index]    
    tvec = forc.index.to_pydatetime()
    time[:]=date2num(tvec, units=time.units, calendar=time.calendar); 
    lon[:]=gis['lon0']; lat[:]=gis['lat0']
    
    Nsto=ncf.createVariable('/nut/nsto','f4',('dtime','dlat','dlon',)); Nsto.units='soil N storage [kgha-1]'
    Nup=ncf.createVariable('/nut/nup','f4',('dtime','dlat','dlon',)); Nup.units='N uptake [kgha-1day-1]'
    Nrel=ncf.createVariable('/nut/nrel','f4',('dtime','dlat','dlon',)); Nrel.units='N release [kgha-1day-1]'
    Nconc=ncf.createVariable('/nut/nconc','f4',('dtime','dlat','dlon',)); Nconc.units='soil N conc [mgl-1]'    

    Psto=ncf.createVariable('/nut/psto','f4',('dtime','dlat','dlon',)); Psto.units='soil P storage [kgha-1]'
    Pup=ncf.createVariable('/nut/pup','f4',('dtime','dlat','dlon',)); Pup.units='P uptake [kgha-1day-1]'
    Prel=ncf.createVariable('/nut/prel','f4',('dtime','dlat','dlon',)); Prel.units='P release [kgha-1day-1]'
    Pconc=ncf.createVariable('/nut/pconc','f4',('dtime','dlat','dlon',)); Pconc.units='soil P conc [mgl-1]'    

    return ncf,ff

def initialize_netCDF_nut(ID, gis, tindex, fpath=None, fname=None):
 
    "netCDF file for NutSpathy outputs"    
    
    from netCDF4 import Dataset 
    from datetime import datetime

    #dimensions
    dlat,dlon=np.shape(gis['cmask'])
    dtime=None 

    #new from here
    print('**** creating NutSpaFHy netCDF4 file: ' + fname + ' ****')
    
    # create dataset & dimensions
    ff=fpath + fname
    ncf = Dataset(ff, 'w')
    ncf.description = 'nutSpatHy results. Catchment : ' + str(ID)
    ncf.history = 'created ' + datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    ncf.source = 'NutSpaFHy v.1.0'

    ncf.createDimension('dtime', dtime)
    ncf.createDimension('dlon', dlon)
    ncf.createDimension('dlat', dlat)


    # create variables into base and groups 'forc','eval','cpy','bu','top'
    # call as createVariable(varname,type,(dimensions))
    time = ncf.createVariable('time', 'f8', ('dtime',))
    time.units = "days since 0001-01-01 00:00:00.0"
    time.calendar = 'standard'
    tvec=tindex.to_pydatetime()
    time[:] = date2num(tvec, units=time.units, calendar=time.calendar)


    lat = ncf.createVariable('lat', 'f4', ('dlat',))
    lat.units = 'ETRS-TM35FIN'
    lon = ncf.createVariable('lon', 'f4', ('dlon',))
    lon.units = 'ETRS-TM35FIN'

    lon[:]=gis['lon0'] 
    lat[:]=gis['lat0']

    Nsto=ncf.createVariable('/nut/nsto','f4',('dtime','dlat','dlon',)); Nsto.units='soil N storage [kgha-1]'
    Nup=ncf.createVariable('/nut/nup','f4',('dtime','dlat','dlon',)); Nup.units='N uptake [kgha-1dt-1]'
    Nup_gv=ncf.createVariable('/nut/nup_gv','f4',('dtime','dlat','dlon',)); Nup_gv.units='ground vegetation N uptake [kgha-1dt-1]'
    Nrel=ncf.createVariable('/nut/nrel','f4',('dtime','dlat','dlon',)); Nrel.units='N release [kgha-1dt-1]'
    Nconc=ncf.createVariable('/nut/nconc','f4',('dtime','dlat','dlon',)); Nconc.units='soil N conc [mgl-1]'    
    Ntogw=ncf.createVariable('/nut/ntogw','f4',('dtime','dlat','dlon',)); Ntogw.units='N down flux to ground water [kg ha-1 dt-1]'
    Nretflow=ncf.createVariable('/nut/nretflow','f4',('dtime','dlat','dlon',)); Nretflow.units='N return flow  [kg ha-1 dt-1]'
    NtoSrun=ncf.createVariable('/nut/ntosrun','f4',('dtime','dlat','dlon',)); NtoSrun.units='N to surface runoff [kg ha-1 dt-1]'
    Psto=ncf.createVariable('/nut/psto','f4',('dtime','dlat','dlon',)); Psto.units='soil P storage [kgha-1]'
    Pup=ncf.createVariable('/nut/pup','f4',('dtime','dlat','dlon',)); Pup.units='P uptake [kgha-1dt-1]'
    Pup_gv=ncf.createVariable('/nut/pup_gv','f4',('dtime','dlat','dlon',)); Pup_gv.units='ground vegetation P uptake [kgha-1dt-1]'
    Prel=ncf.createVariable('/nut/prel','f4',('dtime','dlat','dlon',)); Prel.units='P release [kgha-1dt-1]'
    Pconc=ncf.createVariable('/nut/pconc','f4',('dtime','dlat','dlon',)); Pconc.units='soil P conc [mgl-1]'    
    Ptogw=ncf.createVariable('/nut/ptogw','f4',('dtime','dlat','dlon',)); Ptogw.units='P down flux to ground water [kg ha-1 dt-1]'
    Pretflow=ncf.createVariable('/nut/pretflow','f4',('dtime','dlat','dlon',)); Pretflow.units='P return flow to surface runoff [kg ha-1 dt-1]'
    PtoSrun=ncf.createVariable('/nut/ptosrun','f4',('dtime','dlat','dlon',)); PtoSrun.units='P to surface runoff [kg ha-1 dt-1]'
    
    print ('    + netCDF file created')
    return ncf,ff

def soil_type_pF(spara):
    #soil mask: 1 Coarse textured, 2 Medium textured, 3 Fine textured, 4 Peat
    min_vg ={
    'Humus': [0.96,0.04,0.7827,1.4147],
    'LoamySand':[0.39,0.074,0.035,2.39],
    'SiltyLoam':[0.43,0.061,0.012,1.39],
    'SiltyClay':[0.47,0.163,0.023,1.39]
    }


    lenvp=len(spara['vonP top'])    
    vonP = np.ones(spara['nLyrs'])*spara['vonP bottom']; vonP[0:lenvp] = spara['vonP top']  # degree of  decomposition, von Post scale
    peat_pF, _ = peat_hydrol_properties(vonP, var='H', ptype=spara['peat type'])  # peat hydraulic properties after Päivänen 1973
    #pF shape (nLyrs,4)    
    coarse_pF=np.zeros((spara['nLyrs'],4))
    coarse_pF[0][:] = min_vg['Humus']
    coarse_pF[1:][:] = min_vg['LoamySand']*np.ones((spara['nLyrs']-1,4))

    medium_pF=np.zeros((spara['nLyrs'],4))
    medium_pF[0][:] = min_vg['Humus']
    medium_pF[1:][:] = min_vg['SiltyLoam']*np.ones((spara['nLyrs']-1,4))
    
    fine_pF=np.zeros((spara['nLyrs'],4))
    fine_pF[0][:] = min_vg['Humus']
    fine_pF[1:][:] = min_vg['SiltyClay']*np.ones((spara['nLyrs']-1,4))
    
    return coarse_pF,medium_pF,fine_pF,peat_pF

def local_s_to_gwl(soilmask,local_s,s_to_gwl):

    shape=np.shape(soilmask)
    soilmask=np.ravel(soilmask)
    local_s = np.ravel(local_s)    
    gwl=np.ones(np.shape(soilmask))*-999    
    ixCoarse=np.ravel(np.where(soilmask==1))
    ixMedium=np.ravel(np.where(soilmask==2))
    ixFine=np.ravel(np.where(soilmask==3))
    ixPeat=np.ravel(np.where(soilmask==4))
        
    if np.shape(ixCoarse)[0] >0:    
        gwl[ixCoarse]= s_to_gwl['1']['stoToGwl'](list(s_to_gwl['1']['fullSto']-local_s[ixCoarse]))
    if np.shape(ixMedium)[0] >0:    
        gwl[ixMedium]= s_to_gwl['2']['stoToGwl'](list(s_to_gwl['2']['fullSto']-local_s[ixMedium]))
    if np.shape(ixFine)[0] >0:    
        gwl[ixFine]= s_to_gwl['3']['stoToGwl'](list(s_to_gwl['3']['fullSto']-local_s[ixFine]))
    if np.shape(ixPeat)[0] >0:    
        gwl[ixPeat]= s_to_gwl['4']['stoToGwl'](list(s_to_gwl['4']['fullSto']-local_s[ixPeat]))
    return np.reshape(gwl,shape)

def local_s_to_gwl_functions():
    """
    composes interpolation functions for gwl
    """
    spara = {
    'nLyrs':50, 'dzLyr': 0.05, 'peat type':['S'], 
    'vonP bottom': 8, 'vonP top': [3,4,5,6,7]
    }
    dz=np.ones(spara['nLyrs'])*spara['dzLyr']
    z = np.cumsum(dz)-0.5*spara['dzLyr']
    Ksat = np.ones(spara['nLyrs'])
    coarse_pF,medium_pF,fine_pF,peat_pF=soil_type_pF(spara)    
    coarse_hToSto, coarse_stoToGwl, _, _ = CWTr(spara['nLyrs'], z, dz, coarse_pF, Ksat, direction='negative') # interpolated storage, transmissivity and diff water capacity functions
    medium_hToSto, medium_stoToGwl, _, _ = CWTr(spara['nLyrs'], z, dz, medium_pF, Ksat, direction='negative') # interpolated storage, transmissivity and diff water capacity functions
    fine_hToSto, fine_stoToGwl, _, _ = CWTr(spara['nLyrs'], z, dz, fine_pF, Ksat, direction='negative') # interpolated storage, transmissivity and diff water capacity functions
    peat_hToSto, peat_stoToGwl, _, _ = CWTr(spara['nLyrs'], z, dz, peat_pF, Ksat, direction='negative') # interpolated storage, transmissivity and diff water capacity functions

    s_to_gwl={
    '1':{'name': 'Coarse mineral', 'fullSto': coarse_hToSto(0.0), 'stoToGwl':coarse_stoToGwl},
    '2':{'name': 'Medium mineral', 'fullSto': medium_hToSto(0.0), 'stoToGwl':medium_stoToGwl},
    '3':{'name': 'Fine mineral', 'fullSto': fine_hToSto(0.0), 'stoToGwl':fine_stoToGwl},
    '4':{'name': 'peat', 'fullSto': peat_hToSto(0.0), 'stoToGwl':peat_stoToGwl},   
    }
    return s_to_gwl    

def peat_hydrol_properties(x, unit='g/cm3', var='bd', ptype='A'):
    """
    Peat water retention and saturated hydraulic conductivity as a function of bulk density
    Päivänen 1973. Hydraulic conductivity and water retention in peat soils. Acta forestalia fennica 129.
    see bulk density: page 48, fig 19; degree of humification: page 51 fig 21
    Hydraulic conductivity (cm/s) as a function of bulk density(g/cm3), page 18, as a function of degree of humification see page 51 
    input:
        - x peat inputvariable in: db, bulk density or dgree of humification (von Post)  as array \n
        - bulk density unit 'g/cm3' or 'kg/m3' \n
        - var 'db' if input variable is as bulk density, 'H' if as degree of humification (von Post) \n
        - ptype peat type: 'A': all, 'S': sphagnum, 'C': Carex, 'L': wood, list with length of x 
    output: (ThetaS and ThetaR in m3 m-3)
        van Genuchten water retention parameters as array [ThetaS, ThetaR, alpha, n] \n
        hydraulic conductivity (m/s)
    """
    #paras is dict variable, parameter estimates are stored in tuples, the model is water content = a0 + a1x + a2x2, where x is
    para={}                                                                     #'bd':bulk density in g/ cm3; 'H': von Post degree of humification
    para['bd'] ={'pF0':(97.95, -79.72, 0.0), 'pF1.5':(20.83, 759.69, -2484.3),
            'pF2': (3.81, 705.13, -2036.2), 'pF3':(9.37, 241.69, -364.6),
            'pF4':(-0.06, 249.8, -519.9), 'pF4.2':(0.0, 174.48, -348.9)}
    para['H'] ={'pF0':(95.17, -1.26, 0.0), 'pF1.5':(46.20, 8.32, -0.54),
            'pF2': (27.03, 8.14, -0.43), 'pF3':(17.59, 3.22, -0.07),
            'pF4':(8.81, 3.03, -0.10), 'pF4.2':(5.8, 2.27, -0.08)}
    
    intp_pF1={}                                                                 # interpolation functions for pF1        
    intp_pF1['bd'] = interp1d([0.04,0.08,0.1,0.2],[63.,84.,86.,80.],fill_value='extrapolate')
    intp_pF1['H'] = interp1d([1.,4.,6.,10.],[75.,84.,86.,80.],fill_value='extrapolate')
    
    #Saturatated hydraulic conductivity parameters
    Kpara ={'bd':{'A':(-2.271, -9.80), 'S':(-2.321, -13.22), 'C':(-1.921, -10.702), 'L':(-1.921, -10.702)}, 
            'H':{'A':(-2.261, -0.205), 'S':(-2.471, -0.253), 'C':(-1.850, -0.278), 'L':(-2.399, -0.124)}}
    
    vg_ini=(0.88,	0.09, 0.03, 1.3)                                              # initial van Genuchten parameters (porosity, residual water content, alfa, n)

    x = np.array(x)
    prs = para[var]; pF1=intp_pF1[var]
    if unit=='kg/m3'and var=='db': x=x/1000.
    if  np.shape(x)[0] >1 and len(ptype)==1:
        ptype=np.repeat(ptype, np.shape(x)[0])        
    vgen = np.zeros((np.size(x),4))
    Ksat = np.zeros((np.size(x)))
    wcont = lambda x, a: a[0] + a[1]*x + a[2]*x**2.
    van_g = lambda pot, *p:   p[1] + (p[0] - p[1]) / (1. + (p[2] * pot) **p[3]) **(1. - 1. / p[3])   
    K = lambda x, a: 10.**(a[0] + a[1]*x) / 100.   # to m/s   
    
    potentials =np.array([0.01, 10.,32., 100.,1000.,10000.,15000. ])
    wc = (np.array([wcont(x,prs['pF0']), pF1(x), wcont(x,prs['pF1.5']), wcont(x,prs['pF2']),
               wcont(x,prs['pF3']), wcont(x,prs['pF4']),wcont(x,prs['pF4.2'])]))/100.
        
    for i,s in enumerate(np.transpose(wc)):
        vgen[i],_= curve_fit(van_g,potentials,s, p0=vg_ini)                      # van Genuchten parameters
        
    for i, a, pt in zip(range(len(x)), x, ptype):
        Ksat[i] = K(a, Kpara[var][pt])                                          # hydraulic conductivity (cm/s -> m/s) 
    
    return vgen, Ksat

#print peat_hydrol_properties([5], unit='g/cm3', var='H', ptype='C')

def wrc(pF, x=None, var=None):
    """
    vanGenuchten-Mualem soil water retention curve\n
    IN:
        pF - dict['ThetaS': ,'ThetaR': ,'alpha':, 'n':,] OR
           - list [ThetaS, ThetaR, alpha, n]
        x  - soil water tension [m H2O = 0.1 kPa]
           - volumetric water content [vol/vol]
        var-'Th' is x=vol. wat. cont.
    OUT:
        res - Theta(Psii) or Psii(Theta)
    NOTE:\n
        sole input 'pF' draws water retention curve and returns 'None'. For drawing give only one pF-parameter set. 
        if several pF-curves are given, x can be scalar or len(x)=len(pF). In former case var is pF(x), in latter var[i]=pf[i,x[i]]
               
    Samuli Launiainen, Luke 2/2016
    """
    if type(pF) is dict: #dict input
        #Ts, Tr, alfa, n =pF['ThetaS'], pF['ThetaR'], pF['alpha'], pF['n']
        Ts=np.array(pF['ThetaS'].values()); Tr=np.array( pF['ThetaR'].values()); alfa=np.array( pF['alpha'].values()); n=np.array( pF['n'].values())
        m= 1.0 -np.divide(1.0,n)
    elif type(pF) is list: #list input
        pF=np.array(pF, ndmin=1) #ndmin=1 needed for indexing to work for 0-dim arrays
        Ts=pF[0]; Tr=pF[1]; alfa=pF[2]; n=pF[3] 
        m=1.0 - np.divide(1.0,n)
    elif type(pF) is np.ndarray:
        Ts, Tr, alfa, n = pF.T[0], pF.T[1], pF.T[2], pF.T[3]
        m=1.0 - np.divide(1.0,n)
    else:
        print ('Unknown type in pF')
        
    def theta_psi(x): #'Theta-->Psi'
        x=np.minimum(x,Ts) 
        x=np.maximum(x,Tr) #checks limits
        s= ((Ts - Tr) / (x - Tr))#**(1/m)
        Psi=-1e-2/ alfa*(s**(1/m)-1)**(1/n) # in m
        return Psi
        
    def psi_theta(x): # 'Psi-->Theta'
        x=100*np.minimum(x,0) #cm
        Th = Tr + (Ts-Tr)/(1+abs(alfa*x)**n)**m
        return Th           
 
    if var is 'Th': y=theta_psi(x) #'Theta-->Psi'           
    else: y=psi_theta(x) # 'Psi-->Theta'          
    return y

def CWTr(nLyrs, z, dz, pF, Ksat, direction='positive'):
    """
    Returns interpolation functions 
        sto=f(gwl)  profile water storage as a function ofground water level
        gwl=f(sto)  ground water level
        tra=f(gwl)  transissivity
    Input:
        nLyrs number of soil layers
        d depth of layer midpoint
        dz layer thickness
        pF van Genuchten water retention parameters: ThetaS, ThetaR, alfa, n
        Ksat saturated hydraulic conductivity in m s-1
        direction: positive or negative downwards
    """    
    #-------Parameters ---------------------
    z = np.array(z)   
    dz =np.array(dz)
    #--------- Connection between gwl and water storage------------
    d = 6 if direction == 'positive' else -6   
    gwl=np.linspace(0,d,150)
    if direction == 'positive':
        sto = [sum(wrc(pF, x = np.minimum(z-g, 0.0))*dz) for g in gwl]     #equilibrium head m
    else:
        sto = [sum(wrc(pF, x = np.minimum(z+g, 0.0))*dz) for g in gwl]     #equilibrium head m
    gwlToSto = interp1d(np.array(gwl), np.array(sto), fill_value='extrapolate')
    sto = list(sto); gwl= list(gwl)        
    sto.reverse(); gwl.reverse()
    stoToGwl =interp1d(np.array(sto), np.array(gwl), fill_value='extrapolate')
    C = interp1d(np.array(gwl), np.array(np.gradient(gwlToSto(gwl))/np.gradient(gwl)), fill_value='extrapolate')  #storage coefficient function      
    
    del gwl, sto
        
    #----------Transmissivity-------------------
    K=np.array(Ksat*86400.)   #from m/s to m/day
    tr =[sum(K[t:]*dz[t:]) for t in range(nLyrs)]        
    if direction=='positive':        
        gwlToTra = interS(z, np.array(tr))            
    else:
        z= list(z);  z.reverse(); tr.reverse()
        gwlToTra = interS(-np.array(z), np.array(tr))                    
    del tr
    return gwlToSto, stoToGwl, gwlToTra, C

def daily_to_monthly(sp_res, forc, pnut, balance=True):
    """
    Parameters
    ----------
    sp_res : netCDF4 file
        DESCRIPTION. SpaFHy results in gridded daily timeseries
    forc : pandas dataframe
        DESCRIPTION. daily meteorological input
    pnut : nutrient parameters
        DESCRIPTION. dictionary
    balance : TYPE, boolean
        DESCRIPTION. The default is True.

    Returns
    -------
    Wm: TYPE np-array
        DESCRIPTION. Monthly mean water content, m3 m-3
    lsm: TYPE np-array
        DESCRIPTION. Montly mean Local saturation deficit, m
    Rm: TYPE np-array
        DESCRIPTION. Monthly cumulative runoff, m/month 
    Dm: TYPE np-array
        DESCRIPTION. Monthly cumulative Drainage from root layer to below, m/month
    Dretm: TYPE
        DESCRIPTION. Monthly cumulative Drainage return flow, m/month
    Infiltr: TYPE np-array
        DESCRIPTION. Monthly cumulative infiltration flow, m/month
    s_run TYPE: np-array
        DESCRIPTION. Monthly cumulative surface runoff, m/month
    Tm : TYPE np-array
        DESCRIPTION. Monthly mean air temperature, deg C
    ddsm : TYPE np-array
        DESCRIPTION. Temperature sum, degree-days, threshold 5 deg C
    P : TYPE np-array
        DESCRIPTION. Monthly cumulative precipitation, (m check the unit!)

    """
    
    # this for all
    Runoff = sp_res['top']['Qt'][:]
    dlon = sp_res.variables['lon'][:]
    dlat = sp_res.variables['lat'][:]
    times   =sp_res.variables['time']
    units   =times.units
    calendar=times.calendar
    #firstdate=num2date(times[0],units,calendar=calendar)
    #lastdate =num2date(times[-1],units,calendar=calendar)
    firstdate=num2date(times[0],units,calendar=calendar)
    lastdate =num2date(times[-1],units,calendar=calendar)

    sp_res.close()
    
    #print ('++++++++++++++++++++++++++++')
    #print (sp_res['time'][:].squeeze())
    #print (units, calendar)
    #print (firstdate, lastdate)    
    #print ('++++++++++++++++++++++++++++')
    #import sys; sys.exit()
    # for all:
    #--------- Process xarrays----------------------
    ds0 = xr.open_dataset(pnut['spafhyfile'])
    ds0['dtime'] = ds0['time']
    dtime = ds0['dtime']
    ds0.close()
    
    #----------Monthly mean values-----------------------
    dsbu = xr.open_dataset(pnut['spafhyfile'],group= 'bu', decode_times=True)
    if balance:
        Wliq = dsbu['Wliq_top']
        W = xr.DataArray(Wliq, coords={'dtime':dtime, 'dlon':dlon, 'dlat':dlat})
        Wm = W.resample(dtime='1M').mean()                                          #Monthly mean water content m3 m-3
    
    dstop =xr.open_dataset(pnut['spafhyfile'],group= 'top', decode_times=True)
    local_s= dstop['Sloc']
    ls = xr.DataArray(local_s, coords={'dtime':dtime, 'dlon':dlon, 'dlat':dlat})
    lsm = ls.resample(dtime='1M').mean()                                        #Montly mean Local saturation deficit, m 
    
    #------------Monthly sum values---------------------
    R = xr.DataArray(Runoff, [dtime]) 
    Rm = R.resample(dtime='1M').sum()                                           #Monthly cumulative runoff, m
    
    if balance:
        Drain = dsbu['Drain']
        D = xr.DataArray(Drain, coords={'dtime':dtime, 'dlon':dlon, 'dlat':dlat})
        Dm = D.resample(dtime='1M').sum()                                           #Monthly cumulative Drainage from root layer to below, m/month
        Drainretflow = dsbu['Drainretflow']
        Dret = xr.DataArray(Drainretflow, coords={'dtime':dtime, 'dlon':dlon, 'dlat':dlat})
        Dretm = Dret.resample(dtime='1M').sum()                                     #Monthly cumulative Drainage return flow
        Infi = dsbu['Infil']
        Infil = xr.DataArray(Infi, coords={'dtime':dtime, 'dlon':dlon, 'dlat':dlat})
        Infiltr = Infil.resample(dtime='1M').sum()                                     #Monthly cumulative infiltration flow
        s_runoff = dsbu['Infi_ex']
        s_runo = xr.DataArray(s_runoff, coords={'dtime':dtime, 'dlon':dlon, 'dlat':dlat})
        s_run = s_runo.resample(dtime='1M').sum()                                     #Monthly cumulative surface runoff flow
        
    dsbu.close()
    dstop.close()
    
    #--------------Weather variables------------------
    Tm = forc['T'].resample('M', convention ='end').mean()                    # Monthly mean air temperature 
    base = 5.0    
    dd = forc['T']-base    
    dd[dd<base] = 0.0
    ddsm = dd.resample('M', convention='end').sum()
    ddsm = ddsm.cumsum(axis=0)                                              #cumulative temperature sum degree days
    P = forc['Prec'].resample('M', convention ='end').sum()
    

    if balance is False:
        Wm=0; Dm=0; Dretm=0;Tm=0; P=0; Infiltr=0; s_run =0
        
    return np.array(Wm), np.array(lsm), np.array(Rm), np.array(Dm), np.array(Dretm), np.array(Infiltr), np.array(s_run), Tm, ddsm, P        

def estimate_imm(gisdata):
    cmask = np.ravel(gisdata['cmask'])
    peatm = np.ravel(gisdata['peatm'])
    smc = np.ravel(gisdata['smc'])
    sfc = np.ravel(gisdata['sfc'])
    tot = (cmask==1).sum()
    
    peat_sfc1 = np.logical_and(sfc==1, peatm==1).sum()/tot #peatTot
    peat_sfc2 = np.logical_and(sfc==2, peatm==1).sum()/tot #peatTot
    rich_p = peat_sfc1 + peat_sfc2
    
    min_sfc5  = np.logical_and(sfc==5, peatm==0).sum()/tot #peatTot
    min_sfc6  = np.logical_and(sfc==6, peatm==0).sum()/tot #peatTot
    poor_min = min_sfc5 + min_sfc6
    
    bog_p = (smc==3).sum() / tot 
    vol = np.nanmean(gisdata['vol'])
    conif_vol = np.nanmean(gisdata['p_vol']+gisdata['s_vol'])
    
    #these parameter values with corrected ground vegetation model 
    iN_peat = 0.652 + 0.282*(conif_vol/vol) - 0.15*bog_p        #R2 0.607
    iN_miner = 0.894 + 0.284*poor_min                           #R2 0.301
    iP_peat = 0.757                                             #average over the data
    iP_miner = 0.807                                            #average over the data
    
    #iN_peat = 0.559 - 0.241*bog_p + 0.400*(conif_vol/vol) + 0.024*vol/100.  #R2 0.91
    #iN_miner = 0.671 - 0.199*bog_p + 0.307*(conif_vol/vol)+ 0.017*vol/100.   #R2 0.88
    #iP_peat = 0.954 - 2.673*rich_p    #R2 0.39     
    #iP_miner = 0.839 + 0.071*vol/100.    #0.34
 
    #iN_peat = 0.566 - 0.242*bog_p + 0.395*(conif_vol/vol) + 0.023*vol/100.  #R2 0.92
    #iN_miner = 0.673 - 0.199*bog_p + 0.307*(conif_vol/vol)+ 0.017*vol/100.   #R2 0.88
    #iP_peat = 0.955 - 2.638*rich_p    #R2 0.39     
    #iP_miner = 0.840 + 0.071*vol/100.    #0.34
    
    #iN_peat = 0.549 - 0.256*bog_p + 0.411*(conif_vol/vol) + 0.025*vol/100.  #R2 0.92
    #iN_miner = 0.665 - 0.212*bog_p + 0.313*(conif_vol/vol)+ 0.018*vol/100.   #R2 0.89
    #iP_peat = 0.955 - 2.598*rich_p    #R2 0.40     
    #iP_miner = 0.830 + 0.078*vol/100.    #0.42

    return (iN_peat, iN_miner), (iP_peat, iP_miner)                                                                     

def description_spafhy(sp_res, pgen, cmask, pnut, startd='2023-01-01', endd='2023-12-31'):

    # this for all
    #ET= sp_res['cpy']['ET'][:]
    #Wliq= sp_res['bu']['Wliq'][:]
    #Sloc= sp_res['top']['Sloc'][:]
    
    #Runoff = sp_res['top']['Qt'][:]
    startd=str(startd)
    endd=str(endd)
    
    dlon = sp_res.variables['lon'][:]
    dlat = sp_res.variables['lat'][:]
    times   =sp_res.variables['time']
    units   =times.units
    calendar=times.calendar
    firstdate=num2date(times[0],units,calendar=calendar)
    lastdate =num2date(times[-1],units,calendar=calendar)
    sp_res.close()
    
    # for all:
    #--------- Process xarrays----------------------
    ds0 = xr.open_dataset(pnut['spafhyfile'])#tämähän on sama kuin sp_res
    ds0['dtime'] = ds0['time']
    dtime = ds0['dtime']
    ds0.close()
    
    #----------Monthly mean values-----------------------
    dsbu = xr.open_dataset(pnut['spafhyfile'],group= 'bu', decode_times=True)
    Wliq = dsbu['Wliq_top'] #vai wliq vaan eikä tota top?
    W = xr.DataArray(Wliq, coords={'dtime':dtime, 'dlon':dlon, 'dlat':dlat})
    W=W.loc[startd:endd]
    Wliqm = W.resample(dtime='1M').mean()                                          #Monthly mean water content m3 m-3
    #Mbe = ncf.createVariable('/bu/Mbe', 'f4', ('dtime', 'dlat', 'dlon',))
    #Mbe.units = 'mass-balance error [mm]'
    Mbe=dsbu['Mbe']
    Mbe_x = xr.DataArray(Mbe, coords={'dtime':dtime, 'dlon':dlon, 'dlat':dlat})
    Mbe_x=Mbe_x.loc[startd:endd]
    Mbem = Mbe_x.resample(dtime='1M').mean()


    dstop =xr.open_dataset(pnut['spafhyfile'],group= 'top', decode_times=True)
    local_s= dstop['Sloc']
    ls = xr.DataArray(local_s, coords={'dtime':dtime, 'dlon':dlon, 'dlat':dlat})
    ls=ls.loc[startd:endd]
    Slocm = ls.resample(dtime='1M').mean()                                        #Montly mean Local saturation deficit, m over the whole area
    dscpy =xr.open_dataset(pnut['spafhyfile'],group= 'cpy', decode_times=True)
    ETx = dscpy['ET']
    ET = xr.DataArray(ETx, coords={'dtime':dtime, 'dlon':dlon, 'dlat':dlat})
    ET=ET.loc[startd:endd]
    ETm = ET.resample(dtime='1M').sum()                                     #Monthly cumulative Drainage return flow
    Mbecpy=dscpy['Mbe']
    Mbe_xc = xr.DataArray(Mbe, coords={'dtime':dtime, 'dlon':dlon, 'dlat':dlat})
    Mbe_xc=Mbe_xc.loc[startd:endd]
    Mbecpym = Mbe_xc.resample(dtime='1M').mean()

    scens=[f for f in os.listdir(pgen['gis_scenarios']) if f.endswith('.asc')]
    if scens:
        for s in scens:
            cut, _, _, _, _ = read_AsciiGrid(pgen['gis_scenarios'] + s)
            #ix = np.where(np.isfinite(cmask)) 
            #cut2 = cut[ix].copy()
            ls_cc=ls.where(cut>0)
            Slocmcc = ls_cc.resample(dtime='1M').mean()                                        #Montly mean Local saturation deficit, m over the whole area
            W_cc=W.where(cut>0)
            Wliqmcc = W_cc.resample(dtime='1M').mean()                                        #Montly mean Local saturation deficit, m over the whole area
            
            if Slocmcc.any():
                print("ok")
            else:
                print("did not work")
    else:
        print('No available clear-cut scenario rasters ')
        Slocmcc=None
        Wliqmcc=None
    
    
    
    dsbu.close()
    dstop.close()
    dscpy.close()
    
   
        
       
    return np.array(Wliqm), np.array(Slocm), np.array(ETm), np.array(Slocmcc), np.array(Wliqmcc), np.array(Mbem),np.array(Mbecpym )#, np.array(Runoffm) #,  Tm, ddsm, P        

