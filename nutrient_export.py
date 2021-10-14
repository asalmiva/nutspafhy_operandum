# -*- coding: utf-8 -*-
"""
Created on Mon Feb 25 13:06:43 2019

@author: lauren
"""
import datetime 
import numpy as np
import pandas as pd
import datetime
from scipy.spatial import distance
import matplotlib.pylab as plt
import seaborn as sns

class Export():
    def __init__(self,pgen, local_sm, Runoff, ddsm, resn, gisdata, soildata): #, kn,kp,nout,pout):
        
        """
        Dimensions of netCFD (time, longitude, latitude): ie. time, y, x
        Dimensions in Ascii-gis files: (row, col): ie, x, y
        Why s
      """

        print ('- Initializing nutrient export calculation')
        self.start_date = datetime.datetime.strptime(pgen['start_date'], '%Y-%m-%d')
        self.end_date = datetime.datetime.strptime(pgen['end_date'], '%Y-%m-%d')
        print ('    + Computing from', self.start_date, ' to ', self.end_date)
        sim_yrs = (self.end_date-self.start_date).days / 365.25    
        self.dates = ddsm.index

        """ dimensions of the computation"""
        self.cmask = gisdata['cmask']
        self.gridshape = np.shape(self.cmask)
        self.Nsteps = np.shape(Runoff)[0]                                      # number of months in the smulation
        self.runoff=np.empty(self.Nsteps)
        self.cellsize = gisdata['cellsize']
        stream_ix = np.transpose(np.where(np.isfinite(gisdata['stream'])))     # stream index, coordinates in order y, x
        self.catchment_ix = np.transpose(np.where(np.isfinite(self.cmask)))    # catchment index, coordinates now in order y, x     
        
        
        """ parameters""" 
        stream_depth = 0.5                                                     # strem depth in m        
        self.soild = 0.5                                                       # soil depth below the root zone 
        self.poros = soildata['poros']                                         # porosity
        ksat=soildata['ksat']                                                  # Hydraulic conductivity
        dem = gisdata['dem']                                                   # Digital elevation models
        self.catchment_area = np.nansum(self.cmask)*self.cellsize**2           # in m2
        self.dsdist = gisdata['dsdist']                                        # distance to water body, m        
        ele_at_stream = gisdata['stream']*(dem-stream_depth)                   # elevation at stream m
 
        """initialization of variables"""
        self.runoff[:]=np.array(Runoff)      #res['top']['Qt'][:]
        self.local_sm = local_sm                                                # local monthly mean saturation deficit
        self.gw_storage = self.poros*(self.soild-self.local_sm[0,:,:])*np.square(float(gisdata['cellsize']))        
        self.nsto = resn['nut']['nsto'][0,:,:]                                 # kg/ha in root zone
        self.psto = resn['nut']['psto'][0,:,:]
        self.gwN = np.nansum(self.gw_storage*1e3*self.nsto*1e-6)               # kg N
        self.gwP = np.nansum(self.gw_storage*1e3*self.psto*1e-6)               # kg P
        print ('    + Initial ground water volume', np.round(np.nansum(self.gw_storage), 1), 'm3')
        print ('    + Initial N content', np.round(self.gwN,0), 'kg')            # kg N
        print ('    + Initial P content', np.round(self.gwP,1), 'kg')           # kg P

        
        dist = distance.cdist(stream_ix, self.catchment_ix, 'euclidean')       # distance to stream each node
        min_dist_ix = np.argmin(dist, axis=0)                                  # layer from where the minimum distance cell can be found                
        dist_2 = np.min(dist, axis=0)                                          # finding the minimum distance value, shrinking dimensions to catchmet's
        self.catchment_ix = np.transpose(self.catchment_ix)
        b = np.full( self.gridshape, np.NaN)
        b[np.array(self.catchment_ix[0]),np.array(self.catchment_ix[1])] = min_dist_ix    # tells the index where to fond the stream node in the 2d distance array (index 0)
        
          
        layer_ix = (b[np.where(np.isfinite(b))].astype('int32'))
        receiving_node = np.transpose(stream_ix[layer_ix])
        receiving_ele =dem[receiving_node[0],receiving_node[1]]
        c = np.full(self.gridshape, np.NaN)
        c[np.array(self.catchment_ix[0]),np.array(self.catchment_ix[1])] = receiving_ele - stream_depth  # elevation of the receiving stream node for all  grid points, m 
        #slope = np.maximum((dem-c)/self.dsdist, 0.001)                         # slope using delta h and downslope distance
        
        a = np.full(self.gridshape, np.NaN)
        dist_2 = np.ravel(dist_2)                                              # euclidean distance to nearest stream node, unit here pixels

        a[np.array(self.catchment_ix[0]),np.array(self.catchment_ix[1])] = dist_2*self.cellsize  # distance to nearest stream node, m
        a = np.where( a > self.cellsize, a, 0.0)
        a = np.where( a > 0.0, a, 0.0)
        #slope = np.maximum((dem-c)/a, 0.01)
        slope = np.where(a > 0.0,(dem-c)/a, np.nan)*self.cmask
        #slope = np.where(self.dsdist > 0.0,(dem-c)/self.dsdist, np.nan)*self.cmask

        ksat= self.cmask*1e-4
        vx = (ksat*86400.)/self.poros * slope                                  # m/day
        #self.residence_time = self.cmask*(a/vx).astype(int)
        self.residence_time = np.where(a > 0, (a/vx).astype(int),0)*self.cmask
        #self.residence_time = np.where(self.dsdist > 0, (self.dsdist/vx).astype(int),0)*self.cmask
        self.residence_time = np.where(self.residence_time > 0.0, self.residence_time, 0)*self.cmask
        self.residence_time = (self.residence_time / 30.).astype(int)*self.cmask  #to months
        
        
        #self.residence_time = self.residence_time
        
        outopt = False
        if outopt:
            #print 'catchment', np.where(np.isfinite(gisdata['cmask']))
            dem = gisdata['dem']
            fig = plt.subplot(331)
            plt.imshow(dem)
            plt.colorbar()
            plt.title('dem')        
            plt.imshow(gisdata['stream'])
            plt.subplot(332)
            plt.imshow(ele_at_stream)
            plt.colorbar()
            plt.title('Elevation in stream')
            plt.imshow(gisdata['peatm'], alpha=0.4)
            plt.subplot(333)
            plt.imshow(a)
            plt.colorbar()
            plt.title('Distance to nearest stream node')
            plt.subplot(334)
            plt.imshow(self.dsdist)
            plt.colorbar()
            plt.title('dsdist')
            plt.subplot(335)
            plt.imshow(c)
            plt.colorbar()
            plt.title('receiving elevation')
            plt.subplot(336)
            plt.imshow(slope)
            plt.colorbar()
            plt.title('slope')     
            plt.subplot(337)
            plt.imshow(ksat*86400.)
            plt.colorbar()
            plt.title('ksat')
            plt.subplot(338)
            plt.imshow(vx)
            plt.colorbar()
            plt.title('vx')
    
            plt.subplot(339)
            plt.imshow(self.residence_time)
            plt.colorbar()
            plt.title('residence time')
    
            fig = plt.figure(num='residence')
            plt.plot(111)
            plt.imshow(self.residence_time)
            plt.colorbar()
            plt.title('residence time')
    

       
        #******* Time delay function computed from distance to stream, and water flow velocity  
        
        print ('    + Ksat, m s-1', np.nanmean(ksat)) 
        print ('    + Water velocity, m day-1', np.round(np.nanmean(vx),3)) 
        print ('    + Mean distance to stream, m ', np.round(np.nanmean(dist_2*self.cellsize), 2)) 
        print ('    + Mean slope ', np.round(np.nanmean(slope*self.cmask),3))
        print ('    + Mean time delay, months ', np.round(np.nanmean(self.residence_time),1)) 

        #self.n_cw=np.zeros(self.Nsteps); self.p_cw=np.zeros(self.Nsteps)
        print ('    + Runoff m/yr' , np.round(sum(self.runoff)/sim_yrs,4))

        # calculate retention of N and P as a function of distance to water body: Heikkinen et al. 2018 Ecological Engineering 117: 153-164
        # N: y = 15.41*ln(x)-52.701, restrict between [0...100], where x is distance to water body, m
        # P: y = 19.139*ln(x)-61.972, restrict between [0...100], where x is distance to water body, m

        aa = np.maximum(a-self.cellsize/2.0,  np.zeros(np.shape(a))) # centerpoint
        self.Nretention = np.maximum(np.where(aa > 0.0, 15.41*np.log(aa)-52.701, 0.0), np.zeros(np.shape(aa)))
        self.Pretention = np.maximum(np.where(aa > 0.0, 19.139*np.log(aa)-61.972, 0.0), np.zeros(np.shape(aa)))
        self.Nretention = np.minimum(self.Nretention, np.ones(np.shape(aa))*100.)
        self.Pretention = np.minimum(self.Pretention, np.ones(np.shape(aa))*100.)
        
        cubeshape = np.shape(resn['nut']['ntogw'])    #time, x,y        
        self.ndelayed, self.pdelayed = self.make_cube(resn, cubeshape, self.residence_time, self.Nretention, self.Pretention, figs= False)
        
        self.ntosrun = np.array(resn['nut']['ntosrun'])
        self.ptosrun = np.array(resn['nut']['ptosrun'])
        print ('    -> Nutrient export calculation initialized')

    def make_cube(self, resn, cubeshape, rtime, Nretention,Pretention, figs = False):
        """
        resn - output to gw, kg ha-1 month-1
        cubeshape - shape of the output datacube
        rtime - resindence time, months
        Nretention - retention of N during the transport process % from the released (Heikkinen et al 2018), array in gridshape
        Pretention  - retention of P during the transport process % from the released (Heikkinen et al 2018), array in gridshape
        """
        print ('    + Computing transport delay')
        nmonths = cubeshape[0]                                                  # number of time steps (months)

        nmean = np.mean(resn['nut']['ntogw'][0:int(nmonths/2),:,:], axis=0)
        pmean = np.mean(resn['nut']['ptogw'][0:int(nmonths/2),:,:], axis=0)

        ndelay = np.ones(cubeshape)*nmean                                            # output data cube
        pdelay = np.ones(cubeshape)*pmean  
        cix = np.where(np.isfinite(self.cmask))                                 # catchment indices, x,y
        narr = np.array(resn['nut']['ntogw'][:,:,:])          # input data
        parr = np.array(resn['nut']['ptogw'][:,:,:]) 
        for m in range(nmonths):
            rst = np.array((self.residence_time[np.array(self.catchment_ix[0]),np.array(self.catchment_ix[1])]).astype('int32'))
            rst = rst + m                                                       # shifht by time step
            rst[rst > nmonths-1] = nmonths-1                                    # cut locate to the last time step
            ixs = (rst, np.array(self.catchment_ix[0]),np.array(self.catchment_ix[1]))  # output indices, shifted by the residence time
            n_now = narr[m,:,:] * (1.0-Nretention/100.)                                                  # take the current time step in nutrient balance array
            p_now = parr[m,:,:] * (1.0-Pretention/100.)
            outn = n_now[cix] 
            outp = p_now[cix]                                                  # pick the cells within the catcment
            ndelay[tuple(ixs)] = np.array(outn)
            pdelay[tuple(ixs)] = np.array(outp)
            
        
        if figs:
            fig=plt.figure('ndelayed') 
            plt.subplot(331)
            plt.imshow(ndelay[0,:,:]*self.cmask)
            plt.subplot(332)
            plt.imshow(ndelay[1,:,:]*self.cmask)
            plt.subplot(333)
            plt.imshow(ndelay[2,:,:]*self.cmask)
            plt.subplot(334)
            plt.imshow(ndelay[3,:,:]*self.cmask)
            plt.subplot(335)
            plt.imshow(ndelay[4,:,:]*self.cmask)
            plt.subplot(336)
            plt.imshow(ndelay[5,:,:]*self.cmask)
            plt.subplot(337)
            plt.imshow(ndelay[6,:,:]*self.cmask)
            plt.subplot(338)
            plt.imshow(ndelay[7,:,:]*self.cmask)
            plt.subplot(339)
            plt.imshow(ndelay[8,:,:]*self.cmask)
        print ('    + Transport delay done')
            
        return ndelay, pdelay

    def run(self):

        print ('- Starting nutrient export computation')
        print ('     + Nitrogen & Phosphorus')
        nconclist=[]; pconclist=[]; nexport=[]; pexport=[]; datelist=[]
        ncon = self.gwN/np.nansum(self.gw_storage)    #kg/m3
        pcon = self.gwP/np.nansum(self.gw_storage)    #kg/m3
        for d in range(0, self.Nsteps):
            current_date = datetime.datetime(self.dates[d].year, self.dates[d].month, self.dates[d].day)
            n_out = self.runoff[d]*self.catchment_area*ncon                     # outflux in m3, in kg
            p_out = self.runoff[d]*self.catchment_area*pcon
            local_s = self.local_sm[d,:,:]   #res['top']['Sloc'][d,:,:]                                 # local saturation deficit, m    
            gw_vol = np.nansum((self.poros*self.soild - local_s)*self.cellsize**2) # ground water volume m3

            #ntmp = np.nansum(self.n_delay[d,:,:])
            #ptmp = np.nansum(self.p_delay[d,:,:])
            ntmp = np.nansum(self.ndelayed[d,:,:]) + np.nansum(self.ntosrun[d,:,:])
            ptmp = np.nansum(self.pdelayed[d,:,:]) + np.nansum(self.ptosrun[d,:,:])

            n_in = np.nansum(ntmp)*1e-4*self.cellsize**2   #unit kg/ha/day -> convert to kg/gridpoint
            p_in = np.nansum(ptmp)*1e-4*self.cellsize**2   #unit kg/ha/day -> convert to kg/gridpoint
            
            ncon = max((self.gwN + n_in - n_out)/gw_vol, 0.0)  
            pcon = max((self.gwP + p_in - p_out)/gw_vol, 0.0)  
            self.gwN = max(self.gwN + n_in - n_out, 0.0)
            self.gwP = max(self.gwP + p_in - p_out, 0.0)
            #self.gwN = self.gwN + n_in - n_out
            #self.gwP = self.gwP + p_in - p_out
            
            print ('    ', current_date.year, current_date.month,  'N:',np.round(ncon*1e3,3),'P:', np.round(pcon*1e3,3))
            nconclist.append(ncon*1e3); pconclist.append(pcon*1e3)
            nexport.append(n_out/self.catchment_area*1e4);
            pexport.append(p_out/self.catchment_area*1e4)
            datelist.append(current_date)


        print ('Compiling export output')
        data = {'datetime':datelist, 'runoff[m/day]':self.runoff, 'nexport[kg/ha]':nexport,'nconc[mg/l]':nconclist, \
                'pexport[kg/ha]':pexport,'pconc[mg/l]': pconclist}
        dfOut = pd.DataFrame(data=data)                                          # trimming according to time
        dfOut=dfOut.set_index('datetime')
        return dfOut
 