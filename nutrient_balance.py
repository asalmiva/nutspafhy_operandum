# -*- coding: utf-8 -*-
"""
Created on Sat Feb 23 17:18:08 2019

@author: lauren
"""
import numpy as np
from scipy.interpolate import interp1d
from spafhy_nut_utils import initialize_netCDF_nut, local_s_to_gwl_functions, local_s_to_gwl
import datetime
import spafhy_stand as stand
import matplotlib.pylab as plt

class Grid_nutrient_balance():
    def __init__(self, timeind, pbu, pgen, pnut, gisdata, soildata, motti, Wliq, ddsm, lat, lon, iN=None, iP=None):
        """
        Dimensions of netCFD (time, longitude, latitude) ie. time, y, x
        """
        print ('****** Initializing Grid nutrient balance instance ********')        
        self.tstep = 0        
        self.timeind=timeind        
        self.pbu = pbu
        self.pgen = pgen
        self.gisdata = gisdata
        self.gridshape = np.shape(self.gisdata['cmask'])
        self.soildata = soildata
        self.motti = motti
        self.lat = lat
        self.lon = lon
        self.ddsm = ddsm
        #self.dt = pgen['dt'] / 86400.  # from seconds to days
        self.id = pgen['catchment_id']
        self.ncf,self.outf=initialize_netCDF_nut(self.id,gisdata,self.timeind, 
                                   fpath=pgen['output_folder'], fname=pnut['nutspafhyfile'])                 # create output netCDF-file 

        start_date = datetime.datetime.strptime(self.pgen['start_date'], '%Y-%m-%d')
        end_date = datetime.datetime.strptime(self.pgen['end_date'], '%Y-%m-%d')    
        self.simtime = (end_date-start_date).days / 365.                                 # expressed in yrs    

        self.s_to_gwl_dict =local_s_to_gwl_functions()                                   #interpolation functions (storage->gwl) for the different soil types is soil mask

        self.nUp, self.pUp, expected_yield = stand.uptake(self.gisdata, self.motti, self.simtime)        # total N and P uptake in the simulation time [kg ha-1], and expected yield in m3/ha/simulation time
        print ('    + Max N uptake', np.round(np.nanmax(self.nUp),2) , 'kg in ', np.round(self.simtime,2), 'yrs')
        print ('    + Max P uptake', np.round(np.nanmax(self.pUp),2) , 'kg in ', np.round(self.simtime,2), 'yrs')
        #****************************
        #self.close_ncf()
        self.nUp_gv, self.pUp_gv = stand.understory_uptake(self.lat, self.lon, self.ddsm[-1]/self.simtime, 
                                                           gisdata, expected_yield, self.simtime) #total N and P uptake by ground vegetation in kg ha-1 during the simulation time
        #************************************
        
        self.uprate = stand.ddgrowth(self.ddsm ,self.motti['Tsum'], self.simtime)    # composes daily uptake using temperature sum accumulation

        fconc = interp1d([6673395., 7639733.], [1.0, 0.1], fill_value='extrapolate')    # N concentration in atmospheric deposition as a function of latitude [mg L-1]
        nconcini = fconc(self.lat)*0.05                                                       # 0.2 #mg/l initial concentration in soil water
        self.N_conc_depo = fconc(self.lat)                                                   #mg/l
        self.P_conc_depo = fconc(self.lat)*0.                                               #mg/l
        pconcini = self.P_conc_depo*0.05
        #pconcini = 0.2 /10.                                                            # initial concentration in soil water mg/l
    
        #nconcini = (1.5361*tsum - 1111.8)/1000.                                      # Sakarin ja Mikan aineisto
        #pconcini = (0.06*motti['Tsum'] - 45.745)/1000. *0.25
        print ('    + Initial N concentration ', np.round(nconcini,3), ' mg l-1')
        print ('    + Initial P concentration ', np.round(pconcini,3), ' mg l-1')
        print ('    + Rain N concentration ', np.round(self.N_conc_depo,3), ' mg l-1')
        print ('    + Rain P concentration ', np.round(self.P_conc_depo,3), ' mg l-1')
        print ('    + Temperature sum ', np.round(self.ddsm[-1]/self.simtime), 'degree days')
        print ('    + Motti temperature sum ', motti['Tsum'])
        self.nrel = self.nrelm = self.nrelp = np.zeros(np.shape(gisdata['cmask']))                # initiate N release variables 
        self.prel = self.prelm = self.prelp = np.zeros(np.shape(gisdata['cmask']))                # initiate P release variables 
    
        fimmn_peat = interp1d([700., 1500.],[0.89, 0.84], fill_value= 'extrapolate')
        fimmp_peat = interp1d([700., 1500.],[0.9, 0.9], fill_value= 'extrapolate')
        fimmn_min = interp1d([700., 1500.],[0.89+0.08, 0.84+0.08], fill_value= 'extrapolate')
        fimmp_min = interp1d([700., 1500.],[0.9+0.08, 0.9+0.08], fill_value= 'extrapolate')

        if iN is None:
            self.imm_n_peat = 0.89 #fimmn_peat(ddsm[-1]/self.simtime)     #immobilization as a function temperature sum
            self.imm_n_min = 0.95 #fimmn_min(ddsm[-1]/self.simtime)     #immobilization as a function temperature sum
        else:
            #deltaN = 0.05  #0.06:0.961
            self.imm_n_peat = iN[0]
            self.imm_n_min = iN[1] #+ deltaN
        if iP is None:
            self.imm_p_peat = 0.94 #fimmp_peat(ddsm[-1]/self.simtime)
            self.imm_p_min = 0.98 #fimmp_min(ddsm[-1]/self.simtime)
        else:
            #deltaP = 0.0
            self.imm_p_peat = iP[0]
            self.imm_p_min = iP[1]
        #self.imm_n_peat = 0.89 #fimmn_peat(ddsm[-1]/self.simtime)     #immobilization as a function temperature sum
        #self.imm_n_min = 0.95 #fimmn_min(ddsm[-1]/self.simtime)     #immobilization as a function temperature sum

        #self.imm_p_peat = 0.94 #fimmp_peat(ddsm[-1]/self.simtime)
        #self.imm_p_min = 0.98 #fimmp_min(ddsm[-1]/self.simtime)
        
        print ('    + Temp sum, observed:', int(self.ddsm[-1]/self.simtime))
        print ('    + N, P Immob peat:', np.round(self.imm_n_peat,2), np.round(self.imm_p_peat,2))
        print ('    + N, P Immob min:', np.round(self.imm_n_min,2), np.round(self.imm_p_min,2))
        
        self.sdepth = pbu['depth']                                           # reference depth of computation
        print ('    + Reference depth, m', self.sdepth)
    
        self.nsto = np.array(nconcini*1e-6 * np.ones(np.shape(self.nrelp))*1e4*pbu['depth']*Wliq[0][:][:]*1e3)   # ini storage kg/ha
        self.psto = np.array(pconcini*1e-6 * np.ones(np.shape(self.prelp))*1e4*pbu['depth']*Wliq[0][:][:]*1e3)   # ini storage kg/ha
        #self.nsto = nconcini*1e-6 * np.ones(np.shape(self.nrelp))*1e4*pbu['depth']*res['bu']['Wliq'][0]*1e3   # ini storage kg/ha
        #self.psto = pconcini*1e-6 * np.ones(np.shape(self.prelp))*1e4*pbu['depth']*res['bu']['Wliq'][0]*1e3   # ini storage kg/ha
        #nsto = nconcini*1e-6 * np.ones(np.shape(nrelp))*1e4*self.sdepth*res['bu']['Wliq'][0]*1e3   # ini storage kg/ha
        #psto = pconcini*1e-6 * np.ones(np.shape(prelp))*1e4*self.sdepth*res['bu']['Wliq'][0]*1e3   # ini storage kg/ha
        print ('    + Mean initial N storage ', np.round(np.nanmean(self.nsto),2), ' kg ha-1')
        print ('    + Mean initial P storage ', np.round(np.nanmean(self.psto),2), ' kg ha-1')
        print ('      -> Nutrient balance initialized')

        self.ix=np.where(np.isfinite(self.gisdata['cmask']))
        
        
    def run_timestep(self, Ta, Prec, Infil, Q_s, tamr, Wliq, local_s, draindown, returnflow, dt):
        """
        Ta air temperature
        Prec precipitation mm / month
        Infil infiltration mm/month
        Q_s surface runoff from each pixel, mm/month
        tamr growing season air temperature
        Wliq liquid water in the root zone
        local_s local water deficit
        draindown percolation of water from root zone to ground water, m day-1
        returnflow upward flow of water typical near the brook m day-1 
        dt is timestep in days
        """
        #draindown = draindown + returnflow
        
        gwl=local_s_to_gwl(self.gisdata['soil'],local_s,self.s_to_gwl_dict)                    # ground water level in the area 
        nrelp, prelp, ixp = stand.peatRespiration(self.gisdata,Ta, self.imm_n_peat, self.imm_p_peat, TASmr=tamr, gwl=gwl, dt=dt)  # N release in peat respiration [kg N ha-1 dt-1]
        nrelm, prelm, ixm = stand.soilRespiration(self.gisdata, self.soildata, Ta, Wliq, self.imm_n_min, self.imm_p_min, dt= dt)        # N release in mineral soil respiration [kg N ha-1 dt-1]
        self.nrel[ixp]=nrelp[ixp]; self.nrel[ixm]=nrelm[ixm]                                    # N compile whole catchment
        self.prel[ixp]=prelp[ixp]; self.prel[ixm]=prelm[ixm]                                    # P compile whole catchment

        #self.nrel = self.nrel + (Prec * 1e4 * self.N_conc_depo*1e-6 )   #decomposition + deposition kg/ha/month
        #self.prel = self.prel + (Prec * 1e4 * self.P_conc_depo*1e-6 )

        self.nrel[ixm] = self.nrel[ixm] + (Infil[ixm] * 1e4 * self.N_conc_depo*1e-6 )   #decomposition + deposition kg/ha/month
        self.nrel[ixp] = self.nrel[ixp] + (Infil[ixp] * 1e4 * self.N_conc_depo*1e-6 )   #decomposition + deposition kg/ha/month
        
        self.prel[ixm] = self.prel[ixm] + (Infil[ixm] * 1e4 * self.P_conc_depo*1e-6 )   #decomposition + deposition kg/ha/month
        self.prel[ixp] = self.prel[ixp] + (Infil[ixp] * 1e4 * self.P_conc_depo*1e-6 )   #decomposition + deposition kg/ha/month
   
        self.nsto=stand.nutBalance(self.nsto, self.uprate[self.tstep]*self.nUp, self.uprate[self.tstep]*self.nUp_gv, self.nrel)      # kg ha-1
        self.psto=stand.nutBalance(self.psto, self.uprate[self.tstep]*self.pUp, self.uprate[self.tstep]*self.pUp_gv, self.prel)       # kg ha-1

        nconckgm3 = self.nsto/(1e4*Wliq*self.sdepth)    #kg/m3
        pconckgm3 = self.psto/(1e4*Wliq*self.sdepth)    #kg/m3        
        #nconckgm3 = self.nsto/(1e4*Wliq*self.sdepth + returnflow*1e-3)    #kg/m3, returnflow in mm -> m
        #pconckgm3 = self.psto/(1e4*Wliq*self.sdepth + returnflow*1e-3)    #kg/m3        
        Ntogw = 1e4*draindown*nconckgm3     #kg/ha/day
        Ptogw = 1e4*draindown*pconckgm3     #kg/ha/day
        NtoSrun = 1e4*Q_s*nconckgm3     #kg/ha/day
        PtoSrun = 1e4*Q_s*pconckgm3     #kg/ha/day
        
        #Ntogw = 1e4*(np.fmax(draindown,returnflow))*nconckgm3     #kg/ha/day
        #Ptogw = 1e4*(np.fmax(draindown,returnflow))*pconckgm3     #kg/ha/day
        #Ntogw = 1e4*(draindown +returnflow)*nconckgm3     #kg/ha/day
        #Ptogw = 1e4*(draindown +returnflow)*pconckgm3     #kg/ha/day

        #Nreturn = 1e4*returnflow*nconckgm3     #kg/ha/day
        #Preturn = 1e4*returnflow*pconckgm3     #kg/ha/day
        self.nsto = np.maximum(self.nsto-Ntogw-NtoSrun,0.0)
        self.psto = np.maximum(self.psto-Ptogw-PtoSrun,0.0)

        self.nconc = self.nsto*1e6/(1e4*self.sdepth*Wliq*1e3)                                 # mg/l          
        self.pconc = self.psto*1e6/(1e4*self.sdepth*Wliq*1e3)                                 # mg/l          

        # Ntogw, Ptogw  handle the ground water flow
        # find the volume of the whole catchment scale water volume to get the storage
        # Create variables for the return flow and direct it to instantanously to runoff

        self.ncf['nut']['nsto'][self.tstep,:,:]=self.nsto
        self.ncf['nut']['nrel'][self.tstep,:,:]=self.nrel
        self.ncf['nut']['nup'][self.tstep,:,:]=self.uprate[self.tstep]*self.nUp
        self.ncf['nut']['nup_gv'][self.tstep,:,:]=self.uprate[self.tstep]*self.nUp_gv
        self.ncf['nut']['nconc'][self.tstep,:,:]=self.nconc
        self.ncf['nut']['psto'][self.tstep,:,:]=self.psto
        self.ncf['nut']['prel'][self.tstep,:,:]=self.prel
        self.ncf['nut']['pup'][self.tstep,:,:]=self.uprate[self.tstep]*self.pUp
        self.ncf['nut']['pup_gv'][self.tstep,:,:]=self.uprate[self.tstep]*self.pUp_gv
        self.ncf['nut']['pconc'][self.tstep,:,:]=self.pconc
        self.ncf['nut']['ntogw'][self.tstep,:,:]=Ntogw
        self.ncf['nut']['ptogw'][self.tstep,:,:]=Ptogw
        self.ncf['nut']['ntosrun'][self.tstep,:,:]=NtoSrun
        self.ncf['nut']['ptosrun'][self.tstep,:,:]=PtoSrun

        #self.ncf['nut']['nretflow'][self.tstep,:,:]=Nreturn
        #self.ncf['nut']['pretflow'][self.tstep,:,:]=Preturn
        
        self.tstep+=1

    def _to_grid(self, x):
        """
        converts variable x back to original grid for NetCDF outputs
        """
        if self.ix:
            a = np.full(self.gridshape, np.NaN)
            a[self.ix] = x
        else: # for non-flattened, return
            a = x
        return a

    def set_clear_cut(self, cut):
        cut2 = self._to_grid(cut)
        ixc = np.greater(cut2, np.zeros(np.shape(cut2)))
        self.gisdata['vol'][ixc] = 5.0
        self.gisdata['hc'][ixc] = 2.0
        self.gisdata['LAI_pine'][ixc] = 0.1
        self.gisdata['LAI_spruce'][ixc] = 0.1
        self.gisdata['age'][ixc]=1.0
        #nUp, pUp,expected_yield = stand.uptake(self.gisdata, self.motti, self.simtime-self.tstep/365.)                    # total N and P uptake in the simulation time [kg ha-1]
        nUp, pUp, expected_yield = stand.uptake(self.gisdata, self.motti, self.simtime-self.tstep/12.)                    # total N and P uptake in the simulation time [kg ha-1]
        #here ground vegetation uptake...
        nUp_gv, pUp_gv = stand.understory_uptake(self.lat, self.lon, self.ddsm[-1]/self.simtime, \
                                                 self.gisdata, expected_yield, self.simtime-self.tstep/12.) #total N and P uptake by ground vegetation in kg ha-1 during the simulation time

        down = 1.0
        nUp_gv, pUp_gv = down*nUp_gv, down*pUp_gv 
        self.nUp[ixc]= nUp[ixc]
        self.pUp[ixc] = pUp[ixc]
        self.nUp_gv[ixc]= nUp_gv[ixc]
        self.pUp_gv[ixc] = pUp_gv[ixc]
        
        try:
            print ('    + In clear-cut mean stand N uptake', np.round(np.nanmean(nUp[ixc]),2) , 'kg ha-1 in ', np.round(self.simtime-self.tstep/12,2), 'yrs')
            print ('    + In clear-cut mean gv N uptake', np.round(np.nanmean(nUp_gv[ixc]),2) , 'kg  ha-1')            
            print ('    + In clear-cut mean stand P uptake', np.round(np.nanmean(pUp[ixc]),2) , 'kg ha-1 in ', np.round(self.simtime-self.tstep/12,2), 'yrs')
            print ('    + In clear-cut mean gv P uptake', np.round(np.nanmean(pUp_gv[ixc]),2) , 'kg ha-1 ')

        except:
            print ('Clear-cut array empty')
            
    def close_ncf(self):
        self.ncf.close()
        print ('nutspafhy ncf file closed')