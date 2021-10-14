# -*- coding: utf-8 -*-
"""
Created on Sun Apr 19 11:12:55 2020

@author: alauren
"""
#import matplotlib.pylab as plt
import numpy as np
import pandas as pd
import datetime
import os
import time
import spafhy
from spafhy_nut_utils import get_spafhy_results, get_nutbal_results, get_Motti_results, daily_to_monthly,estimate_imm,description_spafhy
from spafhy_io import write_ncf_res, get_log_area_mp, initialize_netcdf_res, preprocess_soildata, create_catchment, read_FMI_weather, get_clear_cuts, get_clear_cuts_times, get_clear_cuts_area, read_AsciiGrid, get_clear_cuts_pery,write_AsciiGrid
from spafhy_parameters_default import soil_properties, parameters, nutpara, immpara
from nutrient_balance import Grid_nutrient_balance
from nutrient_export import Export
import spafhy_stand as stand
import seaborn as sns
import matplotlib.pylab as plt
import argparse
import json
import xarray_extras
import xarray as xr
from osgeo import gdal, ogr, osr
import rioxarray

sns.set()

parser = argparse.ArgumentParser(description='Nutspafhy for defined catchments & scenarios')
parser.add_argument('--c', type=int,
                   help='catchment to be calculated, use integer')
parser.add_argument('--s',
                   help='logging scenario to be calculated', type=str)
parser.add_argument('--f',
                   help='weather forcing scenario to be used, VKe=vakio10v keski 2001, VKu=vakio10 kuiva 2002, VKo= vakio10kostea, VLKo=vakio10 lammin kostea, Had=HadGEM, MPI=MPI-ESM', type=str)

parser.add_argument('--start',
                   help='startdate e.g. as 2018-01-01', type=str)
parser.add_argument('--end',
                   help='enddate e.g. as 2028-12-31', type=str)

#parser.add_argument('catchment', metavar='N', type=int, nargs='+',
#                   help='catchment to be calculated, use number')


args=parser.parse_args()

hydrology_calc = True
nut_balance_calc = True
export_calc = True
clearcut_area = True
remove_nc = True

folder =r'/content/PURUVESI_input_data/' 
catchments =[f for f in sorted(os.listdir(folder)) if f.startswith("04_184")]
cnum=args.c
catchments = [catchments[cnum]] 
scenarios =[(args.s)] 

forcinglist={'Had':r'/content/PURUVESI_input_data_vmi/Puruvesi_saa2005_2099HadGEM_leap_days.csv',
             'MPI':r'/content/PURUVESI_input_data_vmi/Puruvesi_saa2005_2099MPI_leap_days.csv',
             'normi':r'/content/PURUVESI_input_data_vmi/puruvesi_saa_19810101-2020062.csv',}

forcingfile=forcinglist[args.f]
forci=args.f

startdate=args.start
enddate=args.end

print(catchments)
print(scenarios)
print(type(startdate))
print(startdate)

for i in range(len(catchments)):
    for j in range(len(scenarios)):
        if not os.path.exists(folder+catchments[i]+"/"+scenarios[j]):
            os.makedirs(folder+catchments[i]+"/"+scenarios[j])


def nsy(iN, iP, pgen, pcpy, pbu, ptop, pnut, psoil, gisdata, clear_cuts, forcing, outfold = None):
                       
    """run hydrology"""
    if hydrology_calc: 
        """ create SpatHy object """
                
        spa = spafhy.initialize(pgen, pcpy, pbu, ptop, psoil, gisdata, cpy_outputs=False, 
                         bu_outputs=False, top_outputs=False, flatten=True)
    
        """ create netCDF output file """
       
        dlat, dlon = np.shape(spa.GisData['cmask'])
        end_date = datetime.datetime.strptime(pgen['end_date'], '%Y-%m-%d')
        start_date = datetime.datetime.strptime(pgen['start_date'], '%Y-%m-%d')           
        Nsteps = (end_date-start_date).days
    
        ncf, ncf_file = spafhy.initialize_netCDF(ID=spa.id, fname=spa.ncf_file, lat0=spa.GisData['lat0'], 
                                                 lon0=spa.GisData['lon0'], dlat=dlat, dlon=dlon, dtime=Nsteps)
        
     
        """ read forcing data and catchment runoff file """
        FORC = forcing.copy()
        FORC['Prec'] = FORC['Prec'] / spa.dt                            # mms-1
        FORC['U'] = 2.0                                                 # use constant wind speed ms-1
        
    
        for k in range(0, Nsteps):                                      # keep track on dates
            current_date = datetime.datetime.strptime(pgen['start_date'],'%Y-%m-%d').date() + datetime.timedelta(days=k)
            current_datetime = datetime.datetime.strptime(pgen['start_date'],'%Y-%m-%d') + datetime.timedelta(days=k)
    
            if k%30==0: print('step: ' + str(k), current_date.year, current_date.month, current_date.day)
    
            forc= FORC[['doy', 'Rg', 'Par', 'T', 'Prec', 'VPD', 'CO2','U']].iloc[k]
            spa.run_timestep(forc, current_datetime, ncf=ncf)
    
            if current_date in clear_cuts.keys():
                print ('     +Clear cut ', current_date)            
                spa.set_clear_cut(clear_cuts[current_date])
                
        
        ncf.close()                                                             # close output file
        print (ncf_file, ' closed') 
    
    #*********************************************************************************************************            
    if nut_balance_calc:
        
        soildata = preprocess_soildata(pbu, psoil, gisdata['soilclass'], gisdata['cmask'], pgen['spatial_soil'])
                
        TAmr =stand.TAmr(forcing['T'])                                          # growing season air temperature
        sp_res = get_spafhy_results(pnut['spafhyfile'])                      # retrieve results from netCDF file
        
        
        Wm, lsm, Rm, Dm, Dretm, Infil,s_run,Tm, ddsm, Prec = daily_to_monthly(sp_res, forcing, pnut, balance=True)  # Change hydrology to monthly values
        lat= float(gisdata['loc']['lat'])                                    # coordinates in Spathy
        lon= float(gisdata['loc']['lon'])    
        motti = get_Motti_results(pnut['mottisims'], lat, lon)          # retrieve motti parameters
        Nsteps = np.shape(Rm)[0]                                         # number of months in the smulation
        
        print ('Nsteps', Nsteps)
        print ('***************************')
        """ Create nutrient balance object """
        
        nutSpafhy = Grid_nutrient_balance(Tm.index, pbu, pgen, pnut, gisdata, \
                                          soildata, motti, Wm, ddsm, lat, lon, iN=iN, iP=iP)
        
        for k in range(0, Nsteps):                             
            current_date = datetime.date(ddsm.index[k].year, ddsm.index[k].month,ddsm.index[k].day) # keep track on the date
            Ta = Tm[k]                                                                             # Mean monthly air temperature
            Pm = Prec[k]                                                                    # Precipitation monthly, mm
            Infi = Infil[k,:,:]
            local_s= lsm[k,:,:]                                                          # Mean monthly local saturation deficits in the catchment, m
            Wliq = Wm[k,:,:]                                                             # Mean monthly water content in root layer [m3 m-3]
            draindown = Dm[k,:,:]                                                        # Sum of monthly percolation of water from root zone to down m dt-1
            drainretflow = Dretm[k,:,:]                                                  # Sum of monthly returnflow, vertical upwards [mm dt-1], 
            Q_s = s_run[k,:,:]                                                           # Sum of monthly runoff, [mm dt-1]                     
            dt = (current_date - datetime.date(current_date.year, current_date.month, 1)).days  # length of time step in days (different lenght in different months)

            """ Run time step"""
            nutSpafhy.run_timestep(Ta, Pm, Infi, Q_s, TAmr, Wliq, local_s, draindown, drainretflow, dt)        
        
            print ('step: ' + str(k), current_date.year, current_date.month,current_date.day,) 
            print ('   -> Concentration mg L-1 N:', np.round(np.nanmean(nutSpafhy.nconc),3), 'P:', np.round(np.nanmean(nutSpafhy.pconc),3)) 
            
            print ('date: ', current_date) 
             
            if current_date in clear_cuts.keys():
                print ('     +Clear cut ', current_date)            
                nutSpafhy.set_clear_cut(clear_cuts[current_date])        

        nutSpafhy.close_ncf()
                
    #**************************************************************************************************************            
    if export_calc:
        
        soildata = preprocess_soildata(pbu, psoil, gisdata['soilclass'], gisdata['cmask'], pgen['spatial_soil'])            
        
        sp_res = get_spafhy_results(pnut['spafhyfile'])    
        _, lsm, Rm, _, _, _, _,_, ddsm, _ = daily_to_monthly(sp_res, forcing, pnut, balance=False)  # Change hydrology to monthly values: s_local and runoff 
        
        nsp_res = get_nutbal_results(f=pgen['output_folder']+ pnut['nutspafhyfile'])
    
        """ Initialize and run """                
        ex = Export(pgen, lsm, Rm, ddsm, nsp_res, gisdata,soildata)       
        nsp_res.close()
        dfOut = ex.run()
        print ('Computation done')
        
        df_annual_load = dfOut.resample('Y', convention ='end').sum()
        df_annual_conc = dfOut.resample('Y', convention ='end').mean()
     
        print ('******** Export load, kg ha-1 yr-1 **************') 
        print ( np.round(df_annual_load[['nexport[kg/ha]','pexport[kg/ha]' ]] ,3))  
        
        print ('******** Concentration, mg l-1 *****************')
        print (np.round(df_annual_conc[['nconc[mg/l]','pconc[mg/l]']], 3))
     
                    
        print ('')
        print ('COMPLETED: ', cat  +'  '+ scen)
        return dfOut, ex.Nretention, ex.Pretention
        
    #**************************************************************************************************************            

    if clearcut_area:
        if not os.path.exists(pgen['output_folder']+"cc_area/"):
                    os.makedirs(pgen['output_folder']+"cc_area/")
        area_time_of_cc=get_clear_cuts_times(pgen, gisdata['cmask_cc'])
        area_time_of_cc.to_csv(pgen['output_folder']+"cc_area/"+cat+"_"+scen+"_time_clearcut.csv", sep=',',index=False)

    #**************************************************************************************************************            

    
                



#--------- setup general parameters -------------------------------------
for cat in catchments:  
    for scen in scenarios: 
        print('************************************************')
        print ('****** Catchment: ', cat ,' ', scen, ' *******************')
        print('************************************************')
        pgen,pcpy,pbu,ptop=parameters(cat,scen)
        pgen['forcing_file']=forcingfile
        if startdate is None:
            print('start date is:', pgen['start_date'])
            startdate=pgen['start_date']
        else:
            pgen['start_date']= startdate
        if enddate is None:
            print('end date is:', pgen['end_date'])
            enddate=pgen['end_date']
        else:
            pgen['end_date']= enddate
        pgen['output_folder']=r'/scratch/project_2002470/nutspafhy_puruvesi/output/'+cat+r'/'+ forci + '_' +scen + '_' + startdate + r'/'
        outfold=pgen['output_folder']
        pgen['ncf_file']=r'/scratch/project_2002470/nutspafhy_puruvesi/output/'+cat+r'/'+ forci +'_'+ scen + '_' + startdate + r'/'+'ch.nc'
        psoil = soil_properties()
        pnut = nutpara(pgen, cat, scen)
        gisdata = create_catchment(pgen, fpath=pgen['gis_folder'], plotgrids=False, plotdistr=False)
        iN, iP = estimate_imm(gisdata)
        logfile = pgen['output_folder']+'logfile.txt'
        clear_cuts = get_clear_cuts(pgen, gisdata['cmask_cc'], gisdata['cmask'])    
        if not os.path.exists(pgen['output_folder']):
            os.makedirs(pgen['output_folder'])
        FORC = read_FMI_weather(pgen['catchment_id'],
                                        pgen['start_date'],
                                        pgen['end_date'],
                                        sourcefile=pgen['forcing_file'])
                
        dfOut, nret, pret = nsy(iN, iP, pgen, pcpy, pbu, ptop, pnut, psoil, gisdata, clear_cuts, FORC, outfold=outfold)
        dfOut.to_csv(pgen['output_folder']+cat+"_"+scen+"_export.csv")
        
        if remove_nc:
            if os.path.exists(pgen['ncf_file']):
                sp_res = get_spafhy_results(pnut['spafhyfile'])
                Wliqm,Slocm,ETm, Slocmcc,Wliqmcc, Mbem, Mbecpym = description_spafhy(sp_res, pgen, gisdata['cmask'], pnut, startd=startdate, endd=enddate)
                Nsteps = np.shape(Wliqm)[0]                                        # number of months in the smulation
                wliql=[];slocl=[];etl=[]; slocccl=[]; wliqccl=[]; mbem=[]; mbecpym=[]; #runoffm=[];
                for k in range(0, Nsteps):
                    wliql.append(np.nanmean(Wliqm[k,:,:]))
                    slocl.append(np.nanmean(Slocm[k,:,:]))
                    etl.append(np.nanmean(ETm[k,:,:]))
                    mbem.append(np.nanmean(Mbem[k,:,:]))
                    mbecpym.append(np.nanmean(Mbecpym[k,:,:]))
                    if Slocmcc.any():  
                        slocccl.append(np.nanmean(Slocmcc[k,:,:]))
                        wliqccl.append(np.nanmean(Wliqmcc[k,:,:]))
                    else:
                        slocccl.append(0)
                        wliqccl.append(0)
                        #print("no clearcuts")
                dfdesc = pd.DataFrame(data={"Wliqm": wliql,"Slocm": slocl,"ETm": etl,"Slocmcc":slocccl,"Wliqmcc":wliqccl, "Mbem":mbem, "Mbecpym":mbecpym})
                dfdesc.to_csv(pgen['output_folder']+'spafhy_monthly_desc.csv', sep=',',index=True)
                os.remove(pgen['ncf_file'])
                #print("spafhyfile not removed") 
            else:
                print("The spafhyfile does not exist")
            if os.path.exists(pgen['output_folder']+ 'nutspafhy.nc'):
                ds0 = xr.open_dataset(pgen['output_folder']+ 'nutspafhy.nc')
                ds0['dtime'] = ds0['time']
                dtime = ds0['dtime']
                dlon = ds0.variables['lon'][:]
                dlat = ds0.variables['lat'][:]
                ds0.close()
                dsnut = xr.open_dataset(pgen['output_folder']+ 'nutspafhy.nc', group= 'nut', decode_times=True)
                noutg = (1-nret/100.)*(xr.DataArray(dsnut['ntogw'], coords={'dtime':dtime, 'dlon':dlon, 'dlat':dlat})).sum(axis=(0))*gisdata['cmask']   #'ntogw' 'ntosrun'
                nouts = (xr.DataArray(dsnut['ntosrun'], coords={'dtime':dtime, 'dlon':dlon, 'dlat':dlat})).sum(axis=0)*gisdata['cmask']   #'ntogw' 'ntosrun'
             
                poutg = (1-pret/100.)*(xr.DataArray(dsnut['ptogw'], coords={'dtime':dtime, 'dlon':dlon, 'dlat':dlat})).sum(axis=(0))*gisdata['cmask']   #'ntogw' 'ntosrun'
                pouts = (xr.DataArray(dsnut['ptosrun'], coords={'dtime':dtime, 'dlon':dlon, 'dlat':dlat})).sum(axis=0)*gisdata['cmask']   #'ntogw' 'ntosrun'

                nsum = noutg+nouts
                psum = poutg+pouts
                nsum.rio.set_spatial_dims('dlon', 'dlat')
                psum.rio.set_spatial_dims('dlon', 'dlat')
                nsum.rio.to_raster(pgen['output_folder']+'nsum.tif')
                psum.rio.to_raster(pgen['output_folder']+'psum.tif')
                drv = gdal.GetDriverByName('AAIGrid')
                out_asc = os.path.join(pgen['output_folder'], "nsum.asc")
                intif = gdal.Open(pgen['output_folder']+'nsum.tif')
                outasc=drv.CreateCopy(out_asc,intif)
                intif = None
                outasc = None             
                out_asc = None             
                out_asc = os.path.join(pgen['output_folder'], "psum.asc")
                intif = gdal.Open(pgen['output_folder']+'psum.tif')
                outasc=drv.CreateCopy(out_asc,intif)
                intif = None
                outasc = None             
                out_asc = None             

                #nsp_res = get_nutbal_results(pgen['output_folder']+ 'nutspafhy.nc')  # hakee netcdfän
                #nsum =   np.sum(nsp_res['nut']['ntogw'], axis=0) + np.sum(nsp_res['nut']['ntosrun'], axis=0) # tahan vois lisata et vasta ekan vuoden jalkeen ja lisaks pitais skaalata viim tuloskarttaa tehdessa
                #psum =   np.sum(nsp_res['nut']['ptogw'], axis=0) + np.sum(nsp_res['nut']['ptosrun'], axis=0)
                
                #write_AsciiGrid(pgen['output_folder']+ 'nsum.asc',nsum,gisdata['info'],fmt='%1.4f')
                #write_AsciiGrid(pgen['output_folder']+ 'psum.asc',psum,gisdata['info'],fmt='%1.4f')
                os.remove(pgen['output_folder']+ 'nutspafhy.nc')
                percen = 75
                nhotmask = np.ma.masked_less(nsum, np.nanpercentile(nsum, percen))
                photmask = np.ma.masked_less(psum, np.nanpercentile(psum, percen))
                fig = plt.figure(num='hot', figsize=(12,10))
                plt.subplot(221)
                nsum.plot.hist(bins=50)
                plt.subplot(222)
                plt.title('N export')
                plt.imshow(nhotmask)
                plt.colorbar()
                #plt.imshow(acut, alpha=0.2)

                plt.subplot(223)
                psum.plot.hist(bins=50)
                plt.subplot(224)
                plt.title('P export')
                plt.imshow(photmask)
                plt.colorbar()
                #plt.imshow(acut, alpha=0.2)
                outfigp = pgen['output_folder']+str(percen)+"_hotspots.png"
                plt.savefig(outfigp, dpi=300)
                #print("nutspafhyfile not removed")
            else:
                print("The nutspafhyfile does not exist")
        print ('*******************************') 
        print ('COMPLETED: ', cat  +'  '+ scen)
         


