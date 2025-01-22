CASE='SHAREMED'
#CASE='HARMONIA'
divide_qtity_released_for_simplification = 3.0
(R0,Rf)=(1,360) # 360
(R0,Rf)=(DAY0,DAYF) # 360
NF = "{:02d}".format(int(R0/30.)+1)
NFf = "{:02d}".format(int(Rf/30.)+1)
(iskip,jskip)=(2,2)
apply_weighting_coeffs=True

run_parallel=True
plot_figures=False
show_figures=False

import os,sys
import numpy as np
LTRANS_tools = os.getenv('LTRANS_tools')
sys.path.append(LTRANS_tools+'/pymodules')
if(CASE=='HARMONIA'): grid_fname='/share/scratch/athibaut/LTRANS_Zlev/SIM/input/GridforLTRANS-MANTIS-c1.nc'
elif(CASE=='SHAREMED'): grid_fname='/share/scratch/athibaut/LTRANS_Zlev/SIM/input/GridforLTRANS-CADEAU_760p-c1.nc'
else: 
  print('CASE ',CASE,'not implemented')
  sys.exit()
if(CASE=='HARMONIA'):
  hours=[0,3,6,9,12,15,18,21,24,2*24,3*24,4*24,5*24,6*24,7*24,8*24,9*24,10*24-2]
  input_nc=       ['/share/scratch/athibaut/LTRANS_tools/pyscripts/DATA/DATA_oilspill_inputnc_'+str(release)+'.txt' for release in range(R0,Rf+1) ]   #
  input_OilProps= ['/share/scratch/athibaut/LTRANS_tools/pyscripts/DATA/DATA_oilspill_oilprops_'+str(release)+'.txt' for release in range(R0,Rf+1) ]  #
  iniparloc_files=['/share/scratch/athibaut/LTRANS_tools/pyscripts/DATA/DATA_oilspill_iniparloc_'+str(release)+'.txt' for release in range(R0,Rf+1) ] #
elif(CASE=='SHAREMED'): 
  hours=np.arange(0,222,3)
  input_nc=       ['/share/scratch/claurent/LTRANS_tools_adrien/pyscripts/SHAREMED/LISTFILES/DATA_oilspill_inputnc_'+str(release)+'.txt' for release in range(R0,Rf+1) ]   #
  input_OilProps= ['/share/scratch/claurent/LTRANS_tools_adrien/pyscripts/SHAREMED/LISTFILES/DATA_oilspill_oilprops_'+str(release)+'.txt' for release in range(R0,Rf+1) ]  #
  iniparloc_files=['/share/scratch/claurent/LTRANS_tools_adrien/pyscripts/SHAREMED/LISTFILES/DATA_oilspill_iniparloc_'+str(release)+'.txt' for release in range(R0,Rf+1) ] #
numgroups = len(input_nc)

#-----------------------------------------------------
import lagrangiantools as LT
if(plot_figures):
  if(not show_figures):
    import matplotlib
    matplotlib.use('Agg')
  import matplotlib.pyplot as plt
#-----------------------------------------------------

PARALLEL=LT.Parallel(run_parallel)

GRID=LT.LTRANSgrid(grid_fname,PARALLEL)

sum_weights_and_releases=0.0
for group in range(0,numgroups):
  ######################################################## 
  iniparloc=LT.csv_to_pd_dataframe(iniparloc_files[group],PARALLEL)
  if(apply_weighting_coeffs):
    weighting_coeffs=np.zeros((len(iniparloc.contents)),dtype=float)
    for f in range(len(weighting_coeffs)):
        if(CASE=='HARMONIA'): weighting_coeffs[f]=float(iniparloc.contents[f].columns[6])  # FOR HARMONIA SIMULATIONS ON MANTIS DOMAIN
        elif(CASE=='SHAREMED'): weighting_coeffs[f]=float(iniparloc.contents[f].columns[5])  # FOR SHAREMED SIMULATIONS ON CADEAU DOMAIN
  else:
    weighting_coeffs=np.ones((len(iniparloc.contents)),dtype=float)
  sum_weights_and_releases+=np.sum(weighting_coeffs)
  ######################################################## 

  DATA=LT.LTRANSdata(input_nc[group],PARALLEL,store_all=True)
  if(group==0):
    MPATCH=LT.MeshPatch(GRID,DATA,PARALLEL)
    MPATCH.def_meshpatch_coords_from_grid(iskip,jskip)
    MPATCH.def_meshpatch_properties()
    MPATCH.def_release_patches_from_wet_patches()
  
  DATA.apply_Oil_Properties_on_ijgrid(input_OilProps[group],MPATCH,weighting_coeffs=weighting_coeffs,time=np.int0(np.array(hours)*3600.))
  if(group==0):
    SurfOpnSea_tons_per_km2 = DATA.SurfOpnSea_tons_per_km2  
    DispOpnSea_tons_per_km2 = DATA.DispOpnSea_tons_per_km2
    SurfStrand_tons_per_km2 = DATA.SurfStrand_tons_per_km2
    DispStrand_tons_per_km2 = DATA.DispStrand_tons_per_km2
  else:
    SurfOpnSea_tons_per_km2 += DATA.SurfOpnSea_tons_per_km2  
    DispOpnSea_tons_per_km2 += DATA.DispOpnSea_tons_per_km2
    SurfStrand_tons_per_km2 += DATA.SurfStrand_tons_per_km2
    DispStrand_tons_per_km2 += DATA.DispStrand_tons_per_km2
if(PARALLEL.nranks>1):
  sum_weights_and_releases = PARALLEL.allreducesum_flt64(sum_weights_and_releases)
print('sum_weights_and_releases',sum_weights_and_releases)

DATA.SurfOpnSea_tons_per_km2 = SurfOpnSea_tons_per_km2/sum_weights_and_releases/divide_qtity_released_for_simplification
DATA.DispOpnSea_tons_per_km2 = DispOpnSea_tons_per_km2/sum_weights_and_releases/divide_qtity_released_for_simplification
DATA.SurfStrand_tons_per_km2 = SurfStrand_tons_per_km2/sum_weights_and_releases/divide_qtity_released_for_simplification
DATA.DispStrand_tons_per_km2 = DispStrand_tons_per_km2/sum_weights_and_releases/divide_qtity_released_for_simplification

if(NFf==NF+1):
  DATA.oilspillhazard_to_netcdf("Hazard_reggrid_"+str(iskip)+"x"+str(jskip)+"_"+CASE+"_M"+NF+".nc")
else:
  DATA.oilspillhazard_to_netcdf("Hazard_reggrid_"+str(iskip)+"x"+str(jskip)+"_"+CASE+"_M"+"_"+NF+'-'+NFf+".nc")

if(plot_figures):
     fig,ax=plt.subplots(figsize=(14,5))
     _=ax.plot(np.sum(DATA.SurfOpnSea_tons_per_km2[:,:,:],axis=(1,2)),label='SurfOpnSea_tons/km2',c='b')
     _=ax.plot(np.sum(DATA.DispOpnSea_tons_per_km2[:,:,:],axis=(1,2)),label='DispOpnSea_tons/km2',c='r')
     _=ax.plot(np.sum(DATA.SurfStrand_tons_per_km2[:,:,:],axis=(1,2)),label='SurfStrand_tons/km2',c='b',linestyle='--',alpha=0.8)
     _=ax.plot(np.sum(DATA.DispStrand_tons_per_km2[:,:,:],axis=(1,2)),label='DispStrand_tons/km2',c='r',linestyle='--',alpha=0.8)
     plt.legend()
     if(show_figures):plt.show()
     if(NFf==NF+1):
       fig.savefig('Hazard_reggrid_'+str(iskip)+"x"+str(jskip)+'_'+CASE+"_M"+NF+'.png')
     else:
       fig.savefig('Hazard_reggrid_'+str(iskip)+"x"+str(jskip)+'_'+CASE+"_M"+"_"+NF+'-'+NFf+'.png')
