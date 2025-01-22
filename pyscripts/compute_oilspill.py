run_parallel=True
represent_hazard_on_wet_patches=False
represent_hazard_on_reg_grid=True
plot_figures=False
show_figures=False

hours=[1,3,8,11,14,17,20,23,2*24-1,3*24-1,4*24-1,5*24-1,6*24-1,7*24-1,8*24-1,9*24-1,10*24-2]

grid_fname='/share/scratch/athibaut/LTRANS_Zlev/SIM/input/GridforLTRANS-MANTIS-c1.nc'
input_nc=['/share/scratch/athibaut/LTRANS_Zlev/SIM/output_HARMONIA/ArabLight_Fin_point100_release001.nc',
               '/share/scratch/athibaut/LTRANS_Zlev/SIM/output_HARMONIA/ArabLight_Fin_point500_release001.nc',
               '/share/scratch/athibaut/LTRANS_Zlev/SIM/output_HARMONIA/ArabLight_Fin_point200_release001.nc',
               '/share/scratch/athibaut/LTRANS_Zlev/SIM/output_HARMONIA/ArabLight_Fin_point600_release001.nc',
               '/share/scratch/athibaut/LTRANS_Zlev/SIM/output_HARMONIA/ArabLight_Fin_point300_release001.nc',
               '/share/scratch/athibaut/LTRANS_Zlev/SIM/output_HARMONIA/ArabLight_Fin_point700_release001.nc',
               '/share/scratch/athibaut/LTRANS_Zlev/SIM/output_HARMONIA/ArabLight_Fin_point400_release001.nc',
               '/share/scratch/athibaut/LTRANS_Zlev/SIM/output_HARMONIA/ArabLight_Fin_point800_release001.nc',
               '/share/scratch/athibaut/LTRANS_Zlev/SIM/output_HARMONIA/ArabLight_Fin_point100_release002.nc',
               '/share/scratch/athibaut/LTRANS_Zlev/SIM/output_HARMONIA/ArabLight_Fin_point500_release002.nc',
               '/share/scratch/athibaut/LTRANS_Zlev/SIM/output_HARMONIA/ArabLight_Fin_point200_release002.nc',
               '/share/scratch/athibaut/LTRANS_Zlev/SIM/output_HARMONIA/ArabLight_Fin_point600_release002.nc',
               '/share/scratch/athibaut/LTRANS_Zlev/SIM/output_HARMONIA/ArabLight_Fin_point300_release002.nc',
               '/share/scratch/athibaut/LTRANS_Zlev/SIM/output_HARMONIA/ArabLight_Fin_point700_release002.nc',
               '/share/scratch/athibaut/LTRANS_Zlev/SIM/output_HARMONIA/ArabLight_Fin_point400_release002.nc',
               '/share/scratch/athibaut/LTRANS_Zlev/SIM/output_HARMONIA/ArabLight_Fin_point800_release002.nc',
               '/share/scratch/athibaut/LTRANS_Zlev/SIM/output_HARMONIA/ArabLight_Fin_point100_release003.nc',
               '/share/scratch/athibaut/LTRANS_Zlev/SIM/output_HARMONIA/ArabLight_Fin_point200_release003.nc',
               '/share/scratch/athibaut/LTRANS_Zlev/SIM/output_HARMONIA/ArabLight_Fin_point300_release003.nc',
               '/share/scratch/athibaut/LTRANS_Zlev/SIM/output_HARMONIA/ArabLight_Fin_point400_release003.nc']
input_OilProps=['/share/scratch/athibaut/LTRANS_Zlev/SIM/output_HARMONIA/ArabLight_Fin_point100_release001-OilPropsOut.csv',
                '/share/scratch/athibaut/LTRANS_Zlev/SIM/output_HARMONIA/ArabLight_Fin_point500_release001-OilPropsOut.csv',
                '/share/scratch/athibaut/LTRANS_Zlev/SIM/output_HARMONIA/ArabLight_Fin_point200_release001-OilPropsOut.csv',
                '/share/scratch/athibaut/LTRANS_Zlev/SIM/output_HARMONIA/ArabLight_Fin_point600_release001-OilPropsOut.csv',
                '/share/scratch/athibaut/LTRANS_Zlev/SIM/output_HARMONIA/ArabLight_Fin_point300_release001-OilPropsOut.csv',
                '/share/scratch/athibaut/LTRANS_Zlev/SIM/output_HARMONIA/ArabLight_Fin_point700_release001-OilPropsOut.csv',
                '/share/scratch/athibaut/LTRANS_Zlev/SIM/output_HARMONIA/ArabLight_Fin_point400_release001-OilPropsOut.csv',
                '/share/scratch/athibaut/LTRANS_Zlev/SIM/output_HARMONIA/ArabLight_Fin_point800_release001-OilPropsOut.csv', 
                '/share/scratch/athibaut/LTRANS_Zlev/SIM/output_HARMONIA/ArabLight_Fin_point100_release002-OilPropsOut.csv',
                '/share/scratch/athibaut/LTRANS_Zlev/SIM/output_HARMONIA/ArabLight_Fin_point500_release002-OilPropsOut.csv',
                '/share/scratch/athibaut/LTRANS_Zlev/SIM/output_HARMONIA/ArabLight_Fin_point200_release002-OilPropsOut.csv',
                '/share/scratch/athibaut/LTRANS_Zlev/SIM/output_HARMONIA/ArabLight_Fin_point600_release002-OilPropsOut.csv',
                '/share/scratch/athibaut/LTRANS_Zlev/SIM/output_HARMONIA/ArabLight_Fin_point300_release002-OilPropsOut.csv',
                '/share/scratch/athibaut/LTRANS_Zlev/SIM/output_HARMONIA/ArabLight_Fin_point700_release002-OilPropsOut.csv',
                '/share/scratch/athibaut/LTRANS_Zlev/SIM/output_HARMONIA/ArabLight_Fin_point400_release002-OilPropsOut.csv',
                '/share/scratch/athibaut/LTRANS_Zlev/SIM/output_HARMONIA/ArabLight_Fin_point800_release002-OilPropsOut.csv',
                '/share/scratch/athibaut/LTRANS_Zlev/SIM/output_HARMONIA/ArabLight_Fin_point100_release003-OilPropsOut.csv',             
                '/share/scratch/athibaut/LTRANS_Zlev/SIM/output_HARMONIA/ArabLight_Fin_point500_release003-OilPropsOut.csv',
                '/share/scratch/athibaut/LTRANS_Zlev/SIM/output_HARMONIA/ArabLight_Fin_point200_release003-OilPropsOut.csv',
                '/share/scratch/athibaut/LTRANS_Zlev/SIM/output_HARMONIA/ArabLight_Fin_point600_release003-OilPropsOut.csv']
   
iniparloc_files=['/share/scratch/athibaut/LTRANS_Zlev/SIM/input/Iniparloc_Traffic_Routes/Iniparloc_Traffic_routes_100.csv',
                 '/share/scratch/athibaut/LTRANS_Zlev/SIM/input/Iniparloc_Traffic_Routes/Iniparloc_Traffic_routes_500.csv',
                 '/share/scratch/athibaut/LTRANS_Zlev/SIM/input/Iniparloc_Traffic_Routes/Iniparloc_Traffic_routes_200.csv',
                 '/share/scratch/athibaut/LTRANS_Zlev/SIM/input/Iniparloc_Traffic_Routes/Iniparloc_Traffic_routes_600.csv',
                 '/share/scratch/athibaut/LTRANS_Zlev/SIM/input/Iniparloc_Traffic_Routes/Iniparloc_Traffic_routes_300.csv',
                 '/share/scratch/athibaut/LTRANS_Zlev/SIM/input/Iniparloc_Traffic_Routes/Iniparloc_Traffic_routes_700.csv',
                 '/share/scratch/athibaut/LTRANS_Zlev/SIM/input/Iniparloc_Traffic_Routes/Iniparloc_Traffic_routes_400.csv',
                 '/share/scratch/athibaut/LTRANS_Zlev/SIM/input/Iniparloc_Traffic_Routes/Iniparloc_Traffic_routes_800.csv',
                 '/share/scratch/athibaut/LTRANS_Zlev/SIM/input/Iniparloc_Traffic_Routes/Iniparloc_Traffic_routes_100.csv',
                 '/share/scratch/athibaut/LTRANS_Zlev/SIM/input/Iniparloc_Traffic_Routes/Iniparloc_Traffic_routes_500.csv',
                 '/share/scratch/athibaut/LTRANS_Zlev/SIM/input/Iniparloc_Traffic_Routes/Iniparloc_Traffic_routes_200.csv',
                 '/share/scratch/athibaut/LTRANS_Zlev/SIM/input/Iniparloc_Traffic_Routes/Iniparloc_Traffic_routes_600.csv',
                 '/share/scratch/athibaut/LTRANS_Zlev/SIM/input/Iniparloc_Traffic_Routes/Iniparloc_Traffic_routes_300.csv',
                 '/share/scratch/athibaut/LTRANS_Zlev/SIM/input/Iniparloc_Traffic_Routes/Iniparloc_Traffic_routes_700.csv',
                 '/share/scratch/athibaut/LTRANS_Zlev/SIM/input/Iniparloc_Traffic_Routes/Iniparloc_Traffic_routes_400.csv',
                 '/share/scratch/athibaut/LTRANS_Zlev/SIM/input/Iniparloc_Traffic_Routes/Iniparloc_Traffic_routes_800.csv', 
                 '/share/scratch/athibaut/LTRANS_Zlev/SIM/input/Iniparloc_Traffic_Routes/Iniparloc_Traffic_routes_100.csv',
                 '/share/scratch/athibaut/LTRANS_Zlev/SIM/input/Iniparloc_Traffic_Routes/Iniparloc_Traffic_routes_500.csv',
                 '/share/scratch/athibaut/LTRANS_Zlev/SIM/input/Iniparloc_Traffic_Routes/Iniparloc_Traffic_routes_200.csv',
                 '/share/scratch/athibaut/LTRANS_Zlev/SIM/input/Iniparloc_Traffic_Routes/Iniparloc_Traffic_routes_600.csv',
]

#-----------------------------------------------------
import numpy as np
import os,sys
LTRANS_tools = os.getenv('LTRANS_tools')
sys.path.append(LTRANS_tools+'/pymodules')
import lagrangiantools as LT
if(plot_figures):
  if(not show_figures):
    import matplotlib
    matplotlib.use('Agg')
  import matplotlib.pyplot as plt
#-----------------------------------------------------

PARALLEL=LT.Parallel(run_parallel)
iniparloc=LT.csv_to_pd_dataframe(iniparloc_files,PARALLEL)
weighting_coeffs=np.zeros((len(iniparloc.contents)),dtype=float)
for f in range(len(weighting_coeffs)):
    weighting_coeffs[f]=float(iniparloc.contents[f].columns[6])  # FOR HARMONIA SIMULATIONS ON MANTIS DOMAIN
    #weighting_coeffs[f]=float(iniparloc.contents[f].columns[5])  # FOR SHAREMED SIMULATIONS ON CADEAU DOMAIN
print('rank ',PARALLEL.rank,' has weighting_coeffs computed for every file as iniparloc[6]*iniparloc[7] , coeffs are : \n',weighting_coeffs)

GRID=LT.LTRANSgrid(grid_fname,PARALLEL)
DATA=LT.LTRANSdata(input_nc,PARALLEL,store_all=True)

(iskip,jskip)=(5,5)

if(represent_hazard_on_reg_grid):
   MPATCH=LT.MeshPatch(GRID,DATA,PARALLEL)
   MPATCH.def_meshpatch_coords_from_grid(iskip,jskip)
   MPATCH.def_meshpatch_properties()
   MPATCH.def_release_patches_from_wet_patches()

   DATA.apply_Oil_Properties_on_ijgrid(input_OilProps,MPATCH,weighting_coeffs=weighting_coeffs,time=np.int0(np.array(hours)*3600.))
   DATA.oilspillhazard_to_netcdf("Hazard_reggrid.nc")
   if(plot_figures):
      fig,ax=plt.subplots(figsize=(14,5))
      _=ax.plot(np.sum(DATA.SurfOpnSea_tons_per_km2[:,:,:],axis=(1,2)),label='SurfOpnSea_tons/km2',c='b')
      _=ax.plot(np.sum(DATA.DispOpnSea_tons_per_km2[:,:,:],axis=(1,2)),label='DispOpnSea_tons/km2',c='r')
      _=ax.plot(np.sum(DATA.SurfStrand_tons_per_km2[:,:,:],axis=(1,2)),label='SurfStrand_tons/km2',c='b',linestyle='--',alpha=0.8)
      _=ax.plot(np.sum(DATA.DispStrand_tons_per_km2[:,:,:],axis=(1,2)),label='DispStrand_tons/km2',c='r',linestyle='--',alpha=0.8)
      plt.legend()
      if(show_figures):plt.show()
      fig.savefig('Hazard_reggrid.png')


if(represent_hazard_on_wet_patches):
   MPATCH=LT.MeshPatch(GRID,DATA,PARALLEL)
   MPATCH.def_meshpatch_coords_from_grid(iskip,jskip)
   MPATCH.def_meshpatch_properties()
   MPATCH.def_release_patches_from_wet_patches()
   MPATCH.def_arrival_patches_from_wet_patches()

   DATA.apply_Oil_Properties_on_meshpatch(input_OilProps,MPATCH,weighting_coeffs=weighting_coeffs,time=np.int0(np.array(hours)*3600.))
   DATA.oilspillhazard_to_netcdf("Hazard_bypatch.nc")
   if(plot_figures):
     fig,ax=plt.subplots(figsize=(14,5))
     _=ax.plot(np.sum(DATA.SurfOpnSea_tons_per_km2[:,:],axis=1),label='SurfOpnSea_tons/km2',c='b')
     _=ax.plot(np.sum(DATA.DispOpnSea_tons_per_km2[:,:],axis=1),label='DispOpnSea_tons/km2',c='r')
     _=ax.plot(np.sum(DATA.SurfStrand_tons_per_km2[:,:],axis=1),label='SurfStrand_tons/km2',c='b',linestyle='--',alpha=0.8)
     _=ax.plot(np.sum(DATA.DispStrand_tons_per_km2[:,:],axis=1),label='DispStrand_tons/km2',c='r',linestyle='--',alpha=0.8)
     plt.legend()
     if(show_figures):plt.show()
     fig.savefig('Hazard_bypatch.png')

