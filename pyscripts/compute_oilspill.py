
hours=[0,1,2,3,6,9,12,15,18,21,24,48,72]

grid_fname='/path/to/LTRANS_Zlev/SIM/input/GridforLTRANS.nc'
input_nc=['/path/to/file/output.nc']
input_OilProps=['/path/to/file/output-OilPropsOut.csv']
   
(iskip,jskip)=(1,1)

#-----------------------------------------------------
import numpy as np
import os,sys
LTRANS_tools = os.getenv('LTRANS_tools')
sys.path.append(LTRANS_tools+'/pymodules')
import lagrangiantools as LT
#-----------------------------------------------------

PARALLEL=LT.Parallel(False)

GRID=LT.LTRANSgrid(grid_fname,PARALLEL)
DATA=LT.LTRANSdata(input_nc,PARALLEL,store_all=True)

MPATCH=LT.MeshPatch(GRID,DATA,PARALLEL)
MPATCH.def_meshpatch_coords_from_grid(iskip,jskip)
MPATCH.def_meshpatch_properties()
MPATCH.def_release_patches_from_wet_patches()

DATA.apply_Oil_Properties_on_ijgrid(input_OilProps,MPATCH,time=np.int0(np.array(hours)*3600.))
DATA.oilspillhazard_to_netcdf("Hazard_reggrid.nc",unit='tons/km2')
