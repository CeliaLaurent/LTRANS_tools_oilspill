grid_fname='/path/to/LTRANS_Zlev/SIM/input/GridforLTRANS.nc'
postproc_folder='/path/to/postproc/folder/'
ncfile='/path/to/postproc/folder/hazardfile.nc'
import os,sys
LTRANS_tools = os.getenv('LTRANS_tools')
sys.path.append(LTRANS_tools+'/pymodules')
import lagrangiantools as LT
import netCDF4

import matplotlib.pyplot as plt
import numpy as np
import matplotlib.colors as colors
from matplotlib.colors import ListedColormap
from matplotlib.colors import BoundaryNorm

os.system('mkdir '+postproc_folder)
#-----------------------------------------------------

PARALLEL=LT.Parallel(False)
GRID=LT.LTRANSgrid(grid_fname,PARALLEL)
DATA=netCDF4.Dataset(ncfile,'r')

col_dict={0:(1,1,1),
          1:(1,0.93,0.4),
          2:(1,0.80,0),
          3:(1,0.45,0),
          4:(0.8,0,0.0),
          5:(0,0,0),
         }
levels=[0.001, 0.01, 0.1, 1 ,10, 100,1000]

cm = ListedColormap([col_dict[x] for x in col_dict.keys()])
norm = BoundaryNorm(levels, ncolors=cm.N, clip=True)
#ATA.variables['haz_map_surf']=np.ma.masked_where(DATA.variables['haz_map_surf'][:]<1e-8,DATA.variables['haz_map_surf'])
#ATA.variables['haz_map_disp']=np.ma.masked_where(DATA.variables['haz_map_disp'][:]<1e-8,DATA.variables['haz_map_disp'])
#ATA.variables['haz_map_surf_beach']=np.ma.masked_where(DATA.variables['haz_map_surf_beach'][:]<1e-8,DATA.variables['haz_map_surf_beach'])
#ATA.variables['haz_map_disp_beach']=np.ma.masked_where(DATA.variables['haz_map_disp_beach'][:]<1e-8,DATA.variables['haz_map_disp_beach'])

for t in range(DATA.dimensions['time'].size):
  fig,ax=plt.subplots(nrows=2,ncols=2,figsize=(19,9))
  cbars=()

  ax[0,0].contourf(GRID.lon,GRID.lat,GRID.depth,levels=[-0.1,0,0.1],cmap='Greys',zorder=0)
  mappable=ax[0,0].pcolormesh(DATA.variables['lon'][1:,:],DATA.variables['lat'][1:,:],DATA.variables['haz_map_surf'][t][1:,:]*1000.,cmap=cm,norm=norm,shading='auto',zorder=1)
  cbar=fig.colorbar(mappable, ax=ax[0,0],extend='both') 
  ax[0,0].set_title('OpenSea surface-oil',fontsize=16)                                         
  ax[0,0].set_aspect('equal')       
  cbars+=(cbar,)
                                                                                                                             
  ax[0,1].contourf(GRID.lon,GRID.lat,GRID.depth,levels=[-0.1,0,0.1],cmap='Greys',zorder=0)
  mappable=ax[0,1].pcolormesh(DATA.variables['lon'][1:,:],DATA.variables['lat'][1:,:],DATA.variables['haz_map_disp'][t][1:,:]*1000.,cmap=cm,norm=norm,shading='auto',zorder=1)
  cbar=fig.colorbar(mappable, ax=ax[0,1],extend='both')
  ax[0,1].set_title('OpenSea oil dispersed in the water column',fontsize=16)   
  ax[0,1].set_aspect('equal')                                                                                                                            
  cbars+=(cbar,)

  ax[1,0].contourf(GRID.lon,GRID.lat,GRID.depth,levels=[-0.1,0,0.1],cmap='Greys',zorder=0)
  mappable=ax[1,0].pcolormesh(DATA.variables['lon'][1:,:],DATA.variables['lat'][1:,:],DATA.variables['haz_map_surf_beach'][t][1:,:]*1000.,cmap=cm,norm=norm,shading='auto',zorder=1)
  cbar=fig.colorbar(mappable, ax=ax[1,0],extend='both')  
  ax[1,0].set_title('Stranded surface-oil',fontsize=16)  
  ax[1,0].set_aspect('equal')                                                                                                                            
  cbars+=(cbar,)

  ax[1,1].contourf(GRID.lon,GRID.lat,GRID.depth,levels=[-0.1,0,0.1],cmap='Greys',zorder=0)
  mappable=ax[1,1].pcolormesh(DATA.variables['lon'][1:,:],DATA.variables['lat'][1:,:],DATA.variables['haz_map_disp_beach'][t][1:,:]*1000.,cmap=cm,norm=norm,shading='auto',zorder=1)
  cbar=fig.colorbar(mappable, ax=ax[1,1],extend='both') 
  ax[1,1].set_title('Stranded oil dispersed in the water column',fontsize=16) 
  ax[1,1].set_aspect('equal')
  cbars+=(cbar,)

  for y,label,color in zip(np.linspace(80,925,6),['NULL','VERY-LOW','LOW','MEDIUM','HIGH','VERY-HIGH'],['k','k','k','w','w','w']):
    for num in range(len(cbars)):
       cbars[num].ax.set_title('HAZARD')
       cbars[num].ax.set_yticklabels([' ','0.01','0.1','1','10','100',' '])
       cbars[num].ax.text(2500,500,'kg / km2',size=10,rotation=270, ha='center', va='center',color='k') 
       cbars[num].ax.text(380,y,label, size=10,rotation=270, ha='center', va='center',color=color)
  
  fig.tight_layout()
  plt.subplots_adjust(top=0.92,bottom=0.02,hspace=0.08,wspace=0.05)

  fig.suptitle('H '+str(int(np.round(DATA.variables['time'][t]*24.)))+' after release',fontsize=26)
  fig.savefig(postproc_folder+'/Hazard_{:02d}.png'.format(t))
  print(postproc_folder+'/Hazard_{:02d}.png created'.format(t))
