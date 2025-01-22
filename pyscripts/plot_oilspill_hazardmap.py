#CASE='HARMONIA'
CASE='SHAREMED'
import os,sys
LTRANS_tools = os.getenv('LTRANS_tools')
sys.path.append(LTRANS_tools+'/pymodules')
import netcdftools as nc
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.colors as colors
from matplotlib.colors import ListedColormap
import lagrangiantools as LT
os.system('mkdir '+CASE)
countarg=1
def readstr(countarg,arguments,stringquestion):
 if len(arguments)>countarg : stringworld= str(arguments[countarg]) ; countarg+=1 ;  print(stringquestion+stringworld)
 else: stringworld=str(input(stringquestion))
 return (stringworld,countarg)
#-----------------------------------------------------


PARALLEL=LT.Parallel(False)
if(CASE=='HARMONIA'):   GRID=LT.LTRANSgrid('/share/scratch/athibaut/LTRANS_Zlev/SIM/input/GridforLTRANS-MANTIS-c1.nc',PARALLEL)
elif(CASE=='SHAREMED'): GRID=LT.LTRANSgrid('/share/scratch/athibaut/LTRANS_Zlev/SIM/input/GridforLTRANS-CADEAU_760p-c1.nc',PARALLEL)
else:
    print('CASE',CASE,' not supported')
    sys.exit()

(ncfile,countarg)=readstr(countarg,sys.argv,"netcdf file : ")
DATA=nc.df_from_netcdf(ncfile)

col_dict={0:(1,1,1),
          1:(1,0.93,0.4),
          2:(1,0.80,0),
          3:(1,0.45,0),
          4:(0.8,0,0.0),
          5:(0,0,0),
         }

from matplotlib.colors import BoundaryNorm


levels=[0.001, 0.01, 0.1, 1 ,10, 100,1000]

cm = ListedColormap([col_dict[x] for x in col_dict.keys()])
norm = BoundaryNorm(levels, ncolors=cm.N, clip=True)
DATA.var.haz_map_surf=np.ma.masked_where(DATA.var.haz_map_surf<1e-8,DATA.var.haz_map_surf)
DATA.var.haz_map_disp=np.ma.masked_where(DATA.var.haz_map_disp<1e-8,DATA.var.haz_map_disp)
DATA.var.haz_map_surf_beach=np.ma.masked_where(DATA.var.haz_map_surf_beach<1e-8,DATA.var.haz_map_surf_beach)
DATA.var.haz_map_disp_beach=np.ma.masked_where(DATA.var.haz_map_disp_beach<1e-8,DATA.var.haz_map_disp_beach)

for t in range(DATA.dim.time_dim):
  if(CASE=='HARMONIA'):fig,ax=plt.subplots(nrows=2,ncols=2,figsize=(19,13))
  elif(CASE=='SHAREMED'):fig,ax=plt.subplots(nrows=2,ncols=2,figsize=(22,12))
  else:
      print('CASE',CASE,' not supported')
      sys.exit()

  cbars=()

  ax[0,0].contourf(GRID.lon,GRID.lat,GRID.depth,levels=[-0.1,0,0.1],cmap='Greys',zorder=0)
  mappable=ax[0,0].pcolormesh(DATA.var.lon[1:,:],DATA.var.lat[1:,:],DATA.var.haz_map_surf[t][1:,:]*1000.,cmap=cm,norm=norm,shading='auto',zorder=1)
  cbar=fig.colorbar(mappable, ax=ax[0,0],extend='both') 
  ax[0,0].set_title('OpenSea surface-oil',fontsize=16)                                         
  ax[0,0].set_aspect('equal')       
  cbars+=(cbar,)
                                                                                                                             
  ax[0,1].contourf(GRID.lon,GRID.lat,GRID.depth,levels=[-0.1,0,0.1],cmap='Greys',zorder=0)
  mappable=ax[0,1].pcolormesh(DATA.var.lon[1:,:],DATA.var.lat[1:,:],DATA.var.haz_map_disp[t][1:,:]*1000.,cmap=cm,norm=norm,shading='auto',zorder=1)
  cbar=fig.colorbar(mappable, ax=ax[0,1],extend='both')
  ax[0,1].set_title('OpenSea oil dispersed in the water column',fontsize=16)   
  ax[0,1].set_aspect('equal')                                                                                                                            
  cbars+=(cbar,)

  ax[1,0].contourf(GRID.lon,GRID.lat,GRID.depth,levels=[-0.1,0,0.1],cmap='Greys',zorder=0)
  mappable=ax[1,0].pcolormesh(DATA.var.lon[1:,:],DATA.var.lat[1:,:],DATA.var.haz_map_surf_beach[t][1:,:]*1000.,cmap=cm,norm=norm,shading='auto',zorder=1)
  cbar=fig.colorbar(mappable, ax=ax[1,0],extend='both')  
  ax[1,0].set_title('Stranded surface-oil',fontsize=16)  
  ax[1,0].set_aspect('equal')                                                                                                                            
  cbars+=(cbar,)

  ax[1,1].contourf(GRID.lon,GRID.lat,GRID.depth,levels=[-0.1,0,0.1],cmap='Greys',zorder=0)
  mappable=ax[1,1].pcolormesh(DATA.var.lon[1:,:],DATA.var.lat[1:,:],DATA.var.haz_map_disp_beach[t][1:,:]*1000.,cmap=cm,norm=norm,shading='auto',zorder=1)
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
  plt.subplots_adjust(top=0.92,bottom=0.02,hspace=0.15)

  fig.suptitle('H '+str(int(np.round(DATA.var.time[t]*24.)))+' after release',fontsize=26)
  fig.savefig(CASE+'/Hazard_{:02d}.png'.format(t))
  print(CASE+'/Hazard_{:02d}.png created'.format(t))
