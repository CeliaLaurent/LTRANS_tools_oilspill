#CASE='HARMONIA'
CASE='SHAREMED'
SUBCASE='_singlerelease'
LTRANS_tools = os.getenv('LTRANS_tools')
sys.path.append(LTRANS_tools+'/pymodules')
import netcdftools as nc
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.colors as colors
from matplotlib.colors import LinearSegmentedColormap
import lagrangiantools as LT
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
os.system('mkdir '+ncfile[:-3])
DATA=nc.df_from_netcdf(ncfile)

col_dict={0:(1,1,1),
          #1:(0.84,1.0,0.45),
          1:(1,0.93,0.4),
          2:(1,0.80,0),
          3:(1,0.45,0),
          4:(0.8,0,0.0),
          5:(0,0,0),
         }

from matplotlib.colors import BoundaryNorm

vmin=1e-3
vmax=5e-1
levels=[0.1,0.33,1,3.3,10,33,100]
if(SUBCASE=='_singlerelease'):
  levels=np.hstack((
            np.arange(100,1000,100),
            np.arange(1000,10000,1000),
            np.arange(10000,60000,10000),
            ))
  cm = LinearSegmentedColormap.from_list('whiteyellowredblack', colors=[col_dict[x] for x in col_dict.keys()], N=1000)
  norm = colors.LogNorm(vmin=100,vmax=50000) # BoundaryNorm(levels, ncolors=cm.N, clip=True)
else:
  levels=np.hstack((
            np.arange(0.01,0.1,0.01),
            np.arange(0.1,1,0.1),
            np.arange(1,10,1),
            np.arange(10,100,10),
            np.arange(100,600,100),
            ))
  cm = LinearSegmentedColormap.from_list('whiteyellowredblack', colors=[col_dict[x] for x in col_dict.keys()], N=1000)
  norm = colors.LogNorm(vmin=0.01,vmax=500) # BoundaryNorm(levels, ncolors=cm.N, clip=True)

DATA.var.Hazard_OpnSea_Surf=np.ma.masked_where(DATA.var.Hazard_OpnSea_Surf<1e-15,DATA.var.Hazard_OpnSea_Surf)
DATA.var.Hazard_OpnSea_Disp=np.ma.masked_where(DATA.var.Hazard_OpnSea_Disp<1e-15,DATA.var.Hazard_OpnSea_Disp)
DATA.var.Hazard_Strand_Surf=np.ma.masked_where(DATA.var.Hazard_Strand_Surf<1e-15,DATA.var.Hazard_Strand_Surf)
DATA.var.Hazard_Strand_Disp=np.ma.masked_where(DATA.var.Hazard_Strand_Disp<1e-15,DATA.var.Hazard_Strand_Disp)
for t in range(DATA.dim.timestep_number):
  if(CASE=='HARMONIA'):fig,ax=plt.subplots(nrows=2,ncols=2,figsize=(19,13))
  elif(CASE=='SHAREMED'):fig,ax=plt.subplots(nrows=2,ncols=2,figsize=(22,12))
  else:
      print('CASE',CASE,' not supported')
      sys.exit()

  cbars=()

  ax[0,0].contourf(GRID.lon,GRID.lat,GRID.depth,levels=[-0.1,0,0.1],cmap='Greys',zorder=0)
  mappable=ax[0,0].pcolormesh(DATA.var.lon[1:,:],DATA.var.lat[1:,:],DATA.var.Hazard_OpnSea_Surf[t][1:,:]*1000.,cmap=cm,norm=norm,shading='auto',zorder=1) #,levels=(vmin,vmax,50)
  cbar=fig.colorbar(mappable, ax=ax[0,0],extend='both') 
  ax[0,0].set_title('OpenSea surface-oil',fontsize=16)                                         
  ax[0,0].set_aspect('equal')       
  cbars+=(cbar,)
                                                                                    
  (i,j)=(1,0) 
  ax[i,j].contourf(GRID.lon,GRID.lat,GRID.depth,levels=[-0.1,0,0.1],cmap='Greys',zorder=0)
  mappable=ax[i,j].pcolormesh(DATA.var.lon[1:,:],DATA.var.lat[1:,:],DATA.var.Hazard_OpnSea_Disp[t][1:,:]*1000.,cmap=cm,norm=norm,shading='auto',zorder=1) #,levels=(vmin,vmax,50)
  cbar=fig.colorbar(mappable, ax=ax[i,j],extend='both')
  ax[i,j].set_title('OpenSea oil dispersed in the water column',fontsize=16)   
  ax[i,j].set_aspect('equal')                                                                                                                                      
  cbars+=(cbar,)

  (i,j)=(0,1) 
  ax[i,j].contourf(GRID.lon,GRID.lat,GRID.depth,levels=[-0.1,0,0.1],cmap='Greys',zorder=0)
  mappable=ax[i,j].pcolormesh(DATA.var.lon[1:,:],DATA.var.lat[1:,:],DATA.var.Hazard_Strand_Surf[t][1:,:]*1000.,cmap=cm,norm=norm,shading='auto',zorder=1) #,levels=(vmin,vmax,50)
  cbar=fig.colorbar(mappable, ax=ax[i,j],extend='both')  
  ax[i,j].set_title('Stranded surface-oil',fontsize=16)  
  ax[i,j].set_aspect('equal')                                                                                                                                      
  cbars+=(cbar,)

  ax[1,1].contourf(GRID.lon,GRID.lat,GRID.depth,levels=[-0.1,0,0.1],cmap='Greys',zorder=0)
  mappable=ax[1,1].pcolormesh(DATA.var.lon[1:,:],DATA.var.lat[1:,:],DATA.var.Hazard_Strand_Disp[t][1:,:]*1000.,cmap=cm,norm=norm,shading='auto',zorder=1) #,levels=(vmin,vmax,50)
  cbar=fig.colorbar(mappable, ax=ax[1,1],extend='both') 
  ax[1,1].set_title('Stranded oil dispersed in the water column',fontsize=16) 
  ax[1,1].set_aspect('equal')
  cbars+=(cbar,)

  for num in range(len(cbars)):
       cbars[num].ax.set_title('Density')
       #cbars[num].ax.set_yticklabels([' ','0.5','1','5','10','50',' '])
       cbars[num].ax.text(80000000000,3000,'kg/km2',size=10,rotation=270, ha='center', va='center',color='k') 
  
  fig.tight_layout()
  plt.subplots_adjust(top=0.92,bottom=0.02,hspace=0.15)

  fig.suptitle('H '+str(int(np.round(DATA.var.Time[t]*24.)))+' after release',fontsize=26)
  fig.savefig(ncfile[:-3]+'/Density_{:02d}.png'.format(t))
  print(ncfile[:-3]+'/Density_{:02d}.png created'.format(t))
