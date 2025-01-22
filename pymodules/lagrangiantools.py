import numpy as np
import csv
import math
import argparse
import time
import geopy.distance
import sys
import matplotlib as mpl
import matplotlib.patches as patches
import matplotlib.pylab as pl
try:
  import netCDF4
  print('using netCDF4') 
except:
  import scipy.io.netcdf as NC
import strtools

class Parallel():
  def __init__(self,MPIparallel=True):
    if(MPIparallel):
      from mpi4py import MPI
      self.__MPI=MPI 
      self.comm  = self.__MPI.COMM_WORLD
      self.rank  =  self.comm.Get_rank()
      self.nranks = self.comm.size
      if(self.rank<1):print('Program will run in parallel mode using ',self.nranks,' MPI processes')
    else:
      self.nranks=1
      self.rank=0
      self.comm=None
      print('Program running in sequential mode')

  def rank_range(self,NUM):
    if(self.nranks>1):
      splitindices=np.array([])
      for r in range(1,self.nranks):
        splitindices=np.hstack((splitindices,int(math.ceil(float(NUM)*float(r)/float(self.nranks)))))
      ranksindices=np.hstack(([0],splitindices,[NUM]))
    else:
      ranksindices=[0,NUM]
    return np.int0(np.arange(int(ranksindices[self.rank]),int(ranksindices[self.rank+1]),1))


  def allreducesum_flt32(self,a):  #for np.float32
    rcvBuf = np.copy(a)
    self.comm.Allreduce([a, self.__MPI.FLOAT],
        [rcvBuf,self.__MPI.FLOAT],
        op=self.__MPI.SUM)
    return rcvBuf
  def allreducesum_flt64(self,a):  #  np.float64
    rcvBuf = np.copy(a)
    self.comm.Allreduce([a, self.__MPI.DOUBLE],
        [rcvBuf,self.__MPI.DOUBLE],
        op=self.__MPI.SUM)
    return rcvBuf
  def allreducesum_int(self,a):
    rcvBuf = np.copy(a)
    self.comm.Allreduce([a, self.__MPI.INT],
        [rcvBuf,self.__MPI.INT],
        op=self.__MPI.SUM)
    return rcvBuf



# class Arguments():
#   def __init__(self,*required_args,grid_fname=None,input_nc=None,case=None,out_dir=None,bound_fname=None):
#     parser=Parser()
#     arglist=['grid_fname','input_nc','case','out_dir','bound_fname']
#     for arg in arglist:
#        if(arg in required_args):
#            parser.def_required(arg)
#        else:
#            parser.def_optionals(arg)
#     parser.parse_arguments()
#     if(grid_fname==None):  self.grid_fname = parser.get_value('grid_fname')
#     else:                  self.grid_fname = grid_fname
#     if(input_nc==None):  self.input_nc = parser.get_value('input_nc')
#     else:                self.input_nc = input_nc 
#     if(case==None):  self.case = parser.get_value('case')     
#     else:            self.case = case
#     if(out_dir==None):  self.out_dir = parser.get_value('out_dir') 
#     else:               self.out_dir = out_dir
#     if(bound_fname==None):  self.bound_fname  = parser.get_value('bound_fname') 
#     else:                   self.bound_fname  = bound_fname          
# 
# class Parser():
#   def __init__(self):
#      self.parser = argparse.ArgumentParser()
# 
#   def def_required(self,*args):
#      for arg in args:
#         method = eval("_add_"+arg)
#         method(self,True)
# 
#   def def_optionals(self,*args):
#      for arg in args:
#         method = eval("_add_"+arg)
#         method(self,False)
# 
#   def parse_arguments(self):
#       self.args = self.parser.parse_args()
#   
#   def get_value(self,dest):
#       return eval('self.args.'+str(dest))
# 
#   
# def _add_grid_fname(self,required):
#    self.parser.add_argument('-g', '--grid', required=required, dest='grid_fname', help='grid filename path')
# def _add_input_nc(self,required):
#    self.parser.add_argument('-i', '--input_nc', required=required, dest='input_nc', help='netcdf input file name')
# def _add_case(self,required):
#    self.parser.add_argument('-c', '--case', required=required, dest='case', help='Personalised case name',default='TEST')
# def _add_out_dir(self,required):
#    self.parser.add_argument('-o', '--out_dir', required=required, dest='out_dir', help='Name of output directory',default='TEST')
# def _add_bound_fname(self,required):
#    self.parser.add_argument('-b', '--bound_fname', required=required, dest='bound_fname', help='llbounds.bnl file',default=None)


class LTRANSdata():
  def __init__(self,fname,parallel,store_all):
    if( len(np.shape(fname))==0 ):
      if(fname[-2:]=='nc'):
        namefile=(fname,)
        self.numfiles=1
      else:
        csvfile=open(fname,'r')
        csvreader=csv.reader(csvfile,delimiter=",")
        namefile=()
        for row in csvreader:
          namefile+=(strtools.TRIM(row[0]),)
        self.numfiles=len(namefile)
        csvfile.close()
    else: # case where a list or array of names was given
      namefile=fname
    namefile=np.array(namefile)
    self.my_file_range=parallel.rank_range(len(namefile)) 
    self.my_fnames=namefile[self.my_file_range] 
    self.numfiles_mine=len(self.my_file_range)
    self.nmax = None
    self.tmax = None
    self.dt   = None
    self.lon0   = None
    self.lat0   = None
    self.depth0 = None
    self.__parallel=parallel
    print('rank',self.__parallel.rank,' treats files',self.my_file_range)
    if(self.__parallel.rank<1):
      self.nmax,self.tmax,self.dt,self.lon0,self.lat0,self.depth0=self.__check_common_file_properties(namefile)
      print('\nNetcdf data files contains data for (nmax,tmax,dt)=',(self.nmax,self.tmax,self.dt))
    if(self.__parallel.nranks>1):
      self.nmax= self.__parallel.comm.bcast(self.nmax, root=0)
      self.tmax= self.__parallel.comm.bcast(self.tmax, root=0)
      self.dt  = self.__parallel.comm.bcast(self.dt  , root=0)
      self.lon0    = self.__parallel.comm.bcast(self.lon0    , root=0)
      self.lat0    = self.__parallel.comm.bcast(self.lat0    , root=0)
      self.depth0  = self.__parallel.comm.bcast(self.depth0  , root=0)
    self.my_files_read=0
    if(store_all):
      self.f=()
      self.fnames=()
      for inputfile,outputfilename in zip(self.my_file_range,self.my_fnames): 
        self.f+=(LTRANSnetcdf(outputfilename),)
        self.fnames+=(outputfilename,)
        #print('rank',self.__parallel.rank,'stored contents of file ',outputfilename)
    if(self.__parallel.rank<1):print('')

  def read_next(self):
    outputfilename=self.my_fnames[self.my_files_read]
    print('rank',self.__parallel.rank,' reads',outputfilename)
    self.f=(LTRANSnetcdf(outputfilename),)
    self.my_files_read+=1

  def apply_Oil_Properties_on_meshpatch(self,OilProperties_Fnames,mesh,stranded_oil_was_removed_from_surfoil=True,stranded_oil_was_removed_from_dispoil=False,weighting_coeffs=[None],time=[]):
    self.oil_props_computed_on_patches=True
    OIL=OILTRANS(OilProperties_Fnames,stranded_oil_was_removed_from_surfoil,files_range=self.my_file_range,rank=self.__parallel.rank)
    self.mesh=mesh
    self.OILTRANS=OIL
    for nf in range(len(self.my_file_range)):
      self.f[nf].remove_timesteps_preceeding(OIL.t0)
      self.f[nf].time -= OIL.time0
      if(np.max(np.abs(self.f[nf].time[1:]-OIL.Time[nf,1:]))>1e-5):print('WARNING time shapes=',np.shape(self.f[nf].time[:]),np.shape(OIL.Time[nf,:]),'max time diff=',np.max(np.abs(self.f[nf].time[1:]-OIL.Time[nf,1:])),np.hstack((self.f[nf].time[:2],['...'],self.f[nf].time[-2:])),np.hstack((OIL.Time[nf,:2],['...'],OIL.Time[nf,-2:])))
    if(len(time)>0):
      self.tmax=len(time)
    else:
      self.tmax=len(self.f[0].time[:])
    VolOilSurfOpnSea_byApatch      =np.zeros((self.tmax,self.mesh.nARRpatch),dtype=np.float64)
    VolOilDispOpnSea_byApatch      =np.zeros((self.tmax,self.mesh.nARRpatch),dtype=np.float64)
    VolOilSurfStrand_byApatch=np.zeros((self.tmax,self.mesh.nARRpatch),dtype=np.float64)
    VolOilDispStrand_byApatch=np.zeros((self.tmax,self.mesh.nARRpatch),dtype=np.float64)
    kg_OilSurfOpnSea_byApatch      =np.zeros((self.tmax,self.mesh.nARRpatch),dtype=np.float64)
    kg_OilDispOpnSea_byApatch      =np.zeros((self.tmax,self.mesh.nARRpatch),dtype=np.float64)
    kg_OilSurfStrand_byApatch=np.zeros((self.tmax,self.mesh.nARRpatch),dtype=np.float64)
    kg_OilDispStrand_byApatch=np.zeros((self.tmax,self.mesh.nARRpatch),dtype=np.float64)
    for nf in range(len(self.my_file_range)):
      print('rank',self.__parallel.rank,' applies coeff',weighting_coeffs[min(nf,len(weighting_coeffs)-1)],' and Oil Props ',OIL.fnames[nf],' to particles in',self.fnames[nf])
      if(len(time)>0):
        self.time=time
      else:
        self.time=self.f[nf].time[:]
      VolOilSurf       = np.zeros((self.tmax),dtype=float)
      VolOilSurfOpnSea = np.zeros((self.tmax),dtype=float)
      VolOilSurfStrand = np.zeros((self.tmax),dtype=float)
      VolOilDisp       = np.zeros((self.tmax),dtype=float)
      VolOilDispOpnSea = np.zeros((self.tmax),dtype=float)
      VolOilDispStrand = np.zeros((self.tmax),dtype=float)
      VolOilEvap       = np.zeros((self.tmax),dtype=float)
      VolOilTOT        = np.zeros((self.tmax),dtype=float)
      RhoOil           = np.zeros((self.tmax),dtype=float)
      tmaxoil=len(OIL.Time[nf,:])
      t_oldoil=0
      tdata=()
      for t in range(0,self.tmax):
        found=False
        for toil in range(t_oldoil,tmaxoil-1):
          ta=OIL.Time[nf,toil]
          tb=OIL.Time[nf,toil+1]
          if(tb<self.time[t]):continue
          if(ta-1.<=self.time[t] and tb-1.>self.time[t]):
              Factor1=float(tb-self.time[t])/float(tb-ta)
              Factor2=float(self.time[t]-ta)/float(tb-ta)
              if(stranded_oil_was_removed_from_surfoil and stranded_oil_was_removed_from_dispoil):
                VolOilSurfOpnSea[t:]=(OIL.VolSurfOpnSea[nf,toil  ]*Factor1+OIL.VolSurfOpnSea[nf,toil+1]*Factor2)
                VolOilSurfStrand[t:]=(OIL.VolSurfStrand[nf,toil  ]*Factor1+OIL.VolSurfStrand[nf,toil+1]*Factor2)
              elif(stranded_oil_was_removed_from_surfoil and not stranded_oil_was_removed_from_dispoil):
                VolOilSurfOpnSea[t:]=(OIL.VolSurf[nf,toil  ]*Factor1+OIL.VolSurf[nf,toil+1]*Factor2)
                VolOilSurfStrand[t:]=(OIL.VolStrand[nf,toil  ]*Factor1+OIL.VolStrand[nf,toil+1]*Factor2)
              else:
                 VolOilSurf[t:]=(OIL.VolSurf[nf,toil  ]*Factor1+OIL.VolSurf[nf,toil+1]*Factor2)
              if(stranded_oil_was_removed_from_dispoil):
                 VolOilDispStrand = (OIL.VolDispStrand[nf,toil  ] *Factor1+OIL.VolDispStrand[nf,toil+1]*Factor2)
                 VolOilDispOpnSea = (OIL.VolDispOpnSea[nf,toil  ] *Factor1+OIL.VolDispOpnSea[nf,toil+1]*Factor2)
              else:
                 VolOilDisp[t:]=(OIL.VolDisp[nf,toil  ]*Factor1+OIL.VolDisp[nf,toil+1]*Factor2)
              VolOilEvap[t:]=(OIL.VolEvap[nf,toil  ]*Factor1+OIL.VolEvap[nf,toil+1]*Factor2)
              VolOilTOT[t:]=OIL.VolOilTot[nf,toil  ]*Factor1+OIL.VolOilTot[nf,toil+1]*Factor2
              RhoOil[t:]=OIL.Rho_Oil[nf,toil  ]*Factor1+OIL.Rho_Oil[nf,toil+1]*Factor2
              found=True
          if(found):
            tdata+=(toil,)
            break
        t_oldoil=toil
      tdata=np.array(tdata)
      #-------------------
      numaliveparts=np.zeros(self.tmax,dtype=float)
      numbeachparts=np.zeros(self.tmax,dtype=float)
      for t in range(0,self.tmax):
        numaliveparts[t]=float(len(np.ma.masked_where(self.f[nf].color[tdata[t],:]!=1 ,np.arange(0,self.nmax)).compressed()))
        numbeachparts[t]=float(len(np.ma.masked_where(np.abs(self.f[nf].color[tdata[t],:])!=2 ,np.arange(0,self.nmax)).compressed()))
      #-------------------
      if(not stranded_oil_was_removed_from_surfoil):
        VolOilSurfStrand = np.zeros(np.shape(VolOilSurf),dtype=float)
        VolOilSurfOpnSea = np.copy(VolOilSurf)
        for t in range(0,self.tmax):
          if(t==0):
             beaching_parts=numbeachparts[t]
          else:
             beaching_parts=numbeachparts[t]-numbeachparts[t-1]
          if(beaching_parts>0): 
            VolOilSurfStrand[t:] += VolOilSurfOpnSea[t ]*float(beaching_parts)/float(beaching_parts+numaliveparts[t])
            VolOilSurfOpnSea[t:] -= VolOilSurfOpnSea[t:]*float(beaching_parts)/float(beaching_parts+numaliveparts[t])
            VolOilSurf[t:] -= VolOilSurf[t:]*float(beaching_parts)/float(beaching_parts+numaliveparts[t])
      if(not stranded_oil_was_removed_from_dispoil):
        VolOilDispStrand = np.zeros(np.shape(VolOilDisp),dtype=float)
        VolOilDispOpnSea = np.copy(VolOilDisp)
        for t in range(0,self.tmax):
          if(t==0):
             beaching_parts=numbeachparts[t]
          else:
             beaching_parts=numbeachparts[t]-numbeachparts[t-1]
          if(beaching_parts>0): 
            VolOilDispStrand[t:] += VolOilDispOpnSea[t ]*float(beaching_parts)/float(beaching_parts+numaliveparts[t])
            VolOilDispOpnSea[t:] -= VolOilDispOpnSea[t ]*float(beaching_parts)/float(beaching_parts+numaliveparts[t])
           #VolOilDispOpnSea[t:] -= VolOilDispOpnSea[t:]*float(beaching_parts)/float(beaching_parts+numaliveparts[t])
      if(weighting_coeffs[0]!=None):
        VolOilSurf       *= weighting_coeffs[nf]
        VolOilSurfStrand *= weighting_coeffs[nf]
        VolOilSurfOpnSea *= weighting_coeffs[nf]
        VolOilDisp       *= weighting_coeffs[nf]
        VolOilDispStrand *= weighting_coeffs[nf]
        VolOilDispOpnSea *= weighting_coeffs[nf]
        VolOilEvap       *= weighting_coeffs[nf]
        VolOilTOT        *= weighting_coeffs[nf]
      #-------------------
      for n in range(0,self.nmax):
        tlist=np.ma.masked_where(np.abs(self.f[nf].color[tdata,n])!=2 ,np.arange(0,self.tmax)).compressed()
        if(len(tlist)>0):
          Apatches=self.mesh.coords_to_Apatch(self.f[nf].lon[tdata[tlist],n],self.f[nf].lat[tdata[tlist],n])
          for count,(t,patch) in enumerate(zip(tlist,Apatches)):
           if(count>0 and patch<0): 
               Apatches[count]=Apatches[count-1]
               patch=Apatches[count]
           if(patch>=0):
            if(numbeachparts[t]>0):
              VolOilSurfStrand_byApatch[t,patch]=VolOilSurfStrand_byApatch[t,patch]+1./numbeachparts[t]*VolOilSurfStrand[t]
              VolOilDispStrand_byApatch[t,patch]=VolOilDispStrand_byApatch[t,patch]+1./numbeachparts[t]*VolOilDispStrand[t]
              kg_OilSurfStrand_byApatch[t,patch]=kg_OilSurfStrand_byApatch[t,patch]+1./numbeachparts[t]*(VolOilSurfStrand[t]*RhoOil[t])   # [kg] = [m3] * [kg/m3]
              kg_OilDispStrand_byApatch[t,patch]=kg_OilDispStrand_byApatch[t,patch]+1./numbeachparts[t]*(VolOilDispStrand[t]*RhoOil[t])   # [kg] = [m3] * [kg/m3]
      #-------------------
      for n in range(0,self.nmax):
        tlist=np.ma.masked_where(np.abs(self.f[nf].color[tdata,n])!=1 ,np.arange(0,self.tmax)).compressed()
        if(len(tlist)>0):
          Apatches=self.mesh.coords_to_Apatch(self.f[nf].lon[tdata[tlist],n],self.f[nf].lat[tdata[tlist],n])
          for count,(t,patch) in enumerate(zip(tlist,Apatches)):
           if(patch>=0):
              if(numaliveparts[t]>0):
                  VolOilSurfOpnSea_byApatch[t,patch]=VolOilSurfOpnSea_byApatch[t,patch]+1./numaliveparts[t]*VolOilSurfOpnSea[t] 
                  VolOilDispOpnSea_byApatch[t,patch]=VolOilDispOpnSea_byApatch[t,patch]+1./numaliveparts[t]*VolOilDispOpnSea[t] 
                  kg_OilSurfOpnSea_byApatch[t,patch]=kg_OilSurfOpnSea_byApatch[t,patch]+1./numaliveparts[t]*(VolOilSurfOpnSea[t]*RhoOil[t])   # [kg] = [m3] * [kg/m3]
                  kg_OilDispOpnSea_byApatch[t,patch]=kg_OilDispOpnSea_byApatch[t,patch]+1./numaliveparts[t]*(VolOilDispOpnSea[t]*RhoOil[t])   # [kg] = [m3] * [kg/m3]
      #-------------------
    if(self.__parallel.nranks>1): 
      if(self.__parallel.rank<1):print('REDUCE MATRIX')
      VolOilSurfOpnSea_byApatch=self.__parallel.allreducesum_flt64(VolOilSurfOpnSea_byApatch)
      VolOilDispOpnSea_byApatch=self.__parallel.allreducesum_flt64(VolOilDispOpnSea_byApatch)
      VolOilSurfStrand_byApatch=self.__parallel.allreducesum_flt64(VolOilSurfStrand_byApatch)
      VolOilDispStrand_byApatch=self.__parallel.allreducesum_flt64(VolOilDispStrand_byApatch)
      kg_OilSurfOpnSea_byApatch  =self.__parallel.allreducesum_flt64(kg_OilSurfOpnSea_byApatch)
      kg_OilDispOpnSea_byApatch  =self.__parallel.allreducesum_flt64(kg_OilDispOpnSea_byApatch)
      kg_OilSurfStrand_byApatch=self.__parallel.allreducesum_flt64(kg_OilSurfStrand_byApatch)
      kg_OilDispStrand_byApatch=self.__parallel.allreducesum_flt64(kg_OilDispStrand_byApatch)
      if(self.__parallel.rank<1):print('REDUCE MATRIX DONE')
    
    self.mesh.compute_Area_by_Apatch()
    self.SurfOpnSea_mmThick = VolOilSurfOpnSea_byApatch / self.mesh.Apatch_Area_m2  *1000.  # m3/m2* 1000=mm
    self.DispOpnSea_mmThick = VolOilDispOpnSea_byApatch / self.mesh.Apatch_Area_m2  *1000.  # m3/m2 *1000=mm
    self.SurfStrand_mmThick = VolOilSurfStrand_byApatch / self.mesh.Apatch_Area_m2  *1000.  # m3/m2 *1000=mm
    self.DispStrand_mmThick = VolOilDispStrand_byApatch / self.mesh.Apatch_Area_m2  *1000.  # m3/m2 *1000=mm
    self.SurfOpnSea_tons_per_km2 = kg_OilSurfOpnSea_byApatch * 1e-3 /(self.mesh.Apatch_Area_m2 * 1e-6 ) # [tons=kg/1000] / [km2=m2/1000/1000] = tons/km2
    self.DispOpnSea_tons_per_km2 = kg_OilDispOpnSea_byApatch * 1e-3 /(self.mesh.Apatch_Area_m2 * 1e-6 ) # [tons=kg/1000] / [km2=m2/1000/1000] = tons/km2
    self.SurfStrand_tons_per_km2 = kg_OilSurfStrand_byApatch * 1e-3 /(self.mesh.Apatch_Area_m2 * 1e-6 ) # [tons=kg/1000] / [km2=m2/1000/1000] = tons/km2
    self.DispStrand_tons_per_km2 = kg_OilDispStrand_byApatch * 1e-3 /(self.mesh.Apatch_Area_m2 * 1e-6 ) # [tons=kg/1000] / [km2=m2/1000/1000] = tons/km2
    if(self.__parallel.rank<1):print('OIL Properties applied on selected patches\n')

  def apply_Oil_Properties_on_ijgrid(self,OilProperties_Fnames,mesh,stranded_oil_was_removed_from_surfoil=True,stranded_oil_was_removed_from_dispoil=False,weighting_coeffs=[None],time=[]):
    self.oil_props_computed_on_patches=False
    OIL=OILTRANS(OilProperties_Fnames,stranded_oil_was_removed_from_surfoil,files_range=self.my_file_range,rank=self.__parallel.rank)
    self.mesh=mesh
    self.OILTRANS=OIL
    for nf in range(len(self.my_file_range)):
      self.f[nf].remove_timesteps_preceeding(OIL.t0)
      self.f[nf].time -= OIL.time0
      if(np.max(np.abs(self.f[nf].time[1:]-OIL.Time[nf,1:]))>1e-5):print('WARNING time shapes=',np.shape(self.f[nf].time[:]),np.shape(OIL.Time[nf,:]),'max time diff=',np.max(np.abs(self.f[nf].time[1:]-OIL.Time[nf,1:])),np.hstack((self.f[nf].time[:2],['...'],self.f[nf].time[-2:])),np.hstack((OIL.Time[nf,:2],['...'],OIL.Time[nf,-2:])))
    if(len(time)>0):
      self.tmax=len(time)
    else:
      self.tmax=len(self.f[0].time[:])
    VolOilSurfOpnSea_byJIcell      =np.zeros((self.tmax,len(self.mesh.lon),len(self.mesh.lon[0])),dtype=np.float64)
    VolOilDispOpnSea_byJIcell      =np.zeros((self.tmax,len(self.mesh.lon),len(self.mesh.lon[0])),dtype=np.float64)
    VolOilSurfStrand_byJIcell=np.zeros((self.tmax,len(self.mesh.lon),len(self.mesh.lon[0])),dtype=np.float64)
    VolOilDispStrand_byJIcell=np.zeros((self.tmax,len(self.mesh.lon),len(self.mesh.lon[0])),dtype=np.float64)
    kg_OilSurfOpnSea_byJIcell      =np.zeros((self.tmax,len(self.mesh.lon),len(self.mesh.lon[0])),dtype=np.float64)
    kg_OilDispOpnSea_byJIcell      =np.zeros((self.tmax,len(self.mesh.lon),len(self.mesh.lon[0])),dtype=np.float64)
    kg_OilSurfStrand_byJIcell=np.zeros((self.tmax,len(self.mesh.lon),len(self.mesh.lon[0])),dtype=np.float64)
    kg_OilDispStrand_byJIcell=np.zeros((self.tmax,len(self.mesh.lon),len(self.mesh.lon[0])),dtype=np.float64)
    print('rank',self.__parallel.rank,' applies coeffs',np.shape(weighting_coeffs[:]),' and Oil Props ',np.shape(OIL.fnames[:]),' to particles in',np.shape(self.fnames[:]))
    for nf in range(len(self.my_file_range)):
      print('rank',self.__parallel.rank,' applies coeff',weighting_coeffs[min(nf,len(weighting_coeffs)-1)],' and Oil Props ',OIL.fnames[nf],' to particles in',self.fnames[nf])
      print('A rank ',self.__parallel.rank,' file :',self.fnames[nf],' minmaxsum VolOilSurfOpnSea_byJIcell=',[np.min(VolOilSurfOpnSea_byJIcell),np.max(VolOilSurfOpnSea_byJIcell),np.sum(VolOilSurfOpnSea_byJIcell)])
      if(len(time)>0):
        self.time=time
      else:
        self.time=self.f[nf].time[:]
      VolOilSurf       = np.zeros((self.tmax),dtype=float)
      VolOilSurfOpnSea = np.zeros((self.tmax),dtype=float)
      VolOilSurfStrand = np.zeros((self.tmax),dtype=float)
      VolOilDisp       = np.zeros((self.tmax),dtype=float)
      VolOilDispOpnSea = np.zeros((self.tmax),dtype=float)
      VolOilDispStrand = np.zeros((self.tmax),dtype=float)
      VolOilEvap       = np.zeros((self.tmax),dtype=float)
      VolOilTOT        = np.zeros((self.tmax),dtype=float)
      RhoOil           = np.zeros((self.tmax),dtype=float)
      tmaxoil=len(OIL.Time[nf,:])
      t_oldoil=0
      tdata=()
      for t in range(0,self.tmax):
        found=False
        for toil in range(t_oldoil,tmaxoil-1):
          ta=OIL.Time[nf,toil]
          tb=OIL.Time[nf,toil+1]
          if(tb<self.time[t]):continue
          if(ta-1.<=self.time[t] and tb-1.>self.time[t]):
              Factor1=float(tb-self.time[t])/float(tb-ta)
              Factor2=float(self.time[t]-ta)/float(tb-ta)
              if(stranded_oil_was_removed_from_surfoil and stranded_oil_was_removed_from_dispoil):
                VolOilSurfOpnSea[t:]=(OIL.VolSurfOpnSea[nf,toil  ]*Factor1+OIL.VolSurfOpnSea[nf,toil+1]*Factor2)
                VolOilSurfStrand[t:]=(OIL.VolSurfStrand[nf,toil  ]*Factor1+OIL.VolSurfStrand[nf,toil+1]*Factor2)
              elif(stranded_oil_was_removed_from_surfoil and not stranded_oil_was_removed_from_dispoil):
                VolOilSurfOpnSea[t:]=(OIL.VolSurf[nf,toil  ]*Factor1+OIL.VolSurf[nf,toil+1]*Factor2)
                VolOilSurfStrand[t:]=(OIL.VolStrand[nf,toil  ]*Factor1+OIL.VolStrand[nf,toil+1]*Factor2)
              else:
                 VolOilSurf[t:]=(OIL.VolSurf[nf,toil  ]*Factor1+OIL.VolSurf[nf,toil+1]*Factor2)
              if(stranded_oil_was_removed_from_dispoil):
                 VolOilDispStrand = (OIL.VolDispStrand[nf,toil  ] *Factor1+OIL.VolDispStrand[nf,toil+1]*Factor2)
                 VolOilDispOpnSea = (OIL.VolDispOpnSea[nf,toil  ] *Factor1+OIL.VolDispOpnSea[nf,toil+1]*Factor2)
              else:
                 VolOilDisp[t:]=(OIL.VolDisp[nf,toil  ]*Factor1+OIL.VolDisp[nf,toil+1]*Factor2)
              VolOilEvap[t:]=(OIL.VolEvap[nf,toil  ]*Factor1+OIL.VolEvap[nf,toil+1]*Factor2)
              VolOilTOT[t:]=OIL.VolOilTot[nf,toil  ]*Factor1+OIL.VolOilTot[nf,toil+1]*Factor2
              RhoOil[t:]=OIL.Rho_Oil[nf,toil  ]*Factor1+OIL.Rho_Oil[nf,toil+1]*Factor2
              found=True
          if(found):
            tdata+=(toil,)
            break
        t_oldoil=toil
      tdata=np.array(tdata)
      #-------------------
      numaliveparts=np.zeros(self.tmax,dtype=float)
      numbeachparts=np.zeros(self.tmax,dtype=float)
      for t in range(0,self.tmax):
        numaliveparts[t]=float(len(np.ma.masked_where(self.f[nf].color[tdata[t],:]!=1 ,np.arange(0,self.nmax)).compressed()))
        numbeachparts[t]=float(len(np.ma.masked_where(np.abs(self.f[nf].color[tdata[t],:])!=2 ,np.arange(0,self.nmax)).compressed()))
      #-------------------
      if(not stranded_oil_was_removed_from_surfoil):
        VolOilSurfStrand = np.zeros(np.shape(VolOilSurf),dtype=float)
        VolOilSurfOpnSea = np.copy(VolOilSurf)
        for t in range(0,self.tmax):
          if(t==0):
             beaching_parts=numbeachparts[t]
          else:
             beaching_parts=numbeachparts[t]-numbeachparts[t-1]
          if(beaching_parts>0): 
            VolOilSurfStrand[t:] += VolOilSurfOpnSea[t ]*float(beaching_parts)/float(beaching_parts+numaliveparts[t])
            VolOilSurfOpnSea[t:] -= VolOilSurfOpnSea[t:]*float(beaching_parts)/float(beaching_parts+numaliveparts[t])
            VolOilSurf[t:] -= VolOilSurf[t:]*float(beaching_parts)/float(beaching_parts+numaliveparts[t])
      if(not stranded_oil_was_removed_from_dispoil):
        VolOilDispStrand = np.zeros(np.shape(VolOilDisp),dtype=float)
        VolOilDispOpnSea = np.copy(VolOilDisp)
        for t in range(0,self.tmax):
          if(t==0):
             beaching_parts=numbeachparts[t]
          else:
             beaching_parts=numbeachparts[t]-numbeachparts[t-1]
          if(beaching_parts>0): 
            VolOilDispStrand[t:] += VolOilDispOpnSea[t ]*float(beaching_parts)/float(beaching_parts+numaliveparts[t])
            VolOilDispOpnSea[t:] -= VolOilDispOpnSea[t ]*float(beaching_parts)/float(beaching_parts+numaliveparts[t])
           #VolOilDispOpnSea[t:] -= VolOilDispOpnSea[t:]*float(beaching_parts)/float(beaching_parts+numaliveparts[t])
      if(weighting_coeffs[0]!=None):
        VolOilSurf       *= weighting_coeffs[nf]
        VolOilSurfStrand *= weighting_coeffs[nf]
        VolOilSurfOpnSea *= weighting_coeffs[nf]
        VolOilDisp       *= weighting_coeffs[nf]
        VolOilDispStrand *= weighting_coeffs[nf]
        VolOilDispOpnSea *= weighting_coeffs[nf]
        VolOilEvap       *= weighting_coeffs[nf]
        VolOilTOT        *= weighting_coeffs[nf]
      #-------------------
      for n in range(0,self.nmax):
        tlist=np.ma.masked_where(np.abs(self.f[nf].color[tdata,n])!=2 ,np.arange(0,self.tmax)).compressed()
        if(len(tlist)>0):
          ilist=np.array(np.int0(np.round((self.f[nf].lon[tdata[tlist],n]-self.mesh.lon1)/self.mesh.dlon)))
          jlist=np.array(np.int0(np.round((self.f[nf].lat[tdata[tlist],n]-self.mesh.lat1)/self.mesh.dlat)))
          for t,j,i in zip(tlist,jlist,ilist):
            if(numbeachparts[t]>0):
              VolOilSurfStrand_byJIcell[t,j,i]=VolOilSurfStrand_byJIcell[t,j,i]+1./numbeachparts[t]* VolOilSurfStrand[t]
              VolOilDispStrand_byJIcell[t,j,i]=VolOilDispStrand_byJIcell[t,j,i]+1./numbeachparts[t]* VolOilDispStrand[t]
              kg_OilSurfStrand_byJIcell[t,j,i]=kg_OilSurfStrand_byJIcell[t,j,i]+1./numbeachparts[t]*(VolOilSurfStrand[t]*RhoOil[t])   # [kg] = [m3] * [kg/m3]
              kg_OilDispStrand_byJIcell[t,j,i]=kg_OilDispStrand_byJIcell[t,j,i]+1./numbeachparts[t]*(VolOilDispStrand[t]*RhoOil[t])   # [kg] = [m3] * [kg/m3]
      #-------------------
      for n in range(0,self.nmax):
        tlist=np.ma.masked_where(np.abs(self.f[nf].color[tdata,n])!=1 ,np.arange(0,self.tmax)).compressed()
        if(len(tlist)>0):
          ilist=np.int0(np.round((self.f[nf].lon[tdata[tlist],n]-self.mesh.lon1)/self.mesh.dlon))
          jlist=np.int0(np.round((self.f[nf].lat[tdata[tlist],n]-self.mesh.lat1)/self.mesh.dlat))
          for t,j,i in zip(tlist,jlist,ilist):
              if(numaliveparts[t]>0):
                  VolOilSurfOpnSea_byJIcell[t,j,i]=VolOilSurfOpnSea_byJIcell[t,j,i]+1./numaliveparts[t]* VolOilSurfOpnSea[t]         
                  VolOilDispOpnSea_byJIcell[t,j,i]=VolOilDispOpnSea_byJIcell[t,j,i]+1./numaliveparts[t]* VolOilDispOpnSea[t]         
                  kg_OilSurfOpnSea_byJIcell[t,j,i]=kg_OilSurfOpnSea_byJIcell[t,j,i]+1./numaliveparts[t]*(VolOilSurfOpnSea[t]*RhoOil[t])   # [kg] = [m3] * [kg/m3]
                  kg_OilDispOpnSea_byJIcell[t,j,i]=kg_OilDispOpnSea_byJIcell[t,j,i]+1./numaliveparts[t]*(VolOilDispOpnSea[t]*RhoOil[t])   # [kg] = [m3] * [kg/m3]
      #-------------------
    if(self.__parallel.nranks>1): 
      if(self.__parallel.rank<1):print('REDUCE MATRIX')
      print('B rank ',self.__parallel.rank,' file :',self.fnames[nf],' minmaxsum VolOilSurfOpnSea_byJIcell=',[np.min(VolOilSurfOpnSea_byJIcell),np.max(VolOilSurfOpnSea_byJIcell),np.sum(VolOilSurfOpnSea_byJIcell)])
      VolOilSurfOpnSea_byJIcell=self.__parallel.allreducesum_flt64(VolOilSurfOpnSea_byJIcell)
      VolOilDispOpnSea_byJIcell=self.__parallel.allreducesum_flt64(VolOilDispOpnSea_byJIcell)
      VolOilSurfStrand_byJIcell=self.__parallel.allreducesum_flt64(VolOilSurfStrand_byJIcell)
      VolOilDispStrand_byJIcell=self.__parallel.allreducesum_flt64(VolOilDispStrand_byJIcell)
      kg_OilSurfOpnSea_byJIcell  =self.__parallel.allreducesum_flt64(kg_OilSurfOpnSea_byJIcell)
      kg_OilDispOpnSea_byJIcell  =self.__parallel.allreducesum_flt64(kg_OilDispOpnSea_byJIcell)
      kg_OilSurfStrand_byJIcell=self.__parallel.allreducesum_flt64(kg_OilSurfStrand_byJIcell)
      kg_OilDispStrand_byJIcell=self.__parallel.allreducesum_flt64(kg_OilDispStrand_byJIcell)
      if(self.__parallel.rank<1):print('REDUCE MATRIX DONE')
    if(self.__parallel.rank<1):print('C after allreducesum VolOilSurfOpnSea_byJIcell minmaxsum =',[np.min(VolOilSurfOpnSea_byJIcell),np.max(VolOilSurfOpnSea_byJIcell),np.sum(VolOilSurfOpnSea_byJIcell)])
    self.SurfOpnSea_mmThick = VolOilSurfOpnSea_byJIcell / self.mesh.Area_m2 *1000.  # m3/m2* 1000=mm
    self.DispOpnSea_mmThick = VolOilDispOpnSea_byJIcell / self.mesh.Area_m2 *1000.  # m3/m2 *1000=mm
    self.SurfStrand_mmThick = VolOilSurfStrand_byJIcell / self.mesh.Area_m2 *1000.  # m3/m2 *1000=mm
    self.DispStrand_mmThick = VolOilDispStrand_byJIcell / self.mesh.Area_m2 *1000.  # m3/m2 *1000=mm
    self.SurfOpnSea_tons_per_km2 = kg_OilSurfOpnSea_byJIcell * 1e-3 /(self.mesh.Area_m2 * 1e-6 ) # [tons=kg/1000] / [km2=m2/1000/1000] = tons/km2
    self.DispOpnSea_tons_per_km2 = kg_OilDispOpnSea_byJIcell * 1e-3 /(self.mesh.Area_m2 * 1e-6 ) # [tons=kg/1000] / [km2=m2/1000/1000] = tons/km2
    self.SurfStrand_tons_per_km2 = kg_OilSurfStrand_byJIcell * 1e-3 /(self.mesh.Area_m2 * 1e-6 ) # [tons=kg/1000] / [km2=m2/1000/1000] = tons/km2
    self.DispStrand_tons_per_km2 = kg_OilDispStrand_byJIcell * 1e-3 /(self.mesh.Area_m2 * 1e-6 ) # [tons=kg/1000] / [km2=m2/1000/1000] = tons/km2
    if(self.__parallel.rank<1):print('OIL Properties applied on regular patch-grid\n')


  def __createDimension(self,ncOUT,name,value):
        ncOUT.createDimension(name,value)
        print ('create dimension : '+name+'=',value)
  def __setattr(self,ncOUT,name,value):
        ncOUT.setattr(name,value)
        print ('set attribute : '+name+'=',value)

  def oilspillhazard_to_txt(self,filename,threshold=1e-4):
    try: basename=filename.index('.')
    except: basename=filename

   #f=open(basename+'_LonLatCoordinates_grid.txt',"w")
   #f.write("Lon Lat\n")
   #for j in range(0,len(self.mesh.lon)):
   #  for i in range(0,len(self.mesh.lon[0,:])):
   #    f.write("{0:8F} {1:8F}".format(self.mesh.lon[j,i],self.mesh.lat[j,i] )+"\n")
   #f.close()
    OILTOT= np.sum( self.SurfOpnSea_tons_per_km2[:,:,:]
                   +self.DispOpnSea_tons_per_km2[:,:,:]
                   +self.SurfStrand_tons_per_km2[:,:,:]
                   +self.DispStrand_tons_per_km2[:,:,:],axis=(1,2))
    for t in range(self.tmax):
      tstring=zerosfillednumberstring(int(self.time[t]/3600.),3)
      Surfpercent=self.SurfOpnSea_tons_per_km2[t] / OILTOT[t]*100.0
      Disppercent=self.DispOpnSea_tons_per_km2[t] / OILTOT[t]*100.0
      Beachpercent=(self.SurfStrand_tons_per_km2[t]+self.DispStrand_tons_per_km2[t]) / OILTOT[t]*100.0
      for percentage,oiltype in ( (Surfpercent, "Surface"), (Disppercent, "Dispersed"), (Beachpercent,"Beached") ):
        fsparse=open(basename+'_'+tstring+'_'+oiltype+'.txt',"w")
        fsparse.write("Lon Lat "+oiltype+"_Oil(tons_per_km2)\n")
        for j in range(0,len(GRIDlon)):
          for i in range(0,len(GRIDlon[0,:])):
            if(percentage[t,j,i]>threshold):
               fsparse.write("{0:8F} {1:8F} {2:8F}".format( self.mesh.lon[j,i], self.mesh.lat[j,i], percentage[t,j,i] )+"\n")
        fsparse.close()
      
      
  def oilspillhazard_to_netcdf(self,Binary_filename):
    if(self.__parallel.rank==0):
      if(self.oil_props_computed_on_patches):
        import scipy.io.netcdf as NC
        ncOUT= NC.netcdf_file(Binary_filename,"w")
        print('create netcdf file',Binary_filename)
        self.__createDimension(ncOUT,"patch_number",  self.mesh.nARRpatch)
        self.__createDimension(ncOUT,"timestep_number",  self.tmax )

        self.mesh.compute_patch_latlon()

        ncvar=ncOUT.createVariable("lon",   'f', ('patch_number',))
        setattr(ncvar,'units'              , 'degrees_east')
        setattr(ncvar,'long_name'          , 'longitude of every patch')
        setattr(ncvar,'standard_name'      , 'lon')
        setattr(ncvar,'axis'               , 'X'        )
        setattr(ncvar,'valid_min'          , np.min(self.mesh.arrival_patch_latlon[:,1])      )
        setattr(ncvar,'valid_max'          , np.max(self.mesh.arrival_patch_latlon[:,1])      )
        setattr(ncvar,'_CoordinateAxisType', 'longitude'      )
        ncvar[:]=self.mesh.arrival_patch_latlon[:,1].astype(np.float64)

        ncvar=ncOUT.createVariable("lat",   'f', ('patch_number',))
        setattr(ncvar,'units'              , 'degrees_north')
        setattr(ncvar,'long_name'          , 'latitude of every patch')
        setattr(ncvar,'standard_name'      , 'lat')
        setattr(ncvar,'axis'               , 'Y'        )
        setattr(ncvar,'valid_min'          , np.min(self.mesh.arrival_patch_latlon[:,0])      )
        setattr(ncvar,'valid_max'          , np.max(self.mesh.arrival_patch_latlon[:,0])      )
        setattr(ncvar,'_CoordinateAxisType', 'latitude'      )
        ncvar[:]=self.mesh.arrival_patch_latlon[:,0].astype(np.float64)

        ncvar=ncOUT.createVariable("haz_map_surf",   'f', ('timestep_number','patch_number',))
        setattr(ncvar,'units'              , 'tons/km2')
        setattr(ncvar,'long_name'          , 'surface oil in open sea')
        setattr(ncvar,'valid_min'          , np.min(self.SurfOpnSea_tons_per_km2)      )
        setattr(ncvar,'valid_max'          , np.max(self.SurfOpnSea_tons_per_km2)      )
        ncvar[:,:]=self.SurfOpnSea_tons_per_km2[:,:].astype(np.float64)

        ncvar=ncOUT.createVariable("haz_map_surf_beach",   'f', ('timestep_number','patch_number',))
        setattr(ncvar,'units'              , 'tons/km2')
        setattr(ncvar,'long_name'          , 'surface oil beached on coasts')
        setattr(ncvar,'valid_min'          , np.min(self.SurfStrand_tons_per_km2)      )
        setattr(ncvar,'valid_max'          , np.max(self.SurfStrand_tons_per_km2)      )
        ncvar[:,:]=self.SurfStrand_tons_per_km2[:,:].astype(np.float64)

        ncvar=ncOUT.createVariable("haz_map_disp",   'f', ('timestep_number','patch_number',))
        setattr(ncvar,'units'              , 'tons/km2')
        setattr(ncvar,'long_name'          , 'dispersed oil in open sea')
        setattr(ncvar,'valid_min'          , np.min(self.DispOpnSea_tons_per_km2)      )
        setattr(ncvar,'valid_max'          , np.max(self.DispOpnSea_tons_per_km2)      )
        ncvar[:,:]=self.DispOpnSea_tons_per_km2[:,:].astype(np.float64)

        ncvar=ncOUT.createVariable("haz_map_disp_beach",   'f', ('timestep_number','patch_number',))
        setattr(ncvar,'units'              , 'tons/km2')
        setattr(ncvar,'long_name'          , 'dispersed oil beached on coasts')
        setattr(ncvar,'valid_min'          , np.min(self.DispStrand_tons_per_km2)      )
        setattr(ncvar,'valid_max'          , np.max(self.DispStrand_tons_per_km2)      )
        ncvar[:,:]=self.DispStrand_tons_per_km2[:,:].astype(np.float64)

      else:

        import netCDF4 as nc4
        lonmin=float(np.min(self.mesh.lon[:,:]))
        lonmax=float(np.max(self.mesh.lon[:,:]))
        latmin=float(np.min(self.mesh.lat[:,:]))
        latmax=float(np.max(self.mesh.lat[:,:]))
        ncOUT= nc4.Dataset(Binary_filename,"w",format='NETCDF4')
        print('create netcdf file',Binary_filename)
        ncOUT.createDimension("lat",   len(self.mesh.lon[:,0]))
        ncOUT.createDimension("lon",   len(self.mesh.lon[0,:]))
        ncOUT.createDimension("time",  self.tmax )
        ncOUT.creator_name="Celia Laurent"
        ncOUT.Institution="OGS National Institute of Oceanography and Applied Geophysics https://www.inogs.it/"
        ncOUT.geospatial_bounds_crs="EPSG:4326"
        ncOUT.geospatial_lon_units="degrees_east"
        ncOUT.geospatial_lat_units="degrees_north"
        ncOUT.geospatial_lat_max = np.max(self.mesh.lat[:,:]) ;
        ncOUT.geospatial_lat_min = np.min(self.mesh.lat[:,:]) ;
        ncOUT.geospatial_lat_resolution = (np.max(self.mesh.lat[:,:])-np.min(self.mesh.lat[:,:]))/(len(self.mesh.lon[:,0])-1) ;
        ncOUT.geospatial_lat_units = "degrees_north" ;
        ncOUT.geospatial_lon_max = np.max(self.mesh.lon[:,:]);
        ncOUT.geospatial_lon_min = np.min(self.mesh.lon[:,:]);
        ncOUT.geospatial_lon_resolution = (np.max(self.mesh.lon[:,:])-np.min(self.mesh.lon[:,:]))/(len(self.mesh.lon[:,0])-1) ;
        ncOUT.geospatial_lon_units = "degrees_east" ;
        ncOUT.geospatial_bounds = "POLYGON (("+str(lonmin)+" "+str(latmax)+", "+str(latmin)+" "+str(latmax)+", "+str(latmin)+" "+str(latmin)+", "+str(lonmin)+" "+str(latmin)+", "+str(lonmin)+" "+str(latmax)+"))" 
        ncOUT._CoordSysBuilder = "ucar.nc2.dataset.conv.CF1Convention" ;
        ncOUT.Conventions = "CF-1.5"
        ncOUT.cdm_data_type = "Grid" ;

        #ncOUT.geospatial_bounds = "POLYGON ((14. 43., 16. 43, 16. 40, 14. 40, 14. 43))" 
       #ncvar_crs=ncOUT.createVariable("crs",   'c')
       #ncvar_crs.grid_mapping_name = "latitude_longitude" ;
       #ncvar_crs.long_name = "CRS definition" ;
       #ncvar_crs.longitude_of_prime_meridian = 0. ;
       #ncvar_crs.semi_major_axis = 6378137. ;
       #ncvar_crs.inverse_flattening = 298.257223563 ;
       #ncvar_crs.spatial_ref = "GEOGCS[\"WGS 84\",DATUM[\"WGS_1984\",SPHEROID[\"WGS 84\",6378137,298.257223563,AUTHORITY[\"EPSG\",\"7030\"]],AUTHORITY[\"EPSG\",\"6326\"]],PRIMEM[\"Greenwich\",0,AUTHORITY[\"EPSG\",\"8901\"]],UNIT[\"degree\",0.0174532925199433,AUTHORITY[\"EPSG\",\"9122\"]],AXIS[\"Latitude\",NORTH],AXIS[\"Longitude\",EAST],AUTHORITY[\"EPSG\",\"4326\"]]" ;
       #ncvar_crs.GeoTransform = "0 0 0 0 0 0 " ;

        ncvar_lon=ncOUT.createVariable("lon",   'f', ('lat','lon',))
        ncvar_lon.units               = 'degrees_east'
        ncvar_lon.long_name           = 'longitude'
        ncvar_lon.standard_name       = 'longitude'
        ncvar_lon.axis                = 'X'        
        ncvar_lon.valid_min           = np.min(self.mesh.lon[:,:])      
        ncvar_lon.valid_max           = np.max(self.mesh.lon[:,:])      
        ncvar_lon._CoordinateAxisType = 'lon'      
        ncvar_lon.ioos_category       = 'Location'
        ncvar_lon.actual_range        = np.array([np.min(self.mesh.lon[:,:]),np.max(self.mesh.lon[:,:])])
        ncvar_lon[:,:]=self.mesh.lon[:,:].astype(np.float64)

        ncvar_lat=ncOUT.createVariable("lat",   'f', ('lat','lon',))
        ncvar_lat.units               = 'degrees_north'
        ncvar_lat.long_name           = 'latitude'
        ncvar_lat.standard_name       = 'latitude'
        ncvar_lat.axis                = 'Y'        
        ncvar_lat.valid_min           = np.min(self.mesh.lat[:,:])      
        ncvar_lat.valid_max           = np.max(self.mesh.lat[:,:])      
        ncvar_lat._CoordinateAxisType = 'lat'      
        ncvar_lat.ioos_category       = 'Location'
        ncvar_lat.actual_range        =np.array([np.min(self.mesh.lat[:,:]),np.max(self.mesh.lat[:,:])])
        ncvar_lat[:,:]=self.mesh.lat[:,:].astype(np.float64)

        ncvar_haz_map_surf=ncOUT.createVariable("haz_map_surf",   'f', ('time','lat','lon',))
        ncvar_haz_map_surf.units               = 'tons/km2'
        ncvar_haz_map_surf.ioos_category        = 'Oil'
        ncvar_haz_map_surf.long_name           = 'surface oil in open sea'
        ncvar_haz_map_surf.valid_min           = float(np.min(np.ma.masked_where(self.SurfOpnSea_tons_per_km2==0,self.SurfOpnSea_tons_per_km2)))      
        ncvar_haz_map_surf.valid_max           = np.max(self.SurfOpnSea_tons_per_km2)      
        #ncvar_haz_map_surf.grid_mapping = "crs"
        ncvar_haz_map_surf.missing_value=0
        ncvar_haz_map_surf[:,:,:]=self.SurfOpnSea_tons_per_km2[:,:,:].astype(np.float64)

        ncvar_haz_map_surf_beach=ncOUT.createVariable("haz_map_surf_beach",   'f', ('time','lat','lon',))
        ncvar_haz_map_surf_beach.units               = 'tons/km2'
        ncvar_haz_map_surf_beach.ioos_category        = 'Oil'
        ncvar_haz_map_surf_beach.long_name           = 'surface oil beached on coasts'
        ncvar_haz_map_surf_beach.valid_min           = float(np.min(np.ma.masked_where(self.SurfStrand_tons_per_km2==0,self.SurfStrand_tons_per_km2)))   
        ncvar_haz_map_surf_beach.valid_max           = np.max(self.SurfStrand_tons_per_km2)      
        #ncvar_haz_map_surf_beach.grid_mapping = "crs"
        ncvar_haz_map_surf_beach.missing_value=0
        ncvar_haz_map_surf_beach[:,:,:]=self.SurfStrand_tons_per_km2[:,:,:].astype(np.float64)

        haz_map_surf_sum=np.sum(self.SurfOpnSea_tons_per_km2,axis=0)
        ncvar_haz_map_surf_sum=ncOUT.createVariable("haz_map_surf_sum",   'f', ('lat','lon',))
        ncvar_haz_map_surf_sum.units               = 'tons/km2'
        ncvar_haz_map_surf_sum.ioos_category        = 'Oil'
        ncvar_haz_map_surf_sum.long_name           = 'cumulative sum (in time) of the surface oil in open sea'
        ncvar_haz_map_surf_sum.valid_min           = float(np.min(np.ma.masked_where(haz_map_surf_sum==0,haz_map_surf_sum)))      
        ncvar_haz_map_surf_sum.valid_max           = np.max(haz_map_surf_sum)      
        #ncvar_haz_map_surf_sum.grid_mapping = "crs"
        ncvar_haz_map_surf_sum.missing_value=0
        ncvar_haz_map_surf_sum[:,:]=haz_map_surf_sum[:,:].astype(np.float64)

        ncvar_haz_map_disp=ncOUT.createVariable("haz_map_disp",   'f', ('time','lat','lon',))
        ncvar_haz_map_disp.units               = 'tons/km2'
        ncvar_haz_map_disp.ioos_category        = 'Oil'
        ncvar_haz_map_disp.long_name           = 'dispersed oil in open sea'
        ncvar_haz_map_disp.valid_min           = float(np.min(np.ma.masked_where(self.DispOpnSea_tons_per_km2==0,self.DispOpnSea_tons_per_km2)))      
        ncvar_haz_map_disp.valid_max           = np.max(self.DispOpnSea_tons_per_km2)      
        #ncvar_haz_map_disp.grid_mapping = "crs"
        ncvar_haz_map_disp.missing_value=0
        ncvar_haz_map_disp[:,:,:]=self.DispOpnSea_tons_per_km2[:,:,:].astype(np.float64)

        ncvar_haz_map_disp_beach=ncOUT.createVariable("haz_map_disp_beach",   'f', ('time','lat','lon',))
        ncvar_haz_map_disp_beach.units               = 'tons/km2'
        ncvar_haz_map_disp_beach.ioos_category        = 'Oil'
        ncvar_haz_map_disp_beach.long_name           = 'dispersed oil beached on coasts'
        ncvar_haz_map_disp_beach.valid_min           = float(np.min(np.ma.masked_where(self.DispStrand_tons_per_km2==0,self.DispStrand_tons_per_km2)))      
        ncvar_haz_map_disp_beach.valid_max           = np.max(self.DispStrand_tons_per_km2)      
        #ncvar_haz_map_disp_beach.grid_mapping = "crs"
        ncvar_haz_map_disp_beach.missing_value=0
        ncvar_haz_map_disp_beach[:,:,:]=self.DispStrand_tons_per_km2[:,:,:].astype(np.float64)

        haz_map_disp_sum=np.sum(self.DispOpnSea_tons_per_km2,axis=0)
        ncvar_haz_map_disp_sum=ncOUT.createVariable("haz_map_disp_sum",   'f', ('lat','lon',))
        ncvar_haz_map_disp_sum.units               = 'tons/km2'
        ncvar_haz_map_disp_sum.ioos_category        = 'Oil'
        ncvar_haz_map_disp_sum.long_name           = 'cumulative sum (in time) of the dispersed oil in open sea'
        ncvar_haz_map_disp_sum.valid_min           = float(np.min(np.ma.masked_where(haz_map_disp_sum==0,haz_map_disp_sum)))      
        ncvar_haz_map_disp_sum.valid_max           = np.max(haz_map_disp_sum)      
        #ncvar_haz_map_disp_sum.grid_mapping = "crs"
        ncvar_haz_map_disp_sum.missing_value=0
        ncvar_haz_map_disp_sum[:,:]=haz_map_disp_sum[:,:].astype(np.float64)

      print(self.time/86400.0) 
      ncvar_time=ncOUT.createVariable("time",   'f', ('time',))
      ncvar_time.units               = 'days'
      ncvar_time.long_name           = 'snapshot-time in days after the release'
      ncvar_time.valid_min           = np.min(self.time/86400.0)      
      ncvar_time.valid_max           = np.max(self.time/86400.0)      
      ncvar_time.ioos_category       = 'time'
      ncvar_time[:]=(self.time/86400.0).astype(np.float64)
      
      # ncvar=ncOUT.createVariable("ini_time",   'f', ('time',)
      # ncvar.units                 'days'
      # ncvar.long_name             'starting-date of the time intervals in days after the release'
      # ncvar.valid_min             np.min(self.ini_timeinterval)      
      # ncvar.valid_max             np.max(self.ini_timeinterval)      
      # ncvar[:]=self.ini_timeinterval.astype(np.float64)
        
      # ncvar=ncOUT.createVariable("fin_time",   'f', ('time',)
      # ncvar.units                 'days'
      # ncvar.long_name             'ending-date of the time intervals in days after the release'
      # ncvar.valid_min             np.min(self.fin_timeinterval)      
      # ncvar.valid_max             np.max(self.fin_timeinterval)      
      # ncvar[:]=self.fin_timeinterval.astype(np.float64)
      
      print('netcdf file',Binary_filename,' ready to be written\n')
      ncOUT.close()
      print('netcdf file',Binary_filename,' written and closed\n')

  def __check_common_file_properties(self,namefiles):
      for inputfile,outputfilename in enumerate(namefiles): 
        changes=False
        nmax,tmax,dt,lon0,lat0,depth0=self.__get_file_properties(outputfilename)
        if(inputfile!=0):
          if(nmax!=prev_nmax or tmax!=prev_tmax or dt!=prev_dt): 
            print('ERROR!! file properties changed, dimensions are ',(nmax,tmax,dt),' != ',(prev_nmax,prev_tmax,prev_dt),' in ',outputfilename)
            sys.exit()
          #if( np.max(np.abs(lon0-lon0_prev))>1e-8
          #   or np.max(np.abs(lat0-lat0_prev))>1e-8
          #   or np.max(np.abs(depth0-depth0_prev))>1e-8  ):
          #  changes=True
          #if(changes):print('WARNING, you are opening files containing particles having different initial positions')
        (prev_nmax,prev_tmax,prev_dt,lon0_prev,lat0_prev,depth0_prev)=(nmax,tmax,dt,lon0,lat0,depth0)
      return nmax,tmax,dt,lon0,lat0,depth0

  def __get_file_properties(self,outputfilename):
    try:      f= netCDF4.Dataset(outputfilename,'r')
    except:   f= NC.netcdf_file(outputfilename,'r')
    lon0=np.copy(f.variables['lon'][0,:])
    lat0=np.copy(f.variables['lat'][0,:])
    depth0=np.copy(f.variables['depth'][0,:])
    nmax=len(depth0)
    tmax=len(f.variables['model_time'][:])
    try:dt=f.variables['model_time'][2]-f.variables['model_time'][1]
    except:
      print('ERROR seems that file',outputfilename,' does not contain more than 1 timestep')
      sys.exit()
    f.close()
    return nmax,tmax,dt,lon0,lat0,depth0

class LTRANSnetcdf():
  def __init__(self,outputfilename):
    try:      f= netCDF4.Dataset(outputfilename,'r')
    except:   f= NC.netcdf_file(outputfilename,'r')
    self.time =np.copy(f.variables['model_time'][:])
    self.depth=np.copy(f.variables['depth'][:,:])
    self.lon  =np.copy(f.variables['lon'][:,:])
    self.lat  =np.copy(f.variables['lat'][:,:])
    try: self.color=np.copy(f.variables['color'][:,:])
    except:pass
    try: self.CoastDist=np.copy(f.variables['CoastDist'][:,:])
    except:pass
    try: self.Temp=np.copy(f.variables['Temp'][:,:])
    except:pass
    try: self.Salt=np.copy(f.variables['Salt'][:,:])
    except:pass
    f.close()

  def remove_timesteps_preceeding(self,t0):
    self.time = self.time[t0:]
    self.depth = self.depth[t0:]
    self.lon = self.lon[t0:]
    self.lat = self.lat[t0:]
    try:         self.color     = self.color[t0:]     
    except:pass                                  
    try:         self.CoastDist = self.CoastDist[t0:] 
    except:pass                                  
    try:         self.Temp      = self.Temp[t0:]      
    except:pass                                  
    try:         self.Salt      = self.Salt[t0:]      
    except:pass

class LTRANSgrid():
  def __init__(self,gridfilename,parallel):
        self.__rank=parallel.rank
        self.__nranks=parallel.nranks
        self.__comm=parallel.comm
        if(self.__rank<=0):
                try:
                        import netCDF4
                        g = netCDF4.Dataset(gridfilename,'r')
                except:
                        import scipy.io.netcdf as NC
                        g = NC.netcdf_file(gridfilename,'r')
        else:
                g = None
        if(self.__rank<=0): 
          self.lon=np.copy(g.variables['lon_rho'][:,:])
          self.lat=np.copy(g.variables['lat_rho'][:,:])
          self.lonu=np.copy(g.variables['lon_u'][:,:])
          self.latu=np.copy(g.variables['lat_u'][:,:])
          self.lonv=np.copy(g.variables['lon_v'][:,:])
          self.latv=np.copy(g.variables['lat_v'][:,:])
          try:self.Zp1= np.copy(g.variables['Zp1'][:])
          except:self.Zp1=None
          self.hvar  = np.copy(g.variables['h'][:,:])
          self.mask_rho= np.copy(g.variables['mask_rho'][:,:,:])
          print('\nGrid file',gridfilename,' contains grid-mesh of shape (Nk,Nj,Ni)=',np.shape(self.mask_rho))
        else:
          self.lon     = None
          self.lat     = None
          self.lonu    = None
          self.latu    = None
          self.latu    = None
          self.lonv    = None
          self.hvar    = None
          self.Zp1     =  None
          self.mask_rho= None
        if(self.__nranks>1): 
          self.lon     = self.__comm.bcast(self.lon     , root=0)
          self.lat     = self.__comm.bcast(self.lat     , root=0) 
          self.lonu    = self.__comm.bcast(self.lonu    , root=0) 
          self.latu    = self.__comm.bcast(self.latu    , root=0) 
          self.latu    = self.__comm.bcast(self.latu    , root=0) 
          self.lonv    = self.__comm.bcast(self.lonv    , root=0) 
          self.hvar    = self.__comm.bcast(self.hvar    , root=0) 
          self.Zp1     = self.__comm.bcast(self.Zp1     , root=0) 
          self.mask_rho= self.__comm.bcast(self.mask_rho, root=0) 

        self.depth=-self.hvar
        self.land=np.ma.masked_where(self.depth<0,self.depth)
        if(self.__rank<=0):g.close()
        self.Ni =len(self.lon[0,:])
        self.Nj =len(self.lat[:,0])
        try:self.Nk =len(self.Zp1[:])-1
        except:
                self.Nk =20
                print('Zp1 not found in input file, taking Nk=',self.Nk)
        self.dlon=self.lon[0,1]-self.lon[0,0]
        self.dlat=self.lat[1,0]-self.lat[0,0]
        #self.compute_dX_dY_Area()
        self.lat1            = np.min(self.lat[:,:])
        self.lat2            = np.max(self.lat[:,:])
        self.lon1            = np.min(self.lon[:,:])
        self.lon2            = np.max(self.lon[:,:])

  def compute_dX_dY_Area(self):
    self.dX=np.zeros((self.Nj,self.Ni),dtype=float)
    self.dY=np.zeros((self.Nj,self.Ni),dtype=float)
    for j in range(self.Nj):
      for i in range(self.Ni):
        self.dX[j,i]=geopy.distance.geodesic((self.lon[j,i]-0.5*self.dlon,self.lat[j,i]),
                                             (self.lon[j,i]+0.5*self.dlon,self.lat[j,i])).meters
        self.dY[j,i]=geopy.distance.geodesic((self.lon[j,i],self.lat[j,i]-0.5*self.dlat),
                                             (self.lon[j,i],self.lat[j,i]+0.5*self.dlat)).meters
    self.Area_m2=self.dX*self.dY

  def read_coastseg(self,blnfilename,Xmin=None,Ymin=None,Xmax=None,Ymax=None):
        if(self.__rank<=0): print('reading boundary file '+blnfilename)
        if(Xmin==None):Xmin=np.min(self.lon)
        if(Xmax==None):Xmax=np.max(self.lon)
        if(Ymin==None):Ymin=np.min(self.lat)
        if(Ymax==None):Ymax=np.max(self.lat)
        self.Coast_ncontours=0
        count=0
        revcount=0
        self.Coast_lon=np.array([])
        self.Coast_lat=np.array([])
        num=0 
        if(self.__rank<=0):
                f = open(blnfilename, 'r')        
                self.Coast_ismainbnd=np.array([])
                self.Coast_numcorner=np.array([])
                reader=csv.reader(f,delimiter=",")
                for row in reader:
                        if(revcount==0):
                                self.Coast_numcorner=np.hstack((self.Coast_numcorner,int(row[0])))
                                self.Coast_ismainbnd=np.hstack((self.Coast_ismainbnd,int(row[1])))
                                revcount=self.Coast_numcorner[-1]
                                self.Coast_ncontours += 1
                                writeCoastklvandnum=True
                        else:
                                lonread=float(row[0])
                                latread=float(row[1])
                                if((lonread > Xmin and lonread < Xmax) and (latread > Ymin and latread < Ymax)):
                                        self.Coast_lon=np.hstack((self.Coast_lon,lonread))        
                                        self.Coast_lat=np.hstack((self.Coast_lat,latread))
                                        if(writeCoastklvandnum):        
                                                writeCoastklvandnum=False
                                                klv=int(row[2])
                                                if(num==0):
                                                        self.Coast_klv=np.array([klv])
                                                        self.Coast_num=np.array([0])
                                                else:
                                                        self.Coast_klv=np.hstack((self.Coast_klv,int(klv)))        
                                                        self.Coast_num=np.hstack((self.Coast_num,int(num)))        
                                        num +=1
                                revcount -= 1
                f.close()
        else:
                self.Coast_lon      =None
                self.Coast_lat      =None
                self.Coast_klv      =None
                self.Coast_num      =None
                self.Coast_ncontours=None
        if(self.__nranks>1): 
                self.Coast_lon      =self.__comm.bcast( self.Coast_lon       , root=0)
                self.Coast_lat      =self.__comm.bcast( self.Coast_lat       , root=0)
                self.Coast_klv      =self.__comm.bcast( self.Coast_klv       , root=0)
                self.Coast_num      =self.__comm.bcast( self.Coast_num       , root=0)
                self.Coast_ncontours=self.__comm.bcast( self.Coast_ncontours , root=0)


  def contourf_land(self,ax,contourf_args={'cmap':'Greys','vmin':-0.5,'vmax':0.5}):
    _=ax.contourf(self.lon,self.lat,self.land, **contourf_args)
    try:zorder=contourf_args['zorder']+1
    except:zorder=999
    print('zorder=',zorder)
    _=ax.contour(self.lon,self.lat,self.depth,levels=[-0.1,0,0.1],colors='k',linewidths=2,zorder=zorder,linestyles='solid')

class MeshPatch():
  def __init__(self,grid,data=None,parallel=False):
    try:
     self.__rank=parallel.rank
     self.__nranks=parallel.nranks
     self.__comm=parallel.comm
    except:
     if(parallel==False):  
       self.__rank  =0
       self.__nranks=1
       self.__comm  =None
    try:
      self.part_ini_lon=data.lon0
      self.part_ini_lat=data.lat0
      self.nmax=data.nmax
    except:
      print('DATA object was provided, particles release patches cannot be identified')
    self.__grid=grid
    self.__parallel=parallel
    self.__coords_from_grid=False

  def def_meshpatch_coords_from_grid(self,iskip,jskip):
    self.__coords_from_grid=True
    self.Nlon=int(math.ceil(len(self.__grid.lon[0,:])/iskip))
    self.Nlat=int(math.ceil(len(self.__grid.lat[:,0])/jskip))
    vecdlon=self.__grid.lon[0,1:]-self.__grid.lon[0,:-1]
    vecdlat=self.__grid.lat[1:,0]-self.__grid.lat[:-1,0]
    self.dlon=0.5*(np.min(vecdlon)+np.max(vecdlon))*iskip
    self.dlat=0.5*(np.min(vecdlat)+np.max(vecdlat))*jskip
    self.lon1=np.average(self.__grid.lon[0,:iskip])
    self.lat1=np.average(self.__grid.lat[:jskip,0])
    self.Veclon=np.arange(self.lon1,self.lon1+(self.Nlon-0.9)*self.dlon,self.dlon)  
    self.Veclat=np.arange(self.lat1,self.lat1+(self.Nlat-0.9)*self.dlat,self.dlat)  
    self.lon2=self.Veclon[-1]
    self.lat2=self.Veclat[-1]
    self.iskip=iskip
    self.jskip=jskip
    if(self.__parallel.rank<1):print('Coordinates of the patches defined on a mesh of size:',(self.Nlat,self.Nlon))

  def def_meshpatch_coords_from_scratch(self,delta_lon,delta_lat,lat1=None,lat2=None,lon1=None,lon2=None):
    self.dlon=delta_lon 
    self.dlat=delta_lat 
    self.lat1            = np.min(self.__grid.lat[:,:])
    self.lat2            = np.max(self.__grid.lat[:,:])
    self.lon1            = np.min(self.__grid.lon[:,:])
    self.lon2            = np.max(self.__grid.lon[:,:])
    if(lon1!=None):self.lon1=lon1
    if(lon2!=None):self.lon2=lon2
    if(lat1!=None):self.lat1=lat1
    if(lat2!=None):self.lat2=lat2
    self.Veclon=np.arange(self.lon1,self.lon2+self.dlon,self.dlon)
    self.Veclat=np.arange(self.lat1,self.lat2+self.dlat,self.dlat) 
    self.Nlon=len(self.Veclon)
    self.Nlat=len(self.Veclat)
    print('WARNING : routine mesh_from_scratch has not been checked')

  def def_meshpatch_properties(self):
    self.__mesh2d()
    self.__compute_grid_cells_in_patches()
    self.__compute_water_cells_in_patches()

  def __mesh2d(self):
    self.lon,self.lat=np.meshgrid(self.Veclon,self.Veclat)
    self.CPATCHES_patch_id=np.arange(0,self.Nlon*self.Nlat).reshape(np.shape(self.lon))
    #
    self.__compute_dX_dY_Area()

  def __compute_dX_dY_Area(self):
    self.Vecglon=np.hstack((self.Veclon[0]-0.5*self.dlon,0.5*(self.Veclon[1:]+self.Veclon[0:-1]),self.Veclon[-1]+0.5*self.dlon))
    self.Veclat_interface=np.hstack((self.Veclat[0]-0.5*self.dlat,0.5*(self.Veclat[1:]+self.Veclat[0:-1]),self.Veclat[-1]+0.5*self.dlat))
    self.dX=np.zeros(np.shape(self.lon),dtype=float)
    self.dY=np.zeros(np.shape(self.lat),dtype=float)
    for j in range(self.Nlat):
      for i in range(self.Nlon):
        self.dX[j,i]=geopy.distance.geodesic((self.Vecglon[i],self.Veclat_interface[j]),
                                       (self.Vecglon[i+1],self.Veclat_interface[j])).meters
        self.dY[j,i]=geopy.distance.geodesic((self.Vecglon[i],self.Veclat_interface[j]),
                                       (self.Vecglon[i],self.Veclat_interface[j+1])).meters
    if(self.__parallel.rank<1):print('patch dlon=',self.dlon,' -> dX in range [',np.min(self.dX),np.max(self.dX),']')
    if(self.__parallel.rank<1):print('patch dlat=',self.dlat,' -> dY in range [',np.min(self.dY),np.max(self.dY),']')
    self.Area_km2=self.dX/1000.0*self.dY/1000.0
    self.Area_m2=self.dX*self.dY
    if(self.__parallel.rank<1):print('Built 2d mesh of shape:',np.shape(self.Area_km2),' with an average patch-area = {:8.5f} km2'.format(np.average(self.Area_km2)))
    self.CPATCHES_Elon=self.lon-self.dlon
    self.CPATCHES_Elat=self.lat-self.dlat
    self.CPATCHES_Elon,self.CPATCHES_Elat=np.meshgrid(self.Vecglon,self.Veclat_interface)

  def __compute_grid_cells_in_patches(self):
    self.CPATCHES_NumGridCells=np.zeros(np.shape(self.lon),dtype=int)
    Iplon=np.int0(np.round((self.__grid.lon[0,:]-self.lon1)/self.dlon))
    Iplat=np.int0(np.round((self.__grid.lat[:,0]-self.lat1)/self.dlat))
    for j in Iplat:
      for i in Iplon:
         self.CPATCHES_NumGridCells[j,i]+=1
    if(self.__parallel.rank<1):print('Total number of (surface) grid-cells represented by the mesh:',np.sum(self.CPATCHES_NumGridCells))
        
  def __compute_water_cells_in_patches(self):
    self.CPATCHES_NumWaterCells=np.zeros(np.shape(self.lon),dtype=int)
    Iplon=np.int0(np.round((self.__grid.lon[0,:]-self.lon1)/self.dlon))
    Iplat=np.int0(np.round((self.__grid.lat[:,0]-self.lat1)/self.dlat))
    for jg,j in enumerate(Iplat):
      for ig,i in enumerate(Iplon):
         if(self.__grid.mask_rho[-1,jg,ig]>0):self.CPATCHES_NumWaterCells[j,i]+=1
    if(self.__parallel.rank<1):print('Total number of (surface) grid-water-cells represented by the mesh:',np.sum(self.CPATCHES_NumWaterCells))
        
  def __def_release_patches_by_part(self):
    Iplon=np.round((self.part_ini_lon-self.lon1)/self.dlon)
    Jplat=np.round((self.part_ini_lat-self.lat1)/self.dlat)
    self.release_patch_id_by_part=np.int64(Jplat[:]*float(self.Nlon)+Iplon)
    if(self.__parallel.rank<1):print('self.release_patch_id_by_part min,max=',np.min(self.release_patch_id_by_part),np.max(self.release_patch_id_by_part))

  def def_release_patches_from_patches_with_parts(self):
    self.__def_release_patches_by_part()
    self.release_patch_by_part=np.zeros(self.nmax,dtype=int)
    self.patch_id_per_release_patches=np.unique(np.sort(self.release_patch_id_by_part))
    for n in range(self.nmax):
      self.release_patch_by_part[n]=np.where(self.release_patch_id_by_part[n]==self.patch_id_per_release_patches)[0]
    self.nRELpatch=len(self.patch_id_per_release_patches)
    if(self.__parallel.rank<1):print(self.nRELpatch,' release patches (patches where particles were released)\n')

  def def_release_patches_from_wet_patches(self):
    self.patch_id_per_release_patches=np.where( self.CPATCHES_NumWaterCells.flatten()>0 )[0]
    self.nRELpatch=len(self.patch_id_per_release_patches)
    if(self.__parallel.rank<1):print(self.nRELpatch,' release patches (patches where particles were released)\n')

  def def_arrival_patches_from_wet_patches(self):
    self.patch_id_per_arrival_patches=np.where( self.CPATCHES_NumWaterCells.flatten()>0 )[0]
    self.nARRpatch=len(self.patch_id_per_arrival_patches)
    if(self.__parallel.rank<1):print(self.nARRpatch,' arrival patches (one per water patch of the mesh)')
    self.__compute_ij2Apatch()

  def compute_Area_by_Rpatch(self):
    try:
      self.release_patch_jicoordinates[0,0]
    except:
      self.__compute_patch_jicoordinates()
    self.Rpatch_Area_m2=np.zeros((self.nRELpatch),dtype=int)
    for patch in range(self.nRELpatch):
      (j,i)=self.release_patch_jicoordinates[patch,:]
      self.Rpatch_Area_m2[patch]=self.Area_m2[j,i]
    return self.Rpatch_Area_m2

  def compute_Area_by_Apatch(self):
    try:
      self.arrival_patch_jicoordinates[0,0]
    except:
      self.__compute_patch_jicoordinates()
    self.Apatch_Area_m2=np.zeros((self.nRELpatch),dtype=int)
    for patch in range(self.nARRpatch):
      (j,i)=self.arrival_patch_jicoordinates[patch,:]
      self.Apatch_Area_m2[patch]=self.Area_m2[j,i] 
    return self.Apatch_Area_m2

  def def_parts_release_patches(self):
    self.__def_release_patches_by_part()
    self.release_patch_by_part=np.zeros(self.nmax,dtype=int)
    for n in range(self.nmax):
      self.release_patch_by_part[n]=np.where(self.release_patch_id_by_part[n]==self.patch_id_per_release_patches)[0]
   
  def __compute_ij2Apatch(self):
    self.__ij2Apatch=np.zeros(np.shape(self.lon),dtype=int)
    self.__ij2Apatch[:,:]=-1
    for j in range(self.Nlat):
      for i in range(self.Nlon):
        id_patch=j*self.Nlon+i
        Apatch=np.where(self.patch_id_per_arrival_patches==id_patch)[0]
        if(len(Apatch>0)):
          self.__ij2Apatch[j,i]=Apatch

  def __compute_patch_jicoordinates(self):
    self.release_patch_jicoordinates=np.zeros((self.nRELpatch,2),dtype=int)
    self.arrival_patch_jicoordinates=np.zeros((self.nARRpatch,2),dtype=int)
    for j in range(self.Nlat):
      for i in range(self.Nlon):
        id_patch=j*self.Nlon+i
        Apatch=np.where(self.patch_id_per_arrival_patches==id_patch)[0]
        if(len(Apatch>0)):
          self.arrival_patch_jicoordinates[Apatch,:]=[j,i]
        Rpatch=np.where(self.patch_id_per_release_patches==id_patch)[0]
        if(len(Rpatch>0)):
          self.release_patch_jicoordinates[Rpatch,:]=[j,i]

  def compute_patch_latlon(self):
    self.release_patch_latlon=np.zeros((self.nRELpatch,2),dtype=float)
    self.arrival_patch_latlon=np.zeros((self.nARRpatch,2),dtype=float)
    for j in range(self.Nlat):
      for i in range(self.Nlon):
        lat=self.__grid.lat[j,i]
        lon=self.__grid.lon[j,i]
        id_patch=j*self.Nlon+i
        Apatch=np.where(self.patch_id_per_arrival_patches==id_patch)[0]
        if(len(Apatch>0)):
          self.arrival_patch_latlon[Apatch,:]=(lat,lon)
        Rpatch=np.where(self.patch_id_per_release_patches==id_patch)[0]
        if(len(Rpatch>0)):
          self.release_patch_latlon[Rpatch,:]=(lat,lon)

  def coords_to_Apatch(self,lon,lat):
     jlist=np.rint((lat-self.lat1)/self.dlat).astype(int)
     ilist=np.rint((lon-self.lon1)/self.dlon).astype(int)
     Apatches=self.__ij2Apatch[jlist,ilist] 
     return Apatches 

class OILTRANS():
  def __init__(self,Files_names,stranded_oil_was_removed_from_oiltot=False,files_range=[None],rank=0):
    if( len(np.shape(Files_names))==0 ):
      if('-OilPropsOut.csv' in Files_names):
        NumFiles=1
        Files=[Files_names]
      else:
        csvfile=open(Files_names,'r')
        csvreader=csv.reader(csvfile,delimiter=",")
        Files=()
        for row in csvreader:
          Files+=(strtools.TRIM(row[0]),)
        csvfile.close()
        NumFiles=len(Files)
        Files=np.array(Files)
    else:
      NumFiles=len(Files_names)
      Files=np.array(Files_names)
    if(len(files_range)>1):
        Files=Files[files_range]
        NumFiles=len(Files)
    self.maxlines=0
    minlines=999999999999
    self.fnames=()
    for f in range(0,NumFiles):
        filename=strtools.replacestring(Files[f],' ','')
        self.fnames+=(filename,)
        file = open(filename, 'r')
        reader = csv.reader(file,delimiter=',')
        count=0
        for row in reader:
                count=count+1
        self.maxlines=max(self.maxlines,count)
        minlines=min(minlines,count)
        file.close()
    #print 'self.maxlines=',self.maxlines
    #print 'minlines=',minlines
    if(minlines!=self.maxlines):
       print('error num lines changing in OilProp files',minlines,self.maxlines)
       print(NumFiles,' files')
       print(Files)
       sys.exit('error')
    self.R0oil=np.zeros((NumFiles), dtype=float)
    self.RFoil=np.zeros((NumFiles), dtype=float)
    self.T0oil=np.zeros((NumFiles), dtype=float)
    self.TFoil=np.zeros((NumFiles), dtype=float)
    self.Time=np.zeros((NumFiles,self.maxlines), dtype=float)
    self.VolSurf=np.zeros((NumFiles,self.maxlines), dtype=float)
    self.VolOilTot=np.zeros((NumFiles,self.maxlines), dtype=float)
    self.VolEvap=np.zeros((NumFiles,self.maxlines), dtype=float)
    self.VolDisp=np.zeros((NumFiles,self.maxlines), dtype=float)
    self.VolStrand=np.zeros((NumFiles,self.maxlines), dtype=float)
    self.WAngle=np.zeros((NumFiles,self.maxlines), dtype=float)
    self.Wind=np.zeros((NumFiles,self.maxlines), dtype=float)
    self.Temp=np.zeros((NumFiles,self.maxlines), dtype=float)
    self.Rho_Oil=np.zeros((NumFiles,self.maxlines), dtype=float)
    self.ViscOil=np.zeros((NumFiles,self.maxlines), dtype=float)
    self.AreaOil=np.zeros((NumFiles,self.maxlines), dtype=float)
    self.ThckOil=np.zeros((NumFiles,self.maxlines), dtype=float)
    self.Vol_Slick=np.zeros((NumFiles,self.maxlines), dtype=float)
    self.MassOil=np.zeros((NumFiles,self.maxlines), dtype=float)
    self.Time[:,:]=-99999.0
    self.Time[:,0]=0.0
    self.VolSurf[:,0]=100.0
    old_t0=-1
    self.time0=0
    for f in range(0,NumFiles):
       #print 'OILTRANS file ',f,' is ',Files[f]
        filename=strtools.replacestring(Files[f],' ','')
        #print ('rank',rank,' : library OILTRANS reads file ',filename)
        file = open(filename, 'r')
        reader = csv.reader(file,delimiter=',')
        count=-1
        for row in reader:
                #print 'row:',row
                if(count==-1):
                   count=count+1
                   continue
                self.Time[f,count]=float(row[0])
                self.VolSurf[f,count]=self.__cancelstar(row[1]) # volume[f] 
                if(self.T0oil[f]!=0 and self.TFoil[f]==0 and self.VolSurf[f,count]==0.0):
                       self.RFoil[f]=count
                       self.TFoil[f]=self.Time[f,count] 
                       #print 'file ',f,' tFoil=',self.TFoil[f]*86400,' at rec ',self.RFoil[f]
                self.VolEvap[f,count]=self.__cancelstar(row[2]) # volume
                self.VolDisp[f,count]=self.__cancelstar(row[3]) # volume
                self.VolStrand[f,count]=self.__cancelstar(row[4]) # volume
                if(stranded_oil_was_removed_from_oiltot):
                  self.VolOilTot[f,count]=self.VolSurf[f,count]+self.VolEvap[f,count]+self.VolDisp[f,count]
                else:
                  self.VolOilTot[f,count]=self.VolSurf[f,count]+self.VolEvap[f,count]+self.VolDisp[f,count]+self.VolStrand[f,count]
                if(self.T0oil[f]==0 and self.VolOilTot[f,count]>1.0):
                       self.R0oil[f]=count
                       self.T0oil[f]=self.Time[f,count]
                       #print f,'file ',Files[f],' tOoil=',self.T0oil[f],' at rec ',self.R0oil
                #print count,Surf[f,count],OilTot[f,count]
                #f(self.VolOilTot[f,count]>0.0):
                # self.VolSurf[f,count]=self.VolSurf[f,count]/self.VolOilTot[f,count]     *100.0
                # self.VolEvap[f,count]=self.VolEvap[f,count]/self.VolOilTot[f,count]     *100.0
                # self.VolDisp[f,count]=self.VolDisp[f,count]/self.VolOilTot[f,count]     *100.0
                # if(stranded_oil_was_removed_from_oiltot not True):
                #   self.VolStrand[f,count]=self.VolStrand[f,count]/self.VolOilTot[f,count]     *100.0
                # #print 'read ',count,self.VolSurf[f,count],self.VolEvap[f,count],self.VolDisp[f,count],self.VolStrand_Removed[f,count]
                substract=1
                self.Rho_Oil[f,count]=self.__cancelstar(row[7-substract])                      #  in kg/m3
                self.ViscOil[f,count]=self.__cancelstar(row[8-substract])                      #  
                self.AreaOil[f,count]=self.__cancelstar(row[9-substract])                      # self.AreaOil[f,count]=self.VolSurf[f,count]/self.ThckOil[f,count]
                self.ThckOil[f,count]=self.__cancelstar(row[10-substract])                     # self.ThckOil[f,count]=self.VolSurf[f,count]/self.AreaOil[f,count]
                self.Vol_Slick[f,count]=self.__cancelstar(row[11-substract])                   # self.Vol_Slick[f,count] = self.VolSurf[f,count]*(1-watercontent)
                self.MassOil[f,count]=self.__cancelstar(row[12-substract])                     # self.MassOil[f,count]=self.Rho_Oil[f,count]*self.VolSurf[f,count]
                if(len(row)>=15-substract):
                        self.Wind[f,count]=self.__cancelstar(row[13-substract])
                        self.Temp[f,count]=self.__cancelstar(row[14-substract])
                if(len(row)==16-substract):
                        self.WAngle[f,count]=self.__cancelstar(row[15-substract])
                #print '--------------------------------------------------'
                count=count+1
        file.close()
        self.t0=np.where(self.VolSurf[f,:]>0)[0][0]
        if(old_t0>=0 and self.t0!=old_t0):
          print('ERROR t0 Oil Changed')
          sys.exit()
        old_t0=self.t0
        if(self.t0>0): 
          self.time0=self.Time[f,self.t0]
          self.Time[f,:] -= self.time0
    if(self.t0>0):
            self.Time      =  self.Time[:,self.t0:]      
            self.VolSurf   =  self.VolSurf[:,self.t0:]   
            self.VolOilTot =  self.VolOilTot[:,self.t0:]   
            self.VolEvap   =  self.VolEvap[:,self.t0:]   
            self.VolDisp   =  self.VolDisp[:,self.t0:]   
            self.VolStrand =  self.VolStrand[:,self.t0:]   
            self.WAngle    =  self.WAngle[:,self.t0:]     
            self.Wind      =  self.Wind[:,self.t0:]      
            self.Temp      =  self.Temp[:,self.t0:]      
            self.Rho_Oil   =  self.Rho_Oil[:,self.t0:]   
            self.ViscOil   =  self.ViscOil[:,self.t0:]   
            self.AreaOil   =  self.AreaOil[:,self.t0:]   
            self.ThckOil   =  self.ThckOil[:,self.t0:]   
            self.Vol_Slick =  self.Vol_Slick[:,self.t0:]   
            self.MassOil   =  self.MassOil[:,self.t0:]
    if(self.Time[0,-1]<-99998):
            self.Time[:,-1]      =  self.Time[:,-2]    + (self.Time[:,-2]-self.Time[:,-3]) 
            self.VolSurf[:,-1]   =  self.VolSurf[:,-2]   
            self.VolOilTot[:,-1] =  self.VolOilTot[:,-2]   
            self.VolEvap[:,-1]   =  self.VolEvap[:,-2]   
            self.VolDisp[:,-1]   =  self.VolDisp[:,-2]   
            self.VolStrand[:,-1] =  self.VolStrand[:,-2]   
            self.WAngle[:,-1]    =  self.WAngle[:,-2]     
            self.Wind[:,-1]      =  self.Wind[:,-2]      
            self.Temp[:,-1]      =  self.Temp[:,-2]      
            self.Rho_Oil[:,-1]   =  self.Rho_Oil[:,-2]   
            self.ViscOil[:,-1]   =  self.ViscOil[:,-2]   
            self.AreaOil[:,-1]   =  self.AreaOil[:,-2]   
            self.ThckOil[:,-1]   =  self.ThckOil[:,-2]   
            self.Vol_Slick[:,-1] =  self.Vol_Slick[:,-2]   
            self.MassOil[:,-1]   =  self.MassOil[:,-2]
    if(self.t0>0):
      print('WARNING !!! OILSPILL DATA modified by changing Time : ',np.hstack((self.Time[0,:2],['...'],self.Time[0,-2:])) )

  def __cancelstar(self,string):
                if(string=='***************' or string=='NaN'):
                  return 0.0
                else:
                  return float(string) 

class csv_to_pd_dataframe():
  def __init__(self,Files_names,parallel):
    import pandas as pd
    if( len(np.shape(Files_names))==0 ):
      csvfile=open(Files_names,'r')
      csvreader=csv.reader(csvfile,delimiter=",")
      Files=()
      for row in csvreader:
        Files+=(strtools.TRIM(row[0]),)
      csvfile.close()
    else:
      Files=Files_names
    self.my_file_range=parallel.rank_range(len(Files))
    Files=np.array(Files)[self.my_file_range]
    NumFiles=len(Files)
    self.contents=()
    for f in range(0,NumFiles):
        filename=strtools.replacestring(Files[f],' ','')
        self.contents+=(pd.read_csv(filename),)
