#------------------------------------------------------------
# This program is a useful starting point for generating
# model input files. The function initialize() is used to
# convert arrays of meteorological data into a netcdf file
# formatted to be read as input into the model. Several
# example scripts are provided for testing and as inspiration.
#
#
#------------------------------------------------------------

import traceback
import sys,os
import numpy as np
import math
from netCDF4 import Dataset

VERBOSE = 1


#------------------------------------------------------------
# Global constants
#------------------------------------------------------------
meters_in_degree = 111000.0
p0 = 100000.0
Rd = 287.0
cp = 1004.0
grav = 9.81
meters_per_latitude = meters_in_degree

#------------------------------------------------------------
# Create a zonally uniform baroclinic jet in thermal wind
# balance
#------------------------------------------------------------
def baroclinic_jet(outfile):

	#---------------------------------------------------
	# Control parameters for baroclinic jet
	#----------------------------------------------------
	lower_jet_lat = 35.0	# southern boundary of jet
	upper_jet_lat = 55.0	# northern boundary of jet
	jet_max = 50.0			# maximum jet speed (m/s)
	N02T = 1.3e-4			# tropopause static stability
	N02S = 4.0e-4			# stratosphere static stability
	jet_drop = 3000			# drop in tropopause height across jet (meters)
	tropopause_height = 31	# in grid points
	surface_theta = 280.0	# temperature at southern boundary (Kelvin)

	max_humid = 0.75		# lower-level relative humidity
	min_humid = 0.20		# upper-level relative humidity
	lower_z_humid = 1500.0	# set to max_humid below this level
	upper_z_humid = 4000.0	# set to min_humid above this level

	lat0 = 10				# lower left latitude
	lon0 = 50				# lower left longitude
	#---------------------------------------------------
	# Grid points
	#----------------------------------------------------
	nx = 200
	ny = 80
	nz = 45
	#---------------------------------------------------
	# Grid spacing
	#----------------------------------------------------
	dx = 100000
	dz = 500
	zheight =  np.arange(-250.0,24000.0,dz)
	
	#---------------------------------------------------------
	#----------------------------------------------------------
	# 				CONSTRUCT THE BASIC STATE
	#---------------------------------------------------------
	#----------------------------------------------------------
	lats = np.arange(lat0, lat0 + ny * dx / meters_in_degree, dx / meters_in_degree)
	lons = np.arange(lon0, lon0 + nx * dx / meters_in_degree, dx / meters_in_degree)
	
	#---------------------------------------------------
	# Meteorological fields to initialize
	#----------------------------------------------------
	ubar = np.zeros((nx,ny,nz))
	vbar = np.zeros((nx,ny,nz))
	tbar = np.zeros((nx,ny,nz))
	qbar = np.zeros((nx,ny,nz))
	
	u = np.zeros((nx,ny,nz))
	v = np.zeros((nx,ny,nz))
	w = np.zeros((nx,ny,nz))
	t = np.zeros((nx,ny,nz))
	q = np.zeros((nx,ny,nz))
	
	#---------------------------------------------------
	# Determine height of jet at all latitude points
	#----------------------------------------------------
	jet_height = np.zeros((ny))
		
	jet_height0 = int(get_gridpoint_from_coord(lower_jet_lat,lat0,dx))
	jet_height1 = int(get_gridpoint_from_coord(upper_jet_lat,lat0,dx))
	
	for j in range(0,jet_height0):
		
		jet_height[j] = zheight[tropopause_height]
		
	for j in range(jet_height1,ny):
		
		jet_height[j] = zheight[tropopause_height] - jet_drop

	for j in range(jet_height0,jet_height1):
	
		jet_height[j] = zheight[tropopause_height] - jet_drop * float(j-jet_height0) / float(jet_height1-jet_height0)
	
	#---------------------------------------------------
	# Create the zonal wind profile
	#----------------------------------------------------
	lowlat = int(get_gridpoint_from_coord(lower_jet_lat,lat0,dx))
	higlat = int(get_gridpoint_from_coord(upper_jet_lat,lat0,dx))
	midlat = (higlat+lowlat) // 2
	
	for j in range(lowlat,higlat):
		for k in range(0,nz):
		
			if zheight[k] < jet_height[j]:
				
				ubar[:,j,k] = (zheight[k]/jet_height[j]) * jet_max * 0.5 * cos_profile(-0.5,0.5,float(j-lowlat)/float(higlat-lowlat),1)
			else:
				ubar[:,j,k] = (jet_height[j]/zheight[k]) * jet_max * 0.5 * cos_profile(-0.5,0.5,float(j-lowlat)/float(higlat-lowlat),1)
		
	#---------------------------------------------------
	# Create base state temperature profile
	#----------------------------------------------------
	tb = np.zeros((nz))
	
	tb[0] = surface_theta
	tb[1] = tb[0];

	for k in range(2,tropopause_height):
		
		tb[k] = N02T * tb[k-1] / grav * (zheight[k]-zheight[k-1]) + tb[k-1]
		
	for k in range(tropopause_height,nz-1):
		
		tb[k] = N02S * tb[k-1] / grav * (zheight[k]-zheight[k-1]) + tb[k-1]
		
	tb[nz-1] = tb[nz-2];

	#---------------------------------------------------
	# Create basic state relative humidity profile
	#----------------------------------------------------
	rh_profile = np.zeros((nz))

	for k in range(0,nz):

		if zheight[k] < lower_z_humid:
			
			rh_profile[k] = max_humid
			
		if zheight[k] < upper_z_humid and zheight[k] >= lower_z_humid:
			
			rh_profile[k] = max_humid - (max_humid-min_humid) * (zheight[k]-lower_z_humid) / (upper_z_humid-lower_z_humid)
			
		if zheight[k] >= upper_z_humid:
			
			rh_profile[k] = min_humid
	
	#---------------------------------------------------
	# Create basic state potential temperature field
	#----------------------------------------------------
	fc = 2 * 7.292e-5 * np.sin( lats * np.pi / 180. )
	
	pbar = np.zeros((nx,ny,nz))

	for j in range(1,ny):
	
		pbar[:,j,:] = pbar[:,j-1,:] - 0.5 * (ubar[:,j-1,:]+ubar[:,j,:]) * dx * fc[j]
	
	for k in range(1,nz):

		tbar[:,:,k] = tb[k] * (pbar[:,:,k] - pbar[:,:,k-1]) / (grav * (zheight[k]-zheight[k-1]) ) + tb[k]
	
	tbar[:,:,0] = tbar[:,:,1]
	tbar[:,:,nz-1] = tbar[:,:,nz-2]

	#---------------------------------------------------
	# Create basic state mixing ratio field
	#----------------------------------------------------
	rh3d = np.zeros((nx,ny,nz))
	
	rh3d[:,:,:] = rh_profile[:]
	

	qbar = get_basic_state_mixing_ratio(tbar,pbar,rh3d,dz)

	initialize(u,v,w,t,ubar,vbar,tbar,qbar,dx,dz,lat0,lon0,outfile,qv=q)


#------------------------------------------------------------
# Create an idealized dry, barotropic jet.
#------------------------------------------------------------
def barotropic_jet(outfile):

	lat0 = 10
	lon0 = 50

	nx = 200
	ny = 80
	nz = 45

	dx = 100000
	dz = 500

	ubar = np.zeros((nx,ny,nz))
	vbar = np.zeros((nx,ny,nz))
	tbar = np.zeros((nx,ny,nz))
	qbar = np.zeros((nx,ny,nz))
	
	u = np.zeros((nx,ny,nz))
	v = np.zeros((nx,ny,nz))
	w = np.zeros((nx,ny,nz))
	t = np.zeros((nx,ny,nz))
	q = np.zeros((nx,ny,nz))	

	qr = np.zeros((nx,ny,nz))
	qc = np.zeros((nx,ny,nz))
	qs = np.zeros((nx,ny,nz))
	qi = np.zeros((nx,ny,nz))
	
	jl = ny // 2 - 5 		#int(ny)/4
	jh = ny // 2 + 5	#3*int(ny)/4
	
	for j in range(jl,jh+1):
	
		ubar[:,j,:] = -60 * (np.sin(3.14 * (2*float(j-jl)/float(jh-jl) - 0.5) ) + 1) / 2.0
	
		print j,ubar[0,j,0]
	
	tbar[:,:,:] = 300.0
	
	#v[nx//2-3:nx//2+3,ny//2-3:ny//2+3,:] = 5
	
	#initialize(u,v,w,t,ubar,vbar,tbar,qbar,dx,dz,lat0,lon0)

	initialize(u,v,w,t,ubar,vbar,tbar,qbar,dx,dz,lat0,lon0,outfile,qv=q,qr=qr,qc=qc,qs=qs,qi=qi)

#------------------------------------------------------------
#
# Perturbation values set to 'None' will be excluded from output file
#
#
# u,v,w,t -> perturbation wind components (u,v,w) and potential temperature (t)
# qv	-> perturbation water vapor mixing ratio (optional)
# qc	-> cloud water mixing ratio (optional)
# qr	-> rain water mixing ratio (optional)
# qs	-> snow water mixing ratio (optional)
# qi	-> ice water mixing ratio (optional)
# ubar,vbar -> basic state wind components
# tbar	-> basic state potential temperature
# pbar	-> basic state pressure
# qbar	-> basic state water vapor mixing ratio
# dx	-> horizontal grid spacing
# dz	-> vertical grid spacing
# corner_lat -> latitude of lower-left corner
# corner_lon -> longitude of lower-left corner
#
#
#------------------------------------------------------------
def initialize(u,v,w,t,ubar,vbar,tbar,qbar,dx,dz,corner_lat,corner_lon,filename,qv=None,qr=None,qc=None,qs=None,qi=None,rh=None):

	nx,ny,nz = u.shape
	
	#--------------------------------------------------------
	# Create coordinate arrays
	#--------------------------------------------------------
	levs = np.arange(-dz,(nz-1)*dz,dz)
	lats = np.arange(corner_lat,float(dx*ny)/meters_in_degree+corner_lat,float(dx)/meters_in_degree)
	lons = np.arange(corner_lon,float(dx*nx)/meters_in_degree+corner_lon,float(dx)/meters_in_degree)

	#--------------------------------------------------------
	# Other variables
	#--------------------------------------------------------		
	p = np.zeros((nx,ny,nz))		# perturbation pressure
	wbar = np.zeros((nx,ny,nz))		# basic state vertical velocity
	pbar = np.zeros((nx,ny,nz))		# basic state pressure

	#--------------------------------------------------------
	# Create vertically varying base state arrays
	#--------------------------------------------------------	
	tb = np.zeros((nz))
	qb = np.zeros((nz))
	pb = np.zeros((nz))
	
	tb[:] = tbar[nx//2,ny//2,:]
	qb[:] = qbar[nx//2,ny//2,:]
	pb = get_base_state_pressure(tb,qb,100000.0,dz)
	#print(qb)
	#--------------------------------------------------------
	# Remove 1d base state from 3d basic state
	#--------------------------------------------------------
	tbar[:,:,:] = tbar[:,:,:] - tb[:]
	
	tbv = tb*(1.0+0.61*qb)
	
	pbar[:,:,:] = hydrostatic_balance(tbar,qbar,tbv,dz) + pb[:]

	qbar -= qb

	#--------------------------------------------------------
	# Stagger variables onto grid
	#--------------------------------------------------------	
	u,v,w = stagger_winds(u,v,w)
	ubar,vbar = stagger_winds(ubar,vbar)
	
	wbar = basic_state_vertical_velocity(ubar,vbar,get_density(tb,pb,qb),dx,dx,dz)
	
	create_file(filename,dx,dx,dz,u,v,w,t,p,qv,qr,qc,qs,qi,ubar,vbar,wbar,tbar,pbar,qbar,tb,qb,pb,lons,lats,levs)
	
	
#------------------------------------------------------------
#
#------------------------------------------------------------
def create_file(fileName,dx,dy,dz,u,v,w,t,p,qv,qr,qc,qs,qi,ubar,vbar,wbar,tbar,pbar,qbar,tb,qb,pb,lons0,lats0,levs0):

	nx,ny,nz = u.shape

	dataset = Dataset(fileName,'w',format='NETCDF4_CLASSIC')

	dataset.createDimension('x',nx)
	dataset.createDimension('y',ny)
	dataset.createDimension('z',nz)
	dataset.createDimension('t',None)
	
	dataset.setncattr("dx",float(dx))
	dataset.setncattr("dy",float(dy))
	dataset.setncattr("dz",float(dz))

	dataset.setncattr("latoffset",lats0[0])
	dataset.setncattr("lonoffset",lons0[0])
	

	make_and_set_pert_variable("u-wind",u,dataset)
	make_and_set_pert_variable("v-wind",v,dataset)
	make_and_set_pert_variable("w-wind",w,dataset)
	make_and_set_pert_variable("theta",t,dataset)
	make_and_set_pert_variable("pi",p,dataset)	

	make_and_set_pert_variable("qv",qv,dataset)
	make_and_set_pert_variable("qc",qc,dataset)
	make_and_set_pert_variable("qr",qr,dataset)
	make_and_set_pert_variable("qi",qi,dataset)
	make_and_set_pert_variable("qs",qs,dataset)

	ubf = dataset.createVariable("ubar", np.float32,('x','y','z'))
	vbf = dataset.createVariable("vbar", np.float32,('x','y','z'))
	wbf = dataset.createVariable("wbar", np.float32,('x','y','z'))
	tbf = dataset.createVariable("thbar",np.float32,('x','y','z'))
	pbf = dataset.createVariable("pbar", np.float32,('x','y','z'))
	qbf = dataset.createVariable("qbar", np.float32,('x','y','z'))	

	t0f = dataset.createVariable("tb",  np.float32,('z'))	
	q0f = dataset.createVariable("qb",  np.float32,('z'))
	p0f = dataset.createVariable("pib", np.float32,('z'))
	
	ubf[:,:,:] = ubar[:,:,:]
	vbf[:,:,:] = vbar[:,:,:]
	wbf[:,:,:] = wbar[:,:,:]
	tbf[:,:,:] = tbar[:,:,:]
	pbf[:,:,:] = pbar[:,:,:]		
	qbf[:,:,:] = qbar[:,:,:]
	
	t0f[:] = tb[:]
	q0f[:] = qb[:]
	p0f[:] = pb[:]
	#print(uf[:,:,:,0].shape,u.shape)
	

	levs = dataset.createVariable('zu',  np.float32,('z',))
	lats = dataset.createVariable('lat', np.float32,('y',))
	lons = dataset.createVariable('lon', np.float32,('x',)) 
	time = dataset.createVariable('time', np.float32,('t')) 

	#print(len(lats0),len(lons0),len(levs0))

	lats[:] = lats0
	lons[:] = lons0
	levs[:] = levs0
	time[0] = 0

	dataset.close()

#------------------------------------------------------------
#
#------------------------------------------------------------	
def make_and_set_pert_variable(name,values,dataset):
	
	if values is not None:
	
		varf = dataset.createVariable(name,np.float32,('t','x','y','z'))
		
		varf[0,:,:,:] = values[:,:,:]


#------------------------------------------------------------
# Stagger wind components
#
# u -> zonal wind
# v -> meridional wind
# w -> vertical wind (optional argument)
#------------------------------------------------------------
def stagger_winds(u,v,w=None):
	
	nx,ny,nz = u.shape

	#---------------------------------------------
	# Zonal and meridional components
	#----------------------------------------------	
	us = np.zeros((nx,ny,nz))
	vs = np.zeros((nx,ny,nz))
	
	us[0,:,:] = u[0,:,:]
	vs[:,0,:] = v[:,0,:]
	
	for i in range(1,nx):	
		us[i,:,:] = 0.5*(u[i-1,:,:] + u[i,:,:])

	for j in range(1,ny):	
		vs[:,j,:] = 0.5*(v[:,j-1,:] + v[:,j,:])

	if w is None:
		return u,v

	#---------------------------------------------
	# Vertical component, if argument is not null
	#----------------------------------------------
	ws = np.zeros((nx,ny,nz))

	ws[:,:,0] = w[:,:,0]

	for k in range(1,nz):	
		ws[:,:,k] = 0.5*(w[:,:,k] + w[:,:,k])
		
	return u,v,w

#------------------------------------------------------------
#
#------------------------------------------------------------
def get_density(tb,pib,qb):

	return 100000.0*pib[:]**(717.0/287.0) / (287.0*(tb[:]*(1.0+0.61*qb[:])))

#------------------------------------------------------------
# Create vertical velocity field that satisfies the 
# anelastic continuity field using horizontal wind field
#------------------------------------------------------------
def basic_state_vertical_velocity(ubar,vbar,rhou,dx,dy,dz):
	
	nx,ny,nz = ubar.shape
	
	rhow = np.zeros((nz))
	
	for k in range(1,nz):
		rhow[k] = 0.5*(rhou[k-1]+rhou[k])
	
	rhow[0] = rhow[1]
	
	wbar = np.zeros((nx,ny,nz))

	wbar[:,:,0] = 0		# boundary condition = zero vertical 
	wbar[:,:,1] = 0		# velocity at bottom

	udiv = np.zeros((nx,ny,nz))
	vdiv = np.zeros((nx,ny,nz))

	for i in range(0,nx-1):
		udiv[i,:,:] = (ubar[i+1,:,:] - ubar[i,:,:])/dx
		
	for j in range(1,ny-1):
		vdiv[:,j,:] = (vbar[:,j+1,:] - vbar[:,j,:])/dy

	for k in range(2,nz-2):

		wbar[:,:,k] = (
				(rhow[k-1]/rhow[k]) * wbar[:,:,k-1] 
				+ (rhou[k-1]/rhow[k])*
				(
					-udiv[:,:,k-1] - vdiv[:,:,k-1]
				)*dz
			)
	

	wbar[i,j,nz-1] = wbar[i,j,nz-2]	# velocity at top

	return wbar

#------------------------------------------------------------
# Get the grid point from a given grid point
#------------------------------------------------------------
def get_gridpoint_from_coord(lat,startLat,dx):
	
	return int( (lat - startLat) * meters_in_degree / dx)

#------------------------------------------------------------
# Create a profile of non-dimensional pressure from
# potential temperature and mixing ratio, assuming
# hydrostatic balance
#
#------------------------------------------------------------
def hydrostatic_balance(tbar,qbar,tbv,dz):
	
	nx,ny,nz = tbar.shape
	
	pbar = np.zeros((nx,ny,nz))


	pbar[:,:,nz-1] = 0
	pbar[:,:,nz-2] = -0.5*(grav/cp)*( (tbar[:,:,nz-2]*(1.0+0.61*qbar[:,:,nz-2]))/(tbv[nz-2]*tbv[nz-2]) )*dz

	for k in range(nz-3,0,-1):

		tup = (tbar[:,:,k+1]*(1.+0.61*qbar[:,:,k+1]))/(tbv[k+1]*tbv[k+1])
		tdn = (tbar[:,:,k  ]*(1.+0.61*qbar[:,:,k  ]))/(tbv[k  ]*tbv[k  ])
		
		pbar[:,:,k] = pbar[:,:,k+1] - 0.5*(grav/cp) * (tup+tdn) * dz
		
	return pbar


#------------------------------------------------------------
# Create a profile of non-dimensional pressure from
# potential temperature and mixing ratio, assuming
# hydrostatic balance
#
#------------------------------------------------------------
def get_base_state_pressure(tb,qb,pressfc,dz):

	if len(tb) != len(qb):
		print("Error: array lengths of first two arguments must be equal!")
		return False

	nz = len(tb)
	
	pb = np.zeros((nz))

	# virtual potential temperature
	tbv = tb*(1.0+0.61*qb)

	# surface dimensional pressure
	pisfc = (pressfc/p0)**(Rd/cp)
	
	#------------------------------------------------------
	# Create hydrostatically balanced pressure field
	#------------------------------------------------------	
	pb[1] = pisfc - grav * 0.5 * dz / (cp * tbv[1]);
	pb[0] = pisfc + grav * 0.5 * dz / (cp * tbv[0]);

	for k in range(2,nz):

		pb[k] = pb[k-1] - grav * dz / (cp*(0.5*(tbv[k]+tbv[k-1]) ) )
	
	return pb


#------------------------------------------------------------
# Calculate virtual potential temperature
#
#------------------------------------------------------------
def virtual_potential_temp(tb,qb):
	
	return tb*(1.0+0.61*qb)

#------------------------------------------------------------
# Calculate saturation mixing ratio.
#
# temperature - full, actual temperature field (K, not potential, not perturation)
# pressure - full, dimensional pressure (Pa)
#------------------------------------------------------------
def get_QV_Sat(temperature,pressure):

	esl = 611.2 * np.exp(17.67 * (temperature-273.15) / (temperature - 29.65) )

	return 0.62197 * esl / (pressure-esl)

#------------------------------------------------------------
# Get the mixing ratio field from the potential temperature
# and relative humidity field
#
#------------------------------------------------------------
def get_basic_state_mixing_ratio(tbar,pbar,rh,dz):

	nx,ny,nz = tbar.shape
	
	qbar = np.zeros((nx,ny,nz))

	rhz = rh[nx//2,ny//2,:]
	tb = tbar[nx//2,ny//2,:]
	qb = np.zeros((nz))

	pb = get_base_state_pressure(tb,qb,100000.0,dz)

	for k in range(0,nz):
		
		temp = tb[k] * pb[k]
		pres = p0*pow(pb[k],(cp/Rd))
		
		qb[k] = get_QV_Sat(temp,pres) * rhz[k]
	
	pb = get_base_state_pressure(tb,qb,100000.0,dz)

	tbv = virtual_potential_temp(tb,qb)

	#-------------------------------------------------
	# Convert to non-dimensional pressure
	#--------------------------------------------------
	pbar[:,:,:] = pbar[:,:,:] / (cp*tbv[:]) + pb[:]

	#print(pbar[40,40,:])

	temp = tbar[:,:,:] * pbar[:,:,:]			# full, actual temperature
	pres = p0 * pbar[:,:,:]**(cp/Rd)			# full, dimensional pressure

	#print(pres[40,:,20])

	smr = get_QV_Sat(temp,pres)				# calculate saturation mixing ratio

	qbar[:,:,:] = rh[:,:,:] * smr[:,:,:]

	#print(qbar[40,4,:]*1000)

	qbar[:,:,0] = qbar[:,:,1]
	qbar[:,:,nz-1] = qbar[:,:,nz-2]
	
	
	return qbar
	
#*********************************************************************
# 
#
# @param pl - start of interval in periods
# @param ph - end of interval in periods
# @param frac - fractional distance between pl and ph
# @param shift - shift entire function up or down in units of amplitude
# @return - a number from -one to one
#**********************************************************************
def cos_profile(pl,ph,frac,shift):
	
	if frac >= 0 and frac <= 1.0:
	
		return np.cos(  pl * 2*np.pi + (ph-pl) * 2*np.pi * frac ) + shift;
		
	else:
		return 0;
	



#------------------------------------------------------------
#
#------------------------------------------------------------
def main():

	try:		
		baroclinic_jet("../model_input/baroclinicjet.nc")
		#barotropic_jet("barotropicjet.nc")

	except:

		if VERBOSE:
			traceback.print_exc(file=sys.stderr)
		else:
			print >>sys.stderr,err.msg

		return 1

if __name__ == "__main__":
    sys.exit(main())
