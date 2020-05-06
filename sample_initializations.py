#------------------------------------------------------------
# This program is a useful starting point for generating
# model input files. The function initialize() is used to
# convert arrays of meteorological data into a netcdf file
# formatted to be read as input into the model. A number of
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

from initialize import *


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
	jet_max = 75.0			# maximum jet speed (m/s)
	N02T = 1.3e-4			# tropopause static stability
	N02S = 4.0e-4			# stratosphere static stability
	jet_drop = 3000			# drop in tropopause height across jet (meters)
	tropopause_height = 31	# in grid points
	surface_theta = 290.0	# temperature at southern boundary (Kelvin)

	max_humid = 0.75		# lower-level relative humidity
	min_humid = 0.30		# upper-level relative humidity
	lower_z_humid = 1500.0	# set to max_humid below this level
	upper_z_humid = 5000.0	# set to min_humid above this level

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
	
		#print j,ubar[0,j,0]
	
	tbar[:,:,:] = 300.0
	
	#v[nx//2-3:nx//2+3,ny//2-3:ny//2+3,:] = 5
	
	#initialize(u,v,w,t,ubar,vbar,tbar,qbar,dx,dz,lat0,lon0)

	initialize(u,v,w,t,ubar,vbar,tbar,qbar,dx,dz,lat0,lon0,outfile,qv=q,qr=qr,qc=qc,qs=qs,qi=qi)


#------------------------------------------------------------
# Create a ball of rain for an idealized downburst experiment.
#------------------------------------------------------------
def rain_ball(outfile):

	lat0 = 30
	lon0 = 45

	nx = 200
	ny = 100
	nz = 45

	dx = 5000
	dz = 500

	ubar = np.zeros((nx,ny,nz))
	vbar = np.zeros((nx,ny,nz))
	tbar = np.zeros((nx,ny,nz))
	qbar = np.zeros((nx,ny,nz))
	
	for k in range(0,nz):
	
		tbar[:,:,k] = 300.0 + float(k*2)
	
	u = np.zeros((nx,ny,nz))
	v = np.zeros((nx,ny,nz))
	w = np.zeros((nx,ny,nz))
	t = np.zeros((nx,ny,nz))
	q = np.zeros((nx,ny,nz))	

	qr = np.zeros((nx,ny,nz))
	qc = np.zeros((nx,ny,nz))
	qs = np.zeros((nx,ny,nz))
	qi = np.zeros((nx,ny,nz))
	
	ri_begin = 30+111
	ri_end = 40+111
	rj_begin = 30
	rj_end = 40
	rk_begin = 10+5
	rk_end = 20+5
	r_max = 10.0e-3
		
	for i in range(ri_begin,ri_end):
		for j in range(rj_begin,rj_end):
			for k in range(rk_begin,rk_end):
	
				qr[i,j,k] = (0.5* r_max *
					cos_profile(-0.5,0.5,float(i-ri_begin) / (ri_end-ri_begin), 1)*
					cos_profile(-0.5,0.5,float(j-rj_begin) / (rj_end-rj_begin), 1)*
					cos_profile(-0.5,0.5,float(k-rk_begin) / (rk_end-rk_begin), 1)
					)
	
		
	rhbar = np.zeros((nx,ny,nz))
	rhbar[:,:,:] = 1.0
	
	initialize(u,v,w,t,ubar,vbar,tbar,qbar,dx,dz,lat0,lon0,outfile,qv=q,qr=qr,qc=qc,qs=qs,qi=qi)



#------------------------------------------------------------
# Set up a basic state for a squall line simulation
#------------------------------------------------------------
def thunderstorm_basicstate(outfile):

	#-------------------------------------------------
	# Environmental parameters
	#-------------------------------------------------
	tsurf = 300.0
	qsurf = 0.0161
	q4km = 0.001
	ztr = 12000.0
	temptr = 213.0
	ttr = 343.0
	lev0 = 4000.0
	lev1 = 8000.0

	#-------------------------------------------------
	# Grid parameters
	#-------------------------------------------------	
	lat0 = 30
	lon0 = 50

	nx = 100
	ny = 100

	dx = 5000
	dz = 100

	zu =  np.arange(-dz/2,24000.0,dz)
	nz = len(zu)

	#-------------------------------------------------
	# Meteorological fields
	#-------------------------------------------------
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
	
	#---------------------------------------------------------
	# Set up potential temperature profile
	#---------------------------------------------------------
	for k in range(1,nz):
	
		if zu[k] < ztr:		
			tbar[:,:,k] = tsurf + (ttr-tsurf) * (zu[k]/ztr)**1.25
		else:
			tbar[:,:,k] = ttr * np.exp(grav*(zu[k]-ztr)/(1004.0*temptr))
			
	tbar[:,:,0] = tbar[:,:,1]
	print(zu[:])
	print(tbar[5,5,:])
	
	#---------------------------------------------------------
	# Set up humidity profile
	#---------------------------------------------------------
	for k in range(1,nz):
		
		if zu[k] < lev0:
			
			qbar[:,:,k] = qsurf - (qsurf-q4km) * zu[k] / lev0
			
		elif zu[k] < lev1:
			
			qbar[:,:,k] = q4km - q4km * (zu[k] - lev0) / lev0
			
		else:
			
			qbar[:,:,k] = 0

	qbar[:,:,0] = qbar[:,:,1]


	#---------------------------------------------------------
	# Set up shear profile
	#---------------------------------------------------------
	shear = 20.0
	
	for k in range(0,nz):
		
		if zu[k] > ztr:
		
			ubar[:,:,k] = shear
		else:
			ubar[:,:,k] = -shear/2.0 + shear / ztr * zu[k]

	#---------------------------------------------------------
	# Set up thermal perturbation
	#---------------------------------------------------------
	ri_begin = 30
	ri_end = 40
	rj_begin = 30
	rj_end = 40
	rk_begin = 5
	rk_end = 55
	r_max = 2
	
	for i in range(ri_begin,ri_end):
		for j in range(rj_begin,rj_end):
			for k in range(rk_begin,rk_end):
	
				t[i,j,k] = (0.5* r_max *
					cos_profile(-0.5,0.5,float(i-ri_begin) / (ri_end-ri_begin), 1)*
					cos_profile(-0.5,0.5,float(j-rj_begin) / (rj_end-rj_begin), 1)*
					cos_profile(-0.5,0.5,float(k-rk_begin) / (rk_end-rk_begin), 1)
					)
	
	
	initialize(u,v,w,t,ubar,vbar,tbar,qbar,dx,dz,lat0,lon0,outfile,qv=q,qr=qr,qc=qc,qs=qs,qi=qi)


#------------------------------------------------------------
# Set up a basic state suitable for a hurricane simulation
# using the Jordan (1958) sounding
#------------------------------------------------------------
def jordan_sounding(outfile):

	zu = np.array([-100,0,132,583,1054,1547,2063,2609,3182,3792,4442,5138,5888,6703,7595,8581,9682,10935,12396,13236,14177,15260,16568,17883,19620,20743,22139,23971])
	tb = np.array([298,298,299,300,302,304,307,309,312,315,318,321,324,328,332,335,338,342,345,348,354,364,386,418,468,500,542,597])
	qb = np.array([18.2,18.2,17.6,15.3,13.0,11.0,8.4,7.1,5.8,4.6,3.6,3.2,2.1,1.4,0,0,0,0,0,0,0,0,0,0,0,0,0,0])
	pres = np.array([1015.1,1015.1,1000,950,900,850,800,750,700,650,600,550,500,450,400,350,300,250,200,175,150,125,100,80,60,50,40,30])
	Tc = np.array([26.3,26.3,26.0,23.0,19.8,17.3,14.6,11.8,8.6,5.1,1.4,-2.5,-6.9,-11.9,-17.7,-24.8,-33.2,-43.3,-55.2,-61.5,-67.6,-72.2,-73.5,-69.8,-63.9,-60.6,-57.3,-54.0])
	
	
	#print(len(tb),len(qb),len(zu))
	
	lat0 = -10
	lon0 = -100

	nx = 200
	ny = 150
	nz = len(zu)

	dx = 50000
	dz = 100

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
	
	qb *= 0.001

	tbar[:,:,:] = tb[:]
	qbar[:,:,:] = qb[:]
	
	ri_begin = 3
	ri_end = 5
	rj_begin = 3
	rj_end = 5
	rk_begin = 4
	rk_end = 10
	r_max = 5
	
	#for k in range(0,nz):
		#ubar[:,:,k] = float(k)
	
	#for k in range(0,6):
		
		#qbar[:,:,k] = (17 + (6-k)*2)*1.0e-3
	
	for i in range(ri_begin,ri_end):
		for j in range(rj_begin,rj_end):
			for k in range(rk_begin,rk_end):
	
				t[i,j,k] = (0.5* r_max *
					cos_profile(-0.5,0.5,float(i-ri_begin) / (ri_end-ri_begin), 1)*
					cos_profile(-0.5,0.5,float(j-rj_begin) / (rj_end-rj_begin), 1)*
					cos_profile(-0.5,0.5,float(k-rk_begin) / (rk_end-rk_begin), 1)
					)
	
	
			
	initialize(u,v,w,t,ubar,vbar,tbar,qbar,dx,dz,lat0,lon0,outfile,qv=q,qr=qr,qc=qc,qs=qs,qi=qi,zu=zu,surface_pressure=101510)



#------------------------------------------------------------
#
#------------------------------------------------------------
def main():

	try:
		#rain_ball("../model_input/rainball.nc")
		#thunderstorm_basicstate("../model_input/thunderstorm.nc")
		#jordan_sounding("../model_input/hurricane.nc")
		
		baroclinic_jet("../model_input/baroclinicjet.nc")
		#barotropic_jet("../model_input/barotropicjet.nc")
		
	except:

		if VERBOSE:
			traceback.print_exc(file=sys.stderr)
		else:
			print >>sys.stderr,err.msg

		return 1

if __name__ == "__main__":
    sys.exit(main())
