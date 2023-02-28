#------------------------------------------------------------
# This program plots model output. You may have to change
# the paths of the input and output files at the bottom
# of the file under "main()".
#
# run as : python plot_model_output.py start_time end_time
#------------------------------------------------------------

import traceback
import sys,os
import numpy as np

import matplotlib.pyplot as plt
import math

from netCDF4 import Dataset
from mpl_toolkits.basemap import Basemap
from matplotlib.colors import LogNorm
from matplotlib.ticker import LogFormatter


VERBOSE = 1
	
#------------------------------------------------------------
# Plot wind and temperature
#------------------------------------------------------------
def plot_winds_temperature_pressure(infile,outfile,time):

	height = 10
	dataset = Dataset(infile,mode='r')

	#-----------------------------------
	# Get grid data
	#-----------------------------------
	nx,ny,nz,nt = get_netcdf_dims(dataset)
	levs,lons,lats,times = get_netcdf_coords(dataset)

	if time >= nt:
	    return

	print('time = '+str(time))

	#-----------------------------------
	# Get meteorological data.
	#-----------------------------------	
	u = get_var(dataset,"u-wind",time,height)
	v = get_var(dataset,"v-wind",time,height)
	w = get_var(dataset,"w-wind",time,height)
	p = get_var(dataset,"pi",time,height)
	t = get_var(dataset,"theta",time,height)

	ubar = get_var(dataset,"ubar",-1,height)
	vbar = get_var(dataset,"vbar",-1,height)
	tbar = get_var(dataset,"thbar",-1,height)
	pbar = get_var(dataset,"pbar",-1,height)

	pib = dataset.variables['pib'][:]
	tb = dataset.variables['tb'][:]
	qb = dataset.variables['qb'][:]
	
	rho = get_density(tb,pib,qb)

	temp = t + tbar + tb[height]
	u = u #+ ubar
	v = v #+ vbar
	
	#-----------------------------------
	# Destagger winds.
	#-----------------------------------
	u,v = destagger(u,v,1)

	clevs = np.array([-7,-6,-5,-4,-3,-2,-1,1,2,3,4,5,6,7])

	# create Basemap instance
	m = Basemap(projection='cyl',llcrnrlon=lons[0],llcrnrlat=lats[0],urcrnrlon = lons[nx-1],urcrnrlat = lats[ny-1])

	# compute map projection coordinates for lat/lon grid.
	x, y = m(*np.meshgrid(lons,lats))
	#print(temp.shape)
	#-----------------------------------
	# make filled contour plot
	#-----------------------------------
	cs = m.contourf(x,y,temp,cmap=plt.cm.bwr)
	cb = plt.colorbar(cs, shrink=0.5, extend='both')
	
	#-----------------------------------
	# make line contour plot
	#-----------------------------------	
	cs = m.contour(x,y,p,linewidths=1,colors='k')
	
	#-----------------------------------
	# make vector plot
	#-----------------------------------
	stride = 5

	quiv = m.quiver(x[::stride,::stride], y[::stride,::stride], u[::stride,::stride], v[::stride,::stride],pivot='mid',color='Black',scale=400)

	#-----------------------------------
	# Create map
	#-----------------------------------	
	#m.drawcoastlines(linewidth=1) # draw coastlines
	pd1 = m.drawparallels(np.arange(-90.,120.,15.),labels=[1,0,0,0]) # draw parallels
	pd2 = m.drawmeridians(np.arange(0.,420.,15.),labels=[0,0,0,1]) # draw meridians

	remove_gridlines(pd1,pd2)

	#-----------------------------------
	# Title
	#-----------------------------------
	plt.title(r'Pot. Temp. (shading, K, '
	+str(levs[height])[0:3]+
	r' km), $\vec{v}$ (vectors, '
	+str(levs[height])[0:3]+
	r' km), and Pert. Pres. (contours, '
	+str(levs[height])[0:3]+r' km) at '
	+str(times[time])+' hr')

	#-----------------------------------
	# Create figure
	#-----------------------------------
	figure = plt.gcf() # get current figure
	figure.set_size_inches(14, 10)

	plt.savefig(outfile,bbox_inches='tight')

	plt.close()
	dataset.close()


#------------------------------------------------------------
# Plot vertically integrated condensate and winds
#------------------------------------------------------------
def plot_winds_condensate(infile,outfile,time):

	height = 14
	
	sum1 = 10	# lower level of integration
	sum2 = 35	# upper level of integration

	dataset = Dataset(infile,mode='r')

	nx,ny,nz,nt = get_netcdf_dims(dataset)

	if time >= nt:
	    return

	print('time = '+str(time))

	levs,lons,lats,times = get_netcdf_coords(dataset)

	times = dataset.variables['time']

	u = get_var(dataset,"u-wind",time,height) #+ get_var(dataset,"ubar",-1,height)
	v = get_var(dataset,"v-wind",time,height) #+ get_var(dataset,"vbar",-1,height)
	p = get_var(dataset,"pi",time,height)

	ubar = get_var(dataset,"ubar",-1,height)
	vbar = get_var(dataset,"vbar",-1,height)

	e = dataset.variables['qi'][time,:,:,sum1:sum2]
	s = dataset.variables['qs'][time,:,:,sum1:sum2]
	c = dataset.variables['qc'][time,:,:,sum1:sum2]
	r = dataset.variables['qr'][time,:,:,sum1:sum2]

	pib = dataset.variables['pib'][:]
	tb = dataset.variables['tb'][:]
	qb = dataset.variables['qb'][:]

	pressure = 100000*pib**(1004.0/287.0)

	rho = get_density(tb,pib,qb)

	c2 = mass_integral(r+c+e+s,rho[sum1:],levs[sum1-1:])

	c = np.reshape(np.ravel(c2), (ny,nx), order='F')

	for i in range(0,nx):
		for j in range(0,ny):

			if c[j,i] <= 0:
				c[j,i] = 1e-10

	u,v = destagger(u,v,1)

	#-----------------------------------
	# Mask small vectors
	#-----------------------------------
	vs = np.sqrt(u*u + v*v)
	
	u = np.ma.masked_where(vs < 0.5, u)
	v = np.ma.masked_where(vs < 0.5, v)

	clevs = np.array([-5,-4,-3,-2,-1,1,2,3,4,5])
	clevs2 = np.array([
		0.03125/1024,
		0.03125/512,
		0.03125/256,0.03125/128,
		0.03125/64,0.03125/32,0.03125/16,0.03125/8,
		0.03125/4,0.03125/2,0.03125,0.0625,0.125,0.25,0.5,1,2
	])


	# create Basemap instance
	m = Basemap(projection='cyl',llcrnrlon=lons[0],llcrnrlat=lats[0],urcrnrlon = lons[nx-1],urcrnrlat = lats[ny-1]) 

	x, y = m(*np.meshgrid(lons,lats))

	formatter = LogFormatter(2, labelOnlyBase=False) 

	#-----------------------------------
	# make filled contour plot
	#-----------------------------------
	cs = m.contourf(x,y,c*1000,clevs2,norm=LogNorm(vmin=clevs2[0],vmax=clevs2[len(clevs2)-1]),cmap=plt.cm.gist_ncar)
	cb = plt.colorbar(cs, ticks=clevs2, shrink=0.3, extend='both',format=formatter)

	#-----------------------------------
	# make line contour plot
	#-----------------------------------	
	cs3 = m.contour(x,y,p*0.01,clevs,linewidths=1,colors='k')

	#-----------------------------------
	# make vector plot
	#-----------------------------------
	stride = 8

	quiv = m.quiver(x[::stride,::stride], y[::stride,::stride], u[::stride,::stride], v[::stride,::stride],pivot='mid',color='Black',scale=200)

	#-----------------------------------
	# Draw map
	#-----------------------------------
	#m.drawcoastlines(linewidth=1) # draw coastlines
	pd1 = m.drawparallels(np.arange(-90.,120.,15),labels=[1,0,0,0]) # draw parallels
	pd2 = m.drawmeridians(np.arange(0.,420.,15),labels=[0,0,0,5]) # draw meridians

	remove_gridlines(pd1,pd2)

	#-----------------------------------
	# Title
	#-----------------------------------
	plt.title('Integrated Condensate (shading,'+str(levs[sum1])[0:4]+' to '+
	str(levs[sum2])[0:4]+' km, K),\nPerturbation Winds and Pressure ('
	+str(levs[height])[0:4]+' km ' + str(np.around(pressure[height]*0.01,decimals=0))
	+' hPa) at t = '+str(times[time])[0:4]+' hr')

	#-----------------------------------
	# Create figure
	#-----------------------------------
	figure = plt.gcf() # get current figure
	figure.set_size_inches(14, 10)

	plt.savefig(outfile,bbox_inches='tight')

	plt.close()
	dataset.close()

#------------------------------------------------------------
#
#------------------------------------------------------------
def cross_section(infile,outfile,time):

    west = 0
    east = 500
    north = 95#85#48#75
    top = 42

    dataset = Dataset(infile,mode='r')

    lons = dataset.variables['lon'][west:east]
    lats = dataset.variables['lat'][:]
    times = dataset.variables['time'][:]


    if time >= len(times):
        return

    levs = dataset.variables["zu"][:]

    levs /= 1000

    #north = find_index(31.6,lats)
    #north = find_index(32.25,lats)
    #mid = find_index(51.55,lons)

    north = 95#47#95

    #print(north)

    nx = len(lons)
    nz = len(levs)
    top = nz

#   levs = np.arange(-0.25,nz/2,0.5)

    has_microphysics = False

    tb = dataset.variables['tb'][:]
    pib = dataset.variables['pib'][:]
    qb = dataset.variables['qb'][:]

    pres = 100000.0*pib**(1004.0/287.0)


    print('time = '+str(time))
    #print(pres)

    #print(lats[north])

    
    u2 = dataset.variables['u-wind'][time,west:east,north,0:top] #+ dataset.variables['ubar'][west:east,north,0:top]
    v2 = dataset.variables['v-wind'][time,west:east,north,0:top]   
    p2 = dataset.variables['pi'][    time,west:east,north,0:top]
    t2 = dataset.variables['theta'][ time,west:east,north,0:top]
    w2 = dataset.variables['w-wind'][time,west:east,north,0:top]
    
    u = np.reshape(np.ravel(u2), (nz,nx), order='F')
    w = np.reshape(np.ravel(w2), (nz,nx), order='F')
    v = np.reshape(np.ravel(v2), (nz,nx), order='F')
    p = np.reshape(np.ravel(p2), (nz,nx), order='F')
    t = np.reshape(np.ravel(t2), (nz,nx), order='F')


    if has_microphysics:

        q2 = dataset.variables['qv'][time,west:east,north,0:top]
        c2 = dataset.variables['qc'][time,west:east,north,0:top]
        r2 = dataset.variables['qr'][time,west:east,north,0:top]
        s2 = dataset.variables['qr'][time,west:east,north,0:top]
        e2 = dataset.variables['qr'][time,west:east,north,0:top]

        q = np.reshape(np.ravel(q2), (nz,nx), order='F')
        c = np.reshape(np.ravel(c2), (nz,nx), order='F')
        r = np.reshape(np.ravel(r2), (nz,nx), order='F')
        s = np.reshape(np.ravel(s2), (nz,nx), order='F')
        e = np.reshape(np.ravel(e2), (nz,nx), order='F')

        c = 1000 * (c)
        s = 1000 * (s)
        r = 1000 * (r)
        e = 1000 * (e)

        for i in range(0,nx):
            for j in range(0,nz):

                if c[j,i] <= 0:
                    c[j,i] = 1e-10

                if s[j,i] <= 0:
                    s[j,i] = 1e-10

                if r[j,i] <= 0:
                    r[j,i] = 1e-10


    tmin = dataset.variables['theta'][time,:,:,:]

    ind = np.unravel_index(np.argmax(tmin, axis=None), tmin.shape)
    print(ind)
    print(tmin[ind])


    clevs = np.array([-10,-9,-8,-7,-6,-5,-4,-3,-2,-1,1,2,3,4,5,6,7,8,9,10])

    clevs2 = np.array([0.03125/8.0,0.03125/4.0,0.03125/2.0,0.03125,0.0625,0.125,0.25,0.5,1,2,4,8,16,32])
    clevs3 = np.array([-32,-16,-8,-4,-2,-1,-0.5,-0.25,-0.125,-0.0625,-0.03125,-0.03125/2.0,-0.03125/4,-0.03125/8.0,
                        0.03125/8.0,0.03125/4.0,0.03125/2.0,0.03125,0.0625,0.125,0.25,0.5,1,2,4,8,16,32])

    x, z = np.meshgrid(lons,levs)

    formatter = LogFormatter(2, labelOnlyBase=False) 
    formatter2 = LogFormatter(2, labelOnlyBase=False)


    #-----------------------------------
    # make filled contour plot
    #cs = plt.contourf(x,z,r+c,clevs6,norm=LogNorm(vmin=clevs6[0],vmax=clevs6[len(clevs6)-1]),
    #cmap=plt.cm.gist_ncar)
    #cs = plt.contourf(x,z,r+c,clevs6,norm=LogNorm(vmin=clevs6[0],vmax=clevs6[len(clevs6)-1]),
    #cmap=plt.cm.gist_ncar)
    #cb = plt.colorbar(cs, shrink=0.8, extend='both',format=formatter)

    #cs = plt.contour(x,z,s,clevs6,norm=LogNorm(vmin=clevs6[0],vmax=clevs6[len(clevs6)-1]),
    #cmap=plt.cm.gist_ncar)
    #cb = plt.colorbar(cs, shrink=0.8, extend='both',format=formatter)
    #cs2 = plt.contour(x,z,r,clevs6,norm=LogNorm(vmin=clevs6[0],vmax=clevs6[len(clevs6)-1]),
    #cmap=plt.cm.gist_ncar) 
    #-----------------------------------
    # make line contour plot


    cs2 = plt.contourf(x,z,t,clevs*0.5,cmap=plt.cm.bwr)
    cb = plt.colorbar(cs2, shrink=0.8, extend='both')

    xmin, xmax = plt.xlim()
    ymin, ymax = plt.ylim()

    plt.xlim(45.05,45.3)
    plt.ylim(0,10)

    xstride = 2
    zstride = 2

    quiv = plt.quiver(x[::xstride,::zstride], z[::xstride,::zstride], u[::xstride,::zstride], w[::xstride,::zstride],pivot='mid',color='Black',scale=400)

    index = int(infile[-4])
    plt.title('Rain Water Mixing Ratio (shading, g/kg) at '+
    str(np.round(times[time]*60,decimals=1))+
    " min\n\n")# \nPert. Pot. Temp. (contours, K) at '+str(lats[north])+' N' ) 
    plt.xlabel('Longitude')
    plt.ylabel('Height (km)')

    figure = plt.gcf()
    figure.set_size_inches(12, 6)

    plt.savefig(outfile,bbox_inches='tight')


#   plt.show()
    plt.close()
    dataset.close()

#------------------------------------------------------------
# Get the dimension of a netcdf file.
#------------------------------------------------------------
def get_netcdf_dims(dataset):
	
	nx = dataset.dimensions['x'].size
	ny = dataset.dimensions['y'].size
	nz = dataset.dimensions['z'].size
	nt = dataset.dimensions['t'].size
	
	return nx,ny,nz,nt
	
	
#------------------------------------------------------------
# Get the coordinates of a netcdf file.
#------------------------------------------------------------
def get_netcdf_coords(dataset):
	
	levs = 1.0e-3 * dataset.variables["zu"][:]
	lons = dataset.variables['lon'][:]
	lats = dataset.variables['lat'][:]
	time = dataset.variables['time'][:]
	
	return levs,lons,lats,time


#------------------------------------------------------------
# Retrieve data from a file (2d version).
# Reorders dimensions for plotting.
#------------------------------------------------------------
def get_var(dataFile,name,time,height):

	if time != -1:
		var2 = dataFile.variables[name][time,:,:,height]
	else:
		var2 = dataFile.variables[name][:,:,height]
	
	nx,ny = var2.shape
	
	return np.reshape(np.ravel(var2), (ny,nx), order='F')
	
#------------------------------------------------------------
# Retrieve data from a file (3d version). Reorders dimensions
# for plotting.
#------------------------------------------------------------
def get_var_z(dataFile,name,time,zl,zh):

	if time != -1:
		var2 = dataFile.variables[name][time,:,:,zl:zh]
	else:
		var2 = dataFile.variables[name][:,:,zl:zh]

	nx,ny,nz = var2.shape

	return np.reshape(np.ravel(var2), (nz,ny,nx), order='F')

#------------------------------------------------------------
# Destagger wind components. 
#
# @param order 0 - i,j
#			   1 - j,i
#------------------------------------------------------------
def destagger(u,v,order):
	
	if order==0:
		nx,ny = u.shape
	else:
		ny,nx = u.shape
	
	um = np.zeros((ny,nx))
	vm = np.zeros((ny,nx))
	
	for i in range(0,nx-1):
		
		um[:,i] = 0.5*(u[:,i+1]+u[:,i])
		
	for j in range(0,ny-1):
		
		vm[j,:] = 0.5*(v[j+1,:]+v[j,:])
			
	return um,vm

#------------------------------------------------------------
# Calculate density
#------------------------------------------------------------
def get_density(tb,pib,qb):

	return 100000.0*pib[:]**(717.0/287.0) / (287.0*(tb[:]*(1.0+0.61*qb[:])))


#------------------------------------------------------------
# Remove latitude and longitude lines from a plot
#------------------------------------------------------------
def remove_gridlines(pd1,pd2):	
	
	for key in pd1:
		pd1[key][0][0].remove()
	for key in pd2:
		pd2[key][0][0].remove()

#------------------------------------------------------------
# 
#------------------------------------------------------------
def mass_integral(var,rho,levs):	
	
	nx,ny,nz = var.shape

	c = np.zeros((nx,ny))

	mass = np.zeros((nz))

	for k in range(0,nz):

		mass[k] = (levs[k+2]-levs[k]) * rho[k]

	totalMass = np.sum(mass)

	for k in range(0,nz):

		c[:,:] += var[:,:,k] * mass[k]
        
	c /= totalMass

	return c

#------------------------------------------------------------
#
#------------------------------------------------------------
def main():

	try:		

		start = int(sys.argv[1])
		end = int(sys.argv[2])

		input_file = 'outfile0210_6.nc'
		output_file0 = 'plots/pres_wind0210_6fig'
		output_file1 = '../figures26/rain_wind0122_3fig'

		for i in range(start,end):

			plot_winds_temperature_pressure(input_file,output_file0+str(i)+'.png',i)
			
			#plot_winds_condensate(input_file,output_file1+str(i)+'.png',i)


	except:

		if VERBOSE:
			traceback.print_exc(file=sys.stderr)
		else:
			print >>sys.stderr,err.msg

		return 1

if __name__ == "__main__":
    sys.exit(main())