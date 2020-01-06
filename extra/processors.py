import traceback
import sys,os
import numpy as np

import matplotlib.pyplot as plt


VERBOSE = 1

#------------------------------------------------------------
#
#------------------------------------------------------------
def process_file(fname):

	nx = 140#600#40#600+6
	ny = 35#365

	pr = np.zeros((ny,nx))
	xr = np.zeros((ny,nx))
	yr = np.zeros((ny,nx))

	#-------------------------------------------
	# Open file, read in lines
	#-------------------------------------------
	with open(fname) as f:

		content = f.readlines()

	content = [x.strip() for x in content] 

	line_counter = 0

	while line_counter < len(content):

		procs = content[line_counter].split()

		pr[int(procs[1]),int(procs[0])] = float(procs[2])

		xr[int(procs[1]),int(procs[0])] = float(procs[0])
		yr[int(procs[1]),int(procs[0])] = float(procs[1])

		line_counter += 1

	print(pr)
	

	#cs = plt.contour(xr,yr,pr,cint=1,cmap=plt.cm.flag)

	plt.matshow(pr)

	plt.show()

		#print(procs[2])

#------------------------------------------------------------
#
#------------------------------------------------------------
def process_file2(fname):

	nx = 801
	ny = 801

	pr = np.zeros((ny,nx))
	xr = np.zeros((ny,nx))
	yr = np.zeros((ny,nx))

	#-------------------------------------------
	# Open file, read in lines
	#-------------------------------------------
	with open(fname) as f:

		content = f.readlines()

	content = [x.strip() for x in content] 

	line_counter = 0

	while line_counter < len(content):

		procs = content[line_counter].split()

		pr[int(procs[1]),int(procs[0])] = -float(procs[2])

		xr[int(procs[1]),int(procs[0])] = float(procs[0])
		yr[int(procs[1]),int(procs[0])] = float(procs[1])

		line_counter += 1

	print(pr)
	

	#cs = plt.contour(xr,yr,pr,cint=1,cmap=plt.cm.flag)

	plt.matshow(pr)

	plt.show()

		#print(procs[2])

#------------------------------------------------------------
#
#------------------------------------------------------------
def main():

	try:				
	
		process_file("test.txt")

	except:

		if VERBOSE:
			traceback.print_exc(file=sys.stderr)
		else:
			print >>sys.stderr,err.msg

		return 1

if __name__ == "__main__":
    sys.exit(main())
