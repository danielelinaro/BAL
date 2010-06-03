
from tables import openFile
from pylab import *
from numpy import histogram

def h5read(filename):
	"""
	Reads data from an H5 file.

	Syntax:
		data = h5read(filename)

	where filename is the name of the H5 file and data is a vector of
	dictionaries that stores the data contained in filename.

	Every element of the vector data is a dictionary with the following
	keys:
		par - the parameter vector of the integration
		t - the time vector
		x - the state vector
		labels - the labels associated to every step of the
		integration. The meaning of the labels is as follows:
			-10  integration error
			 -3  equilibrium point
			 -2  initial conditions
			 -1  state at the end of transient evolution
			  0  regular step
			 +i  intersection with the i-th Poincare' section.

	Author:
	Daniele Linaro
	daniele.linaro@unige.it
	November 2008
	"""

	fid = openFile(filename, 'r')
	data = []
	cnt = 0
	for node in fid:
		if cnt > 0:
			data.append({
				'par': node.attrs.parameters,
				't': node[:,0],
				'x': node[:,1:-1],
				'labels': node[:,-1]
			})
		cnt = cnt+1
	fid.close()
	return data


def bif1d(data, coord, ap, coeff = [1.4425, -2.4283,1.9727, -0.0001]) :
	nlevels = 1024
	npars = len(data)
	pmin = data[0]['par'][ap]
	pmax = data[-1]['par'][ap]
	i = 0
	while len(data[i]['t']) < 3:
		i = i+1
	coord_min = min(data[i]['x'][2:,coord])
	coord_max = max(data[i]['x'][2:,coord])
	for k in range(i,npars):
		if len(data[k]['t']) > 2:
			tmp_min = min(data[k]['x'][2:,coord])
			tmp_max = max(data[k]['x'][2:,coord])
			if tmp_min < coord_min:
				coord_min = tmp_min
			if tmp_max > coord_max:
				coord_max = tmp_max
	x = linspace(pmin, pmax, npars)
	y = linspace(coord_min, coord_max, nlevels)

	bif1d = zeros([nlevels, npars])
	for k in range(npars):
		if len(data[k]['t']) == 2:
			continue
		bif1d[:,k] = histogram(data[k]['x'][2:,coord], y)[0]
	
	nth = 50
	bif1d[bif1d > nth] = nth;
	tmp = bif1d / nth;
	tmp = coeff[0]*tmp**3 + coeff[1]*tmp**2 + coeff[2]*tmp + coeff[3];
	tmp = tmp * nth
	imshow(tmp, cmap=cm.winter)
	axis('tight')
	show()

