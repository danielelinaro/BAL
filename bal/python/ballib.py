
from tables import openFile
from pylab import imshow, axis, show, cm
from numpy import histogram, nonzero, zeros, linspace

def readH5file(filename):
	"""
	Reads data from an H5 file.

	Syntax:
		data = readH5file(filename)

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
	June 2010
	"""
	fid = openFile(filename, 'r')
        data = []
	for node in fid:
		if not node._v_name == '/':
			data.append({'par': node.attrs.parameters,'t': node[:,0],'x': node[:,1:-1],'labels': node[:,-1]})
	fid.close()
	return data


def saveH5file(solutions,filename):
	"""
	Saves data to an H5 file.

	Syntax:
		saveH5file(solutions,filename)

	where solutions is an array of balSolution objects and filename is the name 
	of the H5 file where data will be saved.

	Author:
	Daniele Linaro
	daniele.linaro@unige.it
	June 2010
	"""
	import tables as tbl
	from numpy import zeros, array, reshape, shape
	print 'Saving data...'
	# open the file
	fid = tbl.openFile(filename,mode='w',title='BAL data file')
	# create a filter for compression
	filter = tbl.Filters(complevel=5, complib='zlib', shuffle=True, fletcher32=False)
	# the type of data that will be stored
	atom = tbl.Float32Atom()
	# save all solutions
	for k,s in enumerate(solutions):
		# the letter p in the name means that the H5 file has been saved in Python
		name = 'p' + str(k+1).zfill(6)
		m = len(s.data['t'])
		n = len(s.data['x'])/m
		node = fid.createCArray('/',name,atom,(m,n+2),filters=filter)
		node[0:,0] = array(s.data['t'])
		node[0:,1:-1] = reshape(array(s.data['x']),(m,n))
		node[0:,-1] = array(s.data['labels'])
		node.attrs.parameters = array(s.parameters)
	# close the file
	fid.close()


def bif1d(data, coord=0, ap=0, event=1, coeff = [1.4425, -2.4283,1.9727, -0.0001]) :
	nlevels = 1024
	npars = len(data)
	pmin = data[0]['par'][ap]
	pmax = data[-1]['par'][ap]
	i = 0
	while len(data[i]['t']) < 3:
		i = i+1
	idx = nonzero(data[i]['labels'] == event)[0]
	coord_min = min(data[i]['x'][idx,coord])
	coord_max = max(data[i]['x'][idx,coord])
	for k in range(i,npars):
		if len(data[k]['t']) > 2:
			idx = nonzero(data[k]['labels'] == event)[0]
			tmp_min = min(data[k]['x'][idx,coord])
			tmp_max = max(data[k]['x'][idx,coord])
			if tmp_min < coord_min:
				coord_min = tmp_min
			if tmp_max > coord_max:
				coord_max = tmp_max
	x = linspace(pmin, pmax, npars)
	y = linspace(coord_min, coord_max, nlevels+1)

	bifdiag = zeros([nlevels,npars])
	for k,entry in enumerate(data):
		if len(entry['t']) == 2:
			continue
		idx = nonzero(entry['labels'] == event)[0]
		bifdiag[:,k] = histogram(entry['x'][idx,coord], y)[0]
	
	nth = 100
	bifdiag[bifdiag > nth] = nth;
	bifdiag = bifdiag / nth;
	bifdiag = coeff[0]*bifdiag**3 + coeff[1]*bifdiag**2 + coeff[2]*bifdiag + coeff[3];
	imshow(bifdiag, cmap=cm.gray, extent=[pmin,pmax,coord_min,coord_max])
	axis('tight')
	show()

