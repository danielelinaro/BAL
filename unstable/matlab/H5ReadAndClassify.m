function block = H5ReadAndClassify(filename, varargin)
%
%  block = H5ReadAndClassify(filename)
%  block = H5ReadAndClassify(filename,component)
%  block = H5ReadAndClassify(filename,component,tol)
%
%  Reads an h5 file saved with the BAL convention and classifies the
%  data in the file. The whole bunch of trajectories contained in the file
%  is discarded and not returned to the caller.
%
%  Author: Daniele Linaro - daniele.linaro@unige.it
%          March 2009

component = 1;
tol = 1e-3;
if nargin == 2 && isscalar(varargin{1})
    component = varargin{1};
end
if nargin == 3 && isscalar(varargin{2})
    tol = varargin{2};
end

fprintf(1, 'Reading informations from file %s ...', filename);
info = hdf5info(filename);
fprintf(1, ' done.\n');
num_records = numel(info.GroupHierarchy.Datasets);
fprintf(1, 'The number of datasets in the file is %d.\n', num_records);

npar = numel(info.GroupHierarchy.Datasets(1).Attributes.Value);
block = zeros(num_records, npar+2);
npr = floor(num_records / 10);

fprintf(1, 'Reading the datasets...\n')
for ii=1:num_records
    actual = info.GroupHierarchy.Datasets(ii).Name;
	dataset = hdf5read(filename, actual);
    dataset = dataset';
    par = info.GroupHierarchy.Datasets(ii).Attributes.Value;
    entry = struct('parameters',par,'t',dataset(:,1),'x',dataset(:,2:end-1),'labels',dataset(:,end));
	block(ii,:) = ClassifyEntry(entry,component,tol);
	
    if mod(ii,npr) == 0
		fprintf(1, '%6d / %6d\n', ii, num_records);
    end
end


