function data = H5Read(filename)
%
%   H5READ  -   Reads data from an H5 file.
%
%   Syntax:
%           data = H5Read(filename)
%
%   where filename is the name of the H5 file and data is a vector of
%   structures that stores the data contained in filename.
%
%   Every element of the vector data is a structure with the following
%   fields:
%           parameters - the parameter vector of the integration
%           t - the time vector
%           x - the state vector
%           labels - the labels associated to every step of the
%           integration. The meaning of the labels is as follows:
%               -2  initial conditions
%               -1  state at the end of transient evolution
%                0  regular step
%               +i  intersection with the i-th Poincare' section.
%
%   Author:
%       Daniele Linaro
%       daniele.linaro@unige.it
%       November 2008
%

info = hdf5info(filename);
num_records = numel(info.GroupHierarchy.Datasets);

actual = info.GroupHierarchy.Datasets(1).Name;
dataset = hdf5read(filename, actual);
dataset = dataset';
pars = info.GroupHierarchy.Datasets(1).Attributes.Value;
data = struct('parameters',pars,'t',dataset(:,1),'x',dataset(:,2:end-1),'labels',dataset(:,end));

if num_records > 1
    data = repmat(data, [num_records, 1]);

    for ii=2:num_records
        actual = info.GroupHierarchy.Datasets(ii).Name;
        dataset = hdf5read(filename, actual);
        dataset = dataset';
        pars = info.GroupHierarchy.Datasets(ii).Attributes.Value;
        data(ii) = struct('parameters',pars,'t',dataset(:,1),'x',dataset(:,2:end-1),'labels',dataset(:,end));
    end
end