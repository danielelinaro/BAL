function classified = Classifier(data, varargin)
% classified = Classifier(data)
% classified = Classifier(data, component)
% classified = Classifier(data, component, tol)

component = 1;
tol = 1e-3;
if nargin == 2 && isscalar(varargin{1})
    component = varargin{1};
end
if nargin == 3 && isscalar(varargin{2})
    tol = varargin{2};
end

sz = numel(data);
classified = zeros(sz, numel(data(1).parameters)+2);
npr = round(sz / 10);

for ii=1:sz
    classified(ii,:) = ClassifyEntry(data(ii), component, tol);
    if mod(ii,npr) == 0
        fprintf('%d/%d\n',ii,sz);
    end
end