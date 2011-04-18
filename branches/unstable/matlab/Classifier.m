function classified = Classifier(data, event, component, tol)
% classified = Classifier(data, event, component, tol)

if ~ exist('event','var')
    event = 1;
end
if ~ exist('component','var')
    component = 1;
end
if ~ exist('tol','var')
    tol = 1e-3;
end

sz = numel(data);
classified = zeros(sz, numel(data(1).parameters)+2);
npr = round(sz / 10);

for ii=1:sz
    classified(ii,:) = ClassifyEntry(data(ii), event, component, tol);
    if mod(ii,npr) == 0
        fprintf('%d/%d\n',ii,sz);
    end
end