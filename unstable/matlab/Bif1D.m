function Bif1D(data, coord, ap, event, varargin)
%
%   function Bif1D(data, coord, ap, event)
%   function Bif1D(data, coord, ap, event, transient)
%   function Bif1D(data, coord, ap, event, transient, coeff)
%
coeff = [1.4425, -2.4283, 1.9727, -0.0001];
if nargin == 6 && numel(varargin{2}) == 4
    coeff = varargin{1};
end

if nargin == 5 && isscalar(varargin{1})
    transient = varargin{1};
    for ii=1:numel(data)
        ndx = find(data(ii).t >= transient);
        data(ii).t = data(ii).t(ndx);
        data(ii).x = data(ii).x(ndx,:);
        data(ii).labels = data(ii).labels(ndx);
    end
end

nlevels = 1024;
npars = length(data);
pmin = data(1).parameters(ap);
pmax = data(end).parameters(ap);

% find the minimum and maximum of the diagram
ii = 1;
while size(data(ii).x,1) < 3
    ii = ii+1;
end
ind = find(data(ii).labels == event);

coord_min = min(data(ii).x(ind,coord));
coord_max = max(data(ii).x(ind,coord));
for k=ii:npars
    if size(data(k).x,1) > 2
        ind = find(data(k).labels == event);
        if isempty(ind)
            continue;
        end
        tmp_min = min(data(k).x(ind,coord));
        tmp_max = max(data(k).x(ind,coord));
        if tmp_min < coord_min
            coord_min = tmp_min;
        end
        if tmp_max > coord_max
            coord_max = tmp_max;
        end
    end
end

% calculate the histograms
figure('renderer', 'painters'); axes; hold on;
x=linspace(pmin, pmax, npars);
y=linspace(coord_min, coord_max, nlevels);

bif1D = repmat(NaN, [nlevels npars]);

for k=1:npars
    if length(data(k).t) == 2
        continue;
    end
    bif1D(:, k) = hist(data(k).x(data(k).labels == event,coord), y);
    if mod(k,1000) == 0
        fprintf(1, 'k = %d\n', k);
    end
end

% some image processing and finally the plot
nth = 50;
bif1D(bif1D > nth) = nth;
tmp = bif1D/nth;
tmp = coeff(1)*tmp.^3 + coeff(2)*tmp.^2 + coeff(3)*tmp + coeff(4);
tmp = uint8(tmp*nth);
image(x, y, uint8(tmp)); 
axis tight;
colormap(flipud(bone(nth)));


