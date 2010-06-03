function PlotterColormap(data, nlevels, varargin)
% PlotterColormap(data, nlevels, varargin)

[p1steps p2steps options] = guess_size(data);

cmap = flipud(hot(nlevels));
if nargin == 3 && isequal(size(varargin{1}),[nlevels 3])
    cmap = varargin{1};
end

mode = 1;
if data(1,2) == data(2,2)
    mode = 2;
end
switch(mode)
    case 1,
        p1 = data(1:p2steps:end,1);
        p2 = data(1:p2steps,2);
    case 2,
        p1 = data(1:p1steps,1);
        p2 = data(1:p1steps:end,2);
end

asc(1) = p1(1) < p1(2);
asc(2) = p2(1) < p2(2);
p1 = sort(p1);
p2 = sort(p2);

switch(mode)
    case 1,
        dbif = reshape(data(:,end),[p2steps,p1steps]);
    case 2,
        dbif = reshape(data(:,end),[p1steps,p2steps])';
end
dbif(dbif > nlevels) = nlevels;

if ~ asc(1)
    dbif = fliplr(dbif);
end
if ~ asc(2)
    dbif = flipud(dbif);
end

imagesc(p1,p2,dbif);
colormap(cmap);
set(gca,'YDir','normal');

%%
function [par_1_steps par_2_steps options] = guess_size(block)

options = [0 0];
tot = size(block,1);
if diff(block(1:2,1)) == 0
    col = 1;
    tmp = find(block(:,col) == block(1,col));
    par_2_steps = tmp(end);
    par_1_steps = tot/par_2_steps;
else % diff(block(1:2,col)) ~= 0
    col = 2;
    options(1) = 1;
    tmp = find(block(:,col) == block(1,col));
    par_1_steps = tmp(end);
    par_2_steps = tot/par_1_steps;
end
if col == 1 && block(2,2) < block(1,2)
    options(2) = 1;
elseif block(2,1) < block(1,1)
    options(2) = 2;
end
