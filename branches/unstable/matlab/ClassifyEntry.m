function data = ClassifyEntry(entry, event, component, tol)
% data = ClassifyEntry(entry, event component, tol)

if ~ exist('event','var')
    event = 1;
end
if ~ exist('component','var')
    component = 1;
end
if ~ exist('tol','var')
    tol = 1e-3;
end

if entry.labels(end) == -10
    % the integration ended with an error
    turns = -1;
    disp('the integration ended with an error.');
elseif entry.labels(end) == -3 || size(entry.x,1) == 2
    % equilibrium
    turns = 0;
else
    idx = find(entry.labels == event);
    if isempty(idx)
        turns = -1;
    else
        coord = entry.x(idx,component);
        if max(max(diff(coord))) < tol
            turns = 1;
        else
            sz = numel(coord);
            stop = floor(sz / 2);

            turns = -1;
            for ii=2:stop
                num = floor(sz/ii);
                tmp = coord(1:ii*num);
                tmp = reshape(tmp, [ii num]);
                distances = diff(tmp,1,2);
                if max(max(distances')) < tol
                    turns = ii;
                    break
                end
            end
            if ii >= stop
                turns = ii;
            end
        end
    end
end
if turns == -1
    fprintf('Unable to classify trajectory at p = %f\n', entry.parameters');
end
data = [entry.parameters' (entry.t(end)-entry.t(end-turns)) turns];
