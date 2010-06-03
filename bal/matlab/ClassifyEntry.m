function data = ClassifyEntry(entry, component, tol)
% data = ClassifyEntry(entry, component, tol)

if size(entry.x,1) == 2
    turns = 0;
else
    x = entry.x(3:end,:);
    coord = x(:,component);
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
if turns == -1
    keyboard
end
data = [entry.parameters' (entry.t(end)-entry.t(end-turns)) turns];
