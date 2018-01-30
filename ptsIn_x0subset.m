function [idx] = ptsIn_x0subset(x,sPt,sSize,tIdx)
% Find points in subset in x0;
%x0ss

% Define range of points within subset
minPt = sPt;  %min
maxPt = sPt + sSize;    %max

% Find indices of the points within the range
idx = zeros(size(x));
for i = 1:3
     idx(:,i) = x(:,i) >= minPt(i) & x(:,i) < maxPt(i);
end
idx = logical(prod(idx,2));
idx = find(idx);
idx = idx.*(~tIdx(idx));
idx = nonzeros(idx);
end

