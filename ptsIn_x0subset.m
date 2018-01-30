function [idx] = ptsIn_x0subset(x,sPt,sSize,tIdx)
%
% [[idx] = ptsIn_x0subset(x,sPt,sSize,tIdx) Finds untracked particles in
% the subset of the reference image
%
% INPUTS
% -------------------------------------------------------------------------
%   x:          List of paritlces in reference image (nx3)
%   sPT:        Subset starting point (1x3)
%   sSize:      Subset size (1x3)
%   tIdx:       Index of tracked particles in x
%
%  OUTPUTS
%  ------------------------------------------------------------------------
%   idx:        Index of untracked particles in the given subset
%


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
idx = idx.*(~tIdx(idx)); %save only untracked index
idx = nonzeros(idx);

end

