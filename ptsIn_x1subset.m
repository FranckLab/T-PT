function idx = ptsIn_x1subset(x,sPt,sSize,utIdx)
%
% [[idx] = ptsIn_x1subset(x,sPt,sSize,tIdx) Finds untracked particles in
% the subset of the deformed image
%
% INPUTS
% -------------------------------------------------------------------------
%   x:          List of paritlces in deformed image (nx3)
%   sPT:        Subset starting point (1x3)
%   sSize:      Subset size (1x3)
%   tIdx:       Index of tracked particles in x
%
%  OUTPUTS
%  ------------------------------------------------------------------------
%   idx:        Index of untracked particles in the given subset in
%               deforemd image
%

% Define range of points within subset
minPt = sPt - sSize/2;  %min
maxPt = sPt + 3/2*sSize;    %max
% Subset size is twice the size of reference image

% Find indices of the points within the range
idx = zeros(size(x));
for i = 1:3
     idx(:,i) = x(:,i) >= minPt(i) & x(:,i) < maxPt(i);
end
idx = logical(prod(idx,2));
idx = find(idx);

idx = idx.*(~utIdx(idx)); %Only save index on untracked particles
idx = nonzeros(idx);

end