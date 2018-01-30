function [fd,fm] = knnDescriptor(x,IDX,knnFD,knnFM,nSpheres)
%
% [fd,fm] = knnDescriptor(x,IDX,knnFD,knnFM,nSpheres) computes particle
% descriptor from neighboringparticles for particles with IDX index in x
%
% INPUTS
% -------------------------------------------------------------------------
%   x:          Particle position array (n x 3)
%   IDX:        Index of particles in x for which particle descritpor needs
%               to be computed.
%   knnFD:      Number of neighboring particles to compute feature
%               descriptor
%   knnFM:      Number of neighboring particles to be used in similarity of
%               neighborhood test.
%   nSpheres:   Number of shells in feature descriptor
%
%  OUTPUTS
%  ------------------------------------------------------------------------
%   fd:         Feature descriptor for requested particles
%   fm:         Index of neighoring particles for similarity of
%               neighborhood test
%

% Number of neighboring particles to find. Max of knnFD and knnFM
knn = max(knnFD,knnFM);

%Find nearest neighbours
[idx,D] = knnsearch(x,x(IDX,:),'K',knn+1);

% Create nn list for neighbor verification
fm = idx(:,2:knnFM+1);
fm = fm';

% Create knn feature descriptor
idx = idx(:,2:knnFD+1);
D = D(:,2:knnFD+1)';

% Radius feature/ Distance descriptor
D = D.^3;
unitD = D(end,:)/nSpheres;
unitD = 1./unitD;
D = bsxfun(@times,D,unitD)-0.0001;
D = floor(D);

% Vecotrized this section of code to optimize for speed. It is messy but it
% works. Good luck if you want to change it

% Divide knn into quadrants in binary
y = idx';
y = y(:);
y = x(y,:);
y = y';
y = y(:);
y = reshape(y,[3,knn,sum(IDX)]);
z(:,1,:)=x(IDX,:)';
y = bsxfun(@minus,y,z);
y = y>0;

%Combine quadrants and distance descriptors
y(4,:,:) = D;

% Convert binary position of knn to number
z = [1,2,4,8];
zz(:,1,1) = z;
z = bsxfun(@times,y,zz);
z = sum(z,1);
z = squeeze(z+1);

%Create histcounts descriptor for each bin
edgePt = nSpheres*8 +0.5;
edges = 0.5:1:edgePt;
fd = zeros(length(x),length(edges)-1);
fd = zeros(sum(IDX),length(edges)-1);

for i = 1:sum(IDX)
        t1 = histcounts(z(:,i),edges);
        fd(i,:) = t1;
end

end