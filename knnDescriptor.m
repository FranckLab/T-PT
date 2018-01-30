function [fd,fm] = knnDescriptor(x,IDX,knnFD,knnFM,nSpheres)
% Descriptor based on nearest neighbor

% Inputs:
% x (n by 3): n 3 dimensional points
% knn: Number of neighbor used for descriptor.
% Outputs:
% knnD: The descriptor

knn = max(knnFD,knnFM);

%Find nearest neighbours
[idx,D] = knnsearch(x,x(IDX,:),'K',knn+1);

% Create nn list for neighbor verification
fm = idx(:,2:knnFM+1);
fm = fm';

% Create knn feature descriptor
idx = idx(:,2:knnFD+1);
D = D(:,2:knnFD+1)';

% Radius feature
D = D.^3;
unitD = D(end,:)/nSpheres;
unitD = 1./unitD;
D = bsxfun(@times,D,unitD)-0.0001;
D = floor(D);

%%% This is a mess, but it works. Good luck if you want to change it.
% Divide knn into quadrants in binary
y = idx';
y = y(:);
y = x(y,:);
y = y';
y = y(:);
% y = reshape(y,[3,knn,length(x)]);
y = reshape(y,[3,knn,sum(IDX)]);
% z(:,1,:)=x';
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

%Create histcounts descriptor for each point
edgePt = nSpheres*8 +0.5;
edges = 0.5:1:edgePt;
fd = zeros(length(x),length(edges)-1);
fd = zeros(sum(IDX),length(edges)-1);

for i = 1:sum(IDX)
        t1 = histcounts(z(:,i),edges);
        fd(i,:) = t1;
end

end