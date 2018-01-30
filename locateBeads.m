function [x0] = locateBeads(I,beadParameter)

% %Parameters
% thres = beadParameter.thres;    %Threshold value
% d = beadParameter.diameter;     %Bead Diameter
% xy = beadParameter.xy;          %xy um2pixel
% z = beadParameter.z;            %z um2pixel
% 
% % % load('Cell1t0.mat')
% % % tic
% % % I = double(vol{3});
% % % I = I/max(I(:));
% % % d = [2,2,6];
% % % thres = 0.2;
% 
% 
% % Sigma for DoG filtering
% dFilt = [d/xy,d/xy,d/z*3]; dFilt = [3,3,8];
% sigma1 = 1./(1+dFilt.^0.5)./dFilt;
% sigma2 = sqrt(2)*sigma1;
% DoG = imgaussfilt3(I,sigma1)-imgaussfilt3(I,sigma2);
% 
% % Local Maxima
% BW = imregionalmax(DoG,26);
% BW = BW.*(I>thres); %Only consider points whoes value is above the thres
% 
% % Select only one point within the diameter of the bead
% falsePts = 1;       %Intialize
% 
% while sum(falsePts)>0
%         
%     % Find local maxima
%     idx = find(BW);
%     x0 = zeros(length(idx),3);
%     [x0(:,1),x0(:,2),x0(:,3)] = ind2sub(size(BW),idx);
% %     x0 = x0(:,1:2)*xy;  x0 = x0*z;  % Pixel to real space
%     
%     % Nearest neighbor search within distance of diameter
%     [IDX,dist] = knnsearch(x0,x0,'k',2);
%     dist = dist(:,2);   %First point is the same point itself. 
%     dIdx = dist<5;      % Points within the diameter
%     
%     % Intensity of dIdx points in DoG 
%     DoGidx = DoG(IDX(dIdx,:));
%     DoGidx = DoGidx(:,1)-DoGidx(:,2);
%     falsePts = DoGidx<0;
%     dIdx = find(dIdx);
%     BW(idx(IDX(dIdx(falsePts),1))) = 0;
% end
% 
% 
% idx = find(BW);
% x0 = zeros(length(idx),3);
% [x0(:,1),x0(:,2),x0(:,3)] = ind2sub(size(BW),idx);

% 

%Parameters
thres = beadParameter.thres;    %Threshold value
minPixels = beadParameter.minSize;  %Minimum pixel count in blob for bead
maxPixels = beadParameter.maxSize;  %Maximum pixel count in blob for bead

%Binary Image
BW = I>thres;

% Find bead blobs
CC = bwconncomp(BW);
numPixels = cellfun(@numel,CC.PixelIdxList);
beadBlob = numPixels>minPixels & numPixels<maxPixels;

% Find centroid of all the eligible blob;
S = regionprops(CC,'Centroid');
blobPts = round(double(struct2dataset(S)));
blobPts = blobPts(beadBlob,:);
temp = blobPts;

% in m,n,o coordinates
blobPts(:,1) = temp(:,2);
blobPts(:,2) = temp(:,1);
x0 = blobPts;


end

