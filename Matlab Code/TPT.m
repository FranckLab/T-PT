function [track] = TPT(varargin)
%
% [track] = TPT(varargin) is the function which links particles between two
% consecutive image frames
%
% INPUTS
% -------------------------------------------------------------------------
%   x0:             Particle position at t = to
%   x1:             Particle position at t = to+1
%   tptParameter:   TPT algorithm parameters
%   predictor:      Displacement field from multi-attribute particles to
%                   predict displacement field. This can also be extended
%                   to predict displacement computed from previous time
%                   points by extrapolation
%
%  OUTPUTS
%  ------------------------------------------------------------------------
%   track:          Array storing particle link from x0 to x1
%

%%%% Parse Inputs %%%%
%-------------------------------------------------------------------------%

x0 = varargin{1};               % Particle position at t=0
x1 = varargin{2};               % Particle position at t=1
tptParameter = varargin{3};     % TPT parameter
predictor = varargin{4};        % Predictor for displacement field
% In current form it uses bead displacement computed from previous channel
% but this can be extended to predictor displacement computed from previous
% time points by extrapolation.

%Parse tpttParameter
knnFD = tptParameter.knnFD;    % # of NN in feature descriptor
knnFM = tptParameter.knnFM;    % # of NN in neighborhood similarity comparison
fmThres = tptParameter.fmThres;    % Threshold for minimum # of common NN in similiarity comparison
maxIter = tptParameter.maxIter;    % Maximum number of iterations in pspt
nSpheres = tptParameter.nSpheres;  % Number of spheres in feature descriptor
outlrThres = tptParameter.outlrThres;  % Residual threshold in outlier removal
sizeI = tptParameter.sizeI;    % Size of image


%%%% Initalize variables for tracking %%%%
%-------------------------------------------------------------------------%

% For bookeeping operation in tracking
track = zeros(size(x0(:,1)));   % Store track results
x0track = false(size(x0(:,1))); % Binary flag of x0 idx if tracked
x1track = false(size(x1(:,1))); % Binary flag of x1 idx if tracked
x0FD = []; % Initialize for feature descriptor for x0
x0FM = []; % Initialize for feature match for x0


% Defining parameters for iteration convergence criteria
nBeads = size(x0,1);    % Number of beads in x0
nTracked = 0;           % Number of beads tracked
iterCrit = 1;           % Iteration Critera (1: Continue iteration, 0: Stop iteration)
iter = 1;               % Current iteration number


% Create subset for first iteration
sSize(1) = 256;    % Subset Size
% FD and FM for x0. This doesn't change in future iterations
[x0FD(~x0track,:),x0FM(:,~x0track)] = knnDescriptor(x0,~x0track,knnFD,knnFM,nSpheres);
% Save true particle position values. These don't get warped
y1= x1; y0=x0;  
% Timer
tPSPT = tic;


%%%% TPT iterations %%%%
%-------------------------------------------------------------------------%

fprintf(['    TPT based tracking \n']);

while iterCrit == 1
    
    tIter = tic;    %Iterations timer
    
    
    %%%% STEP 2: Iterative deformation warping %%%% -----------------------
    
    % Compute displacement from tracked particles
    x0_ = y0(x0track,:); x1_ = y1(track(x0track),:);
    u = x0_ - x1_;  %reverse displacement
    
    % Utlize predictor information for deformation warping
    if predictor.flag == 1 & size(u,1)<2*size(predictor.u)
        x0_ = [x0_;predictor.x0];
        x1_ = [x1_;predictor.x1]; 
        u = [u;predictor.u];
    end
    
    % Deform x0 and x1 for next iterations
    u1 = zeros(size(x1));   %Initialize to store displacement value
    try
        parfor j =1:3
            %             F = scatteredInterpolant(x0_,u(:,j),'linear','linear');
            %             u0(:,j) = 0.5*F(x0);
            
            F = scatteredInterpolant(x1_,u(:,j),'linear','linear');
            u1(:,j) = F(x1);
        end
        %         x0 = x0+u0;
        x1 = y1+u1;
    end
    
    
    
    %%%% Setting up grids for particle matching %%%% ----------------------
    
    %Create grid for subset
    sGrid = cell(3,1);
    for i = 1:3
        sGrid{i} = 1:sSize(end):sizeI(i);
    end
    [sGrid{1},sGrid{2},sGrid{3}] = ndgrid(sGrid{1},sGrid{2},sGrid{3});
    nSSgrid = length(sGrid{1}(:));
    
    % FD and FM for x1. Computed at each iterations
    x1FD = []; x1FM = [];       %Intialize
    [x1FD(~x1track,:),x1FM(:,~x1track)] = knnDescriptor(x1,~x1track,knnFD,knnFM,nSpheres);    %Compute for all points
    
    % Storing particle at each sGrid point
    match = cell(nSSgrid,1);    %Match
    ia = cell(nSSgrid,1);       %
    
    
    
    %%%% Particle Matching %%%% -------------------------------------------
    
    % Iterative particle matching in each subset
    parfor i = 1 : nSSgrid          % Parallel computing
        
        %Subset defining point
        sPt = [sGrid{1}(i),sGrid{2}(i),sGrid{3}(i)];
        
        %Find index of points within the required subset in x0 and x1, and
        %remove points already in the track and utrack
        idx0 = ptsIn_x0subset(x0,sPt,sSize(end),x0track);
        idx1 = ptsIn_x1subset(x1,sPt,sSize(end),x1track);
        
        % Collect points for FD matching and FM link verification in subset
        x0FDss = x0FD(idx0,:);  %x0fMss = x0fM(:,idx0);
        x1FDss = x1FD(idx1,:);  %x1fMss = x1fM(:,idx1);
        
        % Feature matching by minimizing FD distanceto create temporary
        % particle links
        [ssMatch,~] = knnsearch(x1FDss,x0FDss,'K',1);
        ssMatch = idx1(ssMatch);    %Indexing crapshot: Need to do global indexing
        %ssmatch(D>=2) = 0;      %Selected distance threshold of 2 
        
        %Find uniques in feature match
        [match{i}, ia{i}, ~] = match121(ssMatch,idx0,idx1);
        %ia is utlized in particle link verification
    end
    
    %%%% Link verification %%%% -------------------------------------------
    
    % Find unique match and ia (Because match and ia from x0 particles in
    % different subset can point to same x1 particle. Need to remove them)
    match = cell2mat(match);
    ia = cell2mat(ia);
    
    %remove common match
    [~,iaIdx,~] = unique(match);
    ia = ia(iaIdx);
    match = match(iaIdx);
    
    %Assemble global match for link verification
    tempMatch = track;  % Approved matches from previous iterations
    tempMatch(ia) = match;  % New matches produced from current iteration
    
    % Verify matched feature using neighbors.
    for j = ia'
        % Collect neighboring particles for x0 and x1 particles
        x0n = x0FM(:,j)';   x0n = tempMatch(x0n);
        x1n = x1FM(:,tempMatch(j))';
        
        % Verify particle link if more than fmThres same neighboring
        % particles
        if length(intersect(x0n,x1n)) >= fmThres
            track(j) = tempMatch(j);
        end
    end
    
    % Remove Outliers
    [track] = removeOutlierTPT(x0,x1,track,outlrThres);
    
    % Flag tracked particles in x0 and x1
    x0track = track>0;
    x1track(track(x0track)) = 1;
    
    
    %%%%% Convergence criteria and printing time log %%%% -----------------
    
    nTracked(end + 1) = sum(track>0)/nBeads;    %Fraction of beads tracked
    iterCheck = iter>=maxIter;
    [iterCrit,sSize] = iterCriteria(sSize,nTracked,iterCheck);
    fprintf('    Elapsed time (Iteration %d): %0.4f secs   Partices tracked: %0.2f   Subset Size = %d\n',iter,toc(tIter),sum(track>0)/length(x1)*100,sSize(end-1));
    iter = iter+1;
    
end


    
%%%% Displacement predictor based on neighbors %%%% -----------------------

fprintf(['    Displacement predictor based tracking \n']);
iter = 1;
while iter<5
    tIter = tic;    %Iterations timer
    
    u = y1(track(x0track),:) - y0(x0track,:);  % Displacement between tracked pts
    
    %Interpolate displacement to untracked pts
    uPred = zeros(sum(~x0track),3);
    try
        parfor j =1:3
            F = scatteredInterpolant(y0(x0track,:),u(:,j),'natural','linear');
            uPred(:,j) = F(y0(~x0track,:));
        end
    end
    
    % Find nearest neighbor distance for untracked pts
    idx0 = find(~x0track);
    [~,Dx0] = knnsearch(y0,y0(idx0,:),'k',2);
    Dx0 = Dx0(:,2);
    
    %Find nearest pt in x1 based on x0+uPred
    idx1 = find(~x1track);
    [match,D] = knnsearch(y1(~x1track,:),y0(~x0track,:)+uPred,'k',1);
    D = bsxfun(@minus,D,Dx0/3); D = D<0; %D = sum(D,2); D = D==1;
    
    % Find unique matches
    [match, ia, ~] = match121(idx1(match(D)),idx0(D),idx1);
    
    %Assemble global match
    tempMatch = track;
    tempMatch(ia) = match;
    
    % Verify matched feature using neighbors.
    for j = ia'
        x0n = x0FM(:,j)';   x0n = tempMatch(x0n);
        x1n = x1FM(:,tempMatch(j))';
        if length(intersect(x0n,x1n)) >= fmThres
            track(j) = tempMatch(j);
        end
    end
    
    % Remove Outliers
    [track] = removeOutlierTPT(y0,y1,track,outlrThres);
    
    % Flag tracked particles in x0 and x1
    x0track = track>0;
    x1track(track(x0track)) = 1;
    
    fprintf('    Elapsed time (Iteration %d): %0.4f secs   Partices tracked: %0.2f  \n',iter,toc(tIter),sum(track>0)/length(x1)*100);
    iter = iter+1;    

end

end