function [x, track] = funTPT(varargin)
%
% [x, track] = funTPT(fileInfo, beadParameter, tptParam) is the main
% function that performs the single particle tracking on time increment of
% volumetric images.
%
% VARIABLES OPTIONS
% -------------------------------------------------------------------------
%   fileInfo: string for the filename prefix to load the volumetric images
%             in the current directory.
%             Input options:
%             --- If the image is not within a cell ---
%             1) fileInfo{c}{1} = 'filename*.mat' or 'filename'
%
%             --- If image is within a cell (Recommended approach) ---
%             2) fileInfo{c}{1} = 'filename*.mat' or 'filename'
%                fileInfo{c}{2} = Channel number containing images you want
%                                 to run TPT on. If the channel is not
%                                 provided, then channel = 1.
%
%   beadParam: Parameters to detect and localize particles in the images
%              Input options:
%              1) beadParam{c}.thres = value between 0 & 1. (Default 0.5)
%                                      default value = 0.5
%
%                 The threshold value used to convert input images into
%                 binary images via image thresholding operations. This
%                 operation is used to detect particles in the images.
%
%              2) beadParam{c}.minSize = int value between 0 and Inf
%                                        default value = 3
%
%                 The minimum size of particles in voxels as connected
%                 components in the binary images. The connected components
%                 with size smaller than this parameter are discarded as
%                 noise.
%
%              3) beadParam{c}.maxSize = int value between 1 and Inf
%                                        default value = Inf
%
%                 The maximum size of particles in voxels as connected
%                 components in the binary images. The connected components
%                 with size larger than this parameter are discarded as
%                 multiply-connected particles
%
%              4) beadParam{c}.winSize = size in a three column array
%                                        default value = [7, 7, 7]
%
%                 The image subset size used to localize particle using
%                 radial symmetry method. Select a size such that a single
%                 particle just fits in the image subset.
%
%   tptParam:  Parameters used by T-PT to track particles
%              Input options:
%
%              1) tptParam{c}.knnFM = int value (default = 5)
%
%                 The total number of neighboring particles ('q') in the
%                 that are investigated in similarity of neighborhood test.
%
%              2) tptParam{c}.fmThres = int value (default = 2)
%
%                 The total number of neighboring particles ('p') out of
%                 'q' that need to match between the reference and deformed
%                 images to satisfy simalrity of neighborhood test.
%
%                 Note: Increase p/q ratio for low spatial frequency
%                 displacement field and decrease p/q ratio for high
%                 spatial frequency displacement
%
%              3) tptParam{c}.outlrThres = real value > 0 (default = 5)
%
%                 The threshold value for the normalized residual in the
%                 universal median test. Increase for high spatial
%                 frequency displacement field and decrease for low spatial
%                 frequency displacement field.
%
%              4) tptParam{c}.knnFD = int value > 0 (default = 16)
%
%                  The number of nearest neighboring particles used to
%                  compute particle descriptor
%
%              5) tptParam{c}.nSpheres = int value > 0 (default = 2)
%
%                  The number of concentric shells used to compute particle
%                  descriptor.
%
%   NOTE- 'c' lists the order in which multi-attribute particles are
%         tracked. If only single attribute particles, c = 1. If
%         multi-attribute particles, list all the properties for different
%         c values (particles)
%
%
%  OUTPUTS
%  ------------------------------------------------------------------------
%   x:      Particle positions found for each 'c' particle at each 'time'
%           frame
%
%           x{time}{c} = particle positions at time frame in MxNxO format
%
%   track:  Particles links between two consecutive image frames (t = to &
%           t = t0 + 1) for each 'c' particle at each time frame (t = time)
%
%
%           track{time}{c} = an array of size [length(x{time}{c}, 1]. The
%                            i_th index in the array stores the index of
%                            matched particle in x{time}{c}. If i_th index
%                            value is 0, no particle link is found in the
%                            next image frame


%Parse inputs
[fileInfo, beadParameter, tptParameter] = parseInputs(varargin{:});
nChannel = length(fileInfo); %Number of multi-attribute particles

%% Particle Tracking

% Start Tracking
for t = 2:length(fileInfo{1}.filename) % Loop through all time points
    
    tStart = tic;
    disp(['Current time point: t = ' num2str(t-1)])
    
    for j = 1:nChannel % Loop through different types of multi-attribute particles
        
        % Detect and Localize Particles -----------------------------------
        tPP = tic;
        disp(['  Particle channel: ' num2str(j)])
        disp(['  Current Filename: ' fileInfo{j}.filename{t}])
        
        if t == 2 % timePoint = 1
            I = loadFile(fileInfo{j},1,beadParameter{j}.randNoise);
            x{1}{j} = locateParticles(I,beadParameter{j});
            x{1}{j} = radialcenter3dvec(I,x{1}{j},beadParameter{j});
        end
        
        I = loadFile(fileInfo{j},t,beadParameter{j}.randNoise); %Load image
        x{t}{j} = locateParticles(I,beadParameter{j}); % Detect particles
        x{t}{j} = radialcenter3dvec(I,x{t}{j},beadParameter{j}); % Localize particles
        
        disp(['    Time to localize particles = ', num2str(toc(tPP)),' seconds']);
        disp(['    Number of  particles = ', num2str(size(x{t}{j},1))]);
        
        
        % Prepapre predictor for multi-attribute particles ----------------
        % (nChannel > 1)
        if j>1
            predictor.flag = true;
            predictor.x0 = x{t-1}{j-1};
            temp.x1 = x{t}{j-1};
            temp.track = track{t-1}{j-1};
            predictor.u = temp.x1(temp.track(temp.track>0),:)...
                - predictor.x0(temp.track>0,:);
            predictor.x1 = temp.x1(temp.track(temp.track>0),:);
            
        elseif j == 1
            predictor.flag = false;
            predictor.x0 = [];
            predictor.u = [];
            
        end
        
        % Particle Tracking -----------------------------------------------
        tptParameter{j}.sizeI = size(I);
        track{t-1}{j} = TPT(x{t-1}{j},x{t}{j},tptParameter{j},predictor);
        
    end
    
    disp(['Total Elaspsed Time = ', num2str(toc(tStart)),' seconds']);fprintf('\n'); 
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function varargout = parseInputs(varargin)
% varargout = parseInputs(filename,beadParameter, psptParameter)


%%% Parse filenames
filename = varargin{1};
for i = 1 : length(filename)
    
    % Find channel number within file
    if length(filename{i}) == 1
        fileInfo{i}.datachannel = 1; % Assume channel 1 if not mentioned
    else
        fileInfo{i}.datachannel = filename{i}{2};
    end
    
    % Extract filename to load
    filename{i} = filename{i}{1};
    
    [~,filename{i},~] = fileparts(filename{i});
    filename{i} = dir([filename{i},'.mat']);
    fileInfo{i}.filename = {filename{i}.name};
    
    % Error if file doesn't exist
    if isempty(fileInfo{i})
        error('File name doesn''t exist');
    end
end


%%% Detection and Localization Parameters
beadParameter = varargin{2};

% Define default values
thres = 0.5;
minSize = 5;
maxSize = Inf;
winSize = [7,7,7];
dccd = [1,1,1];
abc = [1,1,1];
forloop = 1;
randNoise = 1/10^7; % Something small
xy = 1; % Assume symmetrical if not given
z = 1;  % Assume symmetrical if not givenf
diameter = 5;   % Assume if not given

for i = 1:length(beadParameter)
    
    p = inputParser;
    addParameter(p,'thres',thres);
    addParameter(p,'minSize',minSize);
    addParameter(p,'maxSize',maxSize);
    addParameter(p,'winSize',winSize);
    addParameter(p,'dccd',dccd);
    addParameter(p,'abc',abc);
    addParameter(p,'forloop',forloop);
    addParameter(p,'randNoise',randNoise);
    addParameter(p,'xy',xy);
    addParameter(p,'z',z);
    addParameter(p,'diameter',diameter);
    
    parse(p,beadParameter{i})
    
    beadParameter{i} = p.Results;
    
end

%%% TPT Parameters
tptParameter = varargin{3};

% Define default values
knnFD =16;
knnFM = 4;
fmThres = 2;
maxIter = 14;
nSpheres = 2;
outlrThres = 5;

for i = 1:length(beadParameter)
    
    p = inputParser;
    addParameter(p,'knnFD',knnFD);
    addParameter(p,'knnFM',knnFM);
    addParameter(p,'fmThres',fmThres);
    addParameter(p,'maxIter',maxIter);
    addParameter(p,'nSpheres',nSpheres);
    addParameter(p,'outlrThres',outlrThres);
    
    parse(p,tptParameter{i})
    
    tptParameter{i} = p.Results;
    
end

%%% Outputs
varargout{1} = fileInfo;
varargout{2} = beadParameter;
varargout{3} = tptParameter;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function I = loadFile(fileInfo,idx,randNoise)
I = load(fileInfo.filename{idx});
fieldName = fieldnames(I);
I = getfield(I,fieldName{1});
if iscell(I),
    if numel(I)==1, I = I{1};
    else
        I = I{fileInfo.datachannel};
    end
end

I = double(I);
I = I/max(I(:));
end
