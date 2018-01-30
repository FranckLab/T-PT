function [x,track] = funPSPT(varargin)
% funPSPT is the main function that performs the single particle tracking
% on time increment of volumetric images.
%
% INPUTS
% -------------------------------------------------------------------------
%   filename: 
%
% OUTPUTS
% -------------------------------------------------------------------------
%   u: 


%Parse inputs
[fileInfo, beadParameter, psptParameter] = parseInputs(varargin{:});
nChannel = length(fileInfo); %Number of channel

%% Particle Tracking

% Find particle position in t = 0
% I = loadFile(fileInfo,1,beadParameter.randNoise);
% x{1} = locateBeads(I,beadParameter);
% x{1} = radialcenter3dvec(I0,x0,winSize,dccd,abc,0);

% Start Tracking
for t = 2:length(fileInfo{1}.filename)
    tStart = tic;
    
    %Start DVC
    disp(['Current time point: t = ' num2str(t-1)])
    
    for j = 1:nChannel        
        %% Find particle position
        tPP = tic;
        disp(['  Particle channel: ' num2str(j)])
        disp(['  Current Filename: ' fileInfo{j}.filename{t}])
        
        if t == 2 % timePoint = 1
            I = loadFile(fileInfo{j},1,beadParameter{j}.randNoise);
            x{1}{j} = locateBeads(I,beadParameter{j});
            x{1}{j} = radialcenter3dvec(I,x{1}{j},beadParameter{j});
        end
        
        % Load file
        I = loadFile(fileInfo{j},t,beadParameter{j}.randNoise); 
        % Detect beads
        x{t}{j} = locateBeads(I,beadParameter{j}); 
        % Localize beads
        x{t}{j} = radialcenter3dvec(I,x{t}{j},beadParameter{j});
        
        disp(['    Time to localize particles = ', num2str(toc(tPP)),' seconds']);
        disp(['    Number of  particles = ', num2str(size(x{t}{j},1))]);
        
        
        %% Prepapre predictor if nChannel>1 and 
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
        
        %% Particle Tracking through P-SPT
        psptParameter{j}.sizeI = size(I);
        [track{t-1}{j}] = PSPT(x{t-1}{j},x{t}{j},psptParameter{j},predictor);
        
    end
    
    
    
    disp(['Total Elaspsed Time = ', num2str(toc(tStart)),' seconds']);fprintf('\n');
    
end

end

function varargout = parseInputs(varargin)
%  = parseInputs(filename,beadParameter, psptParameter)


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

%%% P-SPT Parameters
psptParameter = varargin{3};

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
    
    parse(p,psptParameter{i})
    
    psptParameter{i} = p.Results;
    
end

%%% Outputs
varargout{1} = fileInfo;
varargout{2} = beadParameter;
varargout{3} = psptParameter;

end


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
