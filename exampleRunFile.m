%% Example run file for Topology-based Particle Tracking
%
% The two volumetric images that are going to be run are 'I00.mat' and
% 'I01.mat'. The images are prescribed with a sinusoidal displacement field
% of linearly decreasing wavelength and amplitude with a displacement
% parameter d = 1.5
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
               
clear; close all;

% Multi-attribute particle type 1
fileInfo{1}{1} = 'I*.mat'; %filename
fileInfo{1}{2} = 1; %Channel number

% Bead Parameter
beadParam{1}.thres = 0.5;
beadParam{1}.minSize = 3;
beadParam{1}.maxSize = 1000;
beadParam{1}.winSize = [5, 5, 5];

% TPT Parameter
tptParam{1}.knnFM = 5;
tptParam{1}.fmThres = 2;
tptParam{1}.outlrThres = 5;

% Track Particles
[x, track] = funTPT(fileInfo, beadParam, tptParam);

save('resultsTPT.mat', 'x', 'track');
