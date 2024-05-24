% Running calciumImagingAnalysis (CIAPKG) imaging analysis via the command line
% Biafra Ahanonu
% Below is an example `cacliumImagingAnalysis` pipeline using the command line for users that do not want to use the calciumImagingAnalysis class or want to create their own custom batch analyses.
	% It assumes you have already run `example_downloadTestData` to download the example test data.
	% It will also run cross-day matching at the end.
	% All sections marked "USER INTERFACE" indicate that a GUI will appear to view processed movie, conduct cell sorting, or other interface.
% Changelog
	% 2020.09.15 [19:54:14] - Use ciapkg.getDir() to make sure demo always calls correct path regardless of where user is pointing. Also make playMovie calls have titles to make clearer to new users and allow a GUI-less option.
	% 2020.09.23 [08:35:58] - Updated to add support for cross-session analysis and use ciapkg.demo.runPreprocessing() to process the other imaging sessions.
	% 2020.10.17 [19:30:01] - Update to use ciapkg.signal_extraction.runPcaIca for PCA-ICA to make easier for users to run in the future.
	% 2021.01.17 [21:38:55] - Updated to show detrend example
	% 2021.06.20 [16:04:42] - Added CNMF/CNMF-e and EXTRACT to cell extraction examples.
	% 2021.06.22 [09:20:06] - Updated to make NWB saving, etc. smoother.
	% 2021.08.10 [09:57:36] - Updated to handle CIAtah v4.0 switch to all functions inside ciapkg package.

% =================================================


% Initialize
guiEnabled = 0;
saveAnalysis = 0; % save to NWB
dspace_factor = 2;
dstime_factor = 4;

show_preprocess = 1;

run_cnmfe = 0;

inputDatasetName = '/1';
rawFileRegexp = 'concat';

toptions.cropCoords = [170 120 420 450]/dspace_factor; % speicify the crop_coord beforehead

% Setup folder paths
% analysisFolderPath = [ciapkg.getDir() filesep 'data' filesep '2014_04_01_p203_m19_check01'];
analysisFolderPathAll = {'D:\ChenLab_DATA_Backup\Hemin_MiniScope_Processed\m19\2023_08_18_m19_demo'};
   % 'D:\ChenLab_DATA_Backup\Hemin_MiniScope_Processed\m33\2024_03_23_m33_rtppbaseline03'};
    % 'D:\ChenLab_DATA_Backup\Hemin_MiniScope_Processed\m34\2024_03_21_m34_demo';
    %'D:\ChenLab_DATA_Backup\Hemin_MiniScope_Processed\m35\2024_03_21_m35_demo'};

for i = 1:length(analysisFolderPathAll)

analysisFolderPath = analysisFolderPathAll{i};
[~,folderName,~] = fileparts(analysisFolderPath);

% Setup NWB folder paths
% nwbFilePath = [analysisFolderPath filesep 'nwbFiles' filesep folderName];
% nwbFileFolderPath = [analysisFolderPath filesep 'nwbFiles'];

% Load CIAtah functions
 ciapkg.loadBatchFxns();

% =================================================
% Download test data, only a single session
% ciapkg.api.example_downloadTestData('downloadExtraFiles',0);

% =================================================
% Load movie to analyze
inputMoviePath = ciapkg.api.getFileList(analysisFolderPath,rawFileRegexp,'sortMethod','natural');
% inputMoviePath = [analysisFolderPath filesep 'concat_recording_20140401_180333.h5'];
inputMovie = ciapkg.api.loadMovieList(inputMoviePath,'inputDatasetName',inputDatasetName);


%inputMoviePath_gui_mc = {'D:\ChenLab_DATA_Backup\Hemin_MiniScope_Processed\m33\2024_03_21_m33_demo\gui_mc_2024_03_21_p000_m33_NULL000_dsSpace_spFltBfReg_treg_crop_1.h5'};
% inputMoviePath_gui_mc_dfof = {'D:\ChenLab_DATA_Backup\Hemin_MiniScope_Processed\m33\2024_03_21_m33_demo\2024_03_21_p000_m33_NULL000_dsSpace_spFltBfReg_treg_crop_fxFrms_dfof_1.h5'};
% 
%inputMovie_gui_mc = ciapkg.api.loadMovieList(inputMoviePath_gui_mc,'inputDatasetName',inputDatasetName);
% inputMovie_gui_mc_dfof = ciapkg.api.loadMovieList(inputMoviePath_gui_mc_dfof,'inputDatasetName',inputDatasetName);



% =================================================

% downsample the raw movie

inputMovie_ds = ciapkg.api.downsampleMovie(inputMovie,'downsampleDimension','space','downsampleFactor',dspace_factor);


% USER INTERFACE Visualize slice of the movie
% if guiEnabled==1
% 	ciapkg.api.playMovie(inputMovie(:,:,1:500),'extraTitleText','Raw movie');
% 	% Alternatively, visualize by entering the file path
% 	ciapkg.api.playMovie(inputMoviePath,'extraTitleText','Raw movie directly from file');
% end

% =================================================
% USER INTERFACE Downsample input movie if need to
% if guiEnabled==1
% 	inputMovieD = ciapkg.api.downsampleMovie(inputMovie,'downsampleDimension','space','downsampleFactor',4);
% 	ciapkg.api.playMovie(inputMovie,'extraMovie',inputMovieD,'extraTitleText','Raw movie vs. down-sampled movie');
% end

% =================================================
% Alternatively, if you have Inscopix ISXD files, downsample by reading segments from disk using.
% moviePath = 'PATH_TO_ISXD';
% opts.maxChunkSize = 5000; % Max chunk size in Mb to load into RAM.
% opts.downsampleFactor = 4; % How much to downsample original movie, set to 1 for no downsampling.
% ciapkg.api.convertInscopixIsxdToHdf5(moviePath,'options',opts);

% =================================================
% Remove stripes from movie if needed
% if guiEnabled==1
% 	% Show full filter sequence for one frame
% 	sopts.stripOrientation = 'both';
% 	sopts.meanFilterSize = 1;
% 	sopts.freqLowExclude = 10;
% 	sopts.bandpassType = 'highpass';
% 	ciapkg.api.removeStripsFromMovie(inputMovie(:,:,1),'options',sopts,'showImages',1);
% 	drawnow
% end

% Run on the entire movie
% inputMovie = ciapkg.api.removeStripsFromMovie(inputMovie,'options',sopts,'showImages',0);

% =================================================
% Detrend movie if needed (default linear trend), e.g. to compensate for bleaching over time.
% inputMovie = ciapkg.api.normalizeMovie(inputMovie,'normalizationType','detrend','detrendDegree',1);

% =================================================
% USER INTERFACE Get coordinates to crop from the user separately




% =================================================
% Pre-process movies
% Motion correction
toptions.RegisType = 1;
toptions.turboregRotation = 0;
toptions.removeEdges = 1;
toptions.pxToCrop = 60;
toptions.refFrame = 100;
%toptions.refFrameMatrix = single(inputMovie_ds(:,:,toptions.refFrame));

toptions.SmoothX = 10;
toptions.SmoothY = 10;

% Pre-motion correction
toptions.complementMatrix = 1;
toptions.meanSubtract = 1;
toptions.meanSubtractNormalize = 1;
toptions.normalizeType = 'matlabDisk';
toptions.normalizeFreqLow = 10;
toptions.normalizeFreqHigh = 200;

% Spatial filter
%toptions.normalizeBeforeRegister = 'divideByLowpass';
toptions.normalizeBeforeRegister = 'bandpass';

toptions.bandpassBeforeRegister = 0;
toptions.freqLow = 20;
toptions.freqHigh = 100;

% toptions.returnNormalizedMovie=0;

[inputMovie2, ~] = ciapkg.api.turboregMovie(inputMovie_ds,'options',toptions);

% Compare raw and motion corrected movies
% if guiEnabled==1
%	ciapkg.api.playMovie(inputMovie3,'extraMovie',inputMovie_gui_mc_dfof,'extraTitleText','Raw movie vs. motion-corrected movie');
% end


% Run dF/F
% inputMovie2 = inputMovie2 + 1;
toptions.dfofType = 'dfof';
inputMovie3 = ciapkg.api.dfofMovie(single(inputMovie2),'normalizeBeforeRegister','bandpass','dfofType','dfof');
% inputMovie_gui_mc_cmd_dfof = ciapkg.api.dfofMovie(inputMovie_gui_mc,'dfofType','dfof');


% Run temporal downsampling
inputMovie3 = ciapkg.api.downsampleMovie(inputMovie3,'downsampleDimension','time','downsampleFactor',4);

% USER INTERFACE Final check of movie before cell extraction
if show_preprocess
    ciapkg.api.playMovie(inputMovie3,'extraTitleText','Processed movie for cell extraction');
end

% 
if run_cnmfe
    outputMoviePath3 = [analysisFolderPath filesep folderName '_spFltBfReg_turboreg_crop_dfof_dsTime.h5'];
    ciapkg.api.saveMatrixToFile(inputMovie3,outputMoviePath3);
end

% =================================================
%% Run CNMF or CNMF-e cell extraction
if run_cnmfe


% specify cnmfe options

% ========================
% OVERALL
% turn on parallel
cnmfeOpts.nonCNMF.parallel = 1;
% Binary: 1 = run merging algorithms
cnmfeOpts.runMerge = 1;
% Binary: 1 = remove false positives using CNMF-E algorithm
cnmfeOpts.runRemoveFalsePositives = 1;
% ===COMPUTATION
% Float: GB, memory space you allow to use in MATLAB
cnmfeOpts.memory_size_to_use = 24; %
% Float: GB, memory space you allow to use in MATLAB
cnmfeOpts.memory_size_per_patch = 1.8; % 0.6
% Int vector: patch size in pixels
cnmfeOpts.patch_dims = [64, 64]; % [64, 64]
% Int: number of frames per batch, leave blank to turn off batch
cnmfeOpts.batch_frames = [];
% ===SPATIAL
% Int: pixel, gaussian width of a gaussian kernel for filtering the data. 0 means no filtering
cnmfeOpts.gSig = 4;
% Int: pixel, neuron diameter
cnmfeOpts.gSiz = 13;
% Int: spatial downsampling factor
cnmfeOpts.ssub = 1;
% Binary: movie has dendrites?
cnmfeOpts.with_dendrites = false;
% Int: expand kernel for HALS growing (default: 3) and expansion factor of ellipse (default: 3)
cnmfeOpts.updateA_bSiz = 5;
% Char: hals, hals_thresh, lars, nnls
cnmfeOpts.spatial_algorithm = 'hals_thresh';
% ===TEMPORAL
% Int: temporal downsampling factor
cnmfeOpts.tsub = 1;
% Int: frame rate
cnmfeOpts.Fs = 5;
% Float: minimum spike size. When the value is negative, the actual threshold is abs(smin)*noise level
cnmfeOpts.deconv_smin = -5;
% Int: maximum decay time (unit: frame);
cnmfeOpts.max_tau = 100;
% Int: detrending the slow fluctuation. usually 1 is fine (no detrending)
cnmfeOpts.nk = 3;

% ===BACKGROUND
% Char: model of the background {'ring', 'svd'(default), 'nmf'}
cnmfeOpts.bg_model = 'ring';
% Int: number of background sources for each patch (only be used in SVD and NMF model)
cnmfeOpts.nb = 1;
% Int: when the ring model used, it is the radius of the ring used in the background model. otherwise, it's just the width of the overlapping area
cnmfeOpts.ring_radius = 18;
% Int: downsample background for a faster speed
cnmfeOpts.bg_ssub = 2;

% ===MERGING
% Float: 0 to 1, thresholds for merging neurons; [spatial overlap ratio, temporal correlation of calcium traces, spike correlation]
cnmfeOpts.merge_thr = 0.60;
% Char: method for computing neuron distances {'mean', 'max'}
cnmfeOpts.method_dist = 'max';
% Int: minimum distances between two neurons. it is used together with merge_thr
cnmfeOpts.dmin = 5;
% Int: merge neurons if their distances are smaller than dmin_only.
cnmfeOpts.dmin_only = 2;
% Float vector: merge components with highly correlated spatial shapes (corr=0.8) and small temporal correlations (corr=0.1)
cnmfeOpts.merge_thr_spatial = [0.8, 0.5, -inf];

% ===INITIALIZATION
% Int: maximum number of neurons per patch. when K=[], take as many as possible.
cnmfeOpts.K = [];
% Float: minimum local correlation for a seeding pixel
cnmfeOpts.min_corr = 0.88;
% minimum peak-to-noise ratio for a seeding pixel
cnmfeOpts.min_pnr = 13;
cnmfeOpts.bd = 1;             % number of rows/columns to be ignored in the boundary (mainly for motion corrected data)
cnmfeOpts.use_parallel = true;    % use parallel computation for parallel computing
cnmfeOpts.show_init = true;   % show initialization results
cnmfeOpts.center_psf = true;  % set the value as true when the background fluctuation is large (usually 1p data)
% set the value as false when the background fluctuation is small (2p)
% ===Residual
% Float: 0 to 1, minimum local correlation for initializing a neuron (default: 0.3)
cnmfeOpts.min_corr_res = 0.90;
% Float: stands for minimum peak-to-noise ratio to look for a cell
cnmfeOpts.min_pnr_res = 15;
% Char: method for initializing neurons from the residual. 'auto' or 'manual'
cnmfeOpts.seed_method_res = 'auto';
% Binary: boolean, update noise level for each pixel
cnmfeOpts.update_sn = true;

%

% Run CNMF-e
[success] = ciapkg.api.cnmfVersionDirLoad('cnmfe');
%cnmfeOptions.gSiz = cellWidth;
%cnmfeOptions.gSig = ceil(cellWidth/4);
[cnmfeAnalysisOutput] = ciapkg.api.computeCnmfeSignalExtraction_batch(outputMoviePath3,'options',cnmfeOpts);
cnmfeOutputFileName = [analysisFolderPath filesep folderName '_cnmfeAnalysisOutput.mat'];

save(cnmfeOutputFileName, "cnmfeAnalysisOutput");

% Unload CNMF(-E) directories to keep namespace clean
[success] = ciapkg.api.cnmfVersionDirLoad('none');
end

end