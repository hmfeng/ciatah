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
%% Initialize
guiEnabled = 0;
saveAnalysis = 0; % save to NWB
dspace_factor = 2;
dstime_factor = 4;
run_cnmfe = 1;

inputDatasetName = '/1';
rawFileRegexp = 'concat';

% Setup folder paths
% analysisFolderPath = [ciapkg.getDir() filesep 'data' filesep '2014_04_01_p203_m19_check01'];
analysisFolderPath = ['D:\ChenLab_DATA_Backup\Hemin_MiniScope_Processed\m33\2024_03_21_m33_demo'];
    % 'D:\ChenLab_DATA_Backup\Hemin_MiniScope_Processed\m34\2024_03_21_m34_demo';
    %'D:\ChenLab_DATA_Backup\Hemin_MiniScope_Processed\m35\2024_03_21_m35_demo'};
[~,folderName,~] = fileparts(analysisFolderPath);

% Setup NWB folder paths
% nwbFilePath = [analysisFolderPath filesep 'nwbFiles' filesep folderName];
% nwbFileFolderPath = [analysisFolderPath filesep 'nwbFiles'];

% Load CIAtah functions
 ciapkg.loadBatchFxns();

% =================================================
%% Download test data, only a single session
% ciapkg.api.example_downloadTestData('downloadExtraFiles',0);

% =================================================
%% Load movie to analyze
inputMoviePath = ciapkg.api.getFileList(analysisFolderPath,rawFileRegexp,'sortMethod','natural');
% inputMoviePath = [analysisFolderPath filesep 'concat_recording_20140401_180333.h5'];
inputMovie = ciapkg.api.loadMovieList(inputMoviePath,'inputDatasetName',inputDatasetName);


%inputMoviePath_gui_mc = {'D:\ChenLab_DATA_Backup\Hemin_MiniScope_Processed\m33\2024_03_21_m33_demo\gui_mc_2024_03_21_p000_m33_NULL000_dsSpace_spFltBfReg_treg_crop_1.h5'};
inputMoviePath_gui_mc_dfof = {'D:\ChenLab_DATA_Backup\Hemin_MiniScope_Processed\m33\2024_03_21_m33_demo\2024_03_21_p000_m33_NULL000_dsSpace_spFltBfReg_treg_crop_fxFrms_dfof_1.h5'};
% 
%inputMovie_gui_mc = ciapkg.api.loadMovieList(inputMoviePath_gui_mc,'inputDatasetName',inputDatasetName);
inputMovie_gui_mc_dfof = ciapkg.api.loadMovieList(inputMoviePath_gui_mc_dfof,'inputDatasetName',inputDatasetName);



% =================================================

%% downsample the raw movie

inputMovie_ds = ciapkg.api.downsampleMovie(inputMovie,'downsampleDimension','space','downsampleFactor',dspace_factor);


%% USER INTERFACE Visualize slice of the movie
% if guiEnabled==1
% 	ciapkg.api.playMovie(inputMovie(:,:,1:500),'extraTitleText','Raw movie');
% 	% Alternatively, visualize by entering the file path
% 	ciapkg.api.playMovie(inputMoviePath,'extraTitleText','Raw movie directly from file');
% end

% =================================================
%% USER INTERFACE Downsample input movie if need to
% if guiEnabled==1
% 	inputMovieD = ciapkg.api.downsampleMovie(inputMovie,'downsampleDimension','space','downsampleFactor',4);
% 	ciapkg.api.playMovie(inputMovie,'extraMovie',inputMovieD,'extraTitleText','Raw movie vs. down-sampled movie');
% end

% =================================================
%% Alternatively, if you have Inscopix ISXD files, downsample by reading segments from disk using.
% moviePath = 'PATH_TO_ISXD';
% opts.maxChunkSize = 5000; % Max chunk size in Mb to load into RAM.
% opts.downsampleFactor = 4; % How much to downsample original movie, set to 1 for no downsampling.
% ciapkg.api.convertInscopixIsxdToHdf5(moviePath,'options',opts);

% =================================================
%% Remove stripes from movie if needed
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
%% Detrend movie if needed (default linear trend), e.g. to compensate for bleaching over time.
% inputMovie = ciapkg.api.normalizeMovie(inputMovie,'normalizationType','detrend','detrendDegree',1);

% =================================================
%% USER INTERFACE Get coordinates to crop from the user separately
if guiEnabled==1
	[cropCoords] = ciapkg.api.getCropCoords(squeeze(inputMovie(:,:,1)));
	toptions.cropCoords = cropCoords;
	% Or have turboreg function itself directly ask the user for manual area from which to obtain correction coordinates
	% toptions.cropCoords = 'manual';
else
	%toptions.cropCoords = [26 34 212 188];
    % toptions.cropCoords = [160 160 480 460]/dspace_factor; % [x_left, y_top, x_right, y_bottom]
    toptions.cropCoords = [87 78 259 228];

end

% =================================================
%% Pre-process movies
% Motion correction
toptions.RegisType = 1;
toptions.turboregRotation = 0;
toptions.removeEdges = 1;
toptions.pxToCrop = 4;
toptions.refFrame = 100;
%toptions.refFrameMatrix = single(inputMovie_ds(:,:,toptions.refFrame));

toptions.SmoothX = 10;
toptions.SmoothY = 10;

% Pre-motion correction
toptions.complementMatrix = 1;
toptions.meanSubtract = 1;
toptions.meanSubtractNormalize = 1;
toptions.normalizeType = 'matlabDisk';
toptions.normalizeFreqLow = 40;
toptions.normalizeFreqHigh = 100;

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
	ciapkg.api.playMovie(inputMovie3,'extraMovie',inputMovie_gui_mc_dfof,'extraTitleText','Raw movie vs. motion-corrected movie');
% end



%% Run dF/F
% inputMovie2 = inputMovie2 + 1;
toptions.dfofType = 'dfof';
inputMovie3 = ciapkg.api.dfofMovie(single(inputMovie2),'normalizeBeforeRegister','bandpass','dfofType','dfof');
% inputMovie_gui_mc_cmd_dfof = ciapkg.api.dfofMovie(inputMovie_gui_mc,'dfofType','dfof');


%% Run temporal downsampling
inputMovie3 = ciapkg.api.downsampleMovie(inputMovie3,'downsampleDimension','time','downsampleFactor',4);

%% USER INTERFACE Final check of movie before cell extraction
% if guiEnabled==1
	ciapkg.api.playMovie(inputMovie3,'extraTitleText','Processed movie for cell extraction');
% end

outputMoviePath = [analysisFolderPath filesep folderName '_spatialFiltBfReg_turboreg_crop_dfof_downsampleTime.h5'];

% outputMoviePath2 = [analysisFolderPath filesep folderName '_spatialFiltBfReg_turboreg_crop.h5'];
outputMoviePath3 = [analysisFolderPath filesep folderName '_spatialFiltBfReg_turboreg_crop_dfof_revised.h5'];


% ciapkg.api.saveMatrixToFile(inputMovie3,outputMoviePath3);

% =================================================
%% Run PCA-ICA cell extraction
% nPCs = 300;
% nICs = 225;
% pcaicaStruct = ciapkg.signal_extraction.runPcaIca(inputMovie3,nPCs,nICs,'version',2,'output_units','fl','mu',0.1,'term_tol',5e-6,'max_iter',1e3);

%% Save outputs to NWB format
% mkdir([analysisFolderPath filesep 'nwbFiles']);
% if saveAnalysis==1
% 	ciapkg.api.saveNeurodataWithoutBorders(pcaicaStruct.IcaFilters,{pcaicaStruct.IcaTraces},'pcaica',[nwbFilePath '_pcaicaAnalysis.nwb']);
% end

% =================================================
%% Run CNMF or CNMF-e cell extraction
if run_cnmfe

    numExpectedComponents = 225;
    cellWidth = 10;
    cnmfOptions.otherCNMF.tau = cellWidth/2; % expected width of cells

    % Run CNMF
    % [success] = ciapkg.api.cnmfVersionDirLoad('current');
    % [cnmfAnalysisOutput] = ciapkg.api.computeCnmfSignalExtractionClass(inputMovie3,numExpectedComponents,'options',cnmfOptions);

    % Run CNMF-e
    [success] = ciapkg.api.cnmfVersionDirLoad('cnmfe');
    cnmfeOptions.gSiz = cellWidth;
    cnmfeOptions.gSig = ceil(cellWidth/4);
    [cnmfeAnalysisOutput] = ciapkg.api.computeCnmfeSignalExtraction_batch(outputMoviePath3,'options',cnmfeOptions);

    % Save outputs to NWB format
    % if saveAnalysis==1
    % 	% Save CNMF
    % 	ciapkg.api.saveNeurodataWithoutBorders(cnmfAnalysisOutput.extractedImages,{cnmfAnalysisOutput.extractedSignals,cnmfAnalysisOutput.extractedSignalsEst},'cnmf',[nwbFilePath '_cnmf.nwb']);
    %
    % 	% Save CNMF-E
    % 	ciapkg.api.saveNeurodataWithoutBorders(cnmfeAnalysisOutput.extractedImages,{cnmfeAnalysisOutput.extractedSignals,cnmfeAnalysisOutput.extractedSignalsEst},'cnmfe',[nwbFilePath '_cnmfe.nwb']);
    % end

    % Unload CNMF(-E) directories to keep namespace clean
    [success] = ciapkg.api.cnmfVersionDirLoad('none');

end

% =================================================
%% Run EXTRACT cell extraction. Check each function with "edit" for options.
% Load default configuration
% ciapkg.loadBatchFxns('loadEverything');
% extractConfig = get_defaults([]);
% 
% % See https://github.com/schnitzer-lab/EXTRACT-public#configurations.
% cellWidth = 10;
% extractConfig.avg_cell_radius = cellWidth;
% extractConfig.num_partitions_x = 2;
% extractConfig.num_partitions_y = 2;
% extractConfig.use_sparse_arrays = 0;
% 
% outStruct = extractor(inputMovie3,extractConfig);
% 
% % Grab outputs and put into standard format
% extractAnalysisOutput.filters = outStruct.spatial_weights;
% % permute so it is [nCells frames]
% extractAnalysisOutput.traces = permute(outStruct.temporal_weights, [2 1]);
% 
% % Other run information if saving as a MAT-file.
% extractAnalysisOutput.info = outStruct.info;
% extractAnalysisOutput.config = outStruct.config;
% extractAnalysisOutput.info = outStruct.info;
% extractAnalysisOutput.userInputConfig = extractConfig;
% extractAnalysisOutput.opts = outStruct.config;
% 
% % Save outputs to NWB format
% if saveAnalysis==1
% 	ciapkg.api.saveNeurodataWithoutBorders(extractAnalysisOutput.filters,{extractAnalysisOutput.traces},'extract',[nwbFilePath '_extract.nwb']);
% end
% 
% % Remove EXTRACT from the path.
% ciapkg.loadBatchFxns();
% 
% % =================================================
% %% USER INTERFACE Run cell sorting using matrix outputs from cell extraction.
% if guiEnabled==1
% 	[outImages, outSignals, choices] = ciapkg.api.signalSorter(pcaicaStruct.IcaFilters,pcaicaStruct.IcaTraces,'inputMovie',inputMovie3);
% 
% 	% Plot results of sorting
% 	figure;
% 	subplot(1,2,1);imagesc(max(IcaFilters,[],3));axis equal tight; title('Raw filters')
% 	subplot(1,2,2);imagesc(max(outImages,[],3));axis equal tight; title('Sorted filters')
% end
% 
% % =================================================
% %% Run signal sorting using NWB
% [outImages, outSignals, choices] = ciapkg.api.signalSorter([nwbFilePath '_pcaicaAnalysis.nwb'],[],'inputMovie',inputMovie3);
% 
% %% USER INTERFACE Run signal sorting using multiple NWB files from cell extraction.
% if saveAnalysis==1&guiEnabled==1
% 	disp(repmat('=',1,21));disp('Running signalSorter using NWB file input.')
% 	nwbFileList = ciapkg.api.getFileList(nwbFileFolderPath,'.nwb');
% 	if ~isempty(nwbFileList)
% 		nFiles = length(nwbFileList);
% 		outImages = {};
% 		outSignals = {};
% 		choices = {};
% 		for fileNo = 1:nFiles
% 			[outImages{fileNo}, outSignals{fileNo}, choices{fileNo}] = ciapkg.api.signalSorter(nwbFileList{fileNo},[],'inputMovie',inputMovie3);
% 		end
% 
% 		% Plot results of sorting
% 		for fileNo = 1:nFiles
% 			try
% 				[inputImagesTmp,inputSignalsTmp,infoStructTmp,algorithmStrTmp,inputSignals2Tmp] = ciapkg.api.loadSignalExtraction(nwbFileList{fileNo});
% 				figure;
% 				subplot(1,2,1); 
% 					imagesc(max(inputImagesTmp,[],3));
% 					axis equal tight; 
% 					title([algorithmStrTmp ' | Raw filters'])
% 				subplot(1,2,2); 
% 					imagesc(max(outImages{fileNo},[],3));
% 					axis equal tight; 
% 					title('Sorted filters')
% 			catch err
% 				disp(repmat('@',1,7))
% 				disp(getReport(err,'extended','hyperlinks','on'));
% 				disp(repmat('@',1,7))
% 			end
% 		end
% 	end
% end
% 
% % =================================================
% %% USER INTERFACE Create an overlay of extraction outputs on the movie and signal-based movie
% [inputMovieO] = ciapkg.api.createImageOutlineOnMovie(inputMovie3,IcaFilters,'dilateOutlinesFactor',0);
% if guiEnabled==1
% 	ciapkg.api.playMovie(inputMovieO,'extraMovie',inputMovie3,'extraTitleText','Overlay of cell outlines on processed movie');
% end
% 
% [signalMovie] = ciapkg.api.createSignalBasedMovie(IcaTraces,IcaFilters,'signalType','peak');
% if guiEnabled==1
% 	ciapkg.api.playMovie(signalMovie,'extraMovie',inputMovie3,'extraTitleText','Cell activity-based movie');
% end
% 
% movieM = cellfun(@(x) ciapkg.api.normalizeVector(x,'normRange','zeroToOne'),{inputMovie3,inputMovieO,signalMovie},'UniformOutput',false);
% ciapkg.api.playMovie(cat(2,movieM{:}));
% 
% % =================================================
% %% Run pre-processing on 3 batch movies then do cross-session alignment
% batchMovieList = {...
% [ciapkg.getDir() filesep 'data' filesep 'batch' filesep '2014_08_05_p104_m19_PAV08'],...
% [ciapkg.getDir() filesep 'data' filesep 'batch' filesep '2014_08_06_p104_m19_PAV09'],...
% [ciapkg.getDir() filesep 'data' filesep 'batch' filesep '2014_08_07_p104_m19_PAV10']...
% };
% 
% % USER INTERFACE Get the motion correction crop coordinates
% cropCoordsCell = {};
% nFolders = length(batchMovieList);
% for folderNo = 1:nFolders
% 	analysisFolderPath = batchMovieList{folderNo};
% 	inputMoviePath = ciapkg.api.getFileList(analysisFolderPath,rawFileRegexp,'sortMethod','natural');
% 	% inputMoviePath = [analysisFolderPath filesep 'concat_recording_20140401_180333.h5'];
% 	inputMovie = ciapkg.api.loadMovieList(inputMoviePath,'inputDatasetName',inputDatasetName,'frameList',1:2);
% 
% 	[cropCoords] = ciapkg.api.getCropCoords(squeeze(inputMovie(:,:,1)));
% 	% toptions.cropCoords = cropCoords;
% 	cropCoordsCell{folderNo} = cropCoords;
% end
% 
% %% Run pre-processing on each of the movies.
% procMovieCell = cell([1 nFolders]);
% for folderNo = 1:nFolders
% 	inputMoviePath = ciapkg.api.getFileList(analysisFolderPath,rawFileRegexp,'sortMethod','natural');
% 	inputMovie = ciapkg.api.loadMovieList(inputMoviePath,'inputDatasetName',inputDatasetName,'frameList',[]);
% 	procOpts.motionCorrectionCropCoords = cropCoordsCell{folderNo};
% 	procOpts.dfofMovie = 1;
% 	procOpts.motionCorrectionFlag = 1;
% 	procOpts.normalizeMovieFlag = 1;
% 	procOpts.normalizeType = 'divideByLowpass';
% 	procOpts.freqLow = 0;
% 	procOpts.freqHigh = 7;
% 	procOpts.downsampleTimeFactor = 4;
% 	procMovieCell{folderNo} = ciapkg.demo.runPreprocessing(inputMovie,'options',procOpts);
% end
% disp('Done with pre-processing!')
% 
% %% Run cell-extraction on the movies
% pcaicaStructCell = cell([1 nFolders]);
% nPCs = 300;
% nICs = 225;
% for folderNo = 1:nFolders
% 	inputMoviePath = ciapkg.api.getFileList(analysisFolderPath,rawFileRegexp,'sortMethod','natural');
% 	pcaicaStruct{folderNo} = ciapkg.signal_extraction.runPcaIca(procMovieCell{folderNo},nPCs,nICs,'version',2,'outputUnits','fl','mu',0.1,'term_tol',5e-6,'max_iter',1e3);
% end
% disp('Done with PCA-ICA analysis pre-processing!')
% 
% %% Run cross-session alignment of cells
% % Create input images, cell array of [x y nCells] matrices
% inputImages = cellfun(@(x) x.IcaFilters,pcaicaStruct,'UniformOutput',false);
% 
% % options to change
% opts.maxDistance = 5; % distance in pixels between centroids for them to be grouped
% opts.trialToAlign = 1; % which session to start alignment on
% opts.nCorrections = 1; %number of rounds to register session cell maps.
% opts.RegisTypeFinal = 2; % 3 = rotation/translation and iso scaling; 2 = rotation/translation, no iso scaling
% 
% % Run alignment code
% [alignmentStruct] = ciapkg.api.matchObjBtwnTrials(inputImages,'options',opts);
% 
% % Global IDs is a matrix of [globalID sessionID]
% % Each (globalID, sessionID) pair gives the within session ID for that particular global ID
% globalIDs = alignmentStruct.globalIDs;
% 
% % View the cross-session matched cells, saved to `private\_tmpFiles` sub-folder.
% [success] = ciapkg.api.createMatchObjBtwnTrialsMaps(inputImages,alignmentStruct);
% 
% %% Display cross-session matching movies
% disp('Playing movie frames')
% crossSessionMovie1 = [ciapkg.getDir filesep 'private' filesep '_tmpFiles' filesep 'matchObjColorMap50percentMatchedSession_matchedCells.avi'];
% crossSessionMovie2 = [ciapkg.getDir filesep 'private' filesep '_tmpFiles' filesep 'matchObjColorMapAllMatchedSession_matchedCells.avi'];
% ciapkg.api.playMovie(crossSessionMovie1,'extraMovie',crossSessionMovie2,'rgbDisplay',1);