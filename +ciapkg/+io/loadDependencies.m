function loadDependencies(varargin)
	% Download and load CIAtah dependencies.
	% Biafra Ahanonu
	% started: 2014.07.31
		% 2021.02.01 [15:09:46] branched from CIAtah
		% branch from calciumImagingAnalysis 2020.05.07 [15:47:29]
	% inputs
		%
	% outputs
		%

	% changelog
		% 2020.05.12 [17:40:37] - Updated to enable GUI-less loading of dependencies. In particular for easier unit testing.
		% 2020.06.28 [14:25:04] - Added ability for users to force update.
		% 2021.01.22 [13:42:36] - NWB from specific release to reduce compatibility errors.
		% 2021.02.01 [15:10:41] - Separated into non-class function for use in more functions without needing to load CIAtah class.
		% 2021.02.01 [‏‎15:19:40] - Update `_external_programs` to call ciapkg.getDirExternalPrograms() to standardize call across all functions.
		% 2021.03.20 [18:12:20] - Added EXTRACT support to list of functions to download.
		% 2021.06.19 [23:46:58] - Switched to support for original MIJ calling of ImageJ using just jar files, easier compatibility across MATLAB versions and OSes.
	% TODO
		% Verify all dependencies download and if not ask user to download again.

	%========================
	% DESCRIPTION
	options.externalProgramsDir = ciapkg.getDirExternalPrograms();
	options.guiEnabled = 1;
	options.dependencyStr = {'downloadImageJ','downloadCnmfGithubRepositories','example_downloadTestData','downloadNeuroDataWithoutBorders','downloadEXTRACT','downloadBioFormats','downloadMiji','loadMiji'};

	options.dispStr = {'Download ImageJ','Download CNMF, CNMF-E, and CVX code.','Download test one- and two photon datasets.','Download NWB (NeuroDataWithoutBorders)','Download EXTRACT','Download Bio-Formats','Download Fiji (to run Miji)','Load Fiji/Miji into MATLAB path.'};
	% Int vector: index of options.dependencyStr to run by default with no GUI
	options.depIdxArray = [1 2 3 4 5 6];
	% Binary: 1 = force update even if already downloaded. 0 = skip if already downloaded
	options.forceUpdate = 0;
	% get options
	options = getOptions(options,varargin);
	% display(options)
	% unpack options into current workspace
	% fn=fieldnames(options);
	% for i=1:length(fn)
	% 	eval([fn{i} '=options.' fn{i} ';']);
	% end
	%========================

	scnsize = get(0,'ScreenSize');
	if ischar(options.dispStr)
		options.dispStr = {options.dispStr};
	end
	if ischar(options.dependencyStr)
		options.dependencyStr = {options.dependencyStr};
	end
	dependencyStr = options.dependencyStr;

	dispStr = options.dispStr;
	if options.guiEnabled==1
		[depIdxArray, ~] = listdlg('ListString',dispStr,'ListSize',[scnsize(3)*0.3 scnsize(4)*0.3],'Name','Which dependencies to load? (Can select multiple)','InitialValue',options.depIdxArray);

		forceDownloadVec = [0 1];
		[forceUpdate, ~] = listdlg('ListString',{'No - skip installing dependency if already available.','Yes - force update to most recent version of dependency.'},'ListSize',[scnsize(3)*0.3 scnsize(4)*0.3],'Name','Force download/update? (e.g. "Yes" to update dependencies)','InitialValue',[1]);
		forceUpdate = forceDownloadVec(forceUpdate);
	else
		depIdxArray = options.depIdxArray;
		forceUpdate = 0;
	end
	analysisTypeD = dependencyStr(depIdxArray);
	dispStr = dispStr(depIdxArray);
	for depNo = 1:length(depIdxArray)
		disp([10 repmat('>',1,42)])
		disp(dispStr{depNo})
		switch analysisTypeD{depNo}
			case 'downloadCnmfGithubRepositories'
				[success] = downloadCnmfGithubRepositories('forceUpdate',forceUpdate);
			case 'downloadMiji'
				depStr = {'Save Fiji to default directory','Save Fiji to custom directory'};
				if options.guiEnabled==1
					[depIdxArray, ~] = listdlg('ListString',depStr,'ListSize',[scnsize(3)*0.2 scnsize(4)*0.25],'Name','Where to save Fiji?');
				else
					depIdxArray = 1;
				end
				depStr = depStr{depIdxArray};
				if depIdxArray==1
					downloadMiji();
				else
					downloadMiji('defaultDir','');
				end
				% if exist('pathtoMiji','var')
				% end
			case 'loadMiji'
				modelAddOutsideDependencies('miji');
			case 'example_downloadTestData'
				example_downloadTestData();
			case 'downloadCellExtraction'
				optionsH.forceUpdate = forceUpdate;
				optionsH.signalExtractionDir = options.externalProgramsDir;
				optionsH.gitNameDisp = {'cellmax_clean','extract'};
				optionsH.gitRepos = {'https://github.com/schnitzer-lab/CELLMax_CLEAN','https://github.com/schnitzer-lab/EXTRACT'};
				optionsH.gitRepos = cellfun(@(x) [x '/archive/master.zip'],optionsH.gitRepos,'UniformOutput',false);
				optionsH.outputDir = optionsH.gitNameDisp;
				optionsH.gitName = cellfun(@(x) [x '-master'],optionsH.gitNameDisp,'UniformOutput',false);
				[success] = downloadGithubRepositories('options',optionsH);
			case 'downloadEXTRACT'
				optionsH.forceUpdate = forceUpdate;
				optionsH.signalExtractionDir = options.externalProgramsDir;
				optionsH.gitNameDisp = {'extract'};
				optionsH.gitRepos = {'https://github.com/schnitzer-lab/EXTRACT-public'};
				optionsH.gitRepos = cellfun(@(x) [x '/archive/master.zip'],optionsH.gitRepos,'UniformOutput',false);
				optionsH.outputDir = optionsH.gitNameDisp;
				% optionsH.gitName = cellfun(@(x) [x '-master'],optionsH.gitNameDisp,'UniformOutput',false);
				optionsH.gitName = {'EXTRACT-public-master'};
				[success] = downloadGithubRepositories('options',optionsH);
			case 'downloadNeuroDataWithoutBorders'
				optionsH.forceUpdate = forceUpdate;
				optionsH.signalExtractionDir = options.externalProgramsDir;
				optionsH.gitNameDisp = {'nwb_schnitzer_lab','yamlmatlab','matnwb'};
				optionsH.gitRepos = {'https://github.com/schnitzer-lab/nwb_schnitzer_lab','https://github.com/ewiger/yamlmatlab'};

				% 'https://github.com/NeurodataWithoutBorders/matnwb'
				optionsH.gitRepos = cellfun(@(x) [x '/archive/master.zip'],optionsH.gitRepos,'UniformOutput',false);
				optionsH.gitRepos = [optionsH.gitRepos 'https://github.com/NeurodataWithoutBorders/matnwb/archive/v2.2.5.3.zip'];
				optionsH.outputDir = optionsH.gitNameDisp;
				optionsH.gitName = cellfun(@(x) [x '-master'],optionsH.gitNameDisp,'UniformOutput',false);
				optionsH.gitName{end} = 'matnwb-2.2.5.3';
				[success] = downloadGithubRepositories('options',optionsH);

				% Add NWB folders to path.
				ciapkg.nwb.setupNwb;
				% obj.loadBatchFunctionFolders;
			case 'downloadBioFormats'
				optionsH.forceUpdate = forceUpdate;
				optionsH.signalExtractionDir = options.externalProgramsDir;
				optionsH.gitNameDisp = {'bfmatlab'};
				optionsH.gitRepos = {'https://downloads.openmicroscopy.org/bio-formats/6.6.1/artifacts/bfmatlab.zip'};
				optionsH.outputDir = optionsH.gitNameDisp;
				% optionsH.gitName = cellfun(@(x) [x '-master'],optionsH.gitNameDisp,'UniformOutput',false);
				optionsH.gitName = {'bfmatlab'};
				[success] = downloadGithubRepositories('options',optionsH);
			case 'downloadImageJ'
				% Download mij.jar and ij.ar.
				downloadFiles = {'http://bigwww.epfl.ch/sage/soft/mij/mij.jar','http://rsb.info.nih.gov/ij/upgrade/ij.jar'}
				downloadFileNames = {'mij.jar','ij.jar'};
				imagejPath = [options.externalProgramsDir filesep 'imagej'];
				ciapkg.io.mkdir(imagejPath);
				nFiles = length(downloadFiles);
				for i=1:nFiles
					rawSavePathDownload = [imagejPath filesep downloadFileNames{i}];
					downloadUrl = downloadFiles{i};
					if exist(rawSavePathDownload,'file')~=2|forceUpdate==1
						fprintf('Downloading %s file to %s\n',downloadUrl,rawSavePathDownload)
						websave(rawSavePathDownload,downloadUrl);
					else
						fprintf('Already downloaded %s\n',rawSavePathDownload)
					end
				end
				manageMiji('startStop','setupImageJ');
			otherwise
				% nothing
		end
	end
end