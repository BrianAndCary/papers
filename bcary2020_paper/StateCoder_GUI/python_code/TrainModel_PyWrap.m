function [pyMdl_path,OOBerr,saveMdl_path] = TrainModel_PyWrap(train_set,anim,mlmode,nTrees,python_dir,model_dir,pywrap_path,temp_path)
% This is the function that takes the user-coded training data set and
% trains a ML model with it.
% INPUTS:
% - train_set: this is the training data set
% - anim: animal ID name
% - mlmode: string specifying which algorithm to use. Current options:
%           * 'RF' : Random Forest
%           * 'ECOC' : error-correcting output codes
% - nTrees: number of decision trees to use for RF model training
%
% Author: Alejandro Torrado Pacheco, June 2017, Edited by Brian Cary

%% SETUP;

% set up default arguments
if nargin < 4
    nTrees = 200;
end
if isempty(nTrees)
    nTrees = 200;
end


base_dir = model_dir;
pywrap_dir = pywrap_path;
python_path = python_dir;
tempdir = temp_path;

% make directory for this animal where trained model will be saved
mdl_dir = [base_dir filesep mlmode filesep anim];
if ~exist(mdl_dir,'dir'), mkdir(mdl_dir); end

fprintf('Training %s model.\n',mlmode);

% check which ML algorithm to use
switch mlmode
    
    case 'RF'
        %% RANDOM FOREST
        % first save file to be able to read it into python
        training = train_set(:,1:end-1);
        scoredclass = train_set(:,end);
        homefolder = getuserdir;
        if isempty(homefolder)
            homefolder = 'C:\Users\Admin\';
        end
%         tempdir = [homefolder filesep 'MATLAB' filesep 'tempPyVidcodeFiles'];
        if ~exist(tempdir, 'dir'), mkdir(tempdir); end
        tempfile_code = ['temp_train_file_',anim,'.mat'];
        tempfile = fullfile(tempdir,tempfile_code);
        if exist(tempfile,'file'), delete(tempfile); end
        save(tempfile,'training','scoredclass');
        
        rf_t0 = tic;
        
        
        %% here call Python code instead of Matlab treebagger
        % after call to code, load OOB error from OOB_error.mat in TEMP
        % folder
        
        sys_py_call = [python_path ' ' pywrap_dir filesep ...
            'TrainModel_Python.py -f ' tempfile ' -nt ' num2str(nTrees)...
            ' -sd ' tempdir ' -sn ' anim];

        py_train_status = system(sys_py_call);
        OOB_file_code = 'OOB_error.mat';
        try
        oob_file_load = load([tempdir filesep OOB_file_code]);
        catch
            keyboard;
        end
        OOB_Error = oob_file_load.OOB_error;
        
        rf_t1 = toc(rf_t0);
        %fprintf('Time elapsed: %.2f.\n\n',rf_t1);
        
    case 'ECOC'
        %% ERROR CORRECTING OUTPUT CODE
        fprinf('ECOC not available yet. Sorry! Try RF.\n');
        
end

pyMdl_path = [tempdir filesep 'pyMdl.pkl'];
OOBerr = OOB_Error;
saveMdl_path = [mdl_dir filesep mlmode '_pyMdl.mat'];

temp_dir = dir(tempdir);
savepath = [mdl_dir filesep];
disp(['copying files from temp dir to:',savepath])
for ii = 1:length(temp_dir)
    if temp_dir(ii).bytes > 5
        copyfile([temp_dir(ii).folder filesep temp_dir(ii).name],savepath)
    end
end
        
        