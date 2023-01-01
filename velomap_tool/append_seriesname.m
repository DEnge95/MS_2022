% this script takes a directory full of the output from the dicom_copy_tool
% and appends the series name to the folders. This way there is more of a
% description than just the folder number
%
% Alex Barker, 20181017
% Children's Hospital Colorado
% University of Colorado, Anschutz
% alexander.barker@ucdenver.edu

clear all
TF = true; % turn off dicominfo warnings about VR state (happens on Philips)


% get base folder and read all subfolders and files in base folder
path_orig = pwd;
path_base = uigetdir(...
    'v:\cv_mri\Aorta-4D_Flow\BAV_Calgary\Calgary_002_20180228\');
cd(path_base)
path_listing = dir(path_base);

% open each subfolder and read the first file to find the seriesname & rename subfolder
for i = 3:length(path_listing)
    % test if a numbered folder
    numtest = str2double(path_listing(i).name); %if NAN is not a number 
    if path_listing(i).isdir %&&  ~isnan(numtest)
        cd([path_listing(i).folder filesep path_listing(i).name])
        subdir_listing = dir;
        if ~subdir_listing(3).isdir %check that not folder
            headerinfo = dicominfo(subdir_listing(3).name,'UseDictionaryVR',TF);
            seriesnum  = headerinfo.SeriesNumber;
            if isfield(headerinfo,'SeriesDescription')
                seriesname = headerinfo.SeriesDescription;
            else
                seriesname = headerinfo.SequenceName;
            end
            seriesname  = regexprep(seriesname,'*|<|>|:','');    % get rid of invalid characters '<' '>'        
            seriesname  = regexprep(seriesname,' - | -|- |/|\','-'); % get rid of dash spaces or '/'
            seriesname  = regexprep(seriesname,' ','_');         % get rid of spaces
            subdir_newname = [num2str(seriesnum,'%03.0f') '_' seriesname];
            disp(['new name = ' subdir_newname])
            cd ..
%            movefile(path_listing(i).name,subdir_newname) % this is super slow, going to attempt a system call to rename
            system(['ren ' path_listing(i).name ' ' subdir_newname]); % attempt at system call to rename folder - much faster
        end
    end
end
cd(path_orig)
disp('completed appending series descriptions to subfolders in base folder:')
disp(path_base)