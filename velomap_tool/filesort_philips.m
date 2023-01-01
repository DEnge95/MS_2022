function [outputArg1] = filesort_philips(file_root)
%
% FILESORT_PHILIPS This function sorts files Philips files as exported by the scanner 
% to be loaded into the velomap_tool via the mag/flow folder paradigm
% 
% Note: The files must have already gone through the dicom_copy_tool.m and
% append_seriesname.m
%
% file_root is the folder containing only one FH, RL, and AP folder. If
% there are more than one, then the script will ask you to fix it.
%
% outputArg1 = file_root;
%
% TODO: 
% - modify velomap_internal so that temporal resolution is correct!
%   o current implmentation assumes TR is the temporal resolution, this is
%   incorrect for Philips
%   o am not sure where temporal resolution is stored, one approach is to
%   find the final image and compute reslution based on trigger time and
%   known number of cardiac phases
% - move extra non-image dicom files to quarentine folder?
%   o these often appear as negative image numbers in the filenames
%   o they can also be identified by the drastically different file sizes
%   and the fact that 'Image Type' is not a dicom header field
%
% Alex Barker, 20181017
% Children's Hospital Colorado
% University of Colorado, Anschutz
% alexander.barker@ucdenver.edu
% end

%% Identify files to manipulate

% identify the root folder with the 4D flow images
if nargin < 1
    file_root = uigetdir('c:\_Research_Data\michal_processing\apodaca\_03_filesorting\','Select base folder with Philips 4D flow data (must have ''RL'', ''FH'', or ''AP'' directions identified somewhere in folder names)');
end

% identify folder directions (note must already be sorted with series
% description via dicom_copy and append_series name)
if file_root == 0
    disp('User canceled the program')
    return
else
    dir_temp = dir(file_root);
    cd(file_root)
end

% look for RL FH and AP folder (assumption here is that description has
% RL FH AP directions)
dir_names = {dir_temp(:).name}'; % this assume same char lengths for each folder! dangerous esp if series increment increases, i.e. 9->10
dir_FH    = regexpi(dir_names,'FH');
dir_RL    = regexpi(dir_names,'RL');
dir_AP    = regexpi(dir_names,'AP');

% break if there are more than one FH RL AP folders (does not distinguish)
if (sum(~cellfun(@isempty,dir_RL))+sum(~cellfun(@isempty,dir_FH))+sum(~cellfun(@isempty,dir_AP))) ~= 3
    disp('More than one folder was detected to contain a RL, FH, or AP direction.')
    disp(' ')
    disp([dir_names{~cellfun(@isempty,dir_RL)} ' -> RL direction' ]) % check RL folder
    disp([dir_names{~cellfun(@isempty,dir_FH)} ' -> FH direction' ]) % check FH folder
    disp([dir_names{~cellfun(@isempty,dir_AP)} ' -> AP direction' ]) % check AP folder
    disp(' ')
    disp('More than one encode direction folder detected. Consider ')
    disp('manually moving the FH, AP, &  RL folders of interestinto a')
    disp('dedicated 3dpc processing folder. Make sure to only have only ')
    disp('one of each encode direction in the base folder. Note that a')
    disp('common problem is that mag images are assummed to be in the FH')
    disp('folder, not in a standalone folder.')
    return % break if there are more than one folder in each encode direction
end

% ask if folders are correct?

disp('Check that encode direction for each folder is correct:') 
disp(' ') 
disp([dir_names{~cellfun(@isempty,dir_RL)} ' -> RL direction' ]) % check RL folder
disp([dir_names{~cellfun(@isempty,dir_FH)} ' -> FH direction' ]) % check FH folder
disp([dir_names{~cellfun(@isempty,dir_AP)} ' -> AP direction' ]) % check AP folder
disp(' ') 
prompt = 'Are these correct? [Y]/N: ';
str = input(prompt,'s');
if isempty(str)
    str = 'Y';
end
if regexpi(str,'N')
    return
end

%% break out FH magnitude images into their own folder for reading by velomap_tool
% - TODO: don't assume that FH folder contains magnitude data
% - TODO: add in flexibility to find PCA images via header (which will add
% more time to the process, for now, will only mention in code that the
% easy way around this is to break out the magnitude into a seperate
% folder)
%    * to get fancy, note that dicom header Image Type contains whether mag or phase:
%    * (0008,0008) IMAGE TYPE: ORIGINAL\PRIMARY\PHASE CONTRAST M\P\PCA is for phase iamges
%    * (0008,0008) IMAGE TYPE: CS ORIGINAL\PRIMARY\M_FFE\M\FFE are the mag images

% simply measure number of files in FH folder and split in two, where first
% set of files is assumed to be mag data

dir_name_FH = dir_names{~cellfun(@isempty,dir_FH)};
files       = dir(dir_name_FH);
file_names  = cell2mat({files(3:end).name}');
num_files   = size(file_names,1);

% test that divisible by 2 and separate mag files
if mod(num_files,2) ~= 0
    disp('ERROR: Number of FH files is not divisible by two, check files')
    disp('for validity...watch out for non-image dicoms, they will appear')
    disp('as significantly different size files than the other dcm files')
    return
else
    disp('Processing data...')
    %identify mag files (assume they occur in the first half of the lot)
    num_files_half = num_files./2;
    mkdir([dir_temp(1).folder filesep dir_name_FH '_mag'])
    for i = 1:num_files_half
        movefile([dir_temp(1).folder filesep dir_name_FH filesep file_names(i,:)],...
                 [dir_temp(1).folder filesep dir_name_FH '_mag' filesep file_names(i,:)])
    end
    disp(['  > FH magnitude files moved to: ' dir_name_FH '_mag'])
    num_files = num_files_half; % this is needed for the later code (now that all folders have the same number of files)
end

clear files file_names dir_name_FH num_files_half dir_temp dir_names

%% recursively perform the filesort operation on each folder

dir_temp = dir(file_root);
dir_names = {dir_temp(:).name}'; % this assume same char lengths for each folder! dangerous esp if series increment increases, i.e. 9->10
dir_names = dir_names(~contains(dir_names,'.')); % TF modified for MAC

for i = 1 : size(dir_names) %assumes recursing RL, FH, FH_mag, AP (in no specific order)
    file_dir = dir_names{i,:};
    files      = dir(file_dir);
    file_names = cell2mat({files(3:end).name}');

    % pass to sorting algorithm
    disp(' ')
    disp(['reading and determining reordering sequence for: ' file_dir])
    % see file remap_formula.xlsx for the philips -> siemens decoder ring used
    % in filenum_remap
    [out_dir,out_filelist] = filenum_remap([file_dir],file_names); 
    disp('  > done sort sequencing')

    % to put into mag/flow paradigm, need to change out_dir, out_filelist for flow folder numbering
    % determine direction from dir name
    if regexp(file_dir,'FH') & regexp(file_dir,'mag')
        out_dir = 'mag';
    elseif regexp(file_dir,'FH') % this should not detect the 'mag' FH folder bc of order, may want to test
        out_dir = 'flow';
        % remember to assign filenames in order of FH-AP-RL
        % since already moved the mag files from FH we can simply subtract
        % from the filenumber by the number of files to get 1 to num_files
        ind = strfind(out_filelist(1,:),'.');
        in_filenum  = str2num(out_filelist(:,ind-4:ind-1));
        out_filenum = in_filenum - num_files; % because these occur after the mag numbering
    elseif regexp(file_dir,'AP')
        out_dir = 'flow';
        ind = strfind(out_filelist(1,:),'.');
        in_filenum  = str2num(out_filelist(:,ind-4:ind-1));
        out_filenum = in_filenum + num_files; % because these occur after FH phase numbering
    elseif regexp(file_dir,'RL')
        out_dir = 'flow';
        ind = strfind(out_filelist(1,:),'.');
        in_filenum  = str2num(out_filelist(:,ind-4:ind-1));
        out_filenum = in_filenum + 2.*num_files; % because these occur after AP phase numbering
    end    
    
    % update out_filelist character string
    out_file     = num2str(out_filenum,'%04.0f');
    temp_prefix  = repmat('IM-0001-',num_files,1);
    temp_ext     = repmat('.dcm',num_files,1);
    out_filelist = [temp_prefix out_file temp_ext];
    
    %% Rename files consecutively

    num_files     = size(out_filelist,1);

    % create new folder and filenames
    [status,msg] = mkdir(out_dir);
    disp(['  > writing new files to the folder: ' out_dir]) 

    %repeat file directory paths and create char matrices for copy operation
    file_dirrep    = repmat([file_dir '/'],num_files,1);
    file_dirnewrep = repmat([out_dir  '/'],num_files,1);

    % have to do a loop because i don't think (?) movefile likes multiple files
    file_old = [file_dirrep file_names];
    file_new = [file_dirnewrep out_filelist];

    for j = 1:num_files
        copyfile(file_old(j,:),file_new(j,:)) % perhaps do system call instead?
    end
    disp(['  ! done rewriting file sequence order for: ' file_dir]) 
end
outputArg1 = file_root;

%% debug: load files and view order by filename (note radiant will not show correct order)

% have to add series number to file names
% disp('  > loading data for viewing') 
% D = LoadMRImovie(out_dir);
% disp_cine(D) 
end

function [out_dir, out_filelist] = filenum_remap(in_dir, in_filelist)
% FILENAME_REMAP This function takes a Philips filestructure (as output by
% the sorting function in the super_tool) and remaps it to an order that
% we've used for the velomap_tool. I.e. it takes something that outputs all
% cardiac phases, and then slices, to the opposite
%
% *in_dir is the input directory
% *input_filelist is the expected filelist to be remapped (it expects a
%  'file_num' number format at the end of the file, like this
%  'IM-0001-0001.dcm') and a character matrix along the lines of what would
%  be generated with the following command:
%       files      = dir(file_dir); file_names =
%       cell2mat({files(3:end).name}');
%
% *out_dir is the output directory
% *output_filelist is the remapped filelist in a similar format
%
%   The algorithm is: 
%   filenum_new = filenum_old-(phases-1)*(slice-1)+(slices-1)*(phase-1)
%
% Alex Barker, 20181017
% Children's Hospital Colorado
% University of Colorado, Anschutz
% alexander.barker@ucdenver.edu

%% Get slices/phases

%test values (comment out when in production)
% in_dir      = 'C:\_Research_Data\4Dflow_protocol_optimization\20180829_lorna\20180829_sorted_series\303000_4DFLOW_RL';
% in_filelist = ['1234567891_s41256_i00001.dcm' ; '1234567891_s41256_i00002.dcm'];

num_files = size(in_filelist,1);

info   = dicominfo([in_dir filesep in_filelist(1,:)],'UseDictionaryVR',true);
phases = double(info.Private_2001_1017);
phases = phases(1); % take only the first entry (debug for Fabians's kids data)
slices = double(num_files./phases); % can't find field in Philips data for number of slices (didn't try very hard), so I will calculate this

%confirm I have a round number for slices
if ~floor(slices)==slices || slices == 0 || isinf(slices)
    disp('err: the number of files divided by the number of cardiac phases (as extracted from the dicom header) is zero or non-integer')
    return
end

%% get file numbers and remap

% first find '.' position in file name and then pull number character
ind = strfind(in_filelist(1,:),'.');
in_filenum  = str2num(in_filelist(:,ind-5:ind-1)); % this will cause problems if not 5 digits

% reshape my matrix according to slices and phases
in_filenum = reshape(in_filenum,[phases, slices]); % here row is the phase dimension and column is the slice dimension
phases2 = repmat(phases,phases,slices); % prepare matrices for some matrix math
slices2 = repmat(slices,phases,slices); % prepare matrices for some matrix math
[phase,slice] = ndgrid(1:phases,1:slices); 

% convert files to be in order of slice, then time
out_filenum = in_filenum - (phases2-1).*(slice-1) + (slices2-1).*(phase-1);

% now remap to character array
out_filenum  = reshape(out_filenum,num_files,1); % here row is the phase dimension and column is the slice dimension
out_file     = num2str(out_filenum,'%04.0f');
temp_prefix  = repmat('IM-0001-',num_files,1);
temp_ext     = repmat('.dcm',num_files,1);
out_filelist = [temp_prefix out_file temp_ext];
out_dir      = [in_dir(1:end) '_reshape'];


end

