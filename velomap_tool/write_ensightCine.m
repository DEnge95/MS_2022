function dataStruct = write_ensightCine(dataStruct,mrstruct_cine)

% The write_ensightCine function is included in the velomap_tool and is
% called when a folder containing a 2D time series of magnitude data is
% loaded into the ‘Append’ list window. Standalone functionality also
% exists for the write_ensightCine function under the auspices of the
% "cine2D_tool" button included in all super_tool distributions. This
% standalone implementation is intended to create cine case files without
% the need to process a 4D dataset in the velomap_tool. Both calls to the
% write_ensightCine function will read a folder containing a 2D anatomic
% cine dicom series (for example a trufi_cine of the aortic valve) and
% create an ensight case file to complement the 4D flow case file. After a
% successful dicom to case file conversion, a folder containing the case
% files will be written either the root velomap_tool folder, or the folder
% of your choice (depending on whether write_ensightCine was called from
% the velomap_tool or as a standalone call). Once you confirm that the case
% files exist, the next step will be to add this extra case file to a 4D
% dataset within ensight.
%
% In order to add the cine case file to an existing 4D flow case file:
% 1) Click case>add and hit ‘ok’
% 2) Select the cine case file and hit ‘ok’
% 3) A colored plane should appear. To change the color to the magnitude map, select the
% cine part, then color the part with the Magnitude_Cine variable
% 4) Select a grayscale colormap
%
% Ensight (as tested on v9.03f) can do the temporal interpolation and truncation
% automatically by assigning a 'Master Timeline'. To do this:
% 1) Load in both the 4d flow case files and the trufi case files to Ensight
% 2) Click the 'solution time button' (i.e. the clock in the upper left-hand corner)
% 3) Click the 'Advanced' tab and then click on the button 'Timeset Details'
% 4) Select a timeset (under 'Which Timeset(s)') that you want to use as the master
% timeline (this is usually the 4D flow sensitive case/model defined for 'All Model
% Geometry).
% 5) Click the check box 'Master timeline'
% 6) Click on the newly green 'Update Selected Timeset(s)'.
%
% The remaining cases (such as the 2D anatomic image magnitude data) will be temporally
% interpolated to match that of the 4D case file. You can check that the timelines all
% match by selecting all timesets in the upper left corner of the dialog box ('ctrl' &
% select all timesets).
%
% Alex Barker, 6/2011
% alex.l8arker@gmail.com
% Medizin Physik, Universitätsklinikum Freiburg
%
% Note, please recognize the use of this tool where appropriate. The technique is still
% under development and will be presented in the upcoming paper (currently in preparation):
%
% Barker, A.J., Markl, M., Bürk, J., Lorenz, R., Bock, J., Bauer, S., Schulz-Menger,
% J., and von Knobelsdorff-Brenkenhoff, F. Bicuspid Aortic Valve is Associated with 
% Altered Wall Shear Stress in the Ascending Aorta. Circulation:
% Cardiovascular Imaging. 2012 (In Review)

%% init gui or nogui operation
if nargin == 0
    gui = 1;
elseif nargin == 1
    disp('ERR:write_ensightCine requires 0 or 2 input arguments... beginning operation in gui-mode');
    gui = 1;
elseif nargin == 2
    gui = 0;
end

%% main program execution
if gui == 1
    %query for source location of dataStruct and mrstruct_cine information
    disp('The gui option is still under construction... complete gui is on the way')
    dirRead  = uigetdir('C:\Documents and Settings\barker\My Documents\_temp\writeCineDebug2','Choose a 2D cine series folder to convert to a case file');
    if dirRead == 0
        msgbox('A 2D cine folder was not selected, cine2D_tool canceled')
        return
    else
        dirWrite = uigetdir('C:\Documents and Settings\barker\My Documents\_temp\writeCineDebug2','Choose a place to save the Ensight case file folder');
        if dirWrite == 0
            msgbox('A save location was not selected, cine2D_tool canceled')
            return
        end
    end
    %read the dicoms and get the trigger times
    [dataStruct mrstruct_cine] = readDicomSeries(dirRead);
    dataStruct.dataDirectory = dirWrite;
    if strcmp(dataStruct.structType,'series2D')
        write_ensightFile(dataStruct,mrstruct_cine)
    elseif strcmp(dataStruct.structType,'image')
        write_ensightFile(dataStruct,mrstruct_cine)
    elseif strcmp(dataStruct.structType,'volume')
        write_ensightFile(dataStruct,mrstruct_cine)
    else
        disp('ERR: series type (i.e. series 2D, image, volume, etc...) was unable to be determined')
    end
else
    %run write_ensightFile without user interaction
    write_ensightFile(dataStruct,mrstruct_cine)
end


%% subfunction to read 2D cine dicoms (mostly excerpted from the velomap_tool)
    function [dataStruct mrstruct_anatom] = readDicomSeries(dirRead)
        % get files in the read directory
        [fileNamesMx, dirStr] = local_get_filelist(dirRead,'dcm');
        numFiles  = size(fileNamesMx,1);
        edges = zeros(4,4,numFiles);
        anatomTime = nan(numFiles,1); % need trigger time for SSFP cines
        % read first file
        pathName            = full_filename(dirStr,fileNamesMx(1,:));
        [actMRstruct actMRstruct.user] = dicom_read_singlefile(pathName,1,[],1);
        anatomData          = zeros(size(actMRstruct.dataAy,1),size(actMRstruct.dataAy,2),numFiles);
        edges(:,:,1)        = actMRstruct.edges;
        anatomData(:,:,1)   = actMRstruct.dataAy;
        if isfield(actMRstruct.user,'TriggerTime')
            anatomTime(1)       = actMRstruct.user.TriggerTime;
        else
            anatomTime(1)       = 0;
        end
        % get PatientsName and remove spaces
        if isfield(actMRstruct.user.PatientsName,'FamilyName')
            patientName = sprintf('%s%s',actMRstruct.user.PatientsName.FamilyName,'_',actMRstruct.user.StudyDate);
        else
            patientName = sprintf('%s%s',actMRstruct.user.PatientsName,'_',actMRstruct.user.StudyDate);
        end
        dataStruct.patientName = patientName;
        % get edges, mrstruct, and trigger times for additional files
        if numFiles>1
            for n=2:numFiles
                pathName            = full_filename(dirStr,fileNamesMx(n,:));
                [actMRstruct actMRstruct.user] = dicom_read_singlefile(pathName,1,[],1);
                edges(:,:,n)        = actMRstruct.edges;
                anatomData(:,:,n)   = actMRstruct.dataAy;
                if isfield(actMRstruct.user,'TriggerTime')
                    anatomTime(n)       = actMRstruct.user.TriggerTime;
                else
                    anatomTime(n)       = 0;
                end
            end
            if isequal( edges(:,:,1), edges(:,:,end))
                struct_type = 'series2D';
            else
                struct_type = 'volume';
            end
        else
            struct_type = 'image';
        end
        actMRstruct.dataAy = [];
        actMRstruct.user   = [];
        % calculate edges and write geo file(s)
        actMRstruct.edges  = get_volume_geo(edges);
        %initialize mrstruct
        mrstruct_anatom = mrstruct_init(struct_type,anatomData,actMRstruct);
        if numFiles > 1 && ~all(anatomTime) 
            mrstruct_anatom.user = anatomTime;
        elseif all(anatomTime)
            mrstruct_anatom.user = [];
        end
        clear('anatomData','anatomTime','fileNamesMx','actMRstruct','pathName','edges');
        % scale data to range [0 1]
        maxValue = max(max(max(mrstruct_anatom.dataAy)));
        minValue = min(min(min(mrstruct_anatom.dataAy)));
        mrstruct_anatom.dataAy = (mrstruct_anatom.dataAy - minValue)/(maxValue-minValue);
        % save structure type
        dataStruct.structType = struct_type;
        
        %         if strcmp(struct_type,'series2D') %AJB
        %             write_ensightFile(dataStruct,mrstruct_anatom)
        %         end
    end

%% subfunction to read files in a directory
    function [fileNamesMx, dirStr] = local_get_filelist(dirStr,extensionStr)
        %%% default settings
        fileListStruct = [];
        fileNamesMx    = [];
        %%% get file list
        if isempty(extensionStr)
            askStr = full_filename(dirStr,'*');
        else
            tmpStr = sprintf('*.%s',extensionStr);
            askStr = full_filename(dirStr,tmpStr);
        end
        dirStruct    = dir(askStr);
        noOfEntries  = size(dirStruct,1);
        %%% create fileListStruct
        fileCount = 1;
        for k=1:noOfEntries
            fileStr = dirStruct(k).name;
            if dirStruct(k).isdir==0 && strcmp(fileStr,'..')==0 &&  strcmp(fileStr,'.')==0
                fileListStruct(fileCount).name = fileStr;
                fileCount = fileCount+1;
            end
        end
        %%% convert fileListStruct to matrix of file names
        if any(size(fileListStruct))
            fileNamesMx = sortrows(char(fileListStruct.name));
        end
    end
        
%% subfunction to write ensight case file
    function write_ensightFile(dataStruct,mrstruct_cine)
        if strcmp(dataStruct.structType,'series2D') || strcmp(dataStruct.structType,'image') 
            % get trigger time information for cines
            numPhases      = size(mrstruct_cine.dataAy,3);
            timeStamps     = mrstruct_cine.user;
        elseif strcmp(dataStruct.structType,'volume')
            numPhases      = 1;
            timeStamps     = 0;
        end
        % set up ensight write location
        ensightDirRoot = dataStruct.dataDirectory;
        ensightDirPath = fullfile(ensightDirRoot,['EnSight_' dataStruct.patientName '_cine01']);
        %check to see if additional folders exist, if so, create unique name
        if isdir(ensightDirPath)
            folderNames    = ls(ensightDirRoot);  % need this to get quick character only array of names (and not a structure, as is obtained with the 'dir' command)
            folderInfo     = dir(ensightDirRoot); % need this to get names without character padding
            folderIndex    = strmatch(['EnSight_' dataStruct.patientName '_cine'],folderNames);
            numCineFolders = str2double(folderInfo(folderIndex(end)).name(end-1:end));
            ensightDirPath(end-1:end) = num2str(numCineFolders+1,'%02.0f');
            mkdir(ensightDirPath)
        else
            mkdir(ensightDirPath);
        end
        patientName   = dataStruct.patientName;
        casePathName  = sprintf('%s%s%s%s%s',ensightDirPath,filesep,'EnSight_',patientName,'.case');
        geoFileName   = sprintf('EnSight_%s%s%',patientName,'.geo');
        dataPathName  = sprintf('%s%s%s%s%s',ensightDirPath,filesep,'EnSight_',patientName,'_');
        dataFileName  = sprintf('EnSight_%s%s%',patientName,'_');
        % open text file and write data
        fidCase  = fopen(casePathName, 'wt');
        fprintf(fidCase,'FORMAT\n');
        fprintf(fidCase,'type:	 ensight gold\n');
        fprintf(fidCase,'GEOMETRY\n');
        fprintf(fidCase,'model:	 %s\n',geoFileName);
        fprintf(fidCase,'VARIABLE\n');
        % magnitude data
        fprintf(fidCase,'scalar per node:	 Magnitude_Cine	 %s**.mag\n',dataFileName );
        fprintf(fidCase,'TIME\n');
        fprintf(fidCase,'time set:		 1\n');
        fprintf(fidCase,'number of steps:	 %s\n',num2str(numPhases));
        fprintf(fidCase,'filename start number:	 0\n');
        fprintf(fidCase,'filename increment:	 1\n');
        fprintf(fidCase,'time values:\n');
        fprintf(fidCase,'%.3f\n',timeStamps);
        fclose(fidCase); 
        % set up geometry file
        % set transformation flag, here edges
        transf_flag =[];
        % create geo file
        geoFileName             = sprintf('%s%s%s%s%s',ensightDirPath,filesep,'EnSight_',patientName,'.geo');
        mrstruct_temp           = mrstruct_cine;
        if strcmp(dataStruct.structType,'series2D') || strcmp(dataStruct.structType,'image') 
            mrstruct_temp.dataAy    = mrstruct_cine.dataAy(:,:,1);
        end
        [res, errStr_geo, oArg] = mrstruct_ensight(mrstruct_temp,'geoFile',geoFileName,1,transf_flag,'create');
        for k = 1:numPhases;
            FileName                = sprintf('%s%s%s%',dataPathName,num2str(k-1,'%02d'),'.mag');
            mrstruct_temp           = mrstruct_cine;
            if strcmp(dataStruct.structType,'series2D') || strcmp(dataStruct.structType,'image') 
                mrstruct_temp.dataAy    = mrstruct_cine.dataAy(:,:,k);
            end
            [res, errStr_vel, oArg] = mrstruct_ensight(mrstruct_temp,'dataVolume',FileName,1,'');
            clear mrstruct_temp
        end
    end
end











