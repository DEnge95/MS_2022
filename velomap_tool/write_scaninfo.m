%
% scanInfo
% get selected header information from dicom header of selected file
% write results to excel file
%

function status = write_scaninfo(dataStruct)

% get header information
if dataStruct.eDicom
    patterns = ["FH","AP","RL"];
    dirlist  = dir(dataStruct.dataDirectory);
    dirlist  = dirlist(contains({dirlist.name},patterns));
    fileName        = fullfile(dataStruct.dataDirectory,dirlist(1).name,dataStruct.fileNamesMag);
else    
    fileName        = sprintf('%s%s%s',dataStruct.dataDirectory,'/mag/',dataStruct.fileNamesMag(1,1,:));
end
%dataInfoStruct      = dicominfo(fileName);
[infoStr, dataInfoStruct, msg] = dicom_scan_singlefile(fileName,'all','dicom-dict_philips.txt');
if isfield(dataInfoStruct,'SharedFunctionalGroupsSequence')
    [dataInfoStruct,~,~] = get_necessary_tags_edicom(dataInfoStruct);
end

%% get MR method and generate excel file name
if isfield(dataInfoStruct,'SeriesDescription')
    methodStr = dataInfoStruct.SeriesDescription;
else
    methodStr = '';
end

excelFileStr = sprintf('%s%s%s%s',dataStruct.dataDirectory,'/scanInfo_',dataStruct.userSubjectIDStr,'.xls');


%% get selected header information and store in matrix
%get patient's family name
cellnum = 1;
patientName = sprintf('%s%s%s','subject','_',dataInfoStruct.StudyDate);
outExcel(cellnum,:) = {'Patient Name';patientName};
%if isfield(dataInfoStruct,'PatientsName') 
%    if isfield(dataInfoStruct.PatientsName,'GivenName') 
%        outExcel(cellnum,:) = {'Patient Name';sprintf('%s%s%s',dataInfoStruct.PatientsName.FamilyName,' ',dataInfoStruct.PatientsName.GivenName)};
%    elseif isfield(dataInfoStruct.PatientsName,'FamilyName')
%        outExcel(cellnum,:) = {'Patient Name';sprintf('%s',dataInfoStruct.PatientsName.FamilyName)};    
%    else
%        outExcel(cellnum,:) = {'Patient Name';sprintf('%s',dataInfoStruct.PatientsName)};
%    end
%elseif isfield(dataInfoStruct,'PatientName')
%    if isfield(dataInfoStruct,'PatientName.GivenName') 
%        outExcel(cellnum,:) = {'Patient Name';sprintf('%s%s%s',dataInfoStruct.PatientName.FamilyName,' ',dataInfoStruct.PatientName.GivenName)};
%    else
%        outExcel(cellnum,:) = {'Patient Name';sprintf('%s',dataInfoStruct.PatientName.FamilyName)};
%    end
%else
%    outExcel(cellnum,:) = {'Patient Name';'no entry found'};
%end

% get patient's sex
cellnum = cellnum +1;
if isfield(dataInfoStruct,'PatientsSex')
    outExcel(cellnum,:) = {'Patient Sex';dataInfoStruct.PatientsSex};         
elseif isfield(dataInfoStruct,'PatientSex')        
    outExcel(cellnum,:) = {'Patient Sex';dataInfoStruct.PatientSex};
else
    outExcel(cellnum,:) = {'Patient Sex';'no entry found'};
end
% end of: get patient's sex

% get patient's date of birth
cellnum = cellnum +1;
if isfield(dataInfoStruct,'PatientsBirthDate')        
    outExcel(cellnum,:) = {'Patient DOB';dataInfoStruct.PatientsBirthDate};
elseif isfield(dataInfoStruct,'PatientBirthDate')        
    outExcel(cellnum,:) = {'Patient DOB';dataInfoStruct.PatientBirthDate};
else
    outExcel(cellnum,:) = {'Patient DOB';'no entry found'};
end    
% end of: get patient's date of birth 

% get patient's age
cellnum = cellnum +1;
if isfield(dataInfoStruct,'PatientsAge')        
    outExcel(cellnum,:) = {'Patient Age';dataInfoStruct.PatientsAge};
elseif isfield(dataInfoStruct,'PatientAge')        
    outExcel(cellnum,:) = {'Patient Age';dataInfoStruct.PatientAge};
else
    outExcel(cellnum,:) = {'Patient Age';'no entry found'};
end    
% end of: get patient's age 

% get exam date
cellnum = cellnum +1;
if isfield(dataInfoStruct,'FileModDate')        
    outExcel(cellnum,:) = {'Exam Date / Time';dataInfoStruct.FileModDate};
elseif isfield(dataInfoStruct,'StudyDate')
   outExcel(cellnum,:) = {'Exam Date / Time';dataInfoStruct.StudyDate};
else
    outExcel(cellnum,:) = {'Exam Date / Time';'no entry found'};
end
% end of: get exam date

% get acquisition type
cellnum = cellnum +1;
if isfield(dataInfoStruct,'MRAcquisitionType')        
    outExcel(cellnum,:) = {'Volume Coverage';dataInfoStruct.MRAcquisitionType};
else
    outExcel(cellnum,:) = {'Volume Coverage';'no entry found'};
end    
% end of: get acquisition type

% get image rows
cellnum = cellnum +1;
if isfield(dataInfoStruct,'Width')        
    outExcel(cellnum,:) = {'Image Rows';dataInfoStruct.Width};
elseif isfield(dataInfoStruct,'Rows')
    outExcel(cellnum,:) = {'Image Rows';dataInfoStruct.Rows};
else
    outExcel(cellnum,:) = {'Image Rows';'no entry found'};
end    
 % end of: image rows

% get image columns
cellnum = cellnum +1;
if isfield(dataInfoStruct,'Height')        
   outExcel(cellnum,:) = {'Image Columns';dataInfoStruct.Height};
elseif isfield(dataInfoStruct,'Columns')
   outExcel(cellnum,:) = {'Image Columns';dataInfoStruct.Columns};
else
   outExcel(cellnum,:) = {'Image Columns';'no entry found'};
end    
% end of: image columns

% get number of phase encoding steps
cellnum = cellnum +1;
if isfield(dataInfoStruct,'NumberOfPhaseEncodingSteps')        
   outExcel(cellnum,:) = {'Number PE';dataInfoStruct.NumberOfPhaseEncodingSteps};
else
   outExcel(cellnum,:) = {'Number PE';'no entry found'};
end    
% end of: get number of phase encoding steps

% get phase encoding direction
cellnum = cellnum +1;
if isfield(dataInfoStruct,'PhaseEncodingDirection')        
    outExcel(cellnum,:) = {'PE Direction';dataInfoStruct.PhaseEncodingDirection};
    PE = dataInfoStruct.PhaseEncodingDirection;
elseif isfield(dataInfoStruct,'InPlanePhaseEncodingDirection')        
    outExcel(cellnum,:) = {'PE Direction';dataInfoStruct.InPlanePhaseEncodingDirection};
    PE = dataInfoStruct.InPlanePhaseEncodingDirection;
else
    outExcel(cellnum,:) = {'PE Direction';'no entry found'};
    PE ='';
end    
% end of: get phase encoding direction

% get field of view and resolution
cellnum = cellnum +1;
if isfield(dataInfoStruct,'PixelSpacing')
    pix = dataInfoStruct.PixelSpacing;
%     resCol = pix(2);
%     resRow = pix(1);
    %% get FOV row    
    if isfield(dataInfoStruct,'Width') 
        width = dataInfoStruct.Width;
        FOVrow = round(pix(1) * width);
        outExcel(cellnum,:) = {'FOV Row [mm]';FOVrow};
    elseif isfield(dataInfoStruct,'Rows')
        width = dataInfoStruct.Rows;
        FOVrow = round(pix(1) * width);
        outExcel(cellnum,:) = {'FOV Row [mm]';FOVrow};        
    else
        outExcel(cellnum,:) = {'FOV Row [mm]';'no entry found'};
    end
    %% get FOV column
    cellnum = cellnum +1;
    if isfield(dataInfoStruct,'Height')
        height = dataInfoStruct.Height;
        FOVcol = round(pix(2) * height); 
        outExcel(cellnum,:) = {'FOV Column [mm]';FOVcol};
    elseif isfield (dataInfoStruct,'Columns')
        height = dataInfoStruct.Columns;
        FOVcol = round(pix(2) * height); 
        outExcel(cellnum,:) = {'FOV Column [mm]';FOVcol};
    else
        outExcel(cellnum,:) = {'FOV Column [mm]';'no entry found'};
    end

    %% get resolution
    cellnum = cellnum+ 1;
    if (isfield(dataInfoStruct,'PhaseEncodingDirection')||isfield(dataInfoStruct,'InPhaseEncodingDirection'))...
        && (isfield(dataInfoStruct,'NumberOfPhaseEncodingSteps'))...
        &&(isfield(dataInfoStruct,'Width')||isfield(dataInfoStruct,'Rows'))&&isfield(dataInfoStruct,'SliceThickness')...
        &&(isfield(dataInfoStruct,'Height')||isfield(dataInfoStruct,'Columns')) 
        if(strcmp(PE, 'ROW'))
            resRow = str2double(sprintf('%.2f',pix(1) * height/dataInfoStruct.NumberOfPhaseEncodingSteps));
            resCol = str2double(sprintf('%.2f',pix(2)));
        else
            resRow = str2double(sprintf('%.2f',pix(1)));
            resCol = str2double(sprintf ('%.2f',pix(2) * width/dataInfoStruct.NumberOfPhaseEncodingSteps));
        end  
        resSlice = str2double(sprintf('%.2f',dataInfoStruct.SliceThickness));
        
        outExcel(cellnum,:) = {'Resolution Row [mm]';resRow};
        cellnum = cellnum + 1;
        outExcel(cellnum,:) = {'Resolution Column [mm]';resCol};
        cellnum = cellnum + 1;        
        outExcel(cellnum,:) = {'Resolution Slice [mm]';resSlice};
    else
        outExcel(cellnum,:) = {'Resolution Row [mm]';'no entry found'};
        cellnum = cellnum + 1;
        outExcel(cellnum,:) = {'Resolution Column [mm]';'no entry found'};
        cellnum = cellnum + 1;
        outExcel(cellnum,:) = {'Resolution Slice [mm]';'no entry found'};
    end 
    %% end of: get resolution

else    
   outExcel(cellnum,:) = {'FOV Row [mm]';'no entry found'};
   cellnum = cellnum + 1;
   outExcel(cellnum,:) = {'FOV Column [mm]';'no entry found'};
   cellnum = cellnum + 1;
   outExcel(cellnum,:) = {'Resolution Row [mm]';'no entry found'};
   cellnum = cellnum + 1;
   outExcel(cellnum,:) = {'Resolution Column [mm]';'no entry found'};
   cellnum = cellnum + 1;
   outExcel(cellnum,:) = {'Resolution Slice [mm]';'no entry found'};
end
% end of: get field of view and resolution

% get repetion time
cellnum = cellnum + 1;
if isfield(dataInfoStruct,'RepetitionTime')        
    outExcel(cellnum,:) = {'TR [ms]';dataInfoStruct.RepetitionTime};        
else
   outExcel(cellnum,:) = {'TR [ms]';'no entry found'};       
end    
% end of: get repetion time

% get echo time
cellnum = cellnum + 1;
if isfield(dataInfoStruct,'EchoTime')        
    outExcel(cellnum,:) = {'TE [ms]';dataInfoStruct.EchoTime};
else
   outExcel(cellnum,:) = {'TE [ms]';'no entry found'};
end    
% end of: get echo time

% get number of time frames
cellnum = cellnum + 1;
if isfield(dataInfoStruct,'CardiacNumberOfImages') 
    outExcel(cellnum,:) = {'Time Frames';dataInfoStruct.CardiacNumberOfImages};
else
    outExcel(cellnum,:) = {'Time Frames';'no entry found'};
end

% get flip angle
cellnum = cellnum + 1;
if isfield(dataInfoStruct,'FlipAngle')
    outExcel(cellnum,:) = {'Flip Angle [°]';dataInfoStruct.FlipAngle};
else
    outExcel(cellnum,:) = {'Flip Angle [°]';'no entry found'};
end
% end of: get flip angle

% get pixel band width
cellnum = cellnum + 1;
if isfield (dataInfoStruct,'PixelBandwidth')
    outExcel(cellnum,:) = {'Band Width';dataInfoStruct.PixelBandwidth};
else
    outExcel(:,cellnum ) = {'Band Width';'no entry found'};
end
% end of: get pixel band width

% get manufacturer's model name
cellnum = cellnum + 1;
if isfield(dataInfoStruct,'ManufacturersModelName')
    outExcel(cellnum,:) = {'MR System';dataInfoStruct.ManufacturersModelName};
elseif isfield(dataInfoStruct,'ManufacturerModelName')
    outExcel(cellnum,:) = {'MR System';dataInfoStruct.ManufacturerModelName};
else
    outExcel(cellnum,:) = {'MR System';'no entry found'};
end
% end of: get manufacturer's model name

% get velocity encoding
cellnum = cellnum + 1;
if isfield(dataStruct,'venc')
    if length(dataStruct.venc)==3
    outExcel(cellnum,:) = {'Venc in plane I [cm/s]';dataStruct.venc(1)};
    cellnum = cellnum + 1;
    outExcel(cellnum,:) = {'Venc in plane J [cm/s]';dataStruct.venc(2)};
    cellnum = cellnum + 1;
    outExcel(cellnum,:) = {'Venc through plane K [cm/s]';dataStruct.venc(3)};
    elseif length(dataStruct.venc)==1
        outExcel(cellnum,:) = {'Venc through plane K [cm/s]';dataStruct.venc(1)};
    end
else
    outExcel(cellnum,:) = {'Venc in plane I [cm/s]';'no entry found'};
    cellnum = cellnum + 1;
    outExcel(cellnum,:) = {'Venc in plane J [cm/s]';'no entry found'};
    cellnum = cellnum + 1;
    outExcel(cellnum,:) = {'Venc through plane K [cm/s]';'no entry found'};
end
% end of: get velocity encoding

% get R-R interval (NominalInterval), SuS 06272017
cellnum = cellnum + 1;
if isfield (dataInfoStruct,'NominalInterval')
    outExcel(cellnum,:) = {'R-R Interval [ms]';dataInfoStruct.NominalInterval};
else
    outExcel(cellnum,:) = {'R-R Interval [ms]';'no entry found'};
end
% end of: get R-R interval


% write header information to Excel file
outTable = cell2table(outExcel);
try
    writetable(outTable,excelFileStr,'WriteVariableNames',false,'FileType','spreadsheet','Sheet','scanInfo','Range','A1');
    status = 1;
catch
    status = 0;
end







