% 
% Last edit by Kelly Jarvis and Liliana Ma: August 16, 2017
%       -Modified the pcmra calculations to include an additional sqrt of
%       mean sum of squares option
%       -Did not fix function [volume_pcmra] = proc_pcmra (dataStruct) to
%       be compatible 

% Orientations in Siemens vs NU 4D flow sequences Modified by Liliana Ma,
% July 24, 2017.
% 
%   To do: Figure out transverse acquisition signVijk and signVx/Vy/Vz and
%   peDir ~= i. Store Vx/Vy/Vz within the velomap dataStruct.
% 
%   1) signVijk or dataStruct.signVijk accounts for the way the images are
%   displayed in the velomap tool, and played in the avi, but NOT the
%   coordinate axes directions in the mrStruct or Ensight output. signVijk
%   was tested for phase encoding direction row or 'i', and sagittal and
%   coronal orientations.
% 
%   2) Probably for historical reasons, the Ensight and mrStruct coordinate
%   axes are defined by signVx, signVy, and signVz. The code with these
%   definitions is repeated three times throughout the velomap_internal file.
%   It may be easier to save this as a dataStruct variable in the future. The
%   axes signs have been tested for sagittal and coronal orientations, but
%   the lack of a phase encode direction logic in defining these variables
%   may be why the velocities are sometimes inverted.

%		Jelena Bock
%		2/2006
% 
% 
% 


%		PC
%

function [ret1,ret2] = velomap_internal(handleStruct,dataStruct,commandStr,entryStr,dataValue)


%%%%% init and error check
ret1 = dataStruct;
ret2 = handleStruct;
if nargin~=5
   warndlg(lasterr,'there are not enough arguments');
   return;
end
%%%%% End of: init and error check


%%%%% set dataStruct
if strcmp(commandStr,'database')
   %try
      dataStruct = local_set_dataStruct_entry(handleStruct,dataStruct,entryStr,dataValue);
      ret1       = dataStruct;
      %catch
      %warndlg(lasterr,'Error in template_internal','database');
      %end
%%%%% End of: set dataStruct

%%%%% update GUI
elseif strcmp(commandStr,'gui')
   %try
      [dataStruct, handleStruct] = local_update_GUI(handleStruct,dataStruct,entryStr,dataValue);
      ret1       = dataStruct;
      ret2       = handleStruct;
      %catch
      %warndlg(lasterr,'Error in template_internal','database');
      %end
%%%%% End of: update GUI
 
else 
    warndlg(lasterr,'Sorry, command switch was not recognized in ''velomap_internal''');    
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  main local functions                                                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%% set dataStruct value
function  dataStruct = local_set_dataStruct_entry(handleStruct,dataStruct,entryStr,dataValue)

%%% START set block
if strcmp(entryStr, 'all')
    dataStruct = dataValue; 
   
elseif strcmp(entryStr, 'loadData')
    dataStruct = local_sort_filenames(dataStruct, dataValue);
    
    %% if new data was chosen
    if strcmp(dataStruct.loadStatus, 'ok')
       %% set default values 
        dataStruct.stdFlow          = [];
        dataStruct.stdFlowMax       = [];
        dataStruct.imaFLOW          = [];
        dataStruct.imaMAG           = [];
        dataStruct.pcmraSumSquares  = [];
        dataStruct.pcmraMeanAbsVel  = [];
        dataStruct.pcmraMask        = [];
        dataStruct.eddyMask         = [];
        dataStruct.pcmraPseudoComplDiff  = [];
        dataStruct.pcmraSqrtSumSquares   = []; %% LiliMa
        dataStruct.pcmraVelCorr   = []; %% MBS
        dataStruct.filterMask       = [];
        dataStruct.statMask         = [];
        dataStruct.statRegStatus    = 'off';
        dataStruct.pcmraStatus      = 'off';
        dataStruct.previewStatus    = 'off';
       %% end of: set default values 
    end    
    
elseif strcmp(entryStr, 'outToExcel')
    %% read scan information and write to an excel file
    if ~isempty(dataStruct.fileNamesFlow)        
        status = write_scaninfo(dataStruct);
        if status==1
            dataStruct.status = sprintf('Scan information was successfully saved');
        else
            dataStruct.status = sprintf('Scan information could NOT be saved');
        end
    else
        dataStruct.status = sprintf('There is no scan information to save');
    end
    %% end of: read scan information and write to an excel file
    
elseif strcmp(entryStr, 'preview')
    dataStruct.statRegStatus = 'off';    
    if isempty(get(handleStruct.orderList,'String'))
        dataStruct.previewStatus = 'off';
        dataStruct.status = sprintf('Original data was reloaded!');
    else
        dataStruct.previewStatus='on';                
    end

elseif strcmp(entryStr, 'vencValueEdt')
    dataStruct.venc     = dataValue;
    dataStruct.status   = sprintf('%s\n%s','Values of venc are changed.','Original data is reloaded.');
    
% if preview of static tissue was chosen  
elseif strcmp(entryStr, 'previewStatReg')
    if get(handleStruct.previewStatRegBut,'Value')==1
        set(handleStruct.previewStatRegBut,'String','reset');
        dataStruct.statRegStatus = 'on';
        h = waitbar(1,'Application is busy. Please wait !');
        %% calculate static tissue
        actSlice     = ceil(get(handleStruct.SlicesNumSlide,'value'));
        numSlices    = 1;
        staticFactor = str2double(get(handleStruct.staticValueEdt,'String'));
        %if (isempty(dataStruct.imaFLOW)||isempty(dataStruct.imaMAG));%||isempty(dataStruct.imaFlowOrig))
            dataStruct = local_read_all_phases(dataStruct,actSlice,1);
        %end

        dataStruct = flow_stdev(dataStruct);                      
        dataStruct = create_static_tissue_mask(dataStruct,staticFactor,actSlice);

        %% end of: calculate static tissue 
        dataStruct.status = sprintf('Preview of static regions');                        
        close(h);
    else
        set(handleStruct.previewStatRegBut,'String','preview of static regions');
        dataStruct.statRegStatus = 'off';
    end
% end of: if preview of static tissue was chosen

elseif strcmp(entryStr, 'previewBackground') 
    % get threshold
    tvalue = get(handleStruct.backgroundSlide,'Value');
    maskdata = thresholding_percent(dataStruct.imaMAG,tvalue);
    dataStruct.backgroundMask = maskdata;
    % get magnitude maximum
    
   
%% append a chosen directory with anatomic data to the list
elseif strcmp(entryStr, 'appAnatomData')
    directory_name  = uigetdir(dataStruct.dataDirectory,'select directory with anatomic data');
    %% continue only if a directory name was chosen
    if directory_name~=0
        dir_list = cellstr([get(handleStruct.AnatomDataList,'String');directory_name]);    
        set(handleStruct.AnatomDataList,'String',dir_list);
    end
% end of: append directory with anatomic data to the list

%% delete a chosen directory with anatomic data from the list
elseif strcmp(entryStr, 'deleteAnatomData')
    id_value = get(handleStruct.AnatomDataList,'Value');
    dir_list = get(handleStruct.AnatomDataList,'String');
    if ~isempty(dir_list)
        dir_list(id_value) = [];
        set(handleStruct.AnatomDataList,'String',dir_list,'Value',1);
    end
% end of: delete a chosen directory with anatomic data from the list

%% write preprocess functions settings to an excel file    
elseif strcmp(entryStr, 'saveFuncSetting')  

	% create new directory for results files
	resultsDirStr = ['results_user',dataStruct.userIDStr];
    resultsDirPath = fullfile(dataStruct.dataDirectory,resultsDirStr);
    mkdir(resultsDirPath);
	
    excelfileNameStr = ['preProcessing_settings_',dataStruct.userSubjectIDStr,'.xls'];
    excelFileStr = sprintf('%s%s%s',dataStruct.dataDirectory,filesep,excelfileNameStr);
    func_list    = get(handleStruct.orderList,'String');
    if ~isempty(func_list)
        %% set default values
        rowNum = 1;
        outExcel = cell(25,3);
        %% end of: set default values
        
        %% check if noise masks were combined
        a = strfind(func_list, 'noise-filter');
        b = strfind(func_list, 'stdev-filter');
        if ~isempty(a) && ~isempty(b)
            if get(handleStruct.radioORBut, 'Value')==1
                filt_comb = 'OR';
            else
                filt_comb = 'AND';                
            end
        else
            filt_comb = 'none';
        end
        %% end of: check if noise masks were combined
		
		% check if delete static tissue has been selected
		if get(handleStruct.staticTissueCheck,'Value') == 1
			delStatTissueStr = 'ON';
	    else
			delStatTissueStr = 'OFF';
	    end
				
		% check if median filter for PC-MRA is selected
		if get(handleStruct.pcmraMedianFilterCheck,'Value') == 1
			medFiltStr = 'ON';
	    else
			medFiltStr = 'OFF';
	    end
		
		% get venc threshhold for complex difference
		vencThreshold = str2double(get(handleStruct.vencThresholdEdit, 'String'));
		
		% get user selected min and max time frames for teporally constrained PC-MRA //MM
		tPCMRA_min = ceil( str2double(get(handleStruct.timeFramesMinEdit, 'String')) );
		tPCMRA_max = ceil( str2double(get(handleStruct.timeFramesMaxEdit, 'String')) );
		
	    % get venc threshhold for definition of static tissue for eddy current correction
		eddyThreshold = str2double(get(handleStruct.staticValueEdt, 'String'));
            
        %% write chosen preprocessing functions in a cell array
        for i=1:length(func_list)
            imageFun = func_list{i};
            outExcel(1,:)={'preprocessing functions','threshold value','combination - ON/OFF'};            
            switch(imageFun)
                case('noise-filter')
                    rowNum = rowNum+1;
                    outExcel(rowNum,:)={'noise-filter',get(handleStruct.noiseValueEdt,'String'),filt_comb};
                case('stdev-filter')
                    rowNum = rowNum+1;
                    outExcel(rowNum,:)={'stdev-filter',get(handleStruct.stdevNoiseValueEdt,'String'), filt_comb};
                case('derivative-filter')
                    rowNum = rowNum+1;
                    outExcel(rowNum,:)={'derivative-filter',get(handleStruct.derivatNoiseValueEdt,'String'), filt_comb};
                case('median-filter')
                    rowNum = rowNum+1;
                    outExcel(rowNum,:)={'median-filter','',''};
                case('eddy-current')
                    rowNum = rowNum+1;
                    outExcel(rowNum,:)={'eddy-current',get(handleStruct.staticValueEdt,'String'),''};
                case('2nd-order-corr')
                    rowNum = rowNum+1;
                    outExcel(rowNum,:)={'2nd order correction',get(handleStruct.staticValueEdt,'String'),''};
                case('anti-aliasing')
                    rowNum = rowNum+1;
                    outExcel(rowNum,:)={'anti-aliasing','',''};
                case('manual anti-aliasing')
                    rowNum = rowNum+1;
                    outExcel(rowNum,:)={'manual anti-aliasing','',''};
                otherwise
            end
        end
		
		outExcel(rowNum + 1,:)={'Delete static tissue','', delStatTissueStr};
		outExcel(rowNum + 2,:)={'PC-MRA: median filter','3x3', medFiltStr};
		outExcel(rowNum + 3,:)={'PC-MRA 2: min time frame',num2str(tPCMRA_min), ''};
		outExcel(rowNum + 4,:)={'PC-MRA 2: max time frame',num2str(tPCMRA_max), ''};
		outExcel(rowNum + 5,:)={'PC-MRA 3: venc treshhold',num2str(vencThreshold), ''};
	    outExcel(rowNum + 6,:)={'eddy-current - threshold',num2str(eddyThreshold), ''};
	
        warning off MATLAB:xlswrite:AddSheet;
        outTable = cell2table(outExcel);
        try
            writetable(outTable,excelFileStr,'WriteVariableNames',false,'FileType','spreadsheet','Sheet','functions settings','Range','A1');
            status = 1;
        catch
            status = 0;
        end
        
        if status ==1
            dataStruct.status = sprintf('%s%s', 'Settings are successfully saved in ',excelFileStr); 
        else
            dataStruct.status = 'No settings are saved, please try again.';
        end
    else
        dataStruct.status = 'No settings are saved, please choose preprocessing functions and try again.';
    end
% end of: write preprocess functions setting to an excel file
    
%% if button ''apply to data'' was pressed    
elseif strcmp(entryStr, 'applyTodata')
    	    
    if strcmp(dataValue,'dequeue')
        readKeys    = {'conversion','','ensight','','0'};                
        ensightFlag = floor(str2double(char(inifile(dataStruct.inifilePath,'read',readKeys))));        
    else
        ensightFlag = get(handleStruct.enSightCheck, 'Value');
    end
    
    if ensightFlag
        
        if strcmp(dataValue,'dequeue')
            readKeys    = {'conversion','','pcmra','','0'};                
            pcmraFlag   = floor(str2double(char(inifile(dataStruct.inifilePath,'read',readKeys))));
            readKeys    = {'conversion','','pdmaps','','0'};  
            pdFlag      = floor(str2double(char(inifile(dataStruct.inifilePath,'read',readKeys))));
        else
            pcmraFlag   = get(handleStruct.pcmraCheck,'Value');
            pdFlag      = get(handleStruct.pdCheck,'Value'); %% MM: add flag for pressure difference data
        end
             
        numPhases       = dataStruct.numPhases;
        TR              = dataStruct.TR; % TR in msec
        
        % create new directory for EnSight files
        ensightDirPath              = strcat('EnSight_',dataStruct.userSubjectIDStr);
        ensightDirPath              = fullfile(dataStruct.dataDirectory,ensightDirPath);
        mkdir(ensightDirPath);
        
        % generate case file
        dataPathName  = sprintf('%s%s%s%s%s',ensightDirPath,filesep,'EnSight_',dataStruct.userSubjectIDStr,'_');
        casePathName  = sprintf('%s%s%s%s%s',ensightDirPath,filesep,'EnSight_',dataStruct.userSubjectIDStr,'.case');        
        
        
        geoFileName   = sprintf('EnSight_%s%s%',dataStruct.userSubjectIDStr,'.geo');
        dataFileName  = sprintf('EnSight_%s%s%',dataStruct.userSubjectIDStr,'_');

		%% start TF modified to correctly output timeStamps for Philips data 
        if strcmpi(dataStruct.Manufacturer,'philips')
		    timeStamps = mean(sort(dataStruct.TimeStamps, 2), 1);
		elseif ~isempty(dataStruct.SiemensFlag) || strcmpi(dataStruct.Manufacturer,'siemens') % AJB need to add feature for Philips timestamps 
            timeStamps   = (1:numPhases)*TR - TR/2; %LiliMa: time stamps; may need to change this
        else
            timeStamps = mean(sort(dataStruct.TimeStamps, 2), 1); %AJB: Need to fill this in where it occurs in the rest of the code
        end
		%% end   TF modified
		
        %% open text file and write data
        fidCase  = fopen(casePathName, 'wt');
        fprintf(fidCase,'FORMAT\n');
        fprintf(fidCase,'type:	 ensight gold\n');
        fprintf(fidCase,'GEOMETRY\n');
        fprintf(fidCase,'model:	 %s\n',geoFileName);
        fprintf(fidCase,'VARIABLE\n');
        fprintf(fidCase,'scalar per node:	 Magnitude	 %s**.mag\n',dataFileName );
        
        %% add speed data if selected by user
        if pcmraFlag == 1
            if strcmp(dataValue,'dequeue')
                readKeys                = {'modes pcmra','','modeSumSquares','','0'};                
                modeSumSquares          = floor(str2double(char(inifile(dataStruct.inifilePath,'read',readKeys))));
                readKeys                = {'modes pcmra','','modeMeanAbsVel','','0'};
                modeMeanAbsVel          = floor(str2double(char(inifile(dataStruct.inifilePath,'read',readKeys))));
                readKeys                = {'modes pcmra','','modePseudoComplDiff','','0'};
                modePseudoComplexDiff   = floor(str2double(char(inifile(dataStruct.inifilePath,'read',readKeys))));
				readKeys                = {'modes pcmra','','modeTimeAveMag','','0'};
                modeTimeAveMag          = floor(str2double(char(inifile(dataStruct.inifilePath,'read',readKeys))));
                readKeys                = {'modes pcmra','','modeSqrtSumSquares','','0'}; %% LiliMa
                modeSqrtSumSquare       = floor(str2double(char(inifile(dataStruct.inifilePath,'read',readKeys))));
                readKeys                = {'modes pcmra','','modeVelCorr','','0'}; %% MBS
                modeVelCorr             = floor(str2double(char(inifile(dataStruct.inifilePath,'read',readKeys))));
            else
                modeSumSquares          = get(handleStruct.pcmraSumSquaresCheck,'Value');
                modeMeanAbsVel          = get(handleStruct.pcmraMeanAbsVelCheck,'Value');
                modePseudoComplexDiff   = get(handleStruct.pcmraPseudoComplDiffCheck,'Value');
				modeTimeAveMag          = get(handleStruct.pcmraTimeAveMagCheck,'Value');
                modeSqrtSumSquare       = get(handleStruct.pcmraSqrtSumSquaresCheck, 'Value'); 
                modeVelCorr             = get(handleStruct.pcmraVelCorrCheck, 'Value'); 
            end
            
            if modeSumSquares
                fprintf(fidCase,'scalar per node:	 Speed_SumSquares	     %s**.spd01\n',dataFileName );
            end
            if modeMeanAbsVel
                fprintf(fidCase,'scalar per node:	 Speed_MeanAbsVel	     %s**.spd02\n',dataFileName );
            end
            if modePseudoComplexDiff
                fprintf(fidCase,'scalar per node:	 Speed_PseudoComplDiff	     %s**.spd03\n',dataFileName );
            end
            if modeTimeAveMag
                fprintf(fidCase,'scalar per node:	 Speed_TimeAveMag	     %s**.spd04\n',dataFileName );
            end
            if modeSqrtSumSquare %% LiliMa
                fprintf(fidCase,'scalar per node:	 Speed_SqrtSumSquares	     %s**.spd05\n',dataFileName );
            end
            if modeVelCorr %% MBS
                fprintf(fidCase,'scalar per node:	 Speed_VelCorr	     %s**.spd06\n',dataFileName );
            end
            %             if 1==1
            %                 fprintf(fidCase,'scalar per node:	 Centerline	     %s**.cl\n',dataFileName );
            %             end
            
        end
        %% end of: add speed data if selected by user 
        
        %% add pressure difference data if selected by user
        if pdFlag == 1 
            %fprintf(fidCase,'scalar per node:	 Centerline	     %s**.cl\n',dataFileName ); 
            fprintf(fidCase,'scalar per node:	 pressureDifference	     %s**.pd\n',dataFileName );            
            fprintf(fidCase,'vector per node:	 pressureGradient	 %s**.pgrad\n',dataFileName );   
            fprintf(fidCase,'scalar per node:	 pdmask_smooth	     %s**.pdm\n',dataFileName );            
        end
        %% end of: add pressure difference data if selected by user 
        
        
        fprintf(fidCase,'vector per node:	 Velocity	 %s**.vel\n',dataFileName );
        
        %% add additional anatomic data if selected by user
        if ~isempty(get(handleStruct.AnatomDataList,'String'))            
            numDir = size(get(handleStruct.AnatomDataList,'String'));
            for i=1:numDir(1)
                fprintf(fidCase,'scalar per node:	 Anatomy%d	     %s.an%d\n',i,dataFileName,i);
            end
        end
        %% end of: add additional anatomic data if selected by user
    
        fprintf(fidCase,'TIME\n');
        fprintf(fidCase,'time set:		 1\n');
        fprintf(fidCase,'number of steps:	 %s\n',num2str(numPhases));
        fprintf(fidCase,'filename start number:	 0\n');
        fprintf(fidCase,'filename increment:	 1\n');
        fprintf(fidCase,'time values:\n');
        fprintf(fidCase,'%.3f\n',timeStamps);
        fclose(fidCase); % close file
        %% end of: generate case file
        
        %% write data to dataStruct       
        dataStruct.dataPathName     = dataPathName;
        dataStruct.ensightDirPath   = ensightDirPath;
    end
   
% end of: if button ''apply to data'' was pressed    
    
%%    
elseif (strcmp(entryStr, 'preProcess') || strcmp(entryStr, 'createMRstruct')||strcmp(entryStr, 'dequeueing'))      
    %% set read-flags and number of slices
    %% description of the read flags :
    %% 1 --> only the preview is desired
    %% 2 --> only mrstruct is selected, no ensight or avi conversion
    %% 3 --> conversion to ensight and/or avi-movie
    if strcmp(entryStr, 'preProcess')%% if only the preview is desired
        numSlices = 1;
        lowlimit  = ceil(get(handleStruct.SlicesNumSlide,'value'));
        readFlag  = 1; 
    elseif strcmp(entryStr, 'createMRstruct')
        pcmraFlag     = get(handleStruct.pcmraCheck,'Value');
        
        if (get(handleStruct.enSightCheck,'Value')==1 ||get(handleStruct.aviMovieCheck,'Value')==1||get(handleStruct.pcmraToDICOMCheck,'Value')==1);
            numSlices     = dataStruct.numSlices;
            lowlimit      = 1; 
            readFlag      = 3;
        else %% if only mrstruct selected
            numSlices     = dataStruct.numSlices;
            lowlimit      = 1;
            readFlag      = 2;
        end 
    elseif strcmp(entryStr, 'dequeueing')
        %% get ini info if pcmra desired or not
        readKeys    = {'conversion','','pcmra','','0'};                
        pcmraFlag   = floor(str2double(char(inifile(dataStruct.inifilePath,'read',readKeys))));
        
        
        %% settings for read flags and number of slices
        numSlices     = dataStruct.numSlices;
        lowlimit      = 1;
        readFlag      = 3;
        
    end
    %% end of: set read-flags and number of slices
   
    
    
    %% get list of chosen preprocessing functions
    if strcmp(entryStr, 'dequeueing')
        %% get number of the preprocessing functions
        readKeys ={'order','','number_func','','0'};
        num_func = floor(str2double(char(inifile(dataStruct.inifilePath,'read',readKeys))));
        if num_func > 0
            order_list = cell(num_func,1);
            for j=1:num_func
                %% read the order of the functions from ini file
                readKeys ={'order','',strcat('func',num2str(j)),'','none'};
                order_list{j} = char(inifile(dataStruct.inifilePath,'read',readKeys));                
            end
        else
            order_list ={};
        end
    else        
        order_list = get(handleStruct.orderList,'String');
    end 
    
    %% get order of the preprocessing functions,
    %% it's important for combination of noise filter and for deleting of static tissue
    %% set default values
    indexNoise         = 0;
    indexStdev         = 0;
    indexEddy          = 0;
    indexAliasing      = 0;
    indexCorrInhom     = 0;
    %% end of: set default values
    if ~isempty(order_list)        
        for j=1:length(order_list)            
            imageFun = order_list{j};
            if strcmp(imageFun,'noise-filter')
                indexNoise = j;
            elseif strcmp(imageFun,'stdev-filter')
                indexStdev = j;
            elseif strcmp(imageFun,'eddy-current')
                indexEddy = j;
            elseif strcmp(imageFun,'anti-aliasing')
                indexAliasing = j;
            end
        end  
        order_list{j+1} = 'delStatTissue';
    end
   
    %% load all phases from each slice and modify according to the settings
    waitH       = waitbar(0,'');
    for m = 1:numSlices
       actSlice            = lowlimit+ m-1;
       %% settings for waitbar       
       msgOut = ['Pre-Processing slice #',num2str(actSlice)];
       counter = (m/numSlices);
       waitbar(counter,waitH,msgOut);
        %% end of: settings for waitbar
        
       dataStruct = local_read_all_phases(dataStruct,actSlice,readFlag);
       %% if any preprocessing function has been choosen  
       if ~isempty(order_list)
           %% get threshold for static tissue 
           if strcmp(entryStr,'dequeueing')
               readKeys     = {'func settings','','static_factor','','0.0'};
               staticFactor = str2double(char(inifile(dataStruct.inifilePath,'read',readKeys)));
           else
               staticFactor = str2double(get(handleStruct.staticValueEdt,'String'));
           end
           %% end of: get threshold for static tissue
           
            dataStruct.status   = sprintf('');
        
            for j=1:length(order_list)
            
                switch(order_list{j})
                    case('nonuniformity')
                        if strcmp(entryStr, 'dequeueing')
                            readKeys ={'func settings','','nonuniformity','','0.1'};
                            tvalue  = str2double(char(inifile(dataStruct.inifilePath,'read',readKeys))); 
                        else
                            tvalue = get(handleStruct.backgroundSlide,'Value');
                        end
                        dataValue.backgroundMask = thresholding_percent(dataStruct.imaMAG,tvalue);
                        B0_order = 4;
                        % correct non-uniformity using background mask
                        [biasfield, magcorr] = estimate_bias_field(dataStruct.imaMAG,B0_order,not(dataValue.backgroundMask));
                %%-------------------------------------------------------------
                %% apply noise mask to the data
                    case('noise-filter')
                %____________________________________________________  
                    %% get filter combination settings
                                if ((indexStdev~=0) && (indexStdev < indexNoise))
                                    if strcmp(entryStr, 'dequeueing')
                                        readKeys ={'func settings','','filt_comb','','none'};
                                        combSet = char(inifile(dataStruct.inifilePath,'read',readKeys));
                                        if strcmp(combSet,'AND')
                                            filtComb = 2;
                                        elseif strcmp(combSet,'OR')
                                            filtComb = 1;                                
                                        else
                                            filtComb = 0;
                                        end
                                    else
                                        if get(handleStruct.radioANDBut,'Value')==1
                                            filtComb = 2; %% filter combination with AND
                                        else
                                            filtComb = 1; %% filter combination with OR
                                        end
                                    end
                                else
                                    filtComb = 0; %% no filter combination
                                end
                                %% end of: get filter combination settings
                            %__________________________________________________________

                                %% get noise threshold  
                                if strcmp(entryStr, 'dequeueing')
                                    readKeys ={'func settings','','noise_filter','','0.0'};
                                    noiseFactor = str2double(char(inifile(dataStruct.inifilePath,'read',readKeys)));                        
                                else
                                    noiseFactor = str2double(get(handleStruct.noiseValueEdt,'String'));
                                end
                                %% end of: get noise threshold
                            %__________________________________________________________

                                dataStruct        = flow_noise_filter(dataStruct,noiseFactor,filtComb);
                                dataStruct.status = sprintf('%s%s\n',dataStruct.status,...
                                                        'Data was filtered with noise filter.');
                            %______________________________________________________                    

                            %% end of: apply noise mask to the data
                          %%-----------------------------------------------------------
                            % apply derivative filter
                            case('derivative-filter')
                                %% get derivat noise threshold  
                                if strcmp(entryStr, 'dequeueing')
                                    readKeys ={'func settings','','derivat_filter','','0.0'};
                                    noiseFactor = str2double(char(inifile(dataStruct.inifilePath,'read',readKeys)));                        
                                else
                                    noiseFactor = str2double(get(handleStruct.derivatNoiseValueEdt,'String'));
                                end
                                %% end of: get stdev noise threshold

                                dataStruct = flow_derivat_filter(dataStruct,noiseFactor);                    
                                dataStruct.status = sprintf('%s%s\n',dataStruct.status,...
                                                    'Data was filtered with derivative filter.'); 

                          %%-----------------------------------------------------------
                          %%
                            %% apply standard deviation noise mask to the data
                            case('stdev-filter')
                                %% get stdev noise threshold  
                                if strcmp(entryStr, 'dequeueing')
                                    readKeys ={'func settings','','stdev_filter','','0.0'};
                                    noiseFactor = str2double(char(inifile(dataStruct.inifilePath,'read',readKeys)));                        
                                else
                                    noiseFactor = str2double(get(handleStruct.stdevNoiseValueEdt,'String'));
                                end
                                %% end of: get stdev noise threshold

                                %% get filter combination settings
                                if ((indexNoise~=0) && (indexStdev >indexNoise))
                                    if strcmp(entryStr, 'dequeueing')
                                        readKeys ={'func settings','','filt_comb','','none'};
                                        combSet = char(inifile(dataStruct.inifilePath,'read',readKeys));
                                        if strcmp(combSet,'AND')
                                            filtComb = 2;
                                        elseif strcmp(combSet,'OR')
                                            filtComb = 1;                                
                                        else
                                            filtComb = 0;
                                        end
                                    else
                                        if get(handleStruct.radioANDBut,'Value')==1
                                            filtComb = 2; %% filter combination with AND
                                        else
                                            filtComb = 1; %% filter combination with OR
                                        end
                                    end
                                else
                                    filtComb = 0; %% no filter combination
                                end 
                                %% end of: get filter combination settings

                                if ((indexEddy == 0 )||(indexStdev < indexEddy))
                                    dataStruct = flow_stdev(dataStruct);
                                end
                                dataStruct = flow_stdev_filter(dataStruct,noiseFactor,filtComb);                    
                                dataStruct.status = sprintf('%s%s\n',dataStruct.status,...
                                                    'Data was filtered with standard deviation filter.'); 
                              %_____________________________________________________________________________________________  
                             %% end of: apply standard deviation noise mask to the data
                          %------------------------------------------------------------              

                             %% apply median filter to the data
                             case('median-filter')
                                dataStruct = flow_median_filter(dataStruct);
                                dataStruct.status = sprintf('%s%s\n',dataStruct.status,...
                                                    'Data was filtered with median filter.');
                             %% end of: apply median filter to the data
                           %-----------------------------------------------------------  
                             %% apply eddy current correction to the data
                             case('eddy-current')
                                if ((indexStdev == 0)||(indexStdev>indexEddy))
                                    dataStruct = flow_stdev(dataStruct);
                                end

                                dataStruct = flow_eddy_current_corr(dataStruct,staticFactor,actSlice);

                                dataStruct.status = sprintf('%s%s\n',dataStruct.status,...
                                                                'Eddy current correction was done.');

                            %% end of: apply eddy current correction 
                   
                     
                           %%-----------------------------------------------------------
                             %% apply second order correction    
                            case('2nd-order-corr')
                                % FrancescoChange: 2nd-oder-corr -> 2nd-order-corr
                                %disp('Second order correction')

                                if ((indexStdev == 0) || (indexStdev>indexEddy))
                                    dataStruct = flow_stdev(dataStruct);
                                end
                                dataStruct = flow_secOrder_corr(dataStruct,staticFactor, actSlice);

                                dataStruct.status = sprintf('%s%s\n',dataStruct.status,...
                                                      'Second order correction was done.');

                             %% end of : apply second order correction
                                                   
                             
                           %------------------------------------------------------------
                            %% apply anti-aliasing to the data
                            case('anti-aliasing')
                                if strcmp(entryStr, 'dequeueing')
                                    readKeys ={'func settings','','anti_aliasing','','1'};
                                    numOfIter = floor(str2double(char(inifile(dataStruct.inifilePath,'read',readKeys))));                        
                                else
                                    numOfIter  = floor(get(handleStruct.numberOfIterPop,'Value'));
                                end

                                dataStruct = flow_anti_aliasing(dataStruct,numOfIter);
                                %dataStruct = manual_unwrapping(dataStruct);
                                dataStruct.status = sprintf('%s%s\n',dataStruct.status,...
                                                                'Anti-aliasing was done.');
                            %% end of: apply anti-aliasing to the data
                           %----------------------------------------------------------- 

                           %% apply anti-aliasing to the data
                            case('manual anti-aliasing')
                                % read all necessary parameters
                                if m ==1
                                    filename   = fullfile(dataStruct.dataDirectory,'unwrapData.mat');                        
                                    unwrapData = load(filename,'unwrapData');
                                    %% read all slices in which manual phase unwrapping
                                    %% should be performed
                                    allSlices  = cell2mat(unwrapData.unwrapData(1,:));
                                    dataStruct.unwrapData = unwrapData.unwrapData;
                                end

                                % test if manual phase unwraping should be performed on
                                % the actual slice
                                sliceInd = find(allSlices == actSlice);

                                if ~isempty(sliceInd)
                                    for n=1:size(sliceInd,2)
                                        % get direction and actual phase number
                                        actPhaseNum = cell2mat(unwrapData.unwrapData(2,sliceInd(n)));
                                        dirFE = cell2mat(unwrapData.unwrapData(3,sliceInd(n)));                            
                                        actImage = squeeze(dataStruct.imaFLOW(:,:,actPhaseNum,dirFE));

                                        tolerance = 2*(dataStruct.venc(dirFE))/5;
                                        pointList = cell2mat(unwrapData.unwrapData(4,sliceInd(n))); 
                                        [actImage,unwrapInd] = local_unwrap(dataStruct,actImage,tolerance,pointList);
                                        dataStruct.imaFLOW(:,:,actPhaseNum,dirFE)= actImage;
                                        %dataStruct = phase_unwrap_manual(dataStruct,actImage,actPhaseNum,unwrapInd,dirFE);                            
                                    end
                                end
                                dataStruct.status = sprintf('%s%s\n',dataStruct.status,...
                                                                'Manual anti-aliasing was done.');
                            %%% end of: apply anti-aliasing to the data  
                            
                            case('delStatTissue')                                               
                            % delete static tissue if selected                                                            
                                if strcmp(entryStr, 'dequeueing')
                                    readKeys ={'func settings','','delete_stat_tissue','','0'};
                                    delStaticTissue = floor(str2double(char(inifile(dataStruct.inifilePath,'read',readKeys))));
                                else
                                    delStaticTissue = get(handleStruct.staticTissueCheck,'Value');
                                end

                                if delStaticTissue == 1                     
                                        dataStruct.statRegStatus = 'off'; 
                                        dataStruct = flow_stdev(dataStruct);
                                        dataStruct = create_static_tissue_mask(dataStruct,staticFactor,actSlice);
                                        dataStruct = flow_static_tissue_filter(dataStruct);
                                        dataStruct.status = sprintf('%s%s\n',dataStruct.status,'Static tissue was deleted.');                                                    
                                end
                            %% end of: delete static tissue if selected
                                                 
                            otherwise
                                dataStruct.status = dataStruct.status;   
                                     
                                
                end           
            end          
        end %% end of: load all phases from a slice and modify according to the settings
 
       %% mofify data according to selcted filter and correction option //MM
	   %% only do this if NOT in preview mode - cannot go back to originally loaded data
	   %% note: 'createMRstruct' is flag used to trigger the conversion of the data (not ideal name but historic reasons)  
	   if  (strcmp(entryStr, 'createMRstruct')||strcmp(entryStr, 'dequeueing'))
    
           % magnitude data
		   dataStruct.dataMag(:,:,actSlice,:) = dataStruct.imaMAG;  
	       
		   % flow data
		   if strcmp(dataStruct.encodingType,'velocity') 
               % convert velocity data to m/sec 
               dataStruct.dataFlow(:,:,actSlice,:,:) = dataStruct.imaFLOW * 0.01;
           else
               dataStruct.dataFlow(:,:,actSlice,:,:) = dataStruct.imaFLOW;
           end
	   
           % needed to get edges information, used later for ensight conversion //MM
		   if m ==1 %% get properties for new mrstructs
               dataStruct = local_get_mrstruct_properties(dataStruct);
           end

       end 
    
	end  % end of loop over slice
	
	%% now that all data is corrected and filtered, calculate 3D PC-MRA //MM
	if  (strcmp(entryStr, 'createMRstruct')||strcmp(entryStr, 'dequeueing'))
	  if pcmraFlag == 1
		if dataStruct.nFE == 3

          %% get PC-MRA calculation modes and options
		  if strcmp(entryStr, 'dequeue')
			readKeys                = {'modes pcmra','','modeSumSquares','','0'};                
			modeSumSquares          = floor(str2double(char(inifile(dataStruct.inifilePath,'read',readKeys))));
			readKeys                = {'modes pcmra','','modeMeanAbsVel','','0'};
			modeMeanAbsVel          = floor(str2double(char(inifile(dataStruct.inifilePath,'read',readKeys))));
			readKeys                = {'modes pcmra','','modePseudoComplDiff','','0'};
			modePseudoComplDiff     = floor(str2double(char(inifile(dataStruct.inifilePath,'read',readKeys)))); 
			readKeys                = {'modes pcmra','','modeTimeAveMag','','0'};
			modeTimeAveMag          = floor(str2double(char(inifile(dataStruct.inifilePath,'read',readKeys)))); 
            readKeys                = {'modes pcmra','','modeSqrtSumSquares','','0'}; %% LiliMa
			modeSqrtSumSquares      = floor(str2double(char(inifile(dataStruct.inifilePath,'read',readKeys)))); 
            readKeys                = {'modes pcmra','','modeVelCorr','','0'}; %% MBS
			modeVelCorr             = floor(str2double(char(inifile(dataStruct.inifilePath,'read',readKeys)))); 
			
			% get venc threshhold for complex difference
			readKeys ={'options PCMRA','','vencThreshold','','0.0'};
            vencThreshold = str2double(char(inifile(dataStruct.inifilePath,'read',readKeys)));      
			
			% get user selected min and max time frames for teporally constrained PC-MRA //MM
			readKeys ={'options PCMRA','','timeFramesMin','','0.0'};
            tPCMRA_min = str2double(char(inifile(dataStruct.inifilePath,'read',readKeys)));   
			readKeys ={'options PCMRA','','timeFramesMax','','0.0'};
            tPCMRA_max = str2double(char(inifile(dataStruct.inifilePath,'read',readKeys)));    

			% Check if median filter and/or histogram equalization need to be performed for PC-MRA //MM
			readKeys ={'options PCMRA','','medianFilterPCMRA','','0.0'};
            medFilterPCMRA_Flag = str2double(char(inifile(dataStruct.inifilePath,'read',readKeys)));  
            histEqualPCMRA_Flag = 0;
			
		 else
			modeSumSquares          = get(handleStruct.pcmraSumSquaresCheck,'Value');
			modeMeanAbsVel          = get(handleStruct.pcmraMeanAbsVelCheck,'Value');
			modePseudoComplDiff     = get(handleStruct.pcmraPseudoComplDiffCheck,'Value');
			modeTimeAveMag          = get(handleStruct.pcmraTimeAveMagCheck,'Value');
			modeSqrtSumSquares      = get(handleStruct.pcmraSqrtSumSquaresCheck, 'Value'); %% LiliMa
            modeVelCorr             = get(handleStruct.pcmraVelCorrCheck, 'Value'); %% MBS
		    % get venc threshhold for complex difference
		    vencThreshold = str2double(get(handleStruct.vencThresholdEdit, 'String'));
		
		    % get user selected min and max time frames for temporally constrained PC-MRA //MM
		    tPCMRA_min = ceil( str2double(get(handleStruct.timeFramesMinEdit, 'String')) );
		    tPCMRA_max = ceil( str2double(get(handleStruct.timeFramesMaxEdit, 'String')) );	
			
		    % Check if median filter and/or histogram equalization need to be performed for PC-MRA //MM
		    medFilterPCMRA_Flag = get(handleStruct.pcmraMedianFilterCheck,'Value');
		    %% histEqualPCMRA_Flag = get(handleStruct.pcmraHistEqualCheck,'Value');
            histEqualPCMRA_Flag = 0;
		  end
    
          %% wrtite PC-MRA modes to cell array
		  modPCmra ='';
          if (modeSumSquares == 1)               
              modPCmra = [modPCmra;cellstr('modeSumSquares')];
          end
          if (modeMeanAbsVel == 1)
              modPCmra = [modPCmra;cellstr('modeMeanAbsVel')];               
          end
          if (modePseudoComplDiff == 1)
              modPCmra = [modPCmra;cellstr('modePseudoComplexDiff')];
          end
          if (modeTimeAveMag == 1)
              modPCmra = [modPCmra;cellstr('modeTimeAveMag')];
          end
          
          if (modeSqrtSumSquares == 1) %% LiliMa
              modPCmra = [modPCmra;cellstr('modeSqrtSumSquares')];
          end

          if (modeVelCorr == 1) %% MBS
              modPCmra = [modPCmra;cellstr('modeVelCorr')];
          end

		  % check if min / max time frames are out of range
		  if tPCMRA_min < 1 || tPCMRA_min > dataStruct.numPhases
			tPCMRA_min = 1;
			set(handleStruct.timeFramesMinEdit,'String',1);
		  end
		  if tPCMRA_max < 1 || tPCMRA_max > dataStruct.numPhases
			tPCMRA_max = dataStruct.numPhases;
			set(handleStruct.timeFramesMaxEdit,'String',dataStruct.numPhases);
		  end	
		  if tPCMRA_min > tPCMRA_max
			tPCMRA_max = tPCMRA_min;
			set(handleStruct.timeFramesMaxEdit,'String',tPCMRA_max);
		  end	
		
          %% preallocating the data
          if (modeSumSquares == 1)
                dataStruct.volume_pcmraSumSquares = zeros(dataStruct.imageSize(1),dataStruct.imageSize(2),dataStruct.numSlices,'single');
          end
          if (modeMeanAbsVel == 1)
              dataStruct.volume_pcmraMeanAbsVel = zeros(dataStruct.imageSize(1),dataStruct.imageSize(2),dataStruct.numSlices,'single');
          end
          if (modePseudoComplDiff == 1)
              dataStruct.volume_pcmraPseudoComplDiff = zeros(dataStruct.imageSize(1),dataStruct.imageSize(2),dataStruct.numSlices,'single');
          end
          if (modeTimeAveMag == 1)
              dataStruct.volume_pcmraTimeAveMag = zeros(dataStruct.imageSize(1),dataStruct.imageSize(2),dataStruct.numSlices,'single');
          end
          if (modeSqrtSumSquares == 1) %% LiliMa
              dataStruct.volume_pcmraSqrtSumSquares = zeros(dataStruct.imageSize(1),dataStruct.imageSize(2),dataStruct.numSlices,'single');
          end
          if (modeVelCorr == 1) %% MBS
              dataStruct.volume_pcmraVelCorr = zeros(dataStruct.imageSize(1),dataStruct.imageSize(2),dataStruct.numSlices,'single');
          end
          
          % select region of interest if not aready done //MM
          if isempty(dataStruct.pcmraMask)		  
			%% crop data for PC-MRA calculations
			h_pcmra = figure;
			set(h_pcmra,'Name','Crop 4D flow Data for 3D-PC-MRA');
			% calculate MIP of time-averaged magnitude data
			if dataStruct.numSlices > 1
			    minSlc = floor(1 + dataStruct.numSlices/6);
				maxSlc = ceil(dataStruct.numSlices - dataStruct.numSlices/6);
				ImgMag = squeeze( max( mean(dataStruct.dataMag(:,:,minSlc:maxSlc,:),4),[],3) );
			else 
			    ImgMag = squeeze( max( mean(dataStruct.dataMag,4),[],3) );
			end
			h_im = imshow(ImgMag,[]);
			text(0,-10,'Left click, adjust Rectangular ROI and left click twice, and wait for result','FontSize',[10],'Color', 'b');
			rec = imrect;  % get a region R inside a rectangle, BW is a binary image with 1 and 0 inside or outside the rectangle;
			wait(rec);
			ImgMask = createMask(rec,h_im);
			close(h_pcmra);
			dataStruct.pcmraMask = ImgMask;        		
          end

          %% Compute 3D-PC-MRA volume
		  dataStruct = calculate_pcmra(handleStruct,dataStruct,modPCmra,vencThreshold,tPCMRA_min,tPCMRA_max,medFilterPCMRA_Flag,histEqualPCMRA_Flag);
		
		else
           warndlg(lasterr,'Calculation of PCMRA is not possible - not enough flow directions');
		end
      end		   
	end

    close(waitH);
%%
elseif strcmp(entryStr,'previewPDMask')
   waitPD = waitbar(0,'Calculating PD-Maps-MASK. Please Wait ... ');
   PDmaskData = zeros(dataStruct.imageSize(1),dataStruct.imageSize(2),1,dataStruct.numSlices);
   
   for actSlice=1:dataStruct.numSlices 
        %% load Data 
        dataStruct  = local_read_all_phases(dataStruct,actSlice,1);
        %% apply noise mask
        noiseFactor = str2double(get(handleStruct.noiseValueEdt,'String'));
        dataStruct  = flow_noise_filter(dataStruct,noiseFactor,0);
       
        
        if actSlice==1
            %% get which modes of pc-mra calculations shoud be done
  
		    modPCmra ='';
            if (get(handleStruct.pcmraSumSquaresCheck,'Value'))               
                modPCmra = [modPCmra;cellstr('modeSumSquares')];
            end
            if (get(handleStruct.pcmraMeanAbsVelCheck,'Value'))
                modPCmra = [modPCmra;cellstr('modeMeanAbsVel')];               
            end
            if (get(handleStruct.pcmraPseudoComplDiffCheck,'Value'))
                modPCmra = [modPCmra;cellstr('modePseudoComplexDiff')];                
            end
			if (get(handleStruct.pcmraTimeAveMagCheck,'Value'))
                modPCmra = [modPCmra;cellstr('modeTimeAveMagCheck')];                
            end
            if (get(handleStruct.pcmraSqrtSumSquaresCheck,'Value')) %% LiliMa
                modPCmra = [modPCmra;cellstr('modeSqrtSumSquares')];
            end
            if (get(handleStruct.pcmraVelCorrCheck,'Value')) %% MBS
                modPCmra = [modPCmra;cellstr('modeVelCorr')];
            end
        %% end of: get which modes of pc-mra calculations shoud be done
        end
        %% calculate pc-mra according to the choosen modes
        vencThreshold = str2double(get(handleStruct.vencThresholdEdit, 'String'));    
        dataStruct = calculate_pcmra(handleStruct,dataStruct,modPCmra,vencThreshold);
        %% end of: calculate of pc-mra
        
        threshPDmask = str2double(get(handleStruct.thresholdPDmaskEdit,'String'));
        maskData = dataStruct.pcmraSumSquares> (max(dataStruct.pcmraSumSquares(:))*threshPDmask);
        PDmaskData(:,:,1,actSlice)= squeeze(double(maskData));
        
   end
   close(waitPD);

    h = figure;
    
    scrsz = get(0,'ScreenSize');
    set(h,'Position',scrsz,'Name','Preview of PD-MASK');
    iptsetpref('ImshowInitialMagnification','fit');
    state ='true';
    %%JB setWindowOnTop(h,state);
    
    montage(PDmaskData);

    
    impixelinfo; 
    
%% MM: calculate pressure difference maps
elseif strcmp(entryStr,'createPDmaps')
    
    if strcmp(dataValue, 'dequeue')        
        readKeys     = {'modes pcmra','','modePseudoComplDiff','','0'};
        pcmraFlag    = floor(str2double(char(inifile(dataStruct.inifilePath,'read',readKeys)))); 
        readKeys     = {'conversion','','pdmaps','','0'};
        pdFlag       = floor(str2double(char(inifile(dataStruct.inifilePath,'read',readKeys))));
        
        readKeys     = {'func settings','','pdmaskthreshold','','0.3'};
        threshPDmask = str2double(char(inifile(dataStruct.inifilePath,'read',readKeys)));
        
        readKeys     = {'pdmaps settings','','slice','','12'};
        slice = str2double(char(inifile(dataStruct.inifilePath,'read',readKeys)));
        readKeys     = {'pdmaps settings','','row','','88'};
        row = str2double(char(inifile(dataStruct.inifilePath,'read',readKeys)));
        readKeys     = {'pdmaps settings','','col','','55'};
        col = str2double(char(inifile(dataStruct.inifilePath,'read',readKeys)));
    else
        pcmraFlag    = get(handleStruct.pcmraCheck,'Value'); 
        pdFlag       = get(handleStruct.pdCheck,'Value'); 
        threshPDmask = str2double(get(handleStruct.thresholdPDmaskEdit,'String'));
        fname = [dataStruct.dataDirectory,filesep,'pdpoints.mat'];
        points = load_into_variable(fname);
        row = points(1);
        col = points(2);
        slice = points(3);        
    end
    nFE              = dataStruct.nFE;

    
    if pdFlag == 1 && nFE == 3%&& pcmraFlag == 1
        
        numSlices   = dataStruct.numSlices;
        numPhases   = dataStruct.numPhases;
        TR          = dataStruct.TR;
        
        dataStruct.pcmraSumSquares  = [];
        dataStruct.pcmraMeanAbsVel  = [];
        dataStruct.pcmraPseudoComplDiff =[];
        dataStruct.pcmraSqrtSumSquares = []; 
        dataStruct.pcmraVelCorr = []; 
        dataStruct.filterMask       = [];
        dataStruct.statMask         = [];
        dataStruct.statMaskFull     = [];
        dataStruct.stdFlow          = [];
        dataStruct.stdFlowMax       = [];
        dataStruct.imaMAG           = [];
        dataStruct.imaFLOW          = [];
             
        waitPD = waitbar(0,'Calculating PD-Maps. Please Wait ... ');
        %handles = 1;
        
        %% SETUP MASK HERE
 
        %%% load mask from disk
        pathName  = sprintf('%s%s%s',dataStruct.dataDirectory,filesep,'pdmask.mat');        
        dataStruct.pdStruct.MASK =double(load_into_variable(pathName));

        if dataStruct.numSlices>1
            dataStruct.pdStruct = FloodFill3D(dataStruct.pdStruct,row, col, slice);
            % add border region with zeros 
            [sx sy sz] = size(dataStruct.pdStruct.MASK);
            dataStruct.pdStruct.MASK(1,:,:) = 0;
            dataStruct.pdStruct.MASK(:,1,:) = 0;
            dataStruct.pdStruct.MASK(:,:,1) = 0;
            dataStruct.pdStruct.MASK(sx,:,:) = 0;
            dataStruct.pdStruct.MASK(:,sy,:) = 0;
            dataStruct.pdStruct.MASK(:,:,sz) = 0;
            mrstructType ='volume';
%             smoothMat =[3 3 3];
        else
            dataStruct.pdStruct = FloodFill2D(dataStruct.pdStruct,row, col);            
            mrstructType ='image';
        end

        % save mask to disk
        pdStructPathName = sprintf('%s%s%s',dataStruct.mrstructDirPath,filesep,'pd_mask');
        varTemp = dataStruct.pdStruct.MASK;
        save(pdStructPathName,'varTemp');        
%         clear('varTemp');
        
        % save smoothed mask to disk
%         pdStructPathName = sprintf('%s%s%s',dataStruct.mrstructDirPath,filesep,'pdmask_smooth');
%         if dataStruct.numSlices>1
%             varTemp = smooth3((dataStruct.pdStruct.MASK),'gaussian',smoothMat,1.0);
%         else
%             varTemp = (dataStruct.pdStruct.MASK);
%         end
        mrstruct_pdmask_smooth = mrstruct_init(mrstructType,varTemp,dataStruct.propMRstruct);
        mrStructPathName  = sprintf('%s%s%s%03d',dataStruct.mrstructDirPath,filesep,'mrstruct_pdmask_smooth');
        mrstruct_write(mrstruct_pdmask_smooth,mrStructPathName);
        clear('varTemp','mrstruct_pdmask_smooth');
        
        %%%lots of temp things here
        dataStruct.pdStruct.delX   = dataStruct.resY;   %% Grid spacing in x in mm
        dataStruct.pdStruct.delY   = dataStruct.resX;   %% Grid spacing in y in mm
        dataStruct.pdStruct.delZ   = dataStruct.resZ;   %% Grid spacing in z in mm
        dataStruct.pdStruct.tres   = TR;                %% Time spacing in ms
        
        dataStruct.pdStruct.verbose = false;
        dataStruct.pdStruct.GRAD_OUT = 1;
        dataStruct.pdStruct = setup_pressure_static(dataStruct.pdStruct);
        teval  = 2;
        %% do actual pressure map calculation
        for time = 1:numPhases
            if numPhases>1
                if time == 1 
                    tstart = numPhases;
                    tend   = time + 1;                
                elseif time == numPhases
                    tstart = time - 1;
                    tend   = 1;
                else
                    tstart = time - 1;
                    tend   = time + 1;                
                end
            else
                tstart = 1;
                tend = 1;
                dataStruct.pdStruct.tres   = 1;
            end
            
            % update waitbar
            counter = time / numPhases;
            msgOut = sprintf('calculating PD-Map for time frame %s',num2str(time));
            waitbar(counter,waitPD,msgOut)

            %%%%Initialize Structure
            %% load flow data from disk
            if time ==1
                flowData = single(zeros(dataStruct.imageSize(1),dataStruct.imageSize(2),numSlices,3,nFE));
            end
            for actSlice = 1: numSlices
                mrStructPathName      = sprintf('%s%s%s%03d',dataStruct.mrstructDirPath,filesep,'mrstruct_flow_',actSlice);
                tmpMRstruct           = (mrstruct_read(mrStructPathName));
                flowData(:,:,actSlice,1,:) = single(squeeze(tmpMRstruct.dataAy(:,:,tstart,:)));
                flowData(:,:,actSlice,2,:) = single(squeeze(tmpMRstruct.dataAy(:,:,time,:)));
                flowData(:,:,actSlice,3,:) = single(squeeze(tmpMRstruct.dataAy(:,:,tend,:)));
            end
            clear tmpMRstruct
            %Madison Jelena
            if (~strcmp(dataStruct.swVersion,'VA25')&& ~strcmp(dataStruct.swVersion,'VB12')) && ~strcmpi(dataStruct.Manufacturer,'philips')&&isempty(dataStruct.SiemensFlag)%% if software version is at least VB 13
                signVx =  1;
                signVy = -1;
                signVz = -1;
            elseif isempty(dataStruct.SiemensFlag) && strcmpi(dataStruct.Manufacturer,'siemens') %% 20190927 TF non-NU siemens?
                if strcmp(dataStruct.orientation,'cor')
                    signVx = 1;
                    signVy = -1;
                    signVz = 1;
                elseif strcmp(dataStruct.orientation,'sag')
                    signVx = 1;
                    signVy = -1;
                    signVz = 1;
                elseif strcmp(dataStruct.orientation,'tra')
                    signVx = -1;
                    signVy = -1;
                    signVz = 1;
                end
            elseif ~isempty(dataStruct.SiemensFlag) %LiliMa: for Siemens flag
                if strcmp(dataStruct.orientation,'cor')
                    signVx = 1;
                    signVy = -1;
                    signVz = -1;
                elseif strcmp(dataStruct.orientation,'sag')
                    signVx = 1;
                    signVy = -1;
                    signVz = -1;
                elseif strcmp(dataStruct.orientation,'tra')% need more data.
                    signVx = 1;
                    signVy = 1;
                    signVz = 1;
                end
            elseif isempty(dataStruct.SiemensFlag) && strcmpi(dataStruct.Manufacturer,'philips') %% 20190927 TF set to the same direction as NU siemens 
                if strcmp(dataStruct.orientation,'cor')
                    signVx = 1;
                    signVy = -1;
                    signVz = -1;
                elseif strcmp(dataStruct.orientation,'sag')
                    signVx = 1;
                    signVy = -1;
                    signVz = -1;
                elseif strcmp(dataStruct.orientation,'tra')
                    signVx = 1;
                    signVy = 1;
                    signVz = 1;
                end        
            end
            %end Madison Jelena
            
            
            if dataStruct.numSlices==1
                dataStruct.pdStruct.VELXt(:,:,1,:) = signVy*(squeeze(flowData(:,:,:,:,2)));
                dataStruct.pdStruct.VELYt(:,:,1,:) = signVx*(squeeze(flowData(:,:,:,:,1)));
                dataStruct.pdStruct.VELZt(:,:,1,:) = signVz*(squeeze(flowData(:,:,:,:,3)));
            else
                dataStruct.pdStruct.VELXt = single(signVy*(squeeze(flowData(:,:,:,:,2))));
                dataStruct.pdStruct.VELYt = single(signVx*(squeeze(flowData(:,:,:,:,1))));
                dataStruct.pdStruct.VELZt = single(signVz*(squeeze(flowData(:,:,:,:,3))));
            end

            %%%%% Calulation Routines            
            dataStruct.pdStruct = calc_pd_gradient(dataStruct.pdStruct, teval,dataStruct.encodingType);
            dataStruct.pdStruct = update_pressure(dataStruct.pdStruct);
            dataStruct.pdStruct.errorAy(time,:) = dataStruct.pdStruct.tmp;
            % write results of pressure difference calculations to disk
            mrstruct_pd       = mrstruct_init('series2D',dataStruct.pdStruct.PRESSURE_MAT,dataStruct.propMRstruct);
            mrStructPathName  = sprintf('%s%s%s%03d',dataStruct.mrstructDirPath,filesep,'mrstruct_pdiff_',time);
            mrstruct_write(mrstruct_pd,mrStructPathName);
            clear mrstruct_pd
             
            %% write results of pressure gradient calculations to disk
            mrstruct_grad     = mrstruct_init('series2D',dataStruct.pdStruct.GRADx_MAT,dataStruct.propMRstruct);
            mrStructPathName  = sprintf('%s%s%s%03d',dataStruct.mrstructDirPath,filesep,'mrstruct_pGradX_',time);
            mrstruct_write(mrstruct_grad,mrStructPathName);
            clear mrstruct_grad
                        
            mrstruct_grad     = mrstruct_init('series2D',dataStruct.pdStruct.GRADy_MAT,dataStruct.propMRstruct);
            mrStructPathName  = sprintf('%s%s%s%03d',dataStruct.mrstructDirPath,filesep,'mrstruct_pGradY_',time);
            mrstruct_write(mrstruct_grad,mrStructPathName);
            clear mrstruct_grad
            if dataStruct.numSlices>1
                mrstruct_grad     = mrstruct_init('series2D',dataStruct.pdStruct.GRADz_MAT,dataStruct.propMRstruct);
                mrStructPathName  = sprintf('%s%s%s%03d',dataStruct.mrstructDirPath,filesep,'mrstruct_pGradZ_',time);
                mrstruct_write(mrstruct_grad,mrStructPathName);
            end
            clear mrstruct_grad
            
            
        end
        %%write pressure parameters
        write_pressure_param(dataStruct);

        close(waitPD)
        
    else
        dataStruct.status = sprintf('No PD-Map calculation possible');
    end
    
% %% MM: end of: calculate pressure difference maps   

%% convert pc-mra data to dicom files
elseif strcmp(entryStr,'pcmraToDICOM')

    hWait = waitbar(0,'Converting 3D-PC-MRA to dicom ... ');

    warning off all
    DirPath = fullfile(dataStruct.dataDirectory, strcat('PCmra_dicom_',dataStruct.userSubjectIDStr));
    mkdir(DirPath);
    
    if strcmp(entryStr, 'dequeue')
        readKeys                = {'modes pcmra','','modeSumSquares','','0'};                
        modeSumSquares          = floor(str2double(char(inifile(dataStruct.inifilePath,'read',readKeys))));
        readKeys                = {'modes pcmra','','modeMeanAbsVel','','0'};
        modeMeanAbsVel          = floor(str2double(char(inifile(dataStruct.inifilePath,'read',readKeys))));
        readKeys                = {'modes pcmra','','modePseudoComplDiff','','0'};
        modePseudoComplDiff     = floor(str2double(char(inifile(dataStruct.inifilePath,'read',readKeys)))); 
		readKeys                = {'modes pcmra','','modeSqrtSumSquares','','0'};
        modeSqrtSumSquares      = floor(str2double(char(inifile(dataStruct.inifilePath,'read',readKeys))));
        readKeys                = {'modes pcmra','','modeVelCorr','','0'};
        modeVelCorr             = floor(str2double(char(inifile(dataStruct.inifilePath,'read',readKeys)))); 
    else
        modeSumSquares          = get(handleStruct.pcmraSumSquaresCheck,'Value');
        modeMeanAbsVel          = get(handleStruct.pcmraMeanAbsVelCheck,'Value');
        modePseudoComplDiff     = get(handleStruct.pcmraPseudoComplDiffCheck,'Value');
		modeTimeAveMag          = get(handleStruct.pcmraTimeAveMagCheck,'Value');
        modeSqrtSumSquares      = get(handleStruct.pcmraSqrtSumSquaresCheck,'Value');
        modeVelCorr             = get(handleStruct.pcmraVelCorrCheck,'Value');
    end
    
   
    %% write dicom only if any mode of PC-mra was chosen
    if any ([modeSumSquares,modeMeanAbsVel,modePseudoComplDiff,modeTimeAveMag, modeSqrtSumSquares,modeVelCorr]) %% LiliMa
        if dataStruct.eDicom
            patterns = ["FH","AP","RL"];
            FEdirlist = dir(dataStruct.dataDirectory);
            APdir = FEdirlist(contains({FEdirlist.name},patterns(2)));
            pathFlow  = fullfile(APdir.folder,APdir.name,dataStruct.fileNamesFlow); %LiliMa: replaced this with fileNamesFlow because sometimes bSSFP is used which has diff TR
            [~, dicInfoStruct, ~] = dicom_scan_singlefile(pathFlow,'all','dicom-dict_philips.txt'); %LiliMa:some parts of siemens data are not read from dicom scan but instead from read....
            PerframeInfoStruct = dicInfoStruct.PerframeFunctionalGroupsSequence;
            dicInfoStruct = rmfield(dicInfoStruct,'PerframeFunctionalGroupsSequence');
            for slc = 1 : dataStruct.numSlices
                nitem = (slc-1)*dataStruct.numPhases+1;
                dicInfoStruct.PerframeFunctionalGroupsSequence.(['Item_',num2str(slc)]) = PerframeInfoStruct.(['Item_',num2str(nitem)]);
            end
            if modeSumSquares               
                %% create subdirectory for the chosen mode of PC-mra
                subDirPath = fullfile(DirPath, 'mode_MeanSumSquares');
                mkdir(subDirPath);                   
                %% write data to dicom
                fileName  = sprintf('%s%s%s%s%03d%s',subDirPath,filesep,dataStruct.patientName,'_pcmraMeanSumSquares_',1,'.dcm');
                dicInfoStruct.file = fileName;
                multiVolume = reshape(dataStruct.volume_pcmraSumSquares,[dataStruct.imageSize 1 dataStruct.numSlices]);
                dicomwrite(im2uint16(multiVolume), fileName,dicInfoStruct,'CreateMode','Copy','WritePrivate',true);
                dicomwrite(im2uint16(multiVolume), fileName,dicInfoStruct,'Dictionary','dicom-dict_philips.txt','CreateMode','Copy','WritePrivate',true);
            end              
 
            if modeMeanAbsVel
                %% mode 1
                %% create subdirectory for the chosen mode of PC-mra
                subDirPath = fullfile(DirPath, 'mode_MeanAbsVel');
                mkdir(subDirPath);                   
                %% write data to dicom
                fileName  = sprintf('%s%s%s%s%03d%s',subDirPath,filesep,dataStruct.patientName,'_pcmraMeanSysSumSquares_',1,'.dcm');
                dicInfoStruct.file = fileName;
                multiVolume = reshape(dataStruct.volume_pcmraMeanAbsVel,[dataStruct.imageSize 1 dataStruct.numSlices]);
                dicomwrite(im2uint16(multiVolume), fileName,dicInfoStruct,'CreateMode','Copy','WritePrivate',true);
                dicomwrite(im2uint16(multiVolume), fileName,dicInfoStruct,'Dictionary','dicom-dict_philips.txt','CreateMode','Copy','WritePrivate',true);
            end
 
            if modePseudoComplDiff
                subDirPath = fullfile(DirPath, 'mode_PseudoComplDiff');
                %% create subdirectory for the chosen mode of PC-mra
                mkdir(subDirPath);                   
                %% write data to dicom
                fileName  = sprintf('%s%s%s%s%03d%s',subDirPath,filesep,dataStruct.patientName,'_pcmraPseudoComplDiff_',1,'.dcm');
                dicInfoStruct.file = fileName;
                multiVolume = reshape(dataStruct.volume_pcmraPseudoComplDiff,[dataStruct.imageSize 1 dataStruct.numSlices]);
                dicomwrite(im2uint16(multiVolume), fileName,dicInfoStruct,'CreateMode','Copy','WritePrivate',true);                
                dicomwrite(im2uint16(multiVolume), fileName,dicInfoStruct,'Dictionary','dicom-dict_philips.txt','CreateMode','Copy','WritePrivate',true);                
            end
            
            if modeTimeAveMag
                subDirPath = fullfile(DirPath, 'mode_TimeAveMag');
                %% create subdirectory for the chosen mode of PC-mra
                mkdir(subDirPath);
                %% write data to dicom
                fileName  = sprintf('%s%s%s%s%03d%s',subDirPath,filesep,dataStruct.patientName,'volume_pcmraTimeAveMag_',1,'.dcm');
                dicInfoStruct.file = fileName;
                multiVolume = reshape(dataStruct.volume_pcmraTimeAveMag,[dataStruct.imageSize 1 dataStruct.numSlices]);
                dicomwrite(im2uint16(multiVolume), fileName,dicInfoStruct,'CreateMode','Copy','WritePrivate',true);                
                dicomwrite(im2uint16(multiVolume), fileName,dicInfoStruct,'Dictionary','dicom-dict_philips.txt','CreateMode','Copy','WritePrivate',true);                
            end
            %% LiliMa
            if modeSqrtSumSquares
                %% create subdirectory for the chosen mode of PC-mra
                subDirPath = fullfile(DirPath, 'mode_SqrtMeanSumSquares');
                mkdir(subDirPath);
                %% write data to dicom
                fileName  = sprintf('%s%s%s%s%03d%s',subDirPath,filesep,dataStruct.patientName,'_pcmraSqrtMeanSumSquares_',1,'.dcm');
                dicInfoStruct.file = fileName;
                multiVolume = reshape(dataStruct.volume_pcmraSqrtSumSquares,[dataStruct.imageSize 1 dataStruct.numSlices]);
                dicomwrite(im2uint16(multiVolume), fileName,dicInfoStruct,'CreateMode','Copy','WritePrivate',true);                
                dicomwrite(im2uint16(multiVolume), fileName,dicInfoStruct,'Dictionary','dicom-dict_philips.txt','CreateMode','Copy','WritePrivate',true);                
            end
            
            if modeVelCorr
                %% create subdirectory for the chosen mode of PC-mra
                subDirPath = fullfile(DirPath, 'mode_VelCorr');
                mkdir(subDirPath);
                %% write data to dicom
                fileName  = sprintf('%s%s%s%s%03d%s',subDirPath,filesep,dataStruct.patientName,'_pcmraVelocityCorrelation_',1,'.dcm');
                dicInfoStruct.file = fileName;
                multiVolume = reshape(dataStruct.volume_pcmraVelCorr,[dataStruct.imageSize 1 dataStruct.numSlices]);
                dicomwrite(im2uint16(multiVolume), fileName,dicInfoStruct,'CreateMode','Copy','WritePrivate',true);                
                dicomwrite(im2uint16(multiVolume), fileName,dicInfoStruct,'Dictionary','dicom-dict_philips.txt','CreateMode','Copy','WritePrivate',true);                
            end
            msgOut = 'Converting 3D-PC-MRA to dicom ... Please wait.';
            waitbar(1,hWait,msgOut)
        else        
         
            %% write dicom for each slice
            for slc = 1:dataStruct.numSlices     

                counter =  (slc-1) / dataStruct.numSlices;
                msgOut = 'Converting 3D-PC-MRA to dicom ... Please wait.';
                waitbar(counter,hWait,msgOut)
 
                %% get dicom header
                if isempty(dataStruct.SiemensFlag)
                    dirStrFlow      = fullfile(dataStruct.dataDirectory,'flow');
                else
                    dirStrFlow      = dataStruct.subfoldersFlow(1,:);
                end
                fileNameFlowVz  = squeeze(dataStruct.fileNamesFlow(slc,1,1,:));
                pathFlow        = full_filename(dirStrFlow,fileNameFlowVz);
                [infoStruct, dicInfoStruct, msg] = dicom_scan_singlefile(pathFlow,'all','dicom-dict_philips.txt'); %LiliMa:some parts of siemens data are not read from dicom scan but instead from read....
%               dataInfoStruct = dicominfo(pathFlow); % AJB changed this for Phlilps data trx problems, because the dicom-dict is a custom one!!
                %dicInfoStruct   = dicom(pathFlow);
                %% end of: get dicom header
 
                if modeSumSquares               
                    %% create subdirectory for the chosen mode of PC-mra
                    subDirPath = fullfile(DirPath, 'mode_MeanSumSquares');
                    if slc==1                   
                        mkdir(subDirPath);                   
                    end               
                
                    %% write data to dicom
                    fileName  = sprintf('%s%s%s%s%03d%s',subDirPath,filesep,dataStruct.patientName,'_pcmraMeanSumSquares_slice',slc,'.dcm');
                    dicInfoStruct.file = fileName;
                    dicomwrite(im2uint16(squeeze(dataStruct.volume_pcmraSumSquares(:,:,slc))), fileName,dicInfoStruct)%,'WritePrivate',true);                
                end              
 
                if modeMeanAbsVel
                    %% mode 1
                    %% create subdirectory for the chosen mode of PC-mra
                    subDirPath = fullfile(DirPath, 'mode_MeanAbsVel');
                    if slc==1                       
                        mkdir(subDirPath);                   
                    end
                
                    %% write data to dicom
                    fileName  = sprintf('%s%s%s%s%03d%s',subDirPath,filesep,dataStruct.patientName,'_pcmraMeanSysSumSquares_slice',slc,'.dcm');
                    dicInfoStruct.file = fileName;
                    dicomwrite(im2uint16(squeeze(dataStruct.volume_pcmraMeanAbsVel(:,:,slc))), fileName,dicInfoStruct);
                end
 
                if modePseudoComplDiff
                    subDirPath = fullfile(DirPath, 'mode_PseudoComplDiff');
                    %% create subdirectory for the chosen mode of PC-mra
                    if slc==1                   
                        mkdir(subDirPath);                   
                    end
                
                    %% write data to dicom
                    fileName  = sprintf('%s%s%s%s%03d%s',subDirPath,filesep,dataStruct.patientName,'_pcmraPseudoComplDiff_slice',slc,'.dcm');
                    dicInfoStruct.file = fileName;
                    dicomwrite(im2uint16(squeeze(dataStruct.volume_pcmraPseudoComplDiff(:,:,slc))), fileName,dicInfoStruct);
                end
            
                if modeTimeAveMag
                    subDirPath = fullfile(DirPath, 'mode_TimeAveMag');
                    %% create subdirectory for the chosen mode of PC-mra
                    if slc==1
                        mkdir(subDirPath);
                    end
                
                    %% write data to dicom
                    fileName  = sprintf('%s%s%s%s%03d%s',subDirPath,filesep,dataStruct.patientName,'volume_pcmraTimeAveMag',slc,'.dcm');
                    dicInfoStruct.file = fileName;
                    dicomwrite(im2uint16(squeeze(dataStruct.volume_pcmraTimeAveMag(:,:,slc))), fileName,dicInfoStruct);
                end
                %% LiliMa
                if modeSqrtSumSquares
                    %% create subdirectory for the chosen mode of PC-mra
                    subDirPath = fullfile(DirPath, 'mode_SqrtMeanSumSquares');
                    if slc==1
                        mkdir(subDirPath);
                    end
                
                    %% write data to dicom
                    fileName  = sprintf('%s%s%s%s%03d%s',subDirPath,filesep,dataStruct.patientName,'_pcmraSqrtMeanSumSquares_slice',slc,'.dcm');
                    dicInfoStruct.file = fileName;
                    dicomwrite(im2uint16(squeeze(dataStruct.volume_pcmraSqrtSumSquares(:,:,slc))), fileName,dicInfoStruct);%'WritePrivate',true);
                end
            
                if modeVelCorr
                    %% create subdirectory for the chosen mode of PC-mra
                    subDirPath = fullfile(DirPath, 'mode_VelCorr');
                    if slc==1
                        mkdir(subDirPath);
                    end
                
                    %% write data to dicom
                    fileName  = sprintf('%s%s%s%s%03d%s',subDirPath,filesep,dataStruct.patientName,'_pcmraVelocityCorrelation_slice',slc,'.dcm');
                    dicInfoStruct.file = fileName;
                    dicomwrite(im2uint16(squeeze(dataStruct.volume_pcmraVelCorr(:,:,slc))),fileName,dicInfoStruct);%'WritePrivate',true);
                end
            end
        end
    end                          
    
   close(hWait);
   warning on all
 

%% convert mrstruct to the ensight data
elseif strcmp(entryStr,'enSight')

    numSlices   = dataStruct.numSlices;
    numPhases   = dataStruct.numPhases;
    nFE         = dataStruct.nFE;
    h = waitbar(0,'Converting to EnSight Format ... ');
    errStr = [];
    dataStruct.pdStruct=[];
    
    if strcmp(dataValue, 'dequeue')
        readKeys  = {'conversion','','pcmra','','0'};
        pcmraFlag = floor(str2double(char(inifile(dataStruct.inifilePath,'read',readKeys))));
        readKeys                = {'modes pcmra','','modeSumSquares','','0'};
        modeSumSquares          = floor(str2double(char(inifile(dataStruct.inifilePath,'read',readKeys))));
        readKeys                = {'modes pcmra','','modeMeanAbsVel','','0'};
        modeMeanAbsVel          = floor(str2double(char(inifile(dataStruct.inifilePath,'read',readKeys))));
        readKeys                = {'modes pcmra','','modePseudoComplDiff','','0'};
        modePseudoComplDiff   = floor(str2double(char(inifile(dataStruct.inifilePath,'read',readKeys))));
        readKeys                = {'modes pcmra','','modeSqrtSumSquares','','0'};
        modeSqrtSumSquares      = floor(str2double(char(inifile(dataStruct.inifilePath,'read',readKeys)))); %% LiliMa
        readKeys                = {'modes pcmra','','modeVelCorr','','0'};
        modeVelCorr             = floor(str2double(char(inifile(dataStruct.inifilePath,'read',readKeys)))); %% MBS
        
        
        readKeys  = {'conversion','','pdmaps','','0'};
        pdFlag    = floor(str2double(char(inifile(dataStruct.inifilePath,'read',readKeys)))); 
    else        
        pcmraFlag = get(handleStruct.pcmraCheck,'Value');
        modeSumSquares          = get(handleStruct.pcmraSumSquaresCheck,'Value');
        modeMeanAbsVel          = get(handleStruct.pcmraMeanAbsVelCheck,'Value');
        modePseudoComplDiff     = get(handleStruct.pcmraPseudoComplDiffCheck,'Value');
		modeTimeAveMag          = get(handleStruct.pcmraTimeAveMagCheck,'Value');
        modeSqrtSumSquares      = get(handleStruct.pcmraSqrtSumSquaresCheck,'Value'); %% LiliMa
        modeVelCorr             = get(handleStruct.pcmraVelCorrCheck,'Value'); %% MBS
        pdFlag    = get(handleStruct.pdCheck,'Value'); %% MM: add flag for pressure difference data
    end
    
    %% arrange / mirror velocity data according to the patient coordinate
    %% system
    if (~strcmp(dataStruct.swVersion,'VA25')&& ~strcmp(dataStruct.swVersion,'VB12') && ~strcmpi(dataStruct.Manufacturer,'philips'))%% if software version is at least VB 13
        signVx =  1;
        signVy = -1;
        signVz = -1;
    elseif isempty(dataStruct.SiemensFlag) && strcmpi(dataStruct.Manufacturer,'siemens') %% 20190927 TF NU siemens?
        if strcmp(dataStruct.orientation,'cor')  
            signVx = 1;
            signVy = -1;
            signVz = 1;
        elseif strcmp(dataStruct.orientation,'sag')
            signVx = 1;
            signVy = -1;
            signVz = 1;
        elseif strcmp(dataStruct.orientation,'tra')
            signVx = -1;
            signVy = -1;
            signVz = 1;
        end
    elseif ~isempty(dataStruct.SiemensFlag) %LiliMa: for Siemens flag
        if strcmp(dataStruct.orientation,'cor')
            signVx = 1;
            signVy = -1;
            signVz = -1;
        elseif strcmp(dataStruct.orientation,'sag')
            signVx = 1;
            signVy = -1;
            signVz = -1;
        elseif strcmp(dataStruct.orientation,'tra') %LiliMa: need more data. 
            signVx = 1;
            signVy = 1;
            signVz = 1;
        end
    elseif isempty(dataStruct.SiemensFlag) && strcmpi(dataStruct.Manufacturer,'philips') % 20190927 TF set to the same direction as NU
        if strcmp(dataStruct.orientation,'cor')
            signVx = 1;
            signVy = -1;
            signVz = -1;
        elseif strcmp(dataStruct.orientation,'sag')
            signVx = 1;
            signVy = -1;
            signVz = -1;
        elseif strcmp(dataStruct.orientation,'tra') 
            signVx = 1;
            signVy = 1;
            signVz = 1;
        end
    end
    %% end of: arrange / mirror velocity data 
    
    %% preallocating memory
    magData   = zeros(dataStruct.imageSize(1),dataStruct.imageSize(2),numSlices,'single');
    flowData = zeros(dataStruct.imageSize(1),dataStruct.imageSize(2),numSlices,nFE,'single');
	
    %% convert data frame by frame //MM
    for k = 1:numPhases
	
	  counter = ( ( (k-1) * numPhases ) / numPhases );
      msgOut = 'Converting to EnSight Format ... Please wait.';
      waitbar(counter,h,msgOut)
	
      if isempty(errStr)
	      
		magData = squeeze(dataStruct.dataMag(:,:,:,k));
	    flowData = squeeze(dataStruct.dataFlow(:,:,:,k,:));
	  
        %% adjust sign of velocity data
        if nFE==3            
            flowData(:,:,:,1) = squeeze(flowData(:,:,:,1)* signVx);
            flowData(:,:,:,2) = squeeze(flowData(:,:,:,2)* signVy);
        end
        flowData(:,:,:,nFE) = squeeze(flowData(:,:,:,nFE)* signVz);
        %% end of: adjust sign of velocity data
        
        %% reshape flow data according to matlab coordinate system
        flowData = flowData(:, :, :, [2 1 3]);
        %% end of: reshape flow data according to matlab coordinate system
        
        %% init mrstructs
        mrstruct_mag    = mrstruct_init('volume',single(magData),dataStruct.propMRstruct);
        mrstruct_flow   = mrstruct_init('series3D',single(flowData),dataStruct.propMRstruct);
        %% end of: init mrstructs
        
        %% set transformation flag, here edges
        transf_flag =[];
        %% end of: set transformation flag
        
        %% create geo file
        if k == 1         
         geoFileName            = sprintf('%s%s%s%s%s',dataStruct.ensightDirPath,filesep,'EnSight_',dataStruct.userSubjectIDStr,'.geo');
         [res, errStr_geo, oArg]= mrstruct_ensight(mrstruct_flow,'geoFile',geoFileName,1,transf_flag,'create');
        end
        %% end of: create geo file 
        
        %% create ensight files for magnitude and velocity data
        magFileName             = sprintf('%s%s%s%',dataStruct.dataPathName,num2str(k-1,'%02d'),'.mag');
        flowFileName            = sprintf('%s%s%s%',dataStruct.dataPathName,num2str(k-1,'%02d'),'.vel');
        [res, errStr_mag, oArg] = mrstruct_ensight(mrstruct_mag,'dataVolume',magFileName,1,'');
        clear mrstruct_mag
        [res, errStr_vel, oArg] = mrstruct_ensight(mrstruct_flow,'dataVector',flowFileName,1,'',transf_flag,'none');
        clear mrstruct_flow
        errStr                  = [errStr_geo,errStr_mag,errStr_vel];
        %% end of: create ensight files for magnitude and velocity data
        
        %% create ensight files for angiography data if selected
        if pcmraFlag == 1

                if modePseudoComplDiff
                    mrstruct_pcmra            = mrstruct_init('volume',dataStruct.volume_pcmraPseudoComplDiff,dataStruct.propMRstruct);
                    pcmraFileName             = sprintf('%s%s%s%',dataStruct.dataPathName,num2str(k-1,'%02d'),'.spd03');
                    [res, errStr_pcmra, oArg] = mrstruct_ensight(mrstruct_pcmra,'dataVolume',pcmraFileName,1,'');
                    errStr                    = [errStr_geo,errStr_mag,errStr_vel,errStr_pcmra];
                    clear mrstruct_pcmra
                end
                if modeMeanAbsVel
                    mrstruct_pcmra            = mrstruct_init('volume',dataStruct.volume_pcmraMeanAbsVel,dataStruct.propMRstruct);
                    pcmraFileName             = sprintf('%s%s%s%',dataStruct.dataPathName,num2str(k-1,'%02d'),'.spd02');
                    [res, errStr_pcmra, oArg] = mrstruct_ensight(mrstruct_pcmra,'dataVolume',pcmraFileName,1,'');
                    errStr                    = [errStr_geo,errStr_mag,errStr_vel,errStr_pcmra];
                    clear mrstruct_pcmra
                end
                if modeSumSquares                    
                    mrstruct_pcmra            = mrstruct_init('volume',dataStruct.volume_pcmraSumSquares,dataStruct.propMRstruct);
                    pcmraFileName             = sprintf('%s%s%s%',dataStruct.dataPathName,num2str(k-1,'%02d'),'.spd01');
                    [res, errStr_pcmra, oArg] = mrstruct_ensight(mrstruct_pcmra,'dataVolume',pcmraFileName,1,'');
                    errStr                    = [errStr_geo,errStr_mag,errStr_vel,errStr_pcmra];
                    clear mrstruct_pcmra
                end
                if modeTimeAveMag
                    mrstruct_pcmra            = mrstruct_init('volume',dataStruct.volume_pcmraTimeAveMag,dataStruct.propMRstruct);
                    pcmraFileName             = sprintf('%s%s%s%',dataStruct.dataPathName,num2str(k-1,'%02d'),'.spd04');
                    [res, errStr_pcmra, oArg] = mrstruct_ensight(mrstruct_pcmra,'dataVolume',pcmraFileName,1,'');
                    errStr                    = [errStr_geo,errStr_mag,errStr_vel,errStr_pcmra];
                    clear mrstruct_pcmra
                end 
                if modeSqrtSumSquares %%LiliMa
                    mrstruct_pcmra            = mrstruct_init('volume',dataStruct.volume_pcmraSqrtSumSquares,dataStruct.propMRstruct);
                    pcmraFileName             = sprintf('%s%s%s%',dataStruct.dataPathName,num2str(k-1,'%02d'),'.spd05');
                    [res, errStr_pcmra, oArg] = mrstruct_ensight(mrstruct_pcmra,'dataVolume',pcmraFileName,1,'');
                    errStr                    = [errStr_geo,errStr_mag,errStr_vel,errStr_pcmra];
                    clear mrstruct_pcmra
                end
                if modeVelCorr %%MBS
                    mrstruct_pcmra            = mrstruct_init('volume',dataStruct.volume_pcmraVelCorr,dataStruct.propMRstruct);
                    pcmraFileName             = sprintf('%s%s%s%',dataStruct.dataPathName,num2str(k-1,'%02d'),'.spd06');
                    [res, errStr_pcmra, oArg] = mrstruct_ensight(mrstruct_pcmra,'dataVolume',pcmraFileName,1,'');
                    errStr                    = [errStr_geo,errStr_mag,errStr_vel,errStr_pcmra];
                    clear mrstruct_pcmra
                end
        else
            errStr_pcmra ='';
        end
        %% end of: create ensight files for angiography data if selected
        
        %% MM: create ensight files for pressure difference data if selected  
        if pdFlag == 1% && pcmraFlag == 1
            
            %% write pressure difference maps
            mrStructPathName          = sprintf('%s%s%s%03d',dataStruct.mrstructDirPath,filesep,'mrstruct_pdiff_',k);
            mrstruct_pd               = mrstruct_read(mrStructPathName);
            %%convert Pascal to mmHg
            mrstruct_pd.dataAy        = single((mrstruct_pd.dataAy)/133.32);
            pdFileName                = sprintf('%s%s%s%',dataStruct.dataPathName,num2str(k-1,'%02d'),'.pd');
            [res, errStr_pd, oArg]    = mrstruct_ensight(mrstruct_pd,'dataVolume',pdFileName,1,'');
            clear mrstruct_pd
            errStr                    = [errStr_geo,errStr_mag,errStr_vel,errStr_pcmra,errStr_pd];
            
            %% write pressure gradients to ensight
            if k==1
                %% allocate memory for pressure gradients
               pGrad = zeros(dataStruct.imageSize(1),dataStruct.imageSize(2),numSlices,nFE,'single');
               %% write smoothed pressure difference mask to ensight            
                mrStructPathName          = sprintf('%s%s%s%03d',dataStruct.mrstructDirPath,filesep,'mrstruct_pdmask_smooth');
                mrstruct_pdmask           = mrstruct_read(mrStructPathName);
            end            
            
            pdFileName                = sprintf('%s%s%s%',dataStruct.dataPathName,num2str(k-1,'%02d'),'.pdm');
            [res, errStr_pdm, oArg]   = mrstruct_ensight(mrstruct_pdmask,'dataVolume',pdFileName,1,'');
            errStr                    = [errStr_geo,errStr_mag,errStr_vel,errStr_pcmra,errStr_pd,errStr_pdm];
            
            %% gradient in Xs place according to Matlab coordinate system
            
            mrStructPathName          = sprintf('%s%s%s%03d',dataStruct.mrstructDirPath,filesep,'mrstruct_pGradX_',k);
            mrstruct_pd               = mrstruct_read(mrStructPathName);            
            pGrad(:,:,:,2)            = single((squeeze(mrstruct_pd.dataAy)));
            %% gradient in Y 
            mrStructPathName          = sprintf('%s%s%s%03d',dataStruct.mrstructDirPath,filesep,'mrstruct_pGradY_',k);
            mrstruct_pd               = mrstruct_read(mrStructPathName);            
            pGrad(:,:,:,1)            = single((squeeze(mrstruct_pd.dataAy)));
            if dataStruct.numSlices>1
                %% gradient in Z
                mrStructPathName          = sprintf('%s%s%s%03d',dataStruct.mrstructDirPath,filesep,'mrstruct_pGradZ_',k);
                mrstruct_pd               = mrstruct_read(mrStructPathName);                
                pGrad(:,:,:,3)            = single((squeeze(mrstruct_pd.dataAy)));
            else
                pGrad(:,:,:,3)            = zeros;
            end
            
            mrstruct_pgrad            = mrstruct_init('series3D',pGrad,dataStruct.propMRstruct);
            clear pGrad
            
            pgradFileName              = sprintf('%s%s%s%',dataStruct.dataPathName,num2str(k-1,'%02d'),'.pgrad');
            [res, errStr_pd, oArg]     = mrstruct_ensight(mrstruct_pgrad,'dataVector',pgradFileName,1,'',transf_flag,'none');
            errStr                     = [errStr_geo,errStr_mag,errStr_vel,errStr_pcmra,errStr_pd,errStr_pdm];
            clear mrstruct_pgrad
        end
        %% MM: end of: create ensight files for angiography data if selected
%         clear flow        
        %% return if any error is occured
        if ~isempty(errStr)
          dataStruct.status     = sprintf('%s%s%s','Error ',errStr,' has occurred while converting to ensight');
          close (h);
          return
        end
      end
    end

    %% if no error occured
    if isempty(errStr)
        dataStruct.status = sprintf('Conversion of data to EnSight-files was successful');
        close (h);
    end
%% manual phase unwrapping
elseif strcmp(entryStr,'manualPhaseUnwrap')
    button ='ACCEPT & ADD NEW';       
        
        actSliceNum  = ceil(get(handleStruct.SlicesNumSlide,'value'));
        actPhaseNum  = ceil(get(handleStruct.PhasesNumSlide,'value'));
        %% get flow direction
        dirFE = get(handleStruct.changeViewPop,'value');
        
        if isempty(dataStruct.imaFLOW)%% read all phases for the chosen slice
            numSlices  = 1;
            readFlag   = 1;
            dataStruct = local_read_all_phases(dataStruct,actSliceNum,readFlag);
        end

        
        while strcmp(button,'ACCEPT & ADD NEW')
            actImage = squeeze(dataStruct.imaFLOW(:,:,actPhaseNum,dirFE));

            tolerance = 2*(dataStruct.venc(dirFE))/10;
            [actImage,unwrapInd] = local_unwrap(dataStruct,actImage,tolerance,[]);           
            imshow(actImage,[min(min(actImage)),max(max(actImage))]);
            
             button = questdlg('HOW DO YOU WANT TO CONTINUE?','','ACCEPT & QUIT','ACCEPT & ADD NEW','CANCEL','ACCEPT & QUIT');
             if (strcmp(button, 'ACCEPT & QUIT')||strcmp(button, 'ACCEPT & ADD NEW'))
                 
                 %% save chosen points
                 if ~isempty(unwrapInd)
                    dataStruct.unwrapData=[dataStruct.unwrapData,{actSliceNum;actPhaseNum;dirFE;unwrapInd'}];
                 end
                 dataStruct.imaFLOW (:,:,actPhaseNum,dirFE)= actImage;
                 %% unwrap phase
                 %dataStruct = phase_unwrap_manual(dataStruct, actImage,actPhaseNum, unwrapInd,dirFE);                 
                 
                 if strcmp(button, 'ACCEPT & QUIT')
                     orderList = get(handleStruct.orderList,'String');
                     strFlag = any(strcmp(orderList,'manual anti-aliasing'));
                     if ~strFlag
                        func_list = cellstr([orderList;'manual anti-aliasing']); 
                        set(handleStruct.orderList, 'String',func_list);
                     end                    
                    
                    %% save points which should be unwrapped 
                    unwrapData = dataStruct.unwrapData;
                    filename = fullfile(dataStruct.dataDirectory,'unwrapData.mat');
                    save(filename,'unwrapData');
                 end
                 
                 %close(gcf);
             else % cancel 
                button = questdlg('WHICH SELECTION DO YOU WANT TO CANCEL?','','CANCEL LAST & QUIT','CANCEL LAST & CONTINUE','CANCEL ALL','CANCEL LAST & QUIT');
                if strcmp (button,'CANCEL LAST & CONTINUE')
                    button = 'ACCEPT & ADD NEW';
                else

                    if strcmp (button,'CANCEL ALL')
                        dataStruct.unwrapData = [];
                        %%
                        filename = fullfile(dataStruct.dataDirectory,'unwrapData.mat');
                        delete(filename); 
                    end
                        
                    if  isempty(dataStruct.unwrapData)
                        func_list    = get(handleStruct.orderList,'String');     
                        x = strmatch('manual anti-aliasing',func_list);
                        func_list(x) = [];         
                        set(handleStruct.orderList, 'String',func_list,'Value',1);
                        %%reload
                    end 
                end
             end
             close(gcf);
        end
                
    
%% add further data setting to the queue
elseif strcmp(entryStr,'loadInQueue')
    [filename, pathname] = uigetfile({'*.ini';},'Select ini-file');
    local_path = full_filename(pathname,filename);
    %% update queue list only if necessary
    list   = get(handleStruct.queueList,'String');
    [x]=strcmp(list,local_path);
    if (isempty(list)|| not(any(x)))
        dir_list   = cellstr([get(handleStruct.queueList,'String');local_path]);
        set(handleStruct.queueList,'String',dir_list);
    end
    %% end of: update queue list only if necessary

%% take data in queue
elseif strcmp(entryStr,'takeInQueue')
    local_path = strcat(dataStruct.patientName,'_preprocess_set.ini');
    local_path = fullfile(dataStruct.dataDirectory,local_path);    
    inifile(local_path,'new');
    
    %% update queue list only if necessary
    list   = get(handleStruct.queueList,'String');
    [x]=strcmp(list,local_path);
    if (isempty(list)|| not(any(x)))
        dir_list   = cellstr([get(handleStruct.queueList,'String');local_path]);
        set(handleStruct.queueList,'String',dir_list);
    end
    %% end of: update queue list only if necessary
    
    set(handleStruct.queueList,'Enable','on');
    % write data directory path
    writeKeys={'data path','','dirPath',dataStruct.dataDirectory};
               
    % get order of the preprocessing functions           
    order_list = (get(handleStruct.orderList,'String')); 
    if ~isempty(order_list)
        for k=1:length(order_list)            
            func = char(order_list{k});                        
            writeKeys = [writeKeys; {'order','',strcat('func',num2str(k)),func}];           
            
            noiseFilters = 0;

            % get thresholds for the functions
            if strcmp(func,'noise-filter')                
                noiseFilters = noiseFilters + 1;
                threshold = get(handleStruct.noiseValueEdt,'String');
                writeKeys =[writeKeys;{'func settings','','noise_filter',threshold}];                
            elseif strcmp(func,'stdev-filter')                
                noiseFilters = noiseFilters + 1;
                threshold = get(handleStruct.stdevNoiseValueEdt,'String');
                writeKeys =[writeKeys;{'func settings','','stdev_filter',threshold}];
            elseif strcmp(func,'derivative-filter')                
                noiseFilters = noiseFilters + 1;
                threshold = get(handleStruct.derivatNoiseValueEdt,'String');
                writeKeys =[writeKeys;{'func settings','','derivat_filter',threshold}];
            elseif strcmp(func,'anti-aliasing')
                numIter = num2str(get(handleStruct.numberOfIterPop,'Value'),'%2.0f');
                writeKeys =[writeKeys;{'func settings','','anti_aliasing',numIter}];
            elseif (strcmp(func,'eddy-current')||strcmp(func,'2nd-order-corr'))
                threshold = get(handleStruct.staticValueEdt,'String');
                writeKeys =[writeKeys;{'func settings','','static_factor',threshold}];
            end
        end
        
        if (noiseFilters > 1)
            if get(handleStruct.radioANDBut,'Value')==1
                writeKeys =[writeKeys;{'func settings','','filt_comb','AND'}];                    
            else
                writeKeys =[writeKeys;{'func settings','','filt_comb','OR'}];
            end
        else
           writeKeys =[writeKeys;{'func settings','','filt_comb','none'}]; 
        end
        writeKeys = [writeKeys; {'order','','number_func',length(order_list)}];
    end
    
    % get which format should data be converted to
    ensightCheck = num2str(get(handleStruct.enSightCheck, 'Value'));    
    writeKeys =[writeKeys;{'conversion','','ensight',ensightCheck}];
    % is additional angio chosen
    if ~isempty(get(handleStruct.AnatomDataList,'String'))
      writeKeys =[writeKeys;{'conversion','','angio','1'}];  
      list = get(handleStruct.AnatomDataList,'String');
      numDir = size(list,1);
      for j=1:numDir              
          writeKeys =[writeKeys;{'angio settings','',strcat('angio',num2str(j)),list{j}}];
      end
      writeKeys =[writeKeys;{'angio settings','','numberAngio',num2str(numDir)}];
    end
    
    %% is avi movie desired
    aviCheck  = num2str(get(handleStruct.aviMovieCheck, 'Value'));    
    writeKeys = [writeKeys;{'conversion','','avi_movie',aviCheck}];
    
    %% should be mrstruct be returned
    mrstructCheck = num2str(get(handleStruct.selMRstructCheck, 'Value'));    
    writeKeys =[writeKeys;{'conversion','','mrstruct',mrstructCheck}];
    
    %% should static tissue be deleted
    staticCheck = num2str(get(handleStruct.staticTissueCheck, 'Value'));
    writeKeys = [writeKeys;{'conversion','','delete_static_tissue',staticCheck}];
    
    %% convert pcmra to dicom?
    pcmraToDICOMCheck = num2str(get(handleStruct.pcmraToDICOMCheck, 'Value'));
    writeKeys = [writeKeys;{'conversion','','mra2DICOM',pcmraToDICOMCheck}];
    
    if strcmp(pcmraToDICOMCheck,'1')
        pcmraCheck ='1';
    else
        pcmraCheck = num2str(get(handleStruct.pcmraCheck, 'Value')); 
    end
    
    %% additional pcmra 
    writeKeys = [writeKeys;{'conversion','','pcmra',pcmraCheck}];
    if strcmp(pcmraCheck,'1')
       writeKeys = [writeKeys;{'modes pcmra','','modeSumSquares',get(handleStruct.pcmraSumSquaresCheck,'Value')}];
       writeKeys = [writeKeys;{'modes pcmra','','modeMeanAbsVel',get(handleStruct.pcmraMeanAbsVelCheck,'Value')}];
       writeKeys = [writeKeys;{'modes pcmra','','modePseudoComplDiff',get(handleStruct.pcmraPseudoComplDiffCheck,'Value')}];
	   writeKeys = [writeKeys;{'modes pcmra','','modeTimeAveMag',get(handleStruct.pcmraTimeAveMagCheck,'Value')}];
       writeKeys = [writeKeys;{'modes pcmra','','modeSqrtSumSquares',get(handleStruct.pcmraSqrtSumSquaresCheck,'Value')}]; %% LiliMa
       writeKeys = [writeKeys;{'modes pcmra','','modeVelCorr',get(handleStruct.pcmraVelCorrCheck,'Value')}]; %% MBS
       writeKeys = [writeKeys;{'options PCMRA','','vencThreshold',get(handleStruct.vencThresholdEdit,'String')}]; %% // MM
	   writeKeys = [writeKeys;{'options PCMRA','','medianFilterPCMRA',get(handleStruct.pcmraMedianFilterCheck,'Value')}]; %% // MM
	   writeKeys = [writeKeys;{'options PCMRA','','timeFramesMin',get(handleStruct.timeFramesMinEdit,'String')}]; %% // MM
	   writeKeys = [writeKeys;{'options PCMRA','','timeFramesMax',get(handleStruct.timeFramesMaxEdit,'String')}]; %% // MM
    end

    %%  calculate pressure difference mapping 
    pdMapsCheck_FLAG = num2str(get(handleStruct.pdCheck, 'Value'));
    writeKeys = [writeKeys;{'conversion','','pdmaps',pdMapsCheck_FLAG}];
    if pdMapsCheck_FLAG == 1
        threshold = get(handleStruct.thresholdPDmaskEdit,'String');
        writeKeys =[writeKeys;{'func settings','','pdmaskthreshold',threshold}];
        
        fname = [dataStruct.dataDirectory,filesep,'pdpoints.mat'];
        points = load_into_variable(fname);
        
        writeKeys =[writeKeys;{'pdmaps settings','','row',points(1)}];        
        writeKeys =[writeKeys;{'pdmaps settings','','col',points(2)}];        
        writeKeys =[writeKeys;{'pdmaps settings','','slice',points(3)}];
    end
  
    %% write venc and encoding type
    encType = get(handleStruct.encodingTypePop, 'String');
    v = get(handleStruct.encodingTypePop, 'Value');
    writeKeys = [writeKeys;{'encoding sensitivity','','encodingType',num2str(encType{v})}];
    
    writeKeys = [writeKeys;{'encoding sensitivity','','throughplane',dataStruct.venc(3)}];
    writeKeys = [writeKeys;{'encoding sensitivity','','inplaneI',dataStruct.venc(1)}];
    writeKeys = [writeKeys;{'encoding sensitivity','','inplaneJ',dataStruct.venc(2)}];
    inifile(local_path,'write',writeKeys,'tabbed');
 %% end of: take data in queue

%% delete entry from queue
elseif strcmp(entryStr, 'deleteQueueEntry')
    id_value  = get(handleStruct.queueList,'Value');
    dir_list = get(handleStruct.queueList,'String');
    if ~isempty(dir_list)
        dir_list(id_value) = [];
        set(handleStruct.queueList,'String',dir_list,'Value',1);
    end
 %% end of: delete entry from queue

%% append further anatomy data to ensight files    
elseif strcmp(entryStr,'appendToEnsight')
    if strcmp(dataValue, 'dequeue')
        readKeys  = {'angio settings','','numberAngio','','0'};
        numAngio  = floor(str2double(char(inifile(dataStruct.inifilePath,'read',readKeys))));
        dirList   = cell(numAngio,1);
        for j=1:numAngio            
            %% read paths of additional angio conversion from ini file
            readKeys ={'angio settings','',strcat('angio',num2str(j)),'','none'};
            dirList{j} = char(inifile(dataStruct.inifilePath,'read',readKeys));                
        end        
    else
        dirList = get(handleStruct.AnatomDataList,'String');        
    end
    
    numDir  = size(dirList,1);
	% TF added to resolve error in mrStruct
    if strcmpi(dataStruct.Manufacturer, 'philips')
        dicFileName = 'dicom-dict_philips.txt';
    else
        dicFileName = [];
    end
	
    for k=1:numDir
        
        [fileNamesMx, dirStr] = local_get_filelist(dirList{k},'dcm');
        numFiles  = size(fileNamesMx,1);
        edges = zeros(4,4,numFiles);
        anatomTime = nan(numFiles,1); %AJB, need trigger time for SSFP cines

        pathName            = full_filename(dirStr,fileNamesMx(1,:));
%        actMRstruct         = dicom_read_singlefile(pathName,1); %AJB, need header to get trigger time for SSFP cines
            [actMRstruct actMRstruct.user] = dicom_read_singlefile(pathName,1,dicFileName,1); %AJB, need header to get trigger time for SSFP cines
        anatomData          = zeros(size(actMRstruct.dataAy,1),size(actMRstruct.dataAy,2),numFiles,'single');
        edges(:,:,1)        = actMRstruct.edges;
        anatomData(:,:,1)   = single(actMRstruct.dataAy);
        anatomTime(1)       = actMRstruct.user.TriggerTime; %AJB, need trigger time for SSFP cines
        
        % get edges
        if numFiles>1
            for n=2:numFiles
                pathName            = full_filename(dirStr,fileNamesMx(n,:));
%                 actMRstruct         = dicom_read_singlefile(pathName,1); %AJB, need header to get trigger time for SSFP cines
                [actMRstruct actMRstruct.user] = dicom_read_singlefile(pathName,1,dicFileName,1); %AJB, need header to get trigger time for SSFP cines
                edges(:,:,n)        = actMRstruct.edges;
                anatomData(:,:,n)   = single(actMRstruct.dataAy);
                anatomTime(n)       = actMRstruct.user.TriggerTime; %AJB, need trigger time for SSFP cines
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
        actMRstruct.user   = [];  %AJB
        %% calculate edges and write geo file(s)
        actMRstruct.edges  = get_volume_geo(edges);
		
        %initialize mrstruct
        mrstruct_anatom = mrstruct_init(struct_type,anatomData,actMRstruct);
        if numFiles>1
            mrstruct_anatom.user = anatomTime; %AJB
        end
        
        dataStruct.filterMask = [];
        dataStruct.imaFLOW    = [];
        dataStruct.imaMAG     = [];
%        clear('anatomData','fileNamesMx','actMRstruct','pathName','edges'); %AJB
        clear('anatomData','anatomTime','fileNamesMx','actMRstruct','pathName','edges'); %AJB
                       
        % scale data to range [0 1]
        maxValue = max(max(max(mrstruct_anatom.dataAy)));
        minValue = min(min(min(mrstruct_anatom.dataAy)));        
        mrstruct_anatom.dataAy = (mrstruct_anatom.dataAy - minValue)/(maxValue-minValue);
            
        if strcmp(struct_type,'series2D') %AJB
            %function call for to write a separate case file for 2D cine (time-resolved) data (AJB)
            write_ensightCine(dataStruct,mrstruct_anatom)
        else
        
        
        %% write geo file
        geoFileName            = sprintf('%s%s%s%s%s',dataStruct.ensightDirPath,filesep,'EnSight_',dataStruct.patientName,'.geo');
        [res, errStr_geo, oArg]= mrstruct_ensight(mrstruct_anatom,'geoFile',geoFileName,k+1,[],'append');


        %% write mrstruct to ensight as a many times as number of
        %% phases        
        anatomFileName  = sprintf('%s%s%d',dataStruct.dataPathName,'.an',k);
        [res, errStr, oArg] = mrstruct_ensight(mrstruct_anatom,'dataVolume',anatomFileName,k+1,'');%,[],'none');
		end
    end
    
%% avi movie creation
elseif strcmp(entryStr,'avimovie')
    
    if dataStruct.numPhases>1
        hWait = waitbar(0,'Creating avi Movie ... ');
        
        for j=1:dataStruct.numSlices
            
            %% update waitbar
            msgOut = ['Creating avi Movie, slice #',num2str(j)];
            counter = ((j-1)/dataStruct.numSlices);
            waitbar(counter,hWait,msgOut);
            
            imaMag = squeeze(dataStruct.dataMag(:,:,j,:));
            
            %% scale magnitude images to [0 1]
            minMag     = min(imaMag(:));
            maxMag     = max(imaMag(:));
            imaMag     = (imaMag - minMag)/(maxMag-minMag);
            
            imaFlow = squeeze(dataStruct.dataFlow(:,:,j,:,:));
            
            %% scale velocity images to [0 1]
            minFlow     = min(imaFlow(:));
            maxFlow     = max(imaFlow(:));
            imaFlow     = (imaFlow - minFlow)/(maxFlow-minFlow);
            
            %% initialize avi movie
            if j==1
                h = figure;
                %             scrsz = get(0,'ScreenSize');
                %             set(h,'Position',scrsz);
                %             iptsetpref('ImshowInitialMagnification','fit');
                %             state ='true'; %commented out by AJB because of incompatibility with Matlab 2012a
                %             setWindowOnTop(h,state); %commented out by AJB because of incompatibility with Matlab 2012a
                %% set movie path
                aviStr = sprintf('%s%s%s%s',dataStruct.dataDirectory,filesep,dataStruct.userSubjectIDStr,'.avi');
                
                % open movie file
                %             aviobj = avifile(aviStr); % AJB: win7 64 bit does not support any of the older compression algorithms, and avifile is set for removal in future matlab releases
                aviobj = VideoWriter(aviStr,'Motion JPEG AVI'); % AJB: win7 64 bit does not support any older compression algorithms, and avifile is set for removal in future matlab releases
                %             aviobj.Quality = 100;% AJB
                %             aviobj.Compression = 'Cinepak'; % AJB win7 does not support cinepak compression
                aviobj.Quality   = 92; % AJB, otherwise 100 is much, much bigger
                aviobj.FrameRate = 15; % AJB
                open(aviobj) % AJB
            end
            %% end of:initialize avi movie
            
            %% create array of image frames
            imaOut       = cat(2,squeeze(imaMag(:,:,:)),squeeze(imaFlow(:,:,:,1)));
            if dataStruct.nFE==3
                dataV2V3    = cat(2,squeeze(imaFlow(:,:,:,2)),squeeze(imaFlow(:,:,:,3)));
                imaOut  = cat(1, imaOut, dataV2V3);
            end
            %% end of: create array of image frames
            
            %% show each frame and add it to the avi-movie
            for m = 1:dataStruct.numPhases
                %% clear axes
                cla;
                %% show each frame
                imshow(imaOut(:,:,m),'InitialMagnification','fit');
                axis off;
                colormap(gray);
                drawnow;
                
                % add frame to avi movie
                frame  = getframe(gca);
                %             aviobj = addframe(aviobj,frame); % AJB: win7x64 does not support any compression algorithms for avi addframe
                writeVideo(aviobj,frame); % AJB: win7x64 does not support any compression algorithms for avi addframe
                
            end
            %% end of: show each frame and add it to the avi-movie
        end
        close(h);
        close(aviobj);
        close(hWait);
        dataStruct.status = sprintf('%s\n%s',dataStruct.status,'avi-movie was created.');
        drawnow;
    else
        h = msgbox('Only one timepoint, no movie created.','Error: AVI movie');
    end
%%end of: creating avi-movie 
%%
elseif strcmp(entryStr,'selStruct')
    
	h = waitbar(0,'Exporting MR-Struct Data ...');
 
    %% calculate number of slices and time frames
    numSlices   = dataStruct.numSlices;
    numPhases   = dataStruct.numPhases;

    %% start TF modified to correctly output size_t_step for Philips data 
    if strcmpi(dataStruct.Manufacturer,'philips')
        size_t_step_i      = diff(sort(dataStruct.TimeStamps,2),[],2);
		size_t_step_i(:,1) = [];
		size_t_step        = mean(mean(size_t_step_i,1)); 
    elseif ~isempty(dataStruct.SiemensFlag) || strcmpi(dataStruct.Manufacturer,'siemens')  
        size_t_step = dataStruct.TR; 
    else
        size_t_step_i      = diff(sort(dataStruct.TimeStamps,2),[],2);
		size_t_step_i(:,1) = [];
		size_t_step        = mean(mean(size_t_step_i,1)); 
    end
    %% end   TF modified
    
    %% try to preallocate memory                    
    try
        if dataStruct.nFE==3
                    %%init mrstructs               		
					%dataStruct.propMRstruct.user  = struct('size_t_step', dataStruct.TR,'encoding_type',dataStruct.encodingType,'venc_in_plane', dataStruct.venc(1),'venc_through_plane', dataStruct.venc(3));
                    if isfield(dataStruct,'nominalInterval') %% Had to add test for Philips data //AJB
                        dataStruct.propMRstruct.user  = struct('size_t_step', dataStruct.TR,'encoding_type',dataStruct.encodingType,'venc_in_plane', dataStruct.venc(1),'venc_through_plane', dataStruct.venc(3),'nominal_interval',dataStruct.nominalInterval);
                    else
                        dataStruct.propMRstruct.user  = struct('size_t_step', size_t_step,  'encoding_type',dataStruct.encodingType,'venc_in_plane', dataStruct.venc(1),'venc_through_plane', dataStruct.venc(3));
                    end
					magStruct    = mrstruct_init('volume',single(dataStruct.dataMag),dataStruct.propMRstruct);
                    flowStruct   = mrstruct_init('series3D',single(dataStruct.dataFlow),dataStruct.propMRstruct);
                    %% end of: init mstructs
        elseif dataStruct.nFE==1
                    %%init mrstructs    
                    %dataStruct.propMRstruct.user  = struct('size_t_step', dataStruct.TR,'encoding_type',dataStruct.encodingType,'venc_through_plane', dataStruct.venc(1))
                    dataStruct.propMRstruct.user  = struct('size_t_step', dataStruct.TR,'encoding_type',dataStruct.encodingType,'venc_through_plane', dataStruct.venc(1),'nominal_interval',dataStruct.nominalInterval);
                    magStruct    = mrstruct_init('series3D',single(dataStruct.dataMag),dataStruct.propMRstruct);
                    flowStruct   = mrstruct_init('series3D',single(dataStruct.dataFlow),dataStruct.propMRstruct);
        end
    catch
         errordlg('An error  occured while allocating memory for mrstruct export')
         return
    end
    %% end of: try to preallocate memory 
        
    magStruct.dataAy  = single(dataStruct.dataMag);
    flowStruct.dataAy = single(dataStruct.dataFlow); 
  
    % arrange / mirror velocity data according to the main
    % orientation and in-plane phase encoding direction
    if (~strcmp(dataStruct.swVersion,'VA25')&& ~strcmp(dataStruct.swVersion,'VB12')) && ~strcmpi(dataStruct.Manufacturer,'philips')&&isempty(dataStruct.SiemensFlag)%% if software version is at least VB 13
        signVx =  1;
        signVy = -1;
        signVz = -1;
    elseif isempty(dataStruct.SiemensFlag) && strcmpi(dataStruct.Manufacturer,'siemens') %% 20190927 TF not NU-Siemens?
        if strcmp(dataStruct.orientation,'cor')
            signVx = 1;
            signVy = -1;
            signVz = 1;
        elseif strcmp(dataStruct.orientation,'sag')
            signVx = 1;
            signVy = -1;
            signVz = 1;
        elseif strcmp(dataStruct.orientation,'tra')
            signVx = -1;
            signVy = -1;
            signVz = 1;
        end
    elseif ~isempty(dataStruct.SiemensFlag)
        if strcmp(dataStruct.orientation,'cor')
            signVx = 1;
            signVy = -1;
            signVz = -1;
        elseif strcmp(dataStruct.orientation,'sag')
            signVx = 1;
            signVy = -1;
            signVz = -1;
        elseif strcmp(dataStruct.orientation,'tra') %LiliMa: need more data.
            signVx = 1;
            signVy = 1;
            signVz = 1;
        end
    elseif isempty(dataStruct.SiemensFlag) && strcmpi(dataStruct.Manufacturer,'philips') %% 20190927 TF set to the same direction as NU siemens 
        if strcmp(dataStruct.orientation,'cor')
            signVx = 1;
            signVy = -1;
            signVz = -1;
        elseif strcmp(dataStruct.orientation,'sag')
            signVx = 1;
            signVy = -1;
            signVz = -1;
        elseif strcmp(dataStruct.orientation,'tra')
            signVx = 1;
            signVy = 1;
            signVz = 1;
        end        
    end
    
    % arrange / mirror velocity data
    if dataStruct.nFE==3
        %% permute data in to right order [x,y,slices,nFE,time frames]        
        flowStruct.dataAy = permute(flowStruct.dataAy,[1 2 3 5 4]);
        flowStruct.dataAy(:,:,:,1,:) = signVx*flowStruct.dataAy(:,:,:,1,:);
        flowStruct.dataAy(:,:,:,2,:) = signVy*flowStruct.dataAy(:,:,:,2,:);
        flowStruct.dataAy(:,:,:,3,:) = signVz*flowStruct.dataAy(:,:,:,3,:);
        %% permute data according to matlab directions
        flowStruct.dataAy=flowStruct.dataAy(:,:,:,[2 1 3],:);
    elseif dataStruct.nFE ==1
        flowStruct.dataAy = signVz*(flowStruct.dataAy);
    end
       
    %% convert mrstruct if data is 2D
    if numSlices == 1
        magStruct  = mrstruct_convert(magStruct,'series2D');
        if dataStruct.nFE==3
            flowStruct = mrstruct_convert(flowStruct,'series2DEchos');
        elseif dataStruct.nFE==1
            flowStruct = mrstruct_convert(flowStruct,'series2D');
        end
    end
    %%% end of: convert mrstruct if data is 2D
    
    %% convert mrstruct if data is 2D
    if numPhases == 1
        magStruct  = mrstruct_convert(magStruct,'volume');
        if dataStruct.nFE==3
            flowStruct = mrstruct_convert(flowStruct,'volumeEchos');
        elseif dataStruct.nFE==1
            flowStruct = mrstruct_convert(flowStruct,'volume');
        end
    end
    %%% end of: convert mrstruct if data is 2D
    
	% create new directory for mrstruct files
	mrStrctDirStr = ['mrstruct_',dataStruct.userSubjectIDStr];
    dataStruct.mrstructDirPath = fullfile(dataStruct.dataDirectory,mrStrctDirStr);
    mkdir(dataStruct.mrstructDirPath);
    
    msgOut = 'Writing MR-Struct Data to Disk... ';
    waitbar(1,h,msgOut);
    
    %% 19/3/2018 LiliMa: save time stamps for use for retrospectively gated and CS data
    if numPhases >1
		%% start TF modified to correctly output timeStamps for Philips data 
        if strcmpi(dataStruct.Manufacturer,'philips')
		    timeStamps = mean(sort(dataStruct.TimeStamps, 2), 1);
		elseif ~isempty(dataStruct.SiemensFlag) || strcmpi(dataStruct.Manufacturer,'siemens') % AJB need to add feature for Philips timestamps 
            timeStamps   = (1:numPhases)*dataStruct.TR - dataStruct.TR/2; %LiliMa: time stamps; may need to change this
        else
            timeStamps = mean(sort(dataStruct.TimeStamps, 2), 1); %AJB: Need to fill this in where it occurs in the rest of the code
        end
		%% end   TF modified
        
        flowStruct.user.TimeStamps = timeStamps;
        magStruct.user.TimeStamps = timeStamps;
        flowStruct.user.TD = flowStruct.user.TimeStamps(3)-flowStruct.user.TimeStamps(2); %HH: time difference, for envelop calculation
        magStruct.user.TD = magStruct.user.TimeStamps(3)-magStruct.user.TimeStamps(2);
    end
    %% write mrstructs to the disk and save the path of mrstructs in dataStruct
    mrStructPathName  = sprintf('%s%s%s',dataStruct.mrstructDirPath,filesep,'mag_struct');
    mrstruct_write(magStruct,mrStructPathName);
    dataStruct.returnPath = [mrStructPathName;[]];
    
    mrStructPathName  = sprintf('%s%s%s',dataStruct.mrstructDirPath,filesep,'vel_struct');
    mrstruct_write(flowStruct,mrStructPathName);
    dataStruct.returnPath(2,:) = mrStructPathName;
    
    %if Auto Seg checkbox is flagged, generate and save Auto segmentation and preview; %% //HB
    if get(handleStruct.autoseg,'Value') == 1
        cur_dir = pwd;
        
        sitelist = get(handleStruct.site,'String');
        site = sitelist{get(handleStruct.site,'Value')};
        verlist = get(handleStruct.autosegvers,'String');
        ver = verlist{get(handleStruct.autosegvers,'Value')};
        
        % folder suffix
        %versuffix = 'CHCOversion';
        % search latest version on the server
        path_to_exe = dataStruct.serverPath; 
        %dirlist = dir(path_to_exe);
        %dirlist = dirlist(~ismember({dirlist.name},{'.','..'}));
        %dirlist = dirlist([dirlist.isdir]);
        %dirname = {dirlist.name};
        %CHCOvers = extractAfter(dirname,versuffix);
        %CHCOvernum = cellfun(@str2num,CHCOvers,'un',0);
        %latest = max(cell2mat(CHCOvernum));
        %serverPath = fullfile(path_to_exe,[versuffix, num2str(latest,'%.1f')]);
        %velomappath = which("velomap_internal");
        
        serverPath = fullfile(path_to_exe,[site,'version',ver]);

        
        % See if there is a local copy of the weights
        if exist('C:/temp/auto_seg.cfg','file') == 2
            fileID = fopen('C:/temp/auto_seg.cfg','r');
            autoSegPath = fscanf(fileID,'%s');
            fclose(fileID);
            if exist([autoSegPath filesep() weightsVersion],'dir') ~= 7
                % Transfer the files
                copyfile(serverPath,[autoSegPath filesep() weightsVersion]);
            end
            % Change to the autoseg files path
            cd([autoSegPath filesep() weightsVersion]);
        else
            cd(serverPath);
        end
                
        
        %cd('\\10.254.136.37\data_imaging\cv_mri\MachineLearning\aorta_seg\version2');
        path = sprintf('%s%s%s',dataStruct.mrstructDirPath);
        %str = ['aorta_seg.exe --path="' path '"'];
        arg_pt  = ['_',dataStruct.patientName];
        arg_usr = ['_user',dataStruct.userIDStr];
        arg_ver = ['_',site,'ver',num2str(ver,'%.1f')];
        outdir = ['ML_mrStruct',arg_pt,arg_usr,arg_ver];
        str = ['VesselAutoSeg.exe ' '"' path '" ' '"' arg_pt '" ' '"' arg_usr '" ' '"' arg_ver '"'];
        system(str);
        
        cd(cur_dir)
        colors = jet(2);
        [p,~,~] = fileparts(path);
        data = load(fullfile(p, outdir,'mask_struct.mat'));
        % create slicer dir
        pathslicer = fullfile(p,outdir,'slicer');
        if exist(pathslicer,'dir') ~= 7
            mkdir(pathslicer);
        end
        % prepare temporary struct
        tmpStruct = data.mrStruct;
        tmpStruct.vox = double(tmpStruct.vox);
        tmpStruct.userid = 'ML';
        tmpStruct.Version = 3.0;
        tmpStruct.Date_Created = date;
        
        % The number of parts and names
        npart = max(max(max(data.mrStruct.dataAy)));
        partname = {'aortaml','pulmonaryml'};
        labelfnames = cell(1,npart);
        
        % Create figure
        hb = figure('Name','Segmentation Preview');
        hold on;

        % for each part
        for ipart = 1 : npart
            data2 = (data.mrStruct.dataAy);
            data2(data2 ~= ipart) = 0; 
            datatry = bwareaopen(data2,1000,6);
            % exeption 
            if isempty(datatry(datatry~=0)) % 
                disp("warning: auto-segmentation was swept out by bwareaopen function. use original geometry.");
            else
                data2 = datatry;
            end    
            % if it still contains no mask, skip the following process. 
            if isempty(data2(data2~=0)) 
                disp(['Auto-segmentation of ',num2str(ipart,'%0d') ,'th part was failed. Skip the following process for autoseg.']);
                continue;
            end
            %data2 = imfill(data2, 4, 'holes');
            top = false(size(data2,1),size(data2,2));
            data3 = cat(3,top,data2,top);
            m1patch = isosurface(data3,0.5);
            m1vox = double(data.mrStruct.vox);
            m1patch.vertices(:,1) = m1patch.vertices(:,1)*m1vox(2);
            m1patch.vertices(:,2) = m1patch.vertices(:,2)*m1vox(1);
            m1patch.vertices(:,3) = m1patch.vertices(:,3)*m1vox(3);
            mss = patch(m1patch);
            mss.FaceColor = colors(ipart,:);
            mss.EdgeColor = 'none';
            mss.FaceAlpha = 0.7;
            axis equal ij
            axis off
            view(0,90)
            camlight('right')
            cameratoolbar('Show')
            
            % output nrrd file for slicer labels
            tmpStruct.dataAy = double(data2);
            partcolor = 5 + (ipart-1)*2;
            labelfnames{ipart} = mrstruct_to_nrrd(tmpStruct,partcolor,'short',3,'left-posterior-superior',...
                                                'domain','little','gzip',partname{ipart},pathslicer);
            mrStruct = tmpStruct;
            % save maskStruct
            save(fullfile(p,outdir,[partname{ipart},'_mask_struct.mat']), 'mrStruct');
            % save EnSight Case file
            mrStruct = rmfield(mrStruct,{'userid','Version','Date_Created'});
            if ipart == 1 
                ensightDirPath = dataStruct.ensightDirPath; 
                FileCase   = ['EnSight_',dataStruct.userSubjectIDStr,'.case'];   
                [~,~] = append_data_ensight(mrStruct,magStruct, FileCase,ensightDirPath,[FileCase(1:end-5),'_mlseg.case'],partname{ipart}, '.seg');
            else
                [~,~] = append_data_ensight(mrStruct,magStruct, [FileCase(1:end-5),'_mlseg.case'],ensightDirPath,[FileCase(1:end-5),'_mlseg.case'],partname{ipart}, '.seg');
            end    
        end
        saveas(hb,fullfile(p,outdir,'Seg_Preview.fig'));
        hold off;
        
        % output pcmra for slicer volume
        modeSumSquares      = get(handleStruct.pcmraSumSquaresCheck,'Value');
        modeMeanAbsVel      = get(handleStruct.pcmraMeanAbsVelCheck,'Value');
        modePseudoComplDiff = get(handleStruct.pcmraPseudoComplDiffCheck,'Value');
        modeTimeAveMag      = get(handleStruct.pcmraTimeAveMagCheck,'Value');
        modeSqrtSumSquares  = get(handleStruct.pcmraSqrtSumSquaresCheck, 'Value');
        modeVelCorr         = get(handleStruct.pcmraVelCorrCheck, 'Value');

        nvol = modeSumSquares + modeMeanAbsVel + modePseudoComplDiff + ...
               modeTimeAveMag + modeSqrtSumSquares + modeVelCorr;
        volfnames = cell(1,nvol);

        n = 0;
        if modeSumSquares
            n = n + 1;
            tmpStruct.dataAy = im2uint16(dataStruct.volume_pcmraSumSquares);
            volfnames{n} = mrstruct_to_nrrd(tmpStruct,'none','unsigned short',3,'left-posterior-superior', ...
                                          'domain','little','gzip','pcmraSumSquares',pathslicer);
        end

        if modeMeanAbsVel
            n = n + 1;
            tmpStruct.dataAy = im2uint16(dataStruct.volume_pcmraMeanAbsVel);
            volfnames{n} = mrstruct_to_nrrd(tmpStruct,'none','unsigned short',3,'left-posterior-superior', ...
                                          'domain','little','gzip','pcmraMeanAbsVel',pathslicer);
        end
        
        if modePseudoComplDiff
            n = n + 1;
            tmpStruct.dataAy = im2uint16(dataStruct.volume_pcmraPseudoComplDiff);
            volfnames{n} = mrstruct_to_nrrd(tmpStruct,'none','unsigned short',3,'left-posterior-superior', ...
                                          'domain','little','gzip','pcmraPseudoComplDiff',pathslicer);
        end
        
        if modeTimeAveMag
            n = n + 1;
            tmpStruct.dataAy = im2uint16(dataStruct.volume_pcmraTimeAveMag);
            volfnames{n} = mrstruct_to_nrrd(tmpStruct,'none','unsigned short',3,'left-posterior-superior', ...
                                          'domain','little','gzip','pcmraTimeAveMag',pathslicer);
        end

        if modeSqrtSumSquares
            n = n + 1;
            tmpStruct.dataAy = im2uint16(dataStruct.volume_pcmraSqrtSumSquares);
            volfnames{n} = mrstruct_to_nrrd(tmpStruct,'none','unsigned short',3,'left-posterior-superior', ...
                                          'domain','little','gzip','pcmraSqrtSumSquares',pathslicer);
        end

        if modeVelCorr
            n = n + 1;
            tmpStruct.dataAy = im2uint16(dataStruct.volume_pcmraVelCorr);
            volfnames{n} = mrstruct_to_nrrd(tmpStruct,'none','unsigned short',3,'left-posterior-superior', ...
                                          'domain','little','gzip','pcmraVelCorr',pathslicer);
        end
        create_slicer_mrml(tmpStruct,labelfnames,volfnames,'AutoSeg',pathslicer);
        clear tmpStruct;
     end
    
    close(h);
    
%% detect invalid entry string
else
    msg = sprintf('Sorry, entry switch ''%s'' was not recognized in ''local_set_dataStruct_entry'', no action performed', entryStr);
    warndlg(lasterr,msg);
end
%%% End of: detect invalid entry string

%%%%% End of: work on dataStruct value and update gui
%%%%% 


%% update gui
function  [dataStruct,handleStruct] = local_update_GUI(handleStruct,dataStruct,entryStr,dataValue)
%%% START set block
if strcmp(entryStr, 'loadData')    
    if (~isempty(dataStruct.fileNamesFlow) && strcmp(dataStruct.loadStatus,'ok'))        
        numPhases = dataStruct.numPhases; 
        numSlices = dataStruct.numSlices;        
        patientName =  sprintf('%s%s','Patient : ',dataStruct.patientName);
        set(handleStruct.patientTxt,'String',patientName);
        set(handleStruct.dataInfoEdt,'String',dataStruct.dataInfo);
        
        if dataStruct.nFE ==1
            set(handleStruct.vencThPlaneEdt,'String',num2str(dataStruct.venc(1)),'Enable', 'on');
            set(handleStruct.vencInPlaneXEdt,'String','1','Enable', 'off');
            set(handleStruct.vencInPlaneYEdt,'String','1','Enable', 'off');
            set(handleStruct.pcmraCheck,     'Enable','off','Value',0);
            set(handleStruct.enSightCheck,    'Enable','off','Value',0);            
            set(handleStruct.changeViewPop,    'Enable','off','Value',3);            
        elseif dataStruct.nFE==3
            set(handleStruct.vencInPlaneXEdt,'String',num2str(dataStruct.venc(1)),'Enable', 'on');
            set(handleStruct.vencInPlaneYEdt,'String',num2str(dataStruct.venc(2)),'Enable', 'on');
            set(handleStruct.vencThPlaneEdt,'String',num2str(dataStruct.venc(3)),'Enable', 'on');
            set(handleStruct.pcmraCheck,        'Enable','on','Value',0);
            set(handleStruct.autoseg,        'Enable','on','Value',0); %% //HB
            
            if exist(dataStruct.serverPath,'dir') == 7
                set(handleStruct.autoseg,    'Enable','on');
                set(handleStruct.site,       'Enable','on');
                set(handleStruct.autosegvers,'Enable','on');
                handleStruct = local_set_autosegvers(handleStruct,dataStruct.serverPath);
                set(handleStruct.autosegvers,'Value',1);
            else
                set(handleStruct.autoseg,    'Enable','off');
                set(handleStruct.site,       'Enable','off');
                set(handleStruct.autosegvers,'Enable','off');
            end            
            
            if dataStruct.externflag ~=1
                set(handleStruct.pdCheck,           'Enable','on','Value',0);
            end
            set(handleStruct.enSightCheck,      'Enable','on','Value',0);
            set(handleStruct.changeViewPop,    'Enable','on','Value',3);
        end
        if(numSlices > 1)
            set(handleStruct.SlicesNumSlide,'Enable', 'on','String',1,'Value',1,'Min',1,'Max',numSlices,...
                'SliderStep',[1/(numSlices-1) 1/(numSlices-1)]);
        else
            set(handleStruct.SlicesNumSlide,'Enable', 'off','Value',1);
        end
         if(numPhases > 1)
            set(handleStruct.PhasesNumSlide,'String',1,'Value',1,'Min',1,'Max',numPhases,...
            'SliderStep',[1/(numPhases-1) 1/(numPhases-1)],'Enable','on');
         else
             set(handleStruct.PhasesNumSlide,'Enable', 'off','Value',1);
         end
        set(handleStruct.SlicesNumEdit,'String',1);
        set(handleStruct.PhasesNumEdit,'String',1);
        set(handleStruct.previewBut,        'Enable','on');
        set(handleStruct.outToExcelBut,     'Enable','on');
        set(handleStruct.saveFuncSetBut,    'Enable','on');
		set(handleStruct.timeFramesMinEdit,    'String',1);
		set(handleStruct.timeFramesMaxEdit,    'String',numPhases);
		set(handleStruct.vencThresholdEdit,    'String', ceil( mean(dataStruct.venc)*2/3) );
		set(handleStruct.userID,   'Enable','off');
        if dataStruct.externflag~=1
            set(handleStruct.encodingTypePop,   'Enable','on');
        end
        set(handleStruct.filterCheck,       'Enable','on','Value',0);
        set(handleStruct.eddyCurrentCheck,  'Enable','on','Value',0);
        set(handleStruct.antiAliasCheck,    'Enable','on','Value',0);
        set(handleStruct.aviMovieCheck,     'Enable','on','Value',0);
        set(handleStruct.secOrderCorrCheck, 'Enable','off','Value',0);
        set(handleStruct.pcmraToDICOMCheck, 'Enable','on','Value',0);
        set(handleStruct.selMRstructCheck,  'Enable','on','Value',0);
        dataStruct.status = ('Data was loaded!');
        if ~isempty(dataStruct.SiemensFlag)
            dataStruct.status = sprintf('%s\n%s\n',dataStruct.status,...
                'Siemens sequence selected.');
        end
    elseif  (~isempty(dataStruct.fileNamesFlow) && strcmp(dataStruct.loadStatus,'cancelled'))
        %% no changes
    else
        dataStruct.status = sprintf('Loading data was aborted by user or invalid data!');
    end
    

elseif strcmp(entryStr, 'nonuniformityCheck')
    func_list    = get(handleStruct.orderList,'String');
    if dataValue ==1
        enable_status = 'on';
        func_list = [func_list;cellstr('nonuniformity')];        
    else
        enable_status = 'off';        
        x = strcmp('nonuniformity',func_list);
        func_list(x) = [];        
    end
    set(handleStruct.orderList, 'String',func_list,'Value',1);
    set(handleStruct.backgroundSlide,       'Enable', enable_status);
    set(handleStruct.backgroundEdt,         'Enable', enable_status);
    set(handleStruct.previewBackgroundBut,  'Enable', enable_status);
    set(handleStruct.nonuniformityCheck,    'Value',  dataValue);    
    
elseif strcmp(entryStr, 'updateTabGroupConversion') || strcmp(entryStr, 'updateTabGroupPreprocessing')
    hObject = getfield(handleStruct,dataValue);  %#ok<GFLD>  
    tabid = entryStr(7:end);
    handleStruct = tabpanelfcn('tab_group_handler',hObject, handleStruct, tabid);
       
elseif strcmp(entryStr, 'chooseFilter')  
   func_list =[get(handleStruct.orderList,'String');cellstr(dataValue)];    
   set(handleStruct.orderList, 'String',func_list);
       
elseif strcmp(entryStr, 'filterCheck')
    
   if (dataValue == 1)
     set(handleStruct.noiseValueEdt,'Enable', 'on');
     set(handleStruct.chooseFilterList,'Enable', 'on');
     set(handleStruct.noiseValueSlide,'Enable', 'on');
     set(handleStruct.stdevNoiseValueEdt,'Enable', 'on');
     set(handleStruct.stdevNoiseValueSlide,'Enable', 'on');
     set(handleStruct.derivatNoiseValueEdt,'Enable', 'on');
     set(handleStruct.derivatNoiseValueSlide,'Enable', 'on');
     set(handleStruct.radioANDBut,'Enable', 'on');
     set(handleStruct.radioORBut,'Enable', 'on');
     set(handleStruct.staticTissueCheck, 'Enable','on','Value',0);
   else
     set(handleStruct.filterCheck,'Value', 0);
     set(handleStruct.noiseValueEdt,'Enable', 'off');
     set(handleStruct.chooseFilterList,'Enable', 'off');
     set(handleStruct.noiseValueSlide,'Enable', 'off');
     set(handleStruct.stdevNoiseValueEdt,'Enable', 'off');
     set(handleStruct.stdevNoiseValueSlide,'Enable', 'off');
     set(handleStruct.derivatNoiseValueEdt,'Enable', 'off');
     set(handleStruct.derivatNoiseValueSlide,'Enable', 'off');
     set(handleStruct.radioANDBut,'Enable', 'off');
     set(handleStruct.radioORBut,'Enable', 'off');
     set(handleStruct.staticTissueCheck, 'Enable','off','Value',0);
     if get(handleStruct.eddyCurrentCheck, 'Value')==0
        dataStruct.statRegStatus = 'off';
        set(handleStruct.staticValueEdt,'Enable', 'off');     
        set(handleStruct.staticValueSlide,'Enable', 'off');
        set(handleStruct.previewStatRegBut,'Enable', 'off','Value',0,'String','preview of static region');  
     end
     filter_list  = get(handleStruct.chooseFilterList,'String');
     func_list    = get(handleStruct.orderList,'String');
     for j=1:length(filter_list)
        x = strmatch(filter_list{j},func_list);
        func_list(x) = [];
     end      
     set(handleStruct.orderList, 'String',func_list,'Value',1);
   end

elseif strcmp(entryStr, 'showStaticTissue')
    if (dataValue == 1)
     set(handleStruct.staticValueEdt,'Enable', 'on');     
     set(handleStruct.staticValueSlide,'Enable', 'on');
     set(handleStruct.previewStatRegBut,'Enable', 'on');
    else
     dataStruct.statRegStatus = 'off';
     set(handleStruct.staticValueEdt,'Enable', 'off');     
     set(handleStruct.staticValueSlide,'Enable', 'off');
     set(handleStruct.previewStatRegBut,'Enable', 'off','Value',0,'String','preview of static region');     
    end
elseif strcmp(entryStr, 'eddyCurrentCheck')
    if dataValue == 1
        % check for 2nd order correction and disable it if necessary
        if get(handleStruct.secOrderCorrCheck, 'Value') == 1
            set(handleStruct.secOrderCorrCheck, 'Value', 0)
            func_list    = get(handleStruct.orderList,'String');     
            x = strmatch('2nd-order-corr',func_list);
            func_list(x) = [];         
            set(handleStruct.orderList, 'String',func_list,'Value',1);
        end
            
        func_list = [get(handleStruct.orderList,'String');cellstr('eddy-current')];    
        set(handleStruct.orderList, 'String',func_list);        
    else
       func_list    = get(handleStruct.orderList,'String');     
       x = strmatch('eddy-current',func_list);
       func_list(x) = [];         
       set(handleStruct.orderList, 'String',func_list,'Value',1);
       set(handleStruct.eddyCurrentCheck, 'Value',0);
    end
elseif strcmp(entryStr, 'moveEntry')
    func_id   = get(handleStruct.orderList,'Value');
    func_list = get(handleStruct.orderList,'String');    
    if strcmp(dataValue,'U')% move up
        func_id_new = func_id-1;
        if func_id_new<1
            func_id_new=1;
        end
    elseif strcmp(dataValue,'D')% move down
        func_id_new = func_id+1;
        if func_id_new> size(func_list,1)
            func_id_new = size(func_list,1);
        end
    end
    f1 = func_list{func_id};
    f2 = func_list{func_id_new};
    % if function is manual anti-aliasing
    if ~(strcmp(f1,'manual anti-aliasing')&& strcmp(f2,'anti-aliasing'))...
            && ~(strcmp(f2,'manual anti-aliasing')&& strcmp(f1,'anti-aliasing'))
        func_list{func_id_new}= f1;
        func_list{func_id} = f2;
        set(handleStruct.orderList,'String',func_list);
        set(handleStruct.orderList,'Value',func_id_new);
    end
    
elseif strcmp(entryStr, 'deleteOrderListEntry')
    id_value = get(handleStruct.orderList,'Value');
    func_list = get(handleStruct.orderList,'String');
    if ~isempty(func_list)
     %RPM 
     %find data to be deleted
     match = func_list(id_value(:)); 
     noisefilt={'noise-filter','stdev-flter','derivative-filter','median-filter'};
     for k = 1:numel(match); 
         %reverse index to make last first
         n = (size(id_value,2)+1)-k;
            if strcmp(match{n},'anti-aliasing')        
                % delete manual
                x = strmatch('manual anti-aliasing', match{n});
                set(handleStruct.antiAliasCheck,    'Enable','on','Value',0); 
            elseif  strcmp(match{n},'eddy-current');
                 x = strmatch('eddy-current',match{n});
                 set(handleStruct.eddyCurrentCheck,  'Enable','on','Value',0);
            elseif  strcmp(match{n},'noise-filter')...
                  || strcmp(match{n},'stdev-filter')...
                  || strcmp(match{n},'derivative-filter')...
                  || strcmp(match{n},'median-filter')
                  x = strmatch('noise-filter',match{n}); 
                  set(handleStruct.filterCheck,'Enable','on','Value',0);
            end       
        set(handleStruct.orderList,'String',func_list,'Value',1);
     end
            func_list(id_value(:)) = [] ;
            %give updated func_list to handleStruct 
            set(handleStruct.orderList,'String',func_list(:))

            %Look for remaining noise filtering options
            for k = 1:numel(noisefilt);
               if all(strcmp(func_list(:),noisefilt(1))) == 0
                   set(handleStruct.filterCheck,'Enable','on','Value',1);
               end
            end 

    end
    

elseif strcmp(entryStr, 'secOrderCorrCheck')
    if dataValue == 1
        % check for eddy current correction and disable it if necessary
        if get(handleStruct.eddyCurrentCheck, 'Value') == 1
            set(handleStruct.eddyCurrentCheck, 'Value', 0)
            func_list    = get(handleStruct.orderList,'String');     
            x = strmatch('eddy-current',func_list);
            func_list(x) = [];         
            set(handleStruct.orderList, 'String',func_list,'Value',1);
        end
        
        func_list = [get(handleStruct.orderList,'String');cellstr('2nd-order-corr')];    
        set(handleStruct.orderList, 'String',func_list);        
    else
       func_list    = get(handleStruct.orderList,'String');     
       x = strmatch('2nd-order-corr',func_list);
       func_list(x) = [];         
       set(handleStruct.orderList, 'String',func_list,'Value',1);
       set(handleStruct.secOrderCorrCheck, 'Value',0);
    end
elseif strcmp(entryStr, 'pcmraCheck')
    if dataValue ==1
        enable_status = 'on';
        set(handleStruct.staticTissueCheck,'Value', 0);
    else
        enable_status = 'off';
    end
     set(handleStruct.previewPCmraBut,'Enable', enable_status);
     set(handleStruct.staticValueSlide,'Enable', enable_status);
     set(handleStruct.staticValueEdt,'Enable', enable_status);
     set(handleStruct.previewStatRegBut,'Enable', enable_status);
     set(handleStruct.staticTissueCheck,'Enable', enable_status );
     set(handleStruct.pcmraMeanAbsVelCheck,'Enable', enable_status);  
     set(handleStruct.pcmraSumSquaresCheck,'Enable', enable_status);
     set(handleStruct.pcmraPseudoComplDiffCheck,'Enable', enable_status);
	 set(handleStruct.pcmraTimeAveMagCheck,'Enable', enable_status);
     set(handleStruct.pcmraSqrtSumSquaresCheck, 'Enable', enable_status); %% LiliMa
     set(handleStruct.pcmraVelCorrCheck, 'Enable', enable_status); %% MBS
     set(handleStruct.vencThresholdEdit,'Enable', enable_status);  
	 set(handleStruct.pcmraMedianFilterCheck,'Enable', enable_status); 
     %set(handleStruct.pcmraHistEqualCheck,'Enable', enable_status);  
	 set(handleStruct.timeFramesMinEdit,'Enable', enable_status);  
	 set(handleStruct.timeFramesMaxEdit,'Enable', enable_status);  
	
elseif strcmp(entryStr, 'pdCheck')
    set(handleStruct.previewPDMaskBut,'Enable',dataValue);
    set(handleStruct.thresholdPDmaskEdit,'Enable',dataValue);
    
elseif strcmp(entryStr, 'antiAliasCheck')
    if dataValue ==1
        enable_status = 'on';
        func_list = [get(handleStruct.orderList,'String');cellstr('anti-aliasing')];    
        set(handleStruct.orderList, 'String',func_list);
    else
        enable_status = 'off';
        func_list    = get(handleStruct.orderList,'String');     
         x = strmatch('anti-aliasing',func_list);
         func_list(x) = [];
         x = strmatch('manual anti-aliasing',func_list);
         func_list(x) = [];
         set(handleStruct.orderList, 'String',func_list,'Value',1);
    end
    set(handleStruct.unwrapManualBut, 'Enable',enable_status);
    set(handleStruct.numberOfIterPop, 'Enable',enable_status)
    set(handleStruct.loadUnwrapDataBut, 'Enable',enable_status);
    set(handleStruct.antiAliasCheck,'Value',dataValue) 
    
% load existing unwrap data
elseif strcmp(entryStr,'loadUnwrapData')
    orderList = get(handleStruct.orderList,'String');
    strFlag = any(strcmp(orderList,'manual anti-aliasing'));
    if ~strFlag
        func_list = cellstr([orderList;'manual anti-aliasing']); 
        set(handleStruct.orderList, 'String',func_list);
    end    
elseif strcmp(entryStr, 'applyEnable')
    if dataValue == 1
        enable_status ='on';
    else
        enable_status = 'off';
    end    
    set(handleStruct.applyBut,'Enable', enable_status);
    set(handleStruct.takeQueueBut,'Enable', enable_status);
    
elseif strcmp(entryStr, 'enSightCheck')
    if dataValue == 1
        enable_status ='on';
    else
        enable_status = 'off';
    end
    set(handleStruct.appAnatomDataBut,'Enable', enable_status);
    set(handleStruct.deleteAnatomDataBut,'Enable', enable_status);
    set(handleStruct.AnatomDataList,'Enable', enable_status);
    
elseif strcmp(entryStr,'selMRstruct')
    if dataValue == 1
        enable_status ='on';
    else
        enable_status = 'off';
    end    
    
%% load or exchange images     
elseif strcmp(entryStr, 'image')
   %% initialize
   venc         = dataStruct.venc;
   nFE          = dataStruct.nFE;
   signVi       = dataStruct.signVijk(1);
   signVj       = dataStruct.signVijk(2);
   signVk       = dataStruct.signVijk(3);
   actSliceNum  = ceil(get(handleStruct.SlicesNumSlide,'value'));
   actPhaseNum  = ceil(get(handleStruct.PhasesNumSlide,'value'));
   %% end of: initialize
   
    
   %% load original data and display according to choosen view //MM
   imageMxMAG    = flipud(dataStruct.dataMag(:,:,actSliceNum,actPhaseNum));
   imageMxFlowVk = flipud(dataStruct.dataFlow(:,:,actSliceNum,actPhaseNum,nFE));
   
   if signVk == -1
       imageMxFlowVk = signVk*(imageMxFlowVk)+4095;
   end
   viewID = get(handleStruct.changeViewPop, 'Value');   
   if nFE ==1
       imageMxFlowUp    = [];   
       imageMxFlowDown  = [];
       imageMxFlowBig   = imageMxFlowVk;
       VorigBig    = 'Vk original';
       VorigUp     = '';
       VorigDown   = '';
       viewID = 3;     
   elseif nFE==3
       imageMxFlowVi = flipud(dataStruct.dataFlow(:,:,actSliceNum,actPhaseNum,1));
       if signVi == -1
           imageMxFlowVi = signVi*(imageMxFlowVi)+4095;
       end
       imageMxFlowVj = flipud(dataStruct.dataFlow(:,:,actSliceNum,actPhaseNum,2));
       if signVj == -1
           imageMxFlowVj = signVj*(imageMxFlowVj)+4095;
       end
            
       %% change data according to the choosen view       
       if viewID==1 %% big view is Vi
           imageMxFlowUp    = imageMxFlowVj;   
           imageMxFlowDown  = imageMxFlowVk;
           imageMxFlowBig   = imageMxFlowVi;
           VorigBig    ='Vi original';
           VorigUp     ='Vj original';
           VorigDown   ='Vk original'; 
              
       elseif viewID==2 %% big view is Vj
           imageMxFlowUp    = imageMxFlowVk;   
           imageMxFlowDown  = imageMxFlowVi;
           imageMxFlowBig   = imageMxFlowVj;
           VorigBig    ='Vj original';
           VorigUp     ='Vk original';
           VorigDown   ='Vi original';
           
       elseif viewID==3 %% big view is Vk
           imageMxFlowUp    = imageMxFlowVi;   
           imageMxFlowDown  = imageMxFlowVj;
           imageMxFlowBig   = imageMxFlowVk;
           VorigBig     ='Vk original';
           VorigUp      ='Vi original';
           VorigDown    ='Vj original';
       end 
       %% end of: change data according to the choosen view     
   end
   
   set(handleStruct.VorigBigTxt,'String',VorigBig);
   set(handleStruct.VorigUpTxt,'String',VorigUp);
   set(handleStruct.VorigDownTxt,'String',VorigDown);
   
   VmodUpLabel      = regexprep(VorigUp,'original','modified');
   VmodDownLabel    = regexprep(VorigDown,'original','modified');
   VmodBigLabel     = regexprep(VorigBig,'original','modified');   
   
   set(handleStruct.imageMag, 'CData',imageMxMAG(:,:,1));
   set(handleStruct.imageFlowBig, 'CData',imageMxFlowBig(:,:,1));
   set(handleStruct.imageFlowUp, 'CData',imageMxFlowUp(:,:,1));   
   set(handleStruct.imageFlowDown, 'CData',imageMxFlowDown(:,:,1));    
   %% end of: load original data and display according to choosen view 
   
   %% load and display modified data
   if isempty(dataValue)
      %% if not any preprocessing function has been selected
      %% flow images are just copied 
      set(handleStruct.axesFlowUpMod, 'CLimMode','auto');
      set(handleStruct.axesFlowDownMod, 'CLimMode','auto');
      set(handleStruct.axesFlowBigMod, 'CLimMode','auto');
      
      set(handleStruct.imageFlowUpMod, 'CData',imageMxFlowUp(:,:,1)); 
      set(handleStruct.imageFlowDownMod, 'CData',imageMxFlowDown(:,:,1)); 
      set(handleStruct.imageFlowBigMod, 'CData',imageMxFlowBig(:,:,1));
      
      %% settings for colorbar
      yTick = zeros(1,round((max(venc)/25)));
      yTickLabel = zeros(1,round((max(venc)/25)));
      for k=0:round((max(venc)/25))
           yTick(k+1)       = dataStruct.phaseRange*k/round((max(venc)/25));
           yTickLabel(k+1)  = -max(venc)+k*50;
      end
                               
      v        = ver('matlab');              % get matlab version
      isnewver = str2double(v.Version)>=8.4; % test if pesky version of matlab
      if isnewver % treat colorbar as its own entity
          delete(get(handleStruct.colorbarAx, 'children'));
          set(handleStruct.colorbarAx, 'Tag', '', 'UserData', []);
          % colorbar(handleStruct.colorbarAx, 'peer', handleStruct.axesFlowBigMod);
          set(handleStruct.colorbarAx,'YLimMode','manual','YLim',[0 4096],'YTickMode','manual','YTick',yTick',...
              'YTickLabelMode','manual','YTickLabel',yTickLabel');
      end      %% end of: settings for colorbar
      
   else
       % if background preview is selected
       if dataStruct.previewBackground ==1
           % show mask overlay
           mask = flipud(dataStruct.backgroundMask);
           mask = cat(3,mask,mask,zeros(size(mask),'single'));
           % scale magnitude data
           magdata = single(imageMxMAG/max(imageMxMAG(:)));
           maskedData = single(mask) + repmat(magdata,[1,1,3]);
           maskedData= maskedData/max(maskedData(:));
           set(handleStruct.imageMag, 'CData',maskedData);
       else
           set(handleStruct.imageMag, 'CData',imageMxMAG(:,:,1));
       end
       
       %% if any preprocessing function has been selected
       %% modified flow images are loaded from dataStruct           
           
       %% calculate min and max values to scale the images  
       mV = [min(dataStruct.imaFLOW(:)), max(dataStruct.imaFLOW(:))];      
   
                 
       %% load static regions data instead of Vx-flow-data    
       if strcmp(dataStruct.statRegStatus,'on')
           
           imageMxFlowVkMod = flipud(squeeze(dataStruct.imaFLOW(:,:,actPhaseNum,nFE))); 
           if nFE==3               
               imageMxFlowVjMod = flipud(squeeze(dataStruct.imaFLOW(:,:,actPhaseNum,2)));
               imageMxFlowViMod = flipud(squeeze(dataStruct.imaFLOW(:,:,actPhaseNum,1))); 
               Vmag = sqrt((imageMxFlowViMod/100).^2+(imageMxFlowVjMod/100).^2+(imageMxFlowVkMod/100).^2);           
           elseif nFE==1
               imageMxFlowVjMod =[]; 
               Vmag = sqrt( (imageMxFlowVkMod/100).^2 );  
           end

           Vmag = imageMxMAG.*Vmag;
           Vmag = Vmag/max(Vmag(:));
           Vmag = repmat(Vmag,[1 1 3]);
           mask = flipud(dataStruct.statMask);
           mask = cat(3,mask,mask,zeros(size(mask),'single'));
           imageMxFlowViMod = Vmag + mask;
           imageMxFlowViMod = imageMxFlowViMod/max(imageMxFlowViMod(:));           
       %% load flow images without preview of static regions and pc-mra        
       else           
           if nFE==3
              imageMxFlowViMod = flipud(squeeze(dataStruct.imaFLOW(:,:,actPhaseNum,1)));              
              imageMxFlowVjMod = flipud(squeeze(dataStruct.imaFLOW(:,:,actPhaseNum,2)));              
           elseif nFE==1
              imageMxFlowViMod =[];   
              imageMxFlowVjMod =[]; 
           end
           imageMxFlowVkMod =flipud(squeeze(dataStruct.imaFLOW(:,:,actPhaseNum,nFE)));
       end
       
       %% display images according to the choosen view       
       if viewID==1 %% big view is Vx
           if strcmp(dataStruct.statRegStatus,'on')
              set(handleStruct.axesFlowBigMod, 'CLimMode','auto');
              VmodBigLabel = 'preview of static tissue';              
           else               
               set(handleStruct.axesFlowBigMod, 'CLim',[mV(1) mV(2)]);
               set(handleStruct.axesFlowUpMod, 'CLim',[mV(1) mV(2)]);
               set(handleStruct.axesFlowDownMod, 'CLim',[mV(1) mV(2)]);
           end
           set(handleStruct.imageFlowUpMod, 'CData',imageMxFlowVjMod);
           set(handleStruct.imageFlowDownMod, 'CData',imageMxFlowVkMod); 
           set(handleStruct.imageFlowBigMod, 'CData',imageMxFlowViMod);               
       
       elseif viewID==2 %% big view is Vy
           if strcmp(dataStruct.statRegStatus,'on')
               
               if strcmp(dataStruct.statRegStatus,'on')
                  set(handleStruct.axesFlowDownMod, 'CLimMode','auto');
                  VmodDownLabel = 'preview of static tissue';
               else
                  VmodDownLabel = 'Vi modified'; 
               end
           else
               set(handleStruct.axesFlowBigMod, 'CLim',[mV(1) mV(2)]);
               set(handleStruct.axesFlowUpMod, 'CLim',[mV(1) mV(2)]);
               set(handleStruct.axesFlowDownMod, 'CLim',[mV(1) mV(2)]);
           end
           set(handleStruct.imageFlowUpMod, 'CData',imageMxFlowVkMod);
           set(handleStruct.imageFlowDownMod, 'CData',imageMxFlowViMod); 
           set(handleStruct.imageFlowBigMod, 'CData',imageMxFlowVjMod);
                   
       elseif viewID==3 %% big view is Vz 
           if strcmp(dataStruct.statRegStatus,'on')
               VmodDownLabel = 'Vj modified';
               
               if strcmp(dataStruct.statRegStatus,'on')
                  set(handleStruct.axesFlowUpMod, 'CLimMode','auto');
                  VmodUpLabel = 'preview of static tissue';
               else
                  VmodUpLabel = 'Vi modified';
               end
           else
               set(handleStruct.axesFlowBigMod, 'CLim',[mV(1) mV(2)]);
               set(handleStruct.axesFlowUpMod, 'CLim',[mV(1) mV(2)]);
               set(handleStruct.axesFlowDownMod, 'CLim',[mV(1) mV(2)]);
           end
           set(handleStruct.imageFlowUpMod, 'CData',imageMxFlowViMod);
           set(handleStruct.imageFlowDownMod, 'CData',imageMxFlowVjMod); 
           set(handleStruct.imageFlowBigMod, 'CData',imageMxFlowVkMod);               
       end 
       %% end of: display images according to the choosen view   
      
       %% settings for colorbar 
%        delete(get(handleStruct.colorbarAx, 'children'));
%        set(handleStruct.colorbarAx, 'Tag', '', 'UserData', []);
%        colorbar(handleStruct.colorbarAx, 'peer', handleStruct.axesFlowBigMod);      
   end
   
   %% set labels for modified images
   set(handleStruct.VmodBigTxt,'String',VmodBigLabel);
   set(handleStruct.VmodUpTxt,'String',VmodUpLabel);
   set(handleStruct.VmodDownTxt,'String',VmodDownLabel);
   %% end of: set labels for modified images

elseif strcmp(entryStr,'setVers')
   handleStruct = local_set_autosegvers(handleStruct,dataStruct.serverPath);
   set(handleStruct.autosegvers,'Value',1);
   
% % % % % % % % % % % % % % % % %     
elseif strcmp(entryStr, 'dummy')
%%% detect invalid entry string
else
    msg = sprintf('Sorry, entryswitch ''%s'' was not recognized in ''local_set_dataStruct_entry'', no action performed', entryStr);
    warndlg(lasterr,msg);
end

set(handleStruct.statusEdt, 'String', dataStruct.status);
drawnow

%%% End of: detect invalid entry string

%%%%% End of: work on dataStruct value and update gui
%%%%% 

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  more local functions  (utilities)                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%---> all functions above are strictly NECESSARY for the proposed approach to build
%---> tools. All functions below do not have to be placed in this template_internal
%---> file but may be assigned to files of their own. In both cases make sure
%---> that results from these functions in external files or local function within
%---> this file write back their processing results in dataStruct.
%---> Do not use capital lettes for function names (due to portability reasons).



%% function to choose and sort the image files 
function dataStruct = local_sort_filenames(dataStruct, modstr)
    %% set default values
    fileNamesMagMx  = [];
    fileNamesFlowMx = [];
    dataStruct.PhilipsFlag=[]; % this hasn't been used yet - but may be needed for trigger time
    dataStruct.SiemensFlag=[];
    dataStruct.TimeStamps = [];
    %% end of : set default values
    
    %% read data according to the modus
    if strcmp(modstr, 'queue')
        dirStrMag  = fullfile(dataStruct.dataDirectory,'mag');
        extensionStr = 'dcm';
        [fileNamesMagMx, dirStrMag]  = local_get_filelist(dirStrMag,extensionStr);
    else
        %choose magnitude files
        %[fileNamesMagMx, dirStrMag]  = file_chooser_tool;
        dirStrMag = uigetdir(pwd,'Select ''mag'' Data Directory');
        extensionStr = '';
        if dirStrMag
            [fileNamesMagMx, dirStrMag]  = local_get_filelist(dirStrMag,extensionStr);
            filePathStrMAG  = full_filename(dirStrMag,fileNamesMagMx(1,:)); %LiliMa: if sequence name has fl3, then it is a siemens sequence.
            dirStrGen=fileparts(dirStrMag);
            files_gen=dir(char(dirStrGen));
            dirflags_flow=[files_gen.isdir];
            subfolders=files_gen(dirflags_flow);
            pathstr = fileparts(dirStrMag);
            fileCount = 1;
            fileListStruct=[];
            for k=1:length(subfolders)
                fileStr = char(subfolders(k).name);
                num=regexp(fileStr,'\d','ONCE');
                words=regexp(fileStr,'[a-z, A-Z]','ONCE');
                if subfolders(k).isdir==1 && strcmp(fileStr,'..')==0 &&  strcmp(fileStr,'.')==0 && strcmp(fileStr,'mag')==0&&~isempty(num)&&isempty(words)
                    fileListStruct(fileCount,:) = subfolders(k).name;
                    fileCount=fileCount+1;
                end
            end
            if length(fileListStruct)>2
                dataStruct.SiemensFlag=1;
                dataStruct  = local_get_scan_info(dataStruct,filePathStrMAG); %LiliMa: this implementation checks if there are enough folders, then double checks
                                                                              %in the header if the sequence name contains fl3d; Oct 27, 2017: commented out second check
            end
        end
   end
    
    %% continue if at least one file was chosen

    if ~isempty(fileNamesMagMx)&& isempty(dataStruct.SiemensFlag)
        if size(fileNamesMagMx,1) == 1 % Enhanced DICOM
            patterns = ["FH","AP","RL"];
            FEdirlist = files_gen(contains({files_gen.name},patterns));
            nFE = size(FEdirlist,1);
            fileNamesFlowMx = fileNamesMagMx;
            dataStruct.nFE = nFE;
            dataStruct.dataDirectory    = pathstr;            
            % repeat to get trigger times (reversed to get numOfFliesMag at the last iteration)
          	waitH = waitbar(0,'Reading Enhanced DICOM tag ... ');
            for idir = nFE:-1:1
                count = 1-(idir-1)/nFE;
                waitbar(count,waitH,'Reading Enhanced DICOM tag ... ');
                dirnamei = FEdirlist(contains({FEdirlist.name},patterns(idir)));
                filePathStrFLOW  = fullfile(dirnamei.folder,dirnamei.name,fileNamesFlowMx); %LiliMa: replaced this with fileNamesFlow because sometimes bSSFP is used which has diff TR
                %% get information from dicom header
                [dataStruct,dataInfoStruct,errStr]  = local_get_dcm_header_info(dataStruct,filePathStrFLOW,nFE, 0, modstr);
            end
            close(waitH);
            
            if strcmp(modstr,'queue')
                readKeys    = {'encoding sensitivity','','encodingType','','velocity'};
                dataStruct.encodingType = (char(inifile(dataStruct.inifilePath,'read',readKeys)));
                  
                readKeys    = {'encoding sensitivity','','throughplane','','150'};
                dataStruct.venc(3) =str2double(char(inifile(dataStruct.inifilePath,'read',readKeys)));
                readKeys    = {'encoding sensitivity','','inplaneI','','150'};
                dataStruct.venc(1) =str2double(char(inifile(dataStruct.inifilePath,'read',readKeys)));
                readKeys    = {'encoding sensitivity','','inplaneJ','','150'};
                dataStruct.venc(2) =str2double(char(inifile(dataStruct.inifilePath,'read',readKeys)));                   
            end
            
            if isempty (errStr)
                %% get scan information like number of slices/ phases
                dataStruct  =local_get_scan_info(dataStruct,filePathStrFLOW,dataInfoStruct);
       
                %% sort mag und flow files according to slice and phase
                %% number
                dataStruct.fileNamesMag  = fileNamesMagMx;
                dataStruct.fileNamesFlow = fileNamesFlowMx;
                    
                %% end of: write information to dataStruct
                  
                dataStruct = local_load_allData_edicom_philips(dataStruct,FEdirlist);                 
                %% set loading status
                dataStruct.loadStatus = 'ok';
            else
                errordlg(errStr);
            end

        else
            %% automatical choice of flow files in directory "flow"
            %[pathstr, name, ext, vers]= fileparts(dirStrMag);
        
            dirStrFlow   = fullfile(pathstr,'flow');
            [fileNamesFlowMx, dirStrFlow]  = local_get_filelist(dirStrFlow,extensionStr);
            %% end of: automatic choice of flow files in directory "flow"

            %% continue if there are files with flow information
            if ~isempty(fileNamesFlowMx)
                %% calculate number of flow to encode
                numOfFilesMag  = size(fileNamesMagMx,1);
                numOfFilesFlow = size(fileNamesFlowMx,1);
                nFE = numOfFilesFlow / numOfFilesMag;
                %% end of : calculate number of flow to encode
            
                if nFE==1 || nFE==2 || nFE==3
                    dataStruct.nFE  = nFE;
                    dataStruct.dataDirectory    = pathstr;
                    filePathStrFLOW  = full_filename(dirStrFlow,fileNamesFlowMx(1,:)); %LiliMa: replaced this with fileNamesFlow because sometimes bSSFP is used which has diff TR
                
                    %% get information from dicom header
                    [dataStruct,~,errStr]  = local_get_dcm_header_info(dataStruct,filePathStrFLOW,nFE,numOfFilesMag, modstr);
                    if strcmp(modstr,'queue')
                        readKeys    = {'encoding sensitivity','','encodingType','','velocity'};
                        dataStruct.encodingType = (char(inifile(dataStruct.inifilePath,'read',readKeys)));
                    
                        readKeys    = {'encoding sensitivity','','throughplane','','150'};
                        dataStruct.venc(3) =str2double(char(inifile(dataStruct.inifilePath,'read',readKeys)));
                        readKeys    = {'encoding sensitivity','','inplaneI','','150'};
                        dataStruct.venc(1) =str2double(char(inifile(dataStruct.inifilePath,'read',readKeys)));
                        readKeys    = {'encoding sensitivity','','inplaneJ','','150'};
                        dataStruct.venc(2) =str2double(char(inifile(dataStruct.inifilePath,'read',readKeys)));

                    end                    
                
                    if isempty (errStr)
                        %% get scan information like number of slices/ phases
                        dataStruct  =local_get_scan_info(dataStruct,filePathStrFLOW);
       
                        %% sort mag und flow files according to slice and phase
                        %% number
                        dataStruct.fileNamesMag  = reshape(fileNamesMagMx,dataStruct.numSlices,dataStruct.numPhases,size(fileNamesMagMx,2));
                        dataStruct.fileNamesFlow = reshape(fileNamesFlowMx,dataStruct.numSlices,dataStruct.numPhases,nFE,size(fileNamesFlowMx,2));
                    
                        % if phase encoding direction is 'row' --> permute vx
                        % and vy
                        if (strcmp(dataStruct.peDir,'i')&& nFE >1)
                            helpMX1 = squeeze(dataStruct.fileNamesFlow(:,:,1,:));
                            helpMX2 = squeeze(dataStruct.fileNamesFlow(:,:,2,:));
                            dataStruct.fileNamesFlow(:,:,2,:) = helpMX1;
                            dataStruct.fileNamesFlow(:,:,1,:) = helpMX2;
                        end
                        
                        if (strcmp(dataStruct.orientation,'tra'))% 20190923 TF modifided, but necessary to add other condition, e.g. image orientation
                            helpMX1 = squeeze(dataStruct.fileNamesFlow(:,:,1,:));
                            helpMX2 = squeeze(dataStruct.fileNamesFlow(:,:,3,:));
                            dataStruct.fileNamesFlow(:,:,3,:) = helpMX1;
                            dataStruct.fileNamesFlow(:,:,1,:) = helpMX2;
                        end
                        
                        %% end of: write information to dataStruct
                    
                        % load all data //MM
                        dataStruct = local_load_allData(dataStruct);
                    
                        %% set loading status
                        dataStruct.loadStatus = 'ok';
                    else
                        errordlg(errStr);
                    end
                    %% catch and display error
                else
                    errStr = sprintf('%s\n%s','Number of flows to encode is out of range.',...
                                    'Please check the number of images in the chosen directory!');
                    errordlg(errStr);
                end
            else
                errStr = sprintf('%s\n%s%s%s','The required directory ''flow'' does not exist or is empty',...
                                'Please check the directory  ', pathstr);
                errordlg(errStr);
                dataStruct.loadStatus = 'cancelled';
            end
        end
    elseif ~isempty(fileNamesMagMx) && ~isempty(dataStruct.SiemensFlag)

        folderNamesMx=[];

        if length(fileListStruct)>4
            errStr=sprintf('More than 4 folders have numerical values associated. For Siemens sequences, please keep flow folders only named with numeric values.');
            errordlg(errStr);
        end
        
        if any(size(fileListStruct))
            folderNamesMx = sortrows(char(fileListStruct));
        end
        % gather files of first folder and read header
        folder_path1=fullfile(pathstr,folderNamesMx(1,:));
        [fileNamesFlowtemp, ~]  = local_get_filelist(folder_path1,extensionStr);  %LiliMa: sort as x,y,z (I think?)
        numfiles1=size(fileNamesFlowtemp,1);
        fileNamesFlowMx=fileNamesFlowtemp;
        dataStruct.subfoldersFlow=folder_path1;
        filePathStrFlow1=fullfile(dataStruct.subfoldersFlow(1,:),fileNamesFlowMx(1,:));
        
       [~, dataInfoStruct, msg]=dicom_read_singlefile(filePathStrFlow1);
        if isfield(dataInfoStruct,'Private_0029_1020') || isfield(dataInfoStruct.Private_0029_1120)
            tempStr = lower(dataInfoStruct.Private_0029_1020);
            % Siemens phase contrast
            encoding_dir=[];
            posVenc=[];
            try
                encoding_dir(1)  = strfind(tempStr,lower('sAngio.sFlowArray.asElm[0].nDir')); %LiliMa: Siemens sequences separate flow into different folders. This Siemens header field
            catch
                tempStr = lower(dataInfoStruct.Private_0029_1120);
                encoding_dir(1)  = strfind(tempStr,lower('sAngio.sFlowArray.asElm[0].nDir'));
            end
            %is indexed in folder order (numerical order of flow folders. The value in nDir is the phase encoding direction. 1=phase, 2=read, 4=slice.
            encoding_dir(2)  = strfind(tempStr,lower('sAngio.sFlowArray.asElm[1].nDir'));
            encoding_dir(3)  = strfind(tempStr,lower('sAngio.sFlowArray.asElm[2].nDir'));
            posVenc(1)    = strfind(tempStr,lower('sAngio.sFlowArray.asElm[0].nVelocity')); %This header field includes venc for each flow encoding direction
            posVenc(2)    = strfind(tempStr,lower('sAngio.sFlowArray.asElm[1].nVelocity'));
            posVenc(3)    = strfind(tempStr,lower('sAngio.sFlowArray.asElm[2].nVelocity'));
            posVenc(4)=strfind(tempStr, lower('sprescannormalizefilter.ucmode')); %end of section
        end
        
        folder_dir=zeros(1,size(folderNamesMx,1));
        sort_ind=zeros(1,size(folderNamesMx,1));
        for n=1:size(folderNamesMx,1)
            tempStr_venc = tempStr(posVenc(n):encoding_dir(n)-1);
            tempStr_encoding=tempStr(encoding_dir(n):posVenc(n+1)-1);
            % find equal signs in the string
            pos_eq_venc = strfind(tempStr_venc,'=');
            pos_eq_encoding = strfind(tempStr_encoding,'=');
            % find new lines in the string
            [~, idx_venc] = regexp(tempStr_venc, '\n', 'match', 'start');
            [~, idx_encoding] = regexp(tempStr_encoding, '\n', 'match', 'start');
            venc=str2double(tempStr_venc(pos_eq_venc(1)+1:idx_venc(1)-1));
            folder_dir(n) = str2double(tempStr_encoding(pos_eq_encoding(1)+1:idx_encoding(1)-1)); %1=phase, 2=read, 4=slice          
            if folder_dir(n)==4 %slice
            elseif folder_dir(n) == 2
            end
            if folder_dir(n)==2
                sort_ind(1)=n;
                dataStruct.venc(1:2)=venc;
            elseif folder_dir(n)==1
                sort_ind(2)=n;
            elseif folder_dir(n)==4
                sort_ind(3)=n;
                dataStruct.venc(3)=venc;              
            end
        end
       
        for k=1:size(folderNamesMx,1)
            folder_path=fullfile(pathstr,folderNamesMx(sort_ind(k),:));
            switch k
                case 1
                    [fileNamesFlowtemp, ~]  = local_get_filelist(folder_path,extensionStr);  %LiliMa: sort as x,y,z (I think?)
                    numfiles1=size(fileNamesFlowtemp,1);
                    fileNamesFlowMx=fileNamesFlowtemp;
                    dataStruct.subfoldersFlow=folder_path;
                case 2
                    [fileNamesFlowtemp, ~]  = local_get_filelist(folder_path,extensionStr);
                    %numfiles2=size(fileNamesFlowtemp,1);
                    fileNamesFlowMx((numfiles1+1):(numfiles1*2),:)=fileNamesFlowtemp;
                    dataStruct.subfoldersFlow(k,:)=folder_path;
                case 3
                    [fileNamesFlowtemp, ~]  = local_get_filelist(folder_path,extensionStr);
%                     numfiles3=size(fileNamesFlowtemp,1);
                    fileNamesFlowMx((2*numfiles1+1):(numfiles1*3),:)=fileNamesFlowtemp;
                    dataStruct.subfoldersFlow(k,:)=folder_path;
            end
        end
        

        if ~isempty(fileNamesFlowMx)
            %% calculate number of flow to encode
            numOfFilesMag  = size(fileNamesMagMx,1);
            numOfFilesFlow = size(fileNamesFlowMx,1);
            nFE = numOfFilesFlow / numOfFilesMag;
            %% end of : calculate number of flow to encode
            
            if nFE==1 || nFE==2 || nFE==3
                dataStruct.nFE  = nFE;
                dataStruct.dataDirectory    = pathstr;
                filePathStrMAG  = full_filename(dirStrMag,fileNamesMagMx(1,:));
                if nFE==1
                    if length(fileListStruct)>2
                        errStr=sprintf('More than %d folders have numerical values associated. For Siemens sequences, please keep flow folders only named with numeric values.', nFE);
                        errordlg(errStr);
                    end
                elseif nFE == 2
                    if length(fileListStruct)>3
                        errStr=sprintf('More than %d folders have numerical values associated. For Siemens sequences, please keep flow folders only named with numeric values.',nFE);
                        errordlg(errStr);
                    end
                end
                %% get information from dicom header
                [dataStruct,~,errStr]  = local_get_dcm_header_info(dataStruct,filePathStrMAG,nFE,numOfFilesMag, modstr);
                if strcmp(modstr,'queue')
                    readKeys    = {'encoding sensitivity','','encodingType','','velocity'};
                    dataStruct.encodingType = (char(inifile(dataStruct.inifilePath,'read',readKeys)));
                    
                    readKeys    = {'encoding sensitivity','','throughplane','','150'};
                    dataStruct.venc(3) =str2double(char(inifile(dataStruct.inifilePath,'read',readKeys)));
                    readKeys    = {'encoding sensitivity','','inplaneI','','150'};
                    dataStruct.venc(1) =str2double(char(inifile(dataStruct.inifilePath,'read',readKeys)));
                    readKeys    = {'encoding sensitivity','','inplaneJ','','150'};
                    dataStruct.venc(2) =str2double(char(inifile(dataStruct.inifilePath,'read',readKeys)));
                    
                end
                
%                 if isempty (errStr)
                    %% get scan information like number of slices/ phases
                   % dataStruct  = local_get_scan_info(dataStruct,filePathStrMAG);

                
                
                
                %% sort mag und flow files according to slice and phase
                %% number
                dataStruct.fileNamesMag  = reshape(fileNamesMagMx,dataStruct.numSlices,dataStruct.numPhases,size(fileNamesMagMx,2));
                dataStruct.fileNamesFlow = reshape(fileNamesFlowMx,dataStruct.numSlices,dataStruct.numPhases,nFE,size(fileNamesFlowMx,2));
                % if phase encoding direction is 'row' --> permute vx
                % and vy
                if (strcmp(dataStruct.peDir,'i')&& nFE >1) %LiliMa: not sure about this
                    helpMX1 = squeeze(dataStruct.fileNamesFlow(:,:,1,:));
                    helpMX2 = squeeze(dataStruct.fileNamesFlow(:,:,2,:));
                    dataStruct.fileNamesFlow(:,:,2,:) = helpMX1;
                    dataStruct.fileNamesFlow(:,:,1,:) = helpMX2;
                    temp=dataStruct.subfoldersFlow(1,:);
                    dataStruct.subfoldersFlow(1,:)=dataStruct.subfoldersFlow(2,:);
                    dataStruct.subfoldersFlow(2,:)=temp;
                end
                
                %% end of: write information to dataStruct
                
                % load all data //MM
                dataStruct = local_load_allData(dataStruct);
                
                %% set loading status
                dataStruct.loadStatus = 'ok';
            else
                
                errStr = sprintf('%s\n%s','Number of flows to encode is out of range.',...
                    'Please check the number of images in the chosen directory!');
                errordlg(errStr);
            end
        else
            errStr = sprintf('%s\n%s%s%s','The required flow directories do not exist or is empty',...
                'Please check the directory  ', pathstr);
            errordlg(errStr);
            dataStruct.loadStatus = 'cancelled';
        end
    else
        dataStruct.loadStatus = 'cancelled';
    end
%% scan directory for files
function [fileNamesMx, dirStr] = local_get_filelist(dirStr,extensionStr)

    %%% default settings
    fileListStruct = [];
    fileNamesMx    = [];
    %%% End of: default settings

    %%% get file list
    if isempty(extensionStr)
       askStr = full_filename(dirStr,'*');
    else
       tmpStr = sprintf('*.%s',extensionStr);
       askStr = full_filename(dirStr,tmpStr);
    end
    dirStruct    = dir(askStr);
    if isempty(dirStruct) %LiliMa: some Siemens images are in Ima format
        extensionStr='ima';
        tmpStr = sprintf('*.%s',extensionStr);
        askStr = full_filename(dirStr,tmpStr);
        dirStruct = dir(askStr);
    end
    noOfEntries  = size(dirStruct,1);
    %%% End of: get file list

    %%% create fileListStruct
    fileCount = 1;
    for k=1:noOfEntries
       fileStr = dirStruct(k).name;

       if dirStruct(k).isdir==0 && strcmp(fileStr,'..')==0 &&  strcmp(fileStr,'.')==0
          fileListStruct(fileCount).name = fileStr;
          fileCount = fileCount+1;
       end
    end
    %%% End of: create fileListStruct

    %%% convert fileListStruct to matrix of file names
    if any(size(fileListStruct))
       fileNamesMx = sortrows(char(fileListStruct.name));

    end
    %%% End of: convert fileListStruct to matrix of file names

%% get information from dicom header
function [dataStruct,dataInfoStruct,errStr]  = local_get_dcm_header_info(dataStruct,filePath,nFE,numOfFilesMag,modstr)    

    % scan data for information
     [infoStruct, dataInfoStruct, msg] = dicom_scan_singlefile(filePath,'all','dicom-dict_philips.txt'); % AJB this will affect Phlilps data, because the dicom-dict is a custom one!!
%     dataInfoStruct = dicominfo(filePath); % AJB changed this for Phlilps data with trx problems, because the dicom-dict is a custom one!!
     eDicom = false;
     % Enhanced Dicom headers
     if isfield(dataInfoStruct,'SharedFunctionalGroupsSequence')
         [dataInfoStruct,numOfFilesMag,eDicom] = get_necessary_tags_edicom(dataInfoStruct);
         dataStruct.TimeStamps = vertcat(dataStruct.TimeStamps,dataInfoStruct.TimeStamps);
     end

     errStr='';
          
     if nFE>1
     % get number of phases and slices
        if isfield(dataInfoStruct,'CardiacNumberOfImages')&& nFE==3
            numPhases = dataInfoStruct.CardiacNumberOfImages;
        elseif isfield(dataInfoStruct,'Private_2001_1017')&& nFE==3 % AJB this is for Philips data
            numPhases = dataInfoStruct.Private_2001_1017;           % AJB this is for Philips data
        else
            %errStr    = sprintf('%s\n%s',errStr,'Field ''CardiacNumberOfImages'' does not exist in dicomInfoStruct.');
            numPhases =1;
        end 
     else
        numPhases     = numOfFilesMag; 
     end 
     
     rem_part         = mod(numOfFilesMag,numPhases);
     numSlices        = numOfFilesMag/numPhases;
     
    if rem_part==0
		 %% start TF added to identify the vender to calculate correct TimeStamps
	     %% get vender 
		 if isfield(dataInfoStruct,'Manufacturer') && ~isempty(dataInfoStruct.Manufacturer)
		     vender = strtok(dataInfoStruct.Manufacturer);
         elseif isfield(dataInfoStruct,'LUTLabel')
             vender = dataInfoStruct.LUTLabel;
	     end
		 %% end   TF added
		 
         %% get image orientation
         if isfield(dataInfoStruct,'ImageOrientationPatient')
             % EnSight conversion flags, needed for velocity sign & orientation
             % with respect to anatomical (magnitude) images
             % data is accessed in image (x,y,z) coordinates as loaded by matlab
             % x: coronal & axial L/R,  sagittal A/P
             %    read & PE directions
             % y: coronal & sagittal S/I,  axial A/P
             %    read and PE directions
             % z: through slab, coronal A/P, sagittal R/L, axial S/I 
             %    slice/partition direction (numSlices)
             read_vector  = dataInfoStruct.ImageOrientationPatient(1:3);
             phase_vector = dataInfoStruct.ImageOrientationPatient(4:6);
             slice_vector = cross(read_vector,phase_vector);
             mainOrient   = find(abs(slice_vector) == max(abs(slice_vector)));
             if mainOrient == 1
                orientation = 'sag';  % main slab / slice orientation
             elseif mainOrient == 2
                orientation = 'cor';  % main slab / slice orientation
             else
                orientation = 'tra';  % main slab / slice orientation
             end
         else
             errStr =  sprintf('%s\n%s',errStr,'Field ''ImageOrientationPatient'' does not exist in dicomInfoStruct.');               
         end
         %% end of: get image orientation
         
         %% get image size
         if isfield(dataInfoStruct,'Columns')&& isfield(dataInfoStruct,'Rows')
             imageSize = [dataInfoStruct.Rows,dataInfoStruct.Columns];
         else
             errStr =  sprintf('%s\n%s',errStr,'Field ''Columns & Rows'' does not exist in dicomInfoStruct.');
         end
         %% end of: get image size
         
         %% get R-R interval (nominal interval): SuS 06272017
         if isfield(dataInfoStruct,'NominalInterval')
             nominalInterval = dataInfoStruct.NominalInterval; % in ms
         end
         %% 
         
         %% get phase encoding direction
         if isfield(dataInfoStruct,'PhaseEncodingDirection')
             if strcmp(dataInfoStruct.PhaseEncodingDirection,'ROW')
                peDir = 'i';
             else
                peDir = 'j';
             end
         else
             errStr =  sprintf('%s\n%s',errStr,'Field ''PhaseEncodingDirection'' does not exist in dicomInfoStruct.');               
         end
         %% end of: get phase encoding direction
         % read venc from dicom Header           
           phaseRange  = 4096;
           
           %% get software version
           if isfield(dataInfoStruct,'SoftwareVersions')
                if ~isempty(strfind(dataInfoStruct.SoftwareVersions,'VB13'))||~isempty(strfind(dataInfoStruct.SoftwareVersions,'B15'))...
                        ||~isempty(strfind(dataInfoStruct.SoftwareVersions,'B17'))
                    swVersion = 'VB13';
                elseif ~isempty(strfind(dataInfoStruct.SoftwareVersions,'VA25'))
                    swVersion = 'VA25';
                elseif ~isempty(strfind(dataInfoStruct.SoftwareVersions,'VB12'))
                    swVersion = 'VB12';
                else %% set software version to a default value
                    swVersion = 'VB13';
                end
           else
               swVersion = 'VB13';
           end
           %% end of: get software version
           
           %% set signs for velocity directions
           %LiliMa: These signs are different from the signs used in the
           %outputted data formats (see signVx, Vy, Vz) 
           if (~strcmp(swVersion,'VB12') && ~strcmp(swVersion,'VA25')) && isempty(dataStruct.SiemensFlag) && strcmpi(vender,'siemens')
 %LiliMa: siemens data is not usually flipped I think....
               if strcmp(orientation,'sag')
                   if strcmp(peDir,'i')%% phase encoding direction is row
                       signVijk(1)= -1;
                       signVijk(2)= -1;
                       signVijk(3)= -1;
                   else
                       signVijk(1)= -1;
                       signVijk(2)=  1;
                       signVijk(3)= -1;
                   end
               elseif strcmp(orientation,'cor')%%
                   if strcmp(peDir,'i')%% phase encoding direction is row
                       signVijk(1)= -1;
                       signVijk(2)=  1;
                       signVijk(3)= -1;
                   else
                       signVijk(1)= -1;
                       signVijk(2)= -1;
                       signVijk(3)= -1;
                   end
               elseif strcmp(orientation,'tra') %% transversal or axial direction
                   if strcmp(peDir,'i')%% phase encoding direction is row
                       signVijk(1)= -1;
                       signVijk(2)=  1;
                       signVijk(3)= -1;
                   else
                       signVijk(1)=  1;
                       signVijk(2)=  1;
                       signVijk(3)= -1;
                   end
               end
           elseif ~isempty(dataStruct.SiemensFlag)
               if strcmp(orientation,'sag')
                   if strcmp(peDir,'i')%% phase encoding direction is row
                       signVijk(1)= 1;
                       signVijk(2)= 1;
                       signVijk(3)= -1;
                   else %LiliMa: need to data to figure out for sure what the the signs are for Vijk. I can estimate, but to be safe, I kept it 1,1,1.
                       signVijk(1)=  1;
                       signVijk(2)=  1;
                       signVijk(3)=  1;
                   end
               elseif strcmp(orientation,'cor')%%
                   if strcmp(peDir,'i')%% phase encoding direction is row
                       signVijk(1)= -1;
                       signVijk(2)=  1;
                       signVijk(3)= -1;
                   else %LiliMa: need more data
                       signVijk(1)=  1;
                       signVijk(2)=  1;
                       signVijk(3)=  1;
                   end
               elseif strcmp(orientation,'tra') %% transversal or axial direction
                   %LiliMa: need more data
                   signVijk(1)=  1;
                   signVijk(2)=  1;
                   signVijk(3)=  1;
               end
           elseif strcmpi(vender,'philips') % AJB: test for Philips data, not carefully tested!
               signVijk(1)= 1;
               signVijk(2)= 1;
               signVijk(3)= 1;               
           else
               signVijk(1)= 1;
               signVijk(2)= 1;
               signVijk(3)= 1;
           end
           %% end of: set signs for velocity directions
           
         %% get velocity encoding for all directions
         %if ~isempty(strfind(upper(dataInfoStruct.SequenceName),'fl3d');
         if isfield (dataInfoStruct,'Private_0029_1020')&& (~isempty(dataStruct.SiemensFlag) || strcmpi(vender,'siemens')) %LiliMa: siemens data is taken care of in local_sort_filenames; ask about if this is in all siemens sequences
           % read privat header and find venc(s) according to the software version
           tempStr = lower(dataInfoStruct.Private_0029_1020);
           if strcmp(swVersion,'VB13')               
               posVencInplane    = strfind(tempStr,lower('sWiPMemBlock.adFree[8]'));
               posVencThplane    = strfind(tempStr,lower('sWiPMemBlock.adFree[9]'));
               posVencNextField  = strfind(tempStr,lower('sWiPMemBlock.adFree[10]'));
           elseif (strcmp(swVersion,'VB12')||strcmp(swVersion,'VA25'))               
               posVencInplane    = strfind(tempStr,lower('sWiPMemBlock.adFree[9]'));
               posVencThplane    = strfind(tempStr,lower('sWiPMemBlock.adFree[10]'));
               posVencNextField  = strfind(tempStr,lower('sWiPMemBlock.adFree[11]'));
           else
               posVencInplane    = [];
               posVencThplane    = [];
               posVencNextField  = [];             
           end

           if ~isempty(posVencInplane)&& ~isempty(posVencThplane)&& ~isempty(posVencNextField)
               tempStr = tempStr(posVencInplane:posVencNextField-1);
               % find equal signs in the string
               pos_eq = strfind(tempStr,'=');
               % find new lines in the string
               [~, idx] = regexp(tempStr, '\n', 'match', 'start');
               vencInPlane = str2double(tempStr(pos_eq(1)+1:idx(1)-1));
               vencThPlane = str2double(tempStr(pos_eq(2)+1:idx(2)-1));%      
           
           % read venc from standard header fields in case PC data from product sequence is loaded
           elseif (nFE ==1 && ( isempty(posVencInplane) || isempty(posVencThplane) ))
               posVencInplane    = strfind(tempStr,lower('sAngio.sFlowArray.asElm[0].nVelocity'));
               posVencNextField  = strfind(tempStr,lower('sAngio.sFlowArray.asElm[0].nDir'));
               tempStr = tempStr(posVencInplane:posVencNextField-1);
               % find equal signs in the string
               pos_eq = strfind(tempStr,'=');
               % find new lines in the string
               [~, idx] = regexp(tempStr, '\n', 'match', 'start');
               vencThPlane = str2double(tempStr(pos_eq(1)+1:idx(1)-1)) / 100;
               vencInPlane = vencThPlane / 100;  
               msg = sprintf('%s\n%s','Venc succesfully read from dicom header',...
                                      'Assuming Siemens Sequence'); 
               uiwait(warndlg(msg,'','modal')) 
           elseif (nFE ==3 && ( isempty(posVencInplane) || isempty(posVencThplane) ))
               posVencInplane   = strfind(tempStr,lower('sAngio.sFlowArray.asElm[0].nVelocity')); %This header field includes venc for each flow encoding direction
               posVencNextField = strfind(tempStr,lower('sprescannormalizefilter.ucmode')); %end of section
               % find equal signs in the string
               tempStr = tempStr(posVencInplane:posVencNextField-1);
               pos_eq = strfind(tempStr,'=');
               [~, idx] = regexp(tempStr, '\n', 'match', 'start');
               nDir(1) = str2double(tempStr(pos_eq(2)+1:idx(2)-1));
               nDir(2) = str2double(tempStr(pos_eq(4)+1:idx(4)-1));
               nDir(3) = str2double(tempStr(pos_eq(6)+1:idx(6)-1));
               nVel(1) = str2double(tempStr(pos_eq(1)+1:idx(1)-1));
               nVel(2) = str2double(tempStr(pos_eq(3)+1:idx(3)-1));
               nVel(3) = str2double(tempStr(pos_eq(5)+1:idx(5)-1));
               vencThPlane = nVel(nDir==4) / 100;
               vencInPlane = nVel(nDir==1) / 100;  
           else
               
               if ~strcmp(modstr,'queue') || isempty(modstr) 
                   msg = sprintf('%s\n%s','Venc could not be read from dicom header',...
                       'Venc was set to the default value of 150 cm/s');
                   uiwait(warndlg(msg,'','modal'))
        
               end
           end
           % pass venc(s) to array and convert to cm/sec
           venc(nFE)      = vencThPlane * 100;
           if nFE > 1  
              venc(1:2)      = vencInPlane * 100;              
           end
         % AJB Philips switch for VENC (note only set up for isotropic VENC!!!)
         elseif isfield(dataInfoStruct,'Private_2001_101a')  
             vencInPlane = dataInfoStruct.Private_2001_101a; % note sometimes 101a comes in as 101A!
             %vencThPlane = dataInfoStruct.Private_2001_101A(3);
             venc(1:3) = vencInPlane(3); %assume isoVENC	%AJB Not carefully tested!
         % Philips enhancedDicom
         elseif isfield(dataInfoStruct,'PerframeFunctionalGroupsSequence')
             flds_perframe = fieldnames(dataInfoStruct.PerframeFunctionalGroupsSequence);
             if isfield(dataInfoStruct.PerframeFunctionalGroupsSequence.(flds_perframe{numOfFilesMag+1}).Private_2005_140F.Item_1,'Private_2001_101A')
                vencInPlane = dataInfoStruct.PerframeFunctionalGroupsSequence.(flds_perframe{numOfFilesMag+1}).Private_2005_140F.Item_1.Private_2001_101A;
                venc = vencInPlane;           
             else    
                vencInPlane = dataInfoStruct.PerframeFunctionalGroupsSequence.(flds_perframe{numOfFilesMag+1}).Private_2005_140F.Item_1.Private_2005_100B;
                venc(1:3) = vencInPlane;
             end
         elseif isempty(dataStruct.SiemensFlag)
             vencInPlane = 1.5;
             vencThPlane = 1.5;
             % pass venc(s) to array and convert to cm/sec
             venc(nFE)      = vencThPlane * 100;
             if nFE > 1  
              venc(1:2)      = vencInPlane * 100;              
             end
             msg = sprintf('%s\n%s','Venc could not be read from dicom header',...
                                      'Venc was set to the default value of 150 cm/s'); 
             uiwait(warndlg(msg,'','modal'))             
        end
        %% end of: get velocity encoding for all directions

        %% get patient name from header
        %%if isfield(dataInfoStruct.PatientsName,'FamilyName')
        %%     patientName = sprintf('%s%s',dataInfoStruct.PatientsName.FamilyName,'_',dataInfoStruct.StudyDate);
        %%else
        %%     patientName = sprintf('%s%s',dataInfoStruct.PatientsName,'_',dataInfoStruct.StudyDate);
        %%end	
		
		% use deidentified data in velomap_tool
		patientName = sprintf('%s%s%s','subject','_',dataInfoStruct.StudyDate);
			
		
        %% end of: get patient name from header
    else
        errStr = 'Number of slices is not a proper value. Please check the number of images in the chosen directory! ';
    end   
    
    %% write information to dataStruct 
    if isempty(errStr)
        %% delete whitespaces in patientName
        patientName     = regexprep(patientName,' ','');
        dataStruct.peDir            = peDir;
        if isempty(dataStruct.venc)
            dataStruct.venc             = venc;
        elseif any(dataStruct.venc - venc)
            dataStruct.venc(find(venc))       = venc(find(venc));            
        end
        %% start TF added to identify the vender to calculate correct TimeStamps
        dataStruct.Manufacturer     = vender;
        %% end   TF added
		dataStruct.patientName      = patientName;
        dataStruct.orientation      = orientation;
        dataStruct.numSlices        = numSlices;
        dataStruct.numPhases        = numPhases;
        dataStruct.imageSize        = imageSize;
        dataStruct.swVersion        = swVersion;
        dataStruct.signVijk         = signVijk;
        dataStruct.eDicom           = eDicom; 
%        dataStruct.nominalInterval  = nominalInterval; %SuS 06272017
    end
    %% end of: write information to dataStruct 
%%%% end of: get information from dicom header

%% get and display scan information from dicom header  
function dataStruct = local_get_scan_info(varargin)
    % get selected header information from dicom header of selected file
    if nargin < 2
        error('In Function local_get_scan_info: Input argment must be larger than 1');
    end
    
    if nargin > 1
        dataStruct = varargin{1};
        fileNamePath = varargin{2};
    end
    
    if nargin > 2
        infoStruct = varargin{3};
    end 
    
    if nargin < 3
        %infoStruct      = dicominfo(fileNamePath);
        [str1, infoStruct, msg] = dicom_scan_singlefile(fileNamePath,'all','dicom-dict_philips.txt');%,'dicom-dict_mr.mat');
        %infoStruct = dicominfo(fileNamePath); % AJB changed this for Phlilps data trx problems, because the dicom-dict is a custom one!!
        % enhancedDICOM 
        if isfield(infoStruct,'SharedFunctionalGroupsSequence')
             [infoStruct,~,~] = get_necessary_tags_edicom(infoStruct);
        end
    end
    
    outInfo ='';
    if dataStruct.SiemensFlag==1%LiliMa Siemens data
        %  if (~isfield(infoStruct, 'SequenceName') || ~isempty(strfind(infoStruct.SequenceName,'fl3d'))|| ~isempty(strfind(infoStruct.SequenceName,'fl2d')) ) %LiliMa double check that is siemens data
        outInfo=sprintf('%s%s\n', outInfo, 'Siemens WIP');
        %  else
        %             dataStruct.SiemensFlag=[];
        %         end
    end
    
    % get selected header information and store in matrix
    
    if isfield(infoStruct,'PatientsName')
        outInfo = sprintf('%s%s%s%s%s\n',outInfo,'Patient Name : ','subject','_',infoStruct.StudyDate);
    else
        outInfo = sprintf('%s%s%s\n',outInfo,'Patient Name : ','no entry found');
    end
    % get patient's sex
    if isfield(infoStruct,'PatientsSex')
        outInfo = sprintf('%s%s%s\n',outInfo,'Patient Sex : ',infoStruct.PatientsSex);
    elseif isfield(infoStruct,'PatientSex')
        outInfo = sprintf('%s%s%s\n',outInfo,'Patient Sex : ',infoStruct.PatientSex);
    else
        outInfo = sprintf('%s%s%s\n',outInfo,'Patient Sex : ','no entry found');
    end
% end of: get patient's sex
% get patient's date of birth
    if isfield(infoStruct,'PatientsBirthDate')        
        outInfo = sprintf('%s%s%s\n',outInfo,'Patient DOB : ',infoStruct.PatientsBirthDate);
    elseif isfield(infoStruct,'PatientBirthDate')        
        outInfo = sprintf('%s%s%s\n',outInfo,'Patient DOB : ',infoStruct.PatientBirthDate);
    else
        outInfo = sprintf('%s%s%s\n',outInfo,'Patient DOB : ','no entry found');
    end    
% end of: get patient's date of birth 

% get patient's age
    if isfield(infoStruct,'PatientsAge')        
        outInfo = sprintf('%s%s%s\n',outInfo,'Patient Age : ',infoStruct.PatientsAge);
    elseif isfield(infoStruct,'PatientAge')        
        outInfo = sprintf('%s%s%s\n',outInfo,'Patient Age : ',infoStruct.PatientAge);
    else
        outInfo = sprintf('%s%s%s\n',outInfo,'Patient Age : ','no entry found');
    end    
% end of: get patient's age 
% get exam date
    if isfield(infoStruct,'StudyDate')        
       outInfo = sprintf('%s%s%s\n',outInfo,'Exam Date/Time: ',infoStruct.StudyDate);
    else
        outInfo = sprintf('%s%s%s\n',outInfo,'Exam Date/Time: ','no entry found');
    end
% end of: get exam date
    
% get acquisition type
    if isfield(infoStruct,'mracquisitionType')        
        outInfo = sprintf('%s%s%s\n',outInfo,'Volume Coverage : ',infoStruct.mracquisitionType);
    else
        outInfo = sprintf('%s%s%s\n',outInfo,'Volume Coverage : ','no entry found');
    end    
% end of: get acquisition type
    % get image rows
    if isfield(infoStruct,'Width')        
        outInfo = sprintf('%s%s%3.0f\n',outInfo,'Image Rows : ',infoStruct.Width);
    elseif isfield(infoStruct,'Rows')
        outInfo = sprintf('%s%s%3.0f\n',outInfo,'Image Rows : ',infoStruct.Rows);
    else
        outInfo = sprintf('%s%s%s\n',outInfo,'Image Rows : ','no entry found');
    end    
 % end of: image rows

% get image columns
    if isfield(infoStruct,'Height')        
       outInfo = sprintf('%s%s%3.0f\n',outInfo,'Image Columns : ',infoStruct.Height);
    elseif isfield(infoStruct,'Columns')
        outInfo = sprintf('%s%s%3.0f\n',outInfo,'Image Columns : ',infoStruct.Columns);
    else
       outInfo = sprintf('%s%s%s\n',outInfo,'Image Columns : ','no entry found');
    end    
% end of: image columns
    
% get number of phase encoding steps
    if isfield(infoStruct,'NumberOfPhaseEncodingSteps')        
       outInfo = sprintf('%s%s%2.0f\n',outInfo,'Number PE : ',infoStruct.NumberOfPhaseEncodingSteps);
    else
       outInfo = sprintf('%s%s%s\n',outInfo,'Number PE : ','no entry found');
    end    
% end of: get number of phase encoding steps
% get phase encoding direction
    if isfield(infoStruct,'PhaseEncodingDirection')        
        outInfo = sprintf('%s%s%s\n',outInfo,'PE Direction : ',infoStruct.PhaseEncodingDirection);
        PE = infoStruct.PhaseEncodingDirection;
    else
        outInfo = sprintf('%s%s%s\n',outInfo,'PE Direction :','no entry found');
        PE = '';
    end    
% end of: get phase encoding direction    

% get field of view and resolution
    if isfield(infoStruct,'PixelSpacing')
        pix = infoStruct.PixelSpacing;
        resCol = pix(2);
        resRow = pix(1);
        dataStruct.resX = pix(2);
        dataStruct.resY = pix(1);
        %% get FOV row
        if isfield(infoStruct,'Width')
            width =infoStruct.Width;
            FOVrow = round(pix(1) * width);
            outInfo = sprintf('%s%s%3.1f\n',outInfo,'FOV Row : ',FOVrow);
        elseif isfield(infoStruct,'Rows')
            width  = infoStruct.Rows;
            FOVrow = round(pix(1) * width);
            outInfo = sprintf('%s%s%3.1f\n',outInfo,'FOV Row : ',FOVrow);
        else
            outInfo = sprintf('%s%s%s\n',outInfo,'FOV Row : ','no entry found');
        end
        %% get FOV column
        if isfield(infoStruct,'Height')
            height = infoStruct.Height;
            FOVcol = round(pix(2) * height);
            outInfo = sprintf('%s%s%3.1f\n',outInfo,'FOV Column : ',FOVcol);
        elseif isfield(infoStruct,'Columns')
            height = infoStruct.Columns;
            FOVcol = round(pix(2) * height);
            outInfo = sprintf('%s%s%3.1f\n',outInfo,'FOV Column : ',FOVcol);
            
        else
            outInfo = sprintf('%s%s%s\n',outInfo,'FOV Column : ','no entry found');
        end
        
        %% get resolution        
        if (isfield(infoStruct,'PhaseEncodingDirection')||isfield(infoStruct,'PhaseEncodingDirection'))...
            && (isfield(infoStruct,'NumberOfPhaseEncodingSteps'))...
            &&(isfield(infoStruct,'Width')||isfield(infoStruct,'Rows'))&&isfield(infoStruct,'SliceThickness')...
            &&(isfield(infoStruct,'Height')||isfield(infoStruct,'Columns'))
            if(strcmp(PE, 'ROW'))
                resRow = pix(1) * height/infoStruct.NumberOfPhaseEncodingSteps;
            else 
                resCol = pix(2) * width/infoStruct.NumberOfPhaseEncodingSteps;
            end            
            outInfo = sprintf('%s%s%.2f%s%.2f%s%.2f\n',outInfo,'Resolution : ',resRow,' x ',resCol,' x ',infoStruct.SliceThickness);            
            dataStruct.resZ = infoStruct.SliceThickness;
        else
            outInfo = sprintf('%s%s%s\n',outInfo,'Resolution : ','no entry found');
        end 
        %% end of: get resolution
        
    else
       outInfo = sprintf('%s%s%s\n',outInfo,'FOV Row','no entry found');
       outInfo = sprintf('%s%s%s\n',outInfo,'FOV Column','no entry found');
       outInfo = sprintf('%s%s%s\n',outInfo,'Resolution','no entry found');
    end
% end of: get field of view and resolution
    
% get repetion time
    if isfield(infoStruct,'RepetitionTime')        
        outInfo = sprintf('%s%s%.3f\n',outInfo,'TR : ',infoStruct.RepetitionTime);        
    else
       outInfo = sprintf('%s%s%s\n',outInfo,'TR : ','no entry found');       
    end    
% end of: get repetion time

% get echo time
    if isfield(infoStruct,'EchoTime')        
       outInfo = sprintf('%s%s%.3f\n',outInfo,'TE : ',infoStruct.EchoTime');
    else
       outInfo = sprintf('%s%s%s\n',outInfo,'TE : ','no entry found'); 
    end    
% end of: get echo time

if isfield(infoStruct,'CardiacNumberOfImages') %% AJB - for Phlilps this need to be recorded
   outInfo = sprintf('%s%s%2.0f\n',outInfo,'Time Frames : ',infoStruct.CardiacNumberOfImages);
elseif isfield(infoStruct,'Private_2001_1017')                                              %AJB
   outInfo = sprintf('%s%s%2.0f\n',outInfo,'Time Frames : ',infoStruct.Private_2001_1017);  %AJB
end

% get flip angle
    if isfield(infoStruct,'FlipAngle')
        outInfo = sprintf('%s%s%2.0f\n',outInfo,'Flip Angle : ',infoStruct.FlipAngle);
    else
        outInfo = sprintf('%s%s%s\n',outInfo,'Flip Angle : ','no entry found');
    end
% end of: get flip angle

% get pixel band width
    if isfield (infoStruct,'PixelBandwidth')
        outInfo = sprintf('%s%s%3.1f\n',outInfo,'Band Width : ',infoStruct.PixelBandwidth);
    else
        outInfo = sprintf('%s%s%s\n',outInfo,'Band Width : ','no entry found');
    end
% end of: get pixel band width

% get R-R interval (NominalInterval), SuS 06272017
    if isfield (infoStruct,'NominalInterval')
        outInfo = sprintf('%s%s%4.0f\n',outInfo,'R-R Interval (ms): ',infoStruct.NominalInterval);
    else
        outInfo = sprintf('%s%s%s\n',outInfo,'R-R Interval (NominalInterval) : ','no entry found');
    end
% end of: get pixel band width

% get manufacturer's model name
    if isfield(infoStruct,'ManufacturersModelName')
        outInfo = sprintf('%s%s%s\n',outInfo,'MR System : ',infoStruct.ManufacturersModelName);
    elseif isfield(infoStruct,'ManufacturerModelName')
        outInfo = sprintf('%s%s%s\n',outInfo,'MR System : ',infoStruct.ManufacturerModelName);
    else
        outInfo = sprintf('%s%s%s\n',outInfo,'MR System : ','no entry found');
    end
% end of: get manufacturer's model name
    

    dataStruct.dataInfo = outInfo;
    dataStruct.TR       = (infoStruct.RepetitionTime);
%%%% end of: get and display scan information from dicom header 

%% load all magnitude and flow data //MM 
function dataStruct = local_load_allData(dataStruct)

    numPhases   = dataStruct.numPhases;    
	numSlices   = dataStruct.numSlices;    

    dirStrMag   = fullfile(dataStruct.dataDirectory,'mag');
	
	waitH = waitbar(0,'Loading 4D Flow Data ... ');
    %% settings for waitbar       
    

    %% preallocating arrays for speed
	msgOut = sprintf('Allocating Memory Mag .....');
	counter = 0.1;
	waitbar(counter,waitH,msgOut);
    dataStruct.dataMag  = zeros(dataStruct.imageSize(1),dataStruct.imageSize(2),numSlices,numPhases,'single');
	msgOut = sprintf('Allocating Memory Flow .....');
	counter = 0.1;
	waitbar(counter,waitH,msgOut);
    dataStruct.dataFlow = zeros(dataStruct.imageSize(1),dataStruct.imageSize(2),numSlices,numPhases,dataStruct.nFE,'single');
    
    if isempty(dataStruct.SiemensFlag)
        dataStruct.TimeStamps = zeros(dataStruct.nFE+1, numPhases); % AJB, Philips
        dirStrFlow  = fullfile(dataStruct.dataDirectory,'flow');
        for slc = 1:numSlices
            % waitbar update
            msgOut = ['Loading 4D Flow Data: slice #',num2str(slc),' ...'];
            counter = (slc/numSlices);
            waitbar(counter,waitH,msgOut);
            
            %% get file names for current slice
            fileNamesMag   = squeeze(dataStruct.fileNamesMag(slc,:,:));
            fileNamesFlow  = squeeze(dataStruct.fileNamesFlow(slc,:,:,:));
            
            %% loading data
            for m = 1:numPhases
                
                %% get magnitude data for a single phase
                if numPhases == 1
                    pathMag       = full_filename(dirStrMag, fileNamesMag');
                else %if number of time phases is >1
                    pathMag       = full_filename(dirStrMag, squeeze(fileNamesMag(m,:)));
                end
                
                %%dataStruct.dataMag(:,:,slc,m) = dicom_read_singlefile(pathMag,0);
                dataStruct.dataMag(:,:,slc,m) = single(dicomread(pathMag));   %% // MM: use malab dicomread which is much faster
				%% start TF added to correctly estimate TimeStamps at Ensight output.
                if slc == 1 
                    tempTime = dicominfo(pathMag);
                    if isfield(tempTime, 'TriggerTime')
                        dataStruct.TimeStamps(1,m) = double(tempTime.TriggerTime);
                    else
                        disp('error with trigger time:AJB');
                        return
     				end
                end                
				%% end   TF added
                %% get flow data for single phase
                for n = 1:dataStruct.nFE %LiliMa: for siemens data reading use this as well
                    if dataStruct.nFE > 1 && numPhases > 1
                        pathFlow = full_filename(dirStrFlow, squeeze(fileNamesFlow(m,n,:)));
                    elseif dataStruct.nFE > 1 && numPhases
                        pathFlow = full_filename(dirStrFlow, fileNamesFlow(n,:));
                    else
                        pathFlow = full_filename(dirStrFlow, squeeze(fileNamesFlow(m,:)));
                    end
                    
                    %%dataStruct.dataFlow(:,:,slc,m,n) = dicom_read_singlefile(pathFlow,0);
                    dataStruct.dataFlow(:,:,slc,m,n) = dicomread(pathFlow);  %% // MM: use malab dicomread which is much faster

                    % AJB: trigger times for Philips Data (unproven/untested)
                    if slc ==1
                        tempTime = dicominfo(pathFlow);
                        if isfield(tempTime, 'TriggerTime')
                            dataStruct.TimeStamps(n+1, m) = double(tempTime.TriggerTime);
                        else
                            disp('error with trigger time:AJB');
                            return
                        end
                    end

                
                end
                
            end %% end phases
            
        end  %% end slices
    elseif dataStruct.SiemensFlag==1
        dataStruct.TimeStamps = zeros(dataStruct.nFE+1, numPhases);
        for slc = 1:numSlices
            % waitbar update
            msgOut = ['Loading 4D Flow Data: slice #',num2str(slc),' ...'];
            counter = (slc/numSlices);
            waitbar(counter,waitH,msgOut);
            
            %% get file names for current slice
            fileNamesMag   = squeeze(dataStruct.fileNamesMag(slc,:,:));
            fileNamesFlow  = squeeze(dataStruct.fileNamesFlow(slc,:,:,:));
            
            %% loading data
            for m = 1:numPhases
                
                %% get magnitude data for a single phase
                if numPhases == 1
                    pathMag       = full_filename(dirStrMag, fileNamesMag');
                else %if number of time phases is >1
                    pathMag       = full_filename(dirStrMag, squeeze(fileNamesMag(m,:)));
                end
                
                %%dataStruct.dataMag(:,:,slc,m) = dicom_read_singlefile(pathMag,0);
                dataStruct.dataMag(:,:,slc,m) = single(dicomread(pathMag));   %% // MM: use malab dicomread which is much faster
                if slc == 1 %LiliMa: only do it once so data loads faster
                    tempTime = dicominfo(pathMag);
                    if isfield(tempTime, 'TriggerTime')
                        dataStruct.TimeStamps(1,m) = double(tempTime.TriggerTime);
                    else
                        dataStruct.TimeStamps(1,m) = dataStruct.TR/2 + (m-1)*dataStruct.TR; % if can't find trigger time for some reason, use our default way of calculating time stamps
                    end
                end
                %% get flow data for single phase
                for n = 1:dataStruct.nFE %LiliMa: for siemens data reading use this as well
                    
                    pathstr=dataStruct.subfoldersFlow(n,:);
                    
                    if dataStruct.nFE > 1 && numPhases > 1
                        pathFlow = full_filename(pathstr, squeeze(fileNamesFlow(m,n,:)));
                    elseif dataStruct.nFE > 1 && numPhases
                        pathFlow = full_filename(pathstr, fileNamesFlow(n,:));
                    else
                        pathFlow = full_filename(pathstr, squeeze(fileNamesFlow(m,:)));
                    end
                    dataStruct.dataFlow(:,:,slc,m,n) = dicomread(pathFlow);  %% // MM: use malab dicomread which is much faster
                    
                    % LiliMa: trigger times for Siemens Data
                    if slc ==1
                        tempTime = dicominfo(pathFlow);
                        if isfield(tempTime, 'TriggerTime')
                            dataStruct.TimeStamps(n+1, m) = double(tempTime.TriggerTime);
                        else
                            dataStruct.TimeStamps(n+1,m) = dataStruct.TR/2 + (m-1)*dataStruct.TR;
                        end
                    end
                    
                end
                %%dataStruct.dataFlow(:,:,slc,m,n) = dicom_read_singlefile(pathFlow,0);
                
                
            end %% end phases  
            
        end  %% end slices
    end
        close (waitH);

	
%%% end of: load all magnitude and flow data //MM

%% load all magnitude and flow data //MM 
function dataStruct = local_load_allData_edicom_philips(dataStruct,dirlist)

    numPhases   = dataStruct.numPhases;    
	numSlices   = dataStruct.numSlices;    

    dirStrMag   = fullfile(dataStruct.dataDirectory,'mag');
	
	waitH = waitbar(0,'Loading Enhanced DICOM 4D Flow Data ... ');
    %% settings for waitbar       
    

    %% preallocating arrays for speed
    if isempty(dataStruct.SiemensFlag)
        chardirlist = {'FH','AP','RL'};
        for idir = 1:size(dirlist,1)           
            pathFile  = fullfile(dataStruct.dataDirectory,dirlist(contains({dirlist.name},chardirlist{idir})).name,dataStruct.fileNamesFlow);
            % waitbar update
            msgOut = ['Loading Enhanced DICOM 4D Flow Data:',chardirlist{idir},'...'];
            counter = (idir/size(dirlist,1));
            waitbar(counter,waitH,msgOut);
            %%dataStruct.dataMag(:,:,slc,m) = dicom_read_singlefile(pathMag,0);
            eImages = single(dicomread(pathFile));   %% // MM: use malab dicomread which is much faster
            if strcmp(chardirlist{idir},'FH')
                numFH = size(eImages,4);
                dataMag = eImages(:,:,:,1:end/2);
                dataFlowFH = eImages(:,:,:,end/2+1:end); 
                dataStruct.dataMag = permute(reshape(dataMag,[dataStruct.imageSize(1),dataStruct.imageSize(2),numPhases,numSlices]),[1 2 4 3]);
                dataFlowFH = permute(reshape(dataFlowFH,[dataStruct.imageSize(1),dataStruct.imageSize(2),numPhases,numSlices,1]),[1 2 4 3]); 
            elseif strcmp(chardirlist{idir},'AP')
                if size(eImages,4) == numFH
                    dataFlowAP = eImages(:,:,:,end/2+1:end); 
                else
                    dataFlowAP = eImages;                    
                end    
                dataFlowAP = permute(reshape(dataFlowAP,[dataStruct.imageSize(1),dataStruct.imageSize(2),numPhases,numSlices,1]),[1 2 4 3]); 
            elseif strcmp(chardirlist{idir},'RL')
                if size(eImages,4) == numFH
                    dataFlowRL = eImages(:,:,:,end/2+1:end); 
                else
                    dataFlowRL = eImages;                    
                end    
                dataFlowRL = permute(reshape(dataFlowRL,[dataStruct.imageSize(1),dataStruct.imageSize(2),numPhases,numSlices,1]),[1 2 4 3]);
            end
        end  %% end slices
        dataStruct.dataFlow(:,:,:,:,1) = dataFlowAP;
        dataStruct.dataFlow(:,:,:,:,2) = dataFlowFH;
        dataStruct.dataFlow(:,:,:,:,3) = dataFlowRL;
    elseif dataStruct.SiemensFlag==1
        error('Siemens eDicom data cannot be processed in the current version.');
    end
    close (waitH);
    
%% get properties of mrstruct
function dataStruct = local_get_mrstruct_properties(dataStruct)  
    edges = zeros(4,4,dataStruct.numSlices);
    if dataStruct.eDicom
        patterns = ["FH","AP","RL"];
        dirlist  = dir(dataStruct.dataDirectory);
        dirlist  = dirlist(contains({dirlist.name},patterns(1)));
        actPath  = fullfile(dataStruct.dataDirectory,dirlist.name,dataStruct.fileNamesMag);
        actMRstruct   = dicom_read_singlefile(actPath,1,'dicom-dict_philips.txt');
        edges = actMRstruct.edges;    
    else
        dirStrMag       = fullfile(dataStruct.dataDirectory,'mag');  
        %% get edges for one volume
        for j=1:dataStruct.numSlices
            actPath       = full_filename(dirStrMag,squeeze(dataStruct.fileNamesMag(j,1,:)));
            actMRstruct   = dicom_read_singlefile(actPath,1);
            edges(:,:,j)  = actMRstruct.edges;
        end
        %% end of: get edges for one volume
    end
    dataStruct.propMRstruct        = actMRstruct;
    dataStruct.propMRstruct.dataAy = [];
    %% calculate edges   
    dataStruct.propMRstruct.edges  = get_volume_geo(edges);
%%% end of : get properties of mrstruct 

%% read data for all phases from the choosen slice
function dataStruct = local_read_all_phases(dataStruct,actSlice,flag)

    numPhases   = dataStruct.numPhases;    
    venc        = dataStruct.venc;
    signVijk    = dataStruct.signVijk;
    dirStrFlow  = fullfile(dataStruct.dataDirectory,'flow');
    dirStrMag   = fullfile(dataStruct.dataDirectory,'mag');
    
    %% preallocating arrays for speed
    imaMAG    = zeros(dataStruct.imageSize(1),dataStruct.imageSize(2),numPhases,'single');
    imaMxFlow = zeros(dataStruct.imageSize(1),dataStruct.imageSize(2),numPhases,dataStruct.nFE,'single');
    
    %% get file names for actual slice
    %%fileNamesMag   = squeeze(dataStruct.fileNamesMag(actSlice,:,:));
    %%fileNamesFlow  = squeeze(dataStruct.fileNamesFlow(actSlice,:,:,:)); 
    
    %% start point for the counter
    if flag == 1; %% only for preview
      num = 1;
    elseif flag == 2 %% only mrstruct is selected
      num = actSlice;
    elseif flag == 3 %% avi / EnSight are selected
      num = actSlice;
    end
    
    %% loading mag and flow data for single slice //MM
    for m = 1:numPhases
		imaMAG(:,:,m) = dataStruct.dataMag(:,:,actSlice,m);		
        for n = 1:dataStruct.nFE
			imaMxFlow(:,:,m,n) = dataStruct.dataFlow(:,:,actSlice,m,n);		
        end
    end

    %% scale flow images to [-venc venc]    
    tmpMask1    = ones(size(imaMxFlow));
    for m=1:dataStruct.nFE
        tmpMask1(:,:,:,m)= venc(m)*tmpMask1(:,:,:,m);
    end
    tmpMask2    = tmpMask1;
    tmpMask1    = tmpMask1./(dataStruct.phaseRange/2);
    imaMxFlow   = imaMxFlow.*tmpMask1;
    imaMxFlow   = imaMxFlow - tmpMask2;
    %% end of: scale flow images to [-venc venc]    

    %% write data to dataStruct
    if dataStruct.nFE==1 % velocity through slice
       imaMxFlow  = signVijk(3)*imaMxFlow; 
    elseif dataStruct.nFE==3
        imaMxFlow(:,:,:,1)  = signVijk(1)*imaMxFlow(:,:,:,1);
        imaMxFlow(:,:,:,2)  = signVijk(2)*imaMxFlow(:,:,:,2);
        imaMxFlow(:,:,:,3)  = signVijk(3)*imaMxFlow(:,:,:,3);
    end
    
    dataStruct.imaFLOW  = imaMxFlow;
    dataStruct.imaMAG   = imaMAG;
    %% Read data for all phases from the chosen slice

function handleStruct = local_set_autosegvers(handleStruct,serverPath)     
    sitelist = get(handleStruct.site,'String');
    site     = sitelist{get(handleStruct.site,'Value')};
    filePath   = fullfile(serverPath,[site,'vers.txt']);
    fid = fopen(filePath,'r');
    i = 0;
    while ~feof(fid)
        i = i + 1;
        verlist{i,1} = fgetl(fid);
    end
    fclose(fid);
    set(handleStruct.autosegvers,'String',verlist);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% filter functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%  function to remove noise    
function dataStruct = flow_noise_filter(dataStruct,noiseFactor,filtComb)
    %% using a mask created from the magnitude images,
    %% velocity map pixel is recognized as noise if the
    %% corresponding magnitude pixel value is below the noiseFactor
    
    %%% creating noise mask
    imaMAG      = squeeze( dataStruct.imaMAG);        
    noiseFactor = noiseFactor*(max(imaMAG(:))-min(imaMAG(:)));%% Jelena 20080312
    noiseMask   = (imaMAG > noiseFactor);    
    noiseMask   = any(noiseMask,3);   
    %%% end of: creating noise mask
    
    %% mask combination
    if filtComb == 1
        noiseMask = and(noiseMask,dataStruct.filterMask); 
    elseif filtComb == 2
        noiseMask = or(noiseMask,dataStruct.filterMask); 
    end
    
    %% apply noise mask to velocity images        
    noiseMask4D           = repmat(noiseMask,[1 1 dataStruct.numPhases dataStruct.nFE]);
    dataStruct.imaFLOW    = noiseMask4D.*(dataStruct.imaFLOW);
    %% end of: apply noise mask to velocity images  
    
    %% write data to dataStruct
    dataStruct.filterMask = noiseMask;   
    %% end of: write data to dataStruct
    
%%%% end of : function to remove noise    

%% median filter function
function dataStruct = flow_median_filter(dataStruct)
    imaFLOW   = ones(size(dataStruct.imaFlowOrig));
    
    for k=1:dataStruct.nFE
        for m=1:dataStruct.numPhases
            imaFLOW(:,:,m,k)=medfilt2(dataStruct.imaFLOW(:,:,m,k));
        end
    end
    dataStruct.imaFLOW         = imaFLOW;
%%%% end of : median filter function   

%% derivative filter function (for FRANCESCO)
function dataStruct = flow_derivat_filter(dataStruct, noiseFactor)
    %% Derivative of the time course is calculated. In case of noise this
    %% derivative has got great variability, and pixels above a certain
    %% level are set to zero.
    
    % Setup mask
    %flowSize = size(dataStruct.imaFLOW);
    mask = zeros(dataStruct.imageSize);
    nPhases = dataStruct.numPhases;
    for flowEncDir = 1:dataStruct.nFE
        for phase = 2:nPhases
            mask = mask + abs( dataStruct.imaFLOW(:,:,phase,flowEncDir) - dataStruct.imaFLOW(:,:,phase-1,flowEncDir) );
        end
    end

    mask = mask/max(mask(:));

    threshold = 10^(-noiseFactor*5);
    
    mask = (mask < threshold);
    
    dataStruct.imaFLOW = dataStruct.imaFLOW .* repmat(mask, [1 1 nPhases dataStruct.nFE]);


%% function to calculate standard deviation
function dataStruct = flow_stdev(dataStruct)
    stdFlow = zeros(dataStruct.imageSize(1),dataStruct.imageSize(2),dataStruct.nFE);
    imaMAG = (dataStruct.imaMAG);
    imaMAG = (imaMAG - min(imaMAG(:)))/(max(imaMAG(:))-min(imaMAG(:))) ;
    %% calculate standard deviation along temporal direction
    for n=1:(dataStruct.nFE)
        stdFlow(:,:,n) = squeeze(std(dataStruct.imaFLOW(:,:,:,n).*(1-imaMAG.^2),0,3));
    end
    
    dataStruct.stdFlow     = stdFlow;
    dataStruct.stdFlowMax  = squeeze(max(max(stdFlow)));     
%%%% end of: function to calculate standard deviation

%% function to set noisy pixel to zero
%%% if stdev is above a special value the pixel is set to zero
function dataStruct = flow_stdev_filter(dataStruct,noiseFactor,filtComb)

    nFE         = dataStruct.nFE;
    stdFlow     = dataStruct.stdFlow;
    stdFlowMax  = dataStruct.stdFlowMax;    
    
    %% find indices of noisy pixels
    for m=1:nFE
         indexTemp =find(stdFlow(:,:,m)==0 | stdFlow(:,:,m)>noiseFactor*stdFlowMax(m));
         if m==1
             index = indexTemp;
         else
             index = intersect(index, indexTemp);
         end
    end
    %% end of: find indices of noisy pixels
    
    %% create noise mask    
    noiseMask        = true(dataStruct.imageSize);
    noiseMask(index) = false;
    %% end of: create noise mask
    
    %% or /and combination of noise masks
    if filtComb ==1
        noiseMask = and(noiseMask,dataStruct.filterMask); 
    elseif filtComb==2
        noiseMask = or(noiseMask,dataStruct.filterMask); 
    end
    %% end of: or /and combination of noise masks
    
    %% apply noise mask to data
    noiseMask4D    = repmat(noiseMask,[ 1 1 dataStruct.numPhases nFE]);
    imaFLOW        = noiseMask4D.*(dataStruct.imaFLOW);
    %% end of: apply noise mask to data

    dataStruct.imaFLOW         = imaFLOW ;
    dataStruct.filterMask      = noiseMask;
%%%% end of: function to set noisy pixel to zero


%% function to filter static tissue out
function dataStruct = flow_static_tissue_filter(dataStruct)    
    %% creating noise mask    
    noiseMask = (dataStruct.filterMask)&(~dataStruct.statMaskFull);
    
    %% apply noise mask to velocity images 
    noiseMask4D           = repmat(noiseMask,[ 1 1 dataStruct.numPhases dataStruct.nFE]);
    imaFLOW               = noiseMask4D.*(dataStruct.imaFLOW);
    %% end of: apply noise mask to velocity images 
    
    dataStruct.filterMask = noiseMask;
    dataStruct.imaFLOW    = imaFLOW;
  %%% end of :function to filter static tissue out

%%  function to calculate static tissue mask
function dataStruct = create_static_tissue_mask(dataStruct,staticFactor,actSlice)

    nFE         = dataStruct.nFE;
    stdFlow     = dataStruct.stdFlow;
    stdFlowMax  = dataStruct.stdFlowMax;
    imaMAG      = (dataStruct.imaMAG);
    imaMAG      = (imaMAG - min(imaMAG(:)))/(max(imaMAG(:))-min(imaMAG(:))) ;
    meanMag     = squeeze(mean(imaMAG,3));
    
   % select region of interest for eddy current correction if not aready done //MM
   if strcmp(dataStruct.statRegStatus, 'on')
       if isempty(dataStruct.eddyMask)		  
                %% crop data for PC-MRA calculations
                h_eddy = figure;
                set(h_eddy,'Name','Crop 4D flow data for eddy current correction');
                % calculate MIP of time-averaged magnitude data
                if dataStruct.numSlices > 1
                    minSlc = floor(1 + dataStruct.numSlices/6);
                    maxSlc = ceil(dataStruct.numSlices - dataStruct.numSlices/6);
                    ImgMag = squeeze( max( mean(dataStruct.dataMag(:,:,minSlc:maxSlc,:),4),[],3) );
                else 
                    ImgMag = squeeze( max( mean(dataStruct.dataMag,4),[],3) );
                end
                h_im = imshow(ImgMag,[]);
                text(0,-10,'Left click, adjust Rectangular ROI and left click twice, and wait for result','FontSize',[10],'Color', 'b');
                rec = imrect;  % get a region R inside a rectangle, BW is a binary image with 1 and 0 inside or outside the rectangle;
                wait(rec);
                ImgMask = createMask(rec,h_im);
                close(h_eddy);
                dataStruct.eddyMask = ImgMask;        		
       end
   end
    
        %% Separate static regions from noise and blood-flow

            for k=1:nFE

                 indexTemp =find(stdFlow(:,:,k).*(1-meanMag.^2) < staticFactor*stdFlowMax(k)& stdFlow(:,:,k) >0);
                 if k==1
                     index = indexTemp;
                 else
                     index = intersect(index, indexTemp);
                 end
             end

             [szy szx c]         = size(stdFlow);
             staticMask          = zeros(szy,szx); 
             staticMask(index)   = 1;
             %%% remove single pixels and close holes
             nloops = 100;
             staticMask = bwmorph(staticMask,'majority',nloops);
             
         dataStruct.statMaskFull = staticMask;
         if strcmp(dataStruct.statRegStatus, 'on')
              dataStruct.statMask     = staticMask.*dataStruct.eddyMask;
         else
              dataStruct.statMask     = staticMask; 
         end
         
         %% end of: Separate static regions from noise and blood-flow
%%%% end of : create_static_tissue_mask  
    
  
%%  function to correct eddy currents
function dataStruct = flow_eddy_current_corr(dataStruct,staticFactor,actSlice)

    nFE         = dataStruct.nFE;
    stdFlow     = dataStruct.stdFlow;
    stdFlowMax  = dataStruct.stdFlowMax;
    imaMAG      = (dataStruct.imaMAG);
    imaMAG      = (imaMAG - min(imaMAG(:)))/(max(imaMAG(:))-min(imaMAG(:))) ;
    meanMag     = squeeze(mean(imaMAG,3));
    numPhases   = dataStruct.numPhases;
    imaFLOW     = dataStruct.imaFLOW;
    
   % select region of interest for eddy current correction if not aready done //MM
   if isempty(dataStruct.eddyMask)		  
			%% crop data for PC-MRA calculations
			h_eddy = figure;
			set(h_eddy,'Name','Crop 4D flow data for eddy current correction');
			% calculate MIP of time-averaged magnitude data
			if dataStruct.numSlices > 1
			    minSlc = floor(1 + dataStruct.numSlices/6);
				maxSlc = ceil(dataStruct.numSlices - dataStruct.numSlices/6);
				ImgMag = squeeze( max( mean(dataStruct.dataMag(:,:,minSlc:maxSlc,:),4),[],3) );
			else 
			    ImgMag = squeeze( max( mean(dataStruct.dataMag,4),[],3) );
			end
			h_im = imshow(ImgMag,[]);
            figure(h_eddy);
			text(0,-10,'Left click, adjust Rectangular ROI and left click twice, and wait for result','FontSize',[10],'Color', 'b');
			rec = imrect;  % get a region R inside a rectangle, BW is a binary image with 1 and 0 inside or outside the rectangle;
			wait(rec);
			ImgMask = createMask(rec,h_im);
			close(h_eddy);
			dataStruct.eddyMask = ImgMask;        		
   end
     
    if strcmp(dataStruct.encodingType,'velocity')
        
            %% Separate static regions from noise and blood-flow

            for k=1:nFE

                 indexTemp =find(stdFlow(:,:,k).*(1-meanMag.^2) < staticFactor*stdFlowMax(k)& stdFlow(:,:,k) >0);
                 if k==1
                     index = indexTemp;
                 else
                     index = intersect(index, indexTemp);
                 end
             end

             [szy szx c]         = size(stdFlow);
             staticMask          = zeros(szy,szx); 
             staticMask(index)   = 1;
             %%% remove single pixels and close holes
             nloops = 100;
             staticMask = bwmorph(staticMask,'majority',nloops);
             
         dataStruct.statMaskFull = staticMask;
         dataStruct.statMask     = staticMask.*dataStruct.eddyMask;
         
         %% end of: Separate static regions from noise and blood-flow
         
         if sum(dataStruct.statMask(:))>0
             %% correct eddy-currents if selected 
                % Fit plane to static regions (least squares) as:
                % w = alpha*x + beta*y + phi
                % Weight fit by the (squared) magnitude image
                % use last phase (late diastole) for correction
                % apply correction only to non-zero pixels
                
             % needed for least square fit
                X = double(repmat(1:szx,[szy 1]));     
                Y = double(repmat(fliplr(1:szy)',[ 1 szx]));
                for k=1:nFE
                    %% calculate factors for last square fit
                    [phi alpha beta]        = lsqf2(imaFLOW(:,:,numPhases,k),dataStruct.statMask,X,Y);
                    %% calculate fit plane & mask of non-zero pixels
                    fitPlane                = (phi + alpha*X + beta*Y);
                    
                    
                    %start LiliMa: Consider using fit in the future for
                    %ease of changing between first and higher order fits. 
%                     X2 = reshape(X, szy*szx, 1);
%                     Y2 = reshape(Y,szy*szx, 1);
%                     tempMx =reshape(imaFLOW(:,:,numPhases,k).*dataStruct.statMask, szy*szx,1);
%                     
%                     fit2_LM = fit([X2,Y2], tempMx, 'poly22','exclude', tempMx == 0); %second order correction 
% 
%                     phi2 = fit2_LM.p00; alpha2 = fit2_LM.p10; beta2 = fit2_LM.p01;
%                     alpha_2 = fit2_LM.p20; beta_2 = fit2_LM.p02; gamma = fit2_LM.p11; 
%                     
%                     fitPlane  = (phi2 + alpha2*X + beta2*Y + alpha_2*X.^2 + beta_2*Y.^2 + gamma* X.*Y);
%                     
                    %end LiliMa
                    zeroMask                = (imaFLOW(:,:,numPhases,k) ~= 0);
                    
                    
                    
                    %% reproduce masked fit plane up to number of time frames
                    factor3D                = repmat(fitPlane.*zeroMask,[1 1 numPhases]);
                    %% substract fit planes from velocity data
                    imaFLOW(:,:,:,k)        = imaFLOW(:,:,:,k) - factor3D;
                    
                    ec_offset(k) = phi;
                    ec_alpha(k)  = alpha;
                    ec_beta(k)   = beta;
                    
                end
                                    
         end
    elseif strcmp(dataStruct.encodingType,'acceleration')
        for k=1:nFE
            %% offset using last time frame
            offset = squeeze(imaFLOW(:,:,numPhases-2:numPhases,k));
            offset = squeeze(mean(offset,3));
%             %% offset using last and first time frames
%             offset = squeeze((imaFLOW(:,:,numPhases,k)+imaFLOW(:,:,1,k))/2);
%             %% offset using first time frames
%             offset = squeeze(imaFLOW(:,:,1,k));
            factor3D  = repmat(offset,[1 1 numPhases]);
            %% substract fit planes from velocity data
            imaFLOW(:,:,:,k)        = imaFLOW(:,:,:,k) - factor3D;
        end        
    end
    
    dataStruct.imaFLOW         = imaFLOW ;
     %% end of: correct eddy-currents if selected 
%%%% end of : function to correct eddy current
%%
function [actImage,unwrapInd] = local_unwrap(dataStruct,actImage,range, unwrapInd)
    %% if points to unwrap would not defined yet
    if isempty(unwrapInd)
        bin_mask = magicwand2(actImage, range);        
        unwrapInd = find(bin_mask==true);        
    end
    newframe = zeros(dataStruct.imageSize);    
    value    = actImage(unwrapInd); 
    newframe(unwrapInd) = (sign(value)*(dataStruct.venc(2))*2);
    actImage = actImage - newframe;
    
%% function to calculate 2nd order correction for FRANCESCO
function dataStruct = flow_secOrder_corr(dataStruct,staticFactor, actSlice)
    nFE         = dataStruct.nFE;
    stdFlow     = dataStruct.stdFlow;
    stdFlowMax  = dataStruct.stdFlowMax;
    numPhases   = dataStruct.numPhases;
    imaFLOW     = dataStruct.imaFLOW;
    meanMag     = squeeze(mean(dataStruct.imaMAG,3));
    
    % select region of interest for eddy current correction if not aready done //MM
   if isempty(dataStruct.eddyMask)		  
			%% crop data for PC-MRA calculations
			h_eddy = figure;
			set(h_eddy,'Name','Crop 4D flow data for eddy current correction');
			% calculate MIP of time-averaged magnitude data
			if dataStruct.numSlices > 1
			    minSlc = floor(1 + dataStruct.numSlices/6);
				maxSlc = ceil(dataStruct.numSlices - dataStruct.numSlices/6);
				ImgMag = squeeze( max( mean(dataStruct.dataMag(:,:,minSlc:maxSlc,:),4),[],3) );
			else 
			    ImgMag = squeeze( max( mean(dataStruct.dataMag,4),[],3) );
			end
			h_im = imshow(ImgMag,[]);
			text(0,-10,'Left click, adjust Rectangular ROI and left click twice, and wait for result','FontSize',[10],'Color', 'b');
			rec = imrect;  % get a region R inside a rectangle, BW is a binary image with 1 and 0 inside or outside the rectangle;
			wait(rec);
			ImgMask = createMask(rec,h_im);
			close(h_eddy);
			dataStruct.eddyMask = ImgMask;        		
   end
 
            %% Separate static regions from noise and blood-flow

            for k=1:nFE

                 indexTemp =find(stdFlow(:,:,k).*(1-meanMag.^2) < staticFactor*stdFlowMax(k)& stdFlow(:,:,k) >0);
                 if k==1
                     index = indexTemp;
                 else
                     index = intersect(index, indexTemp);
                 end
             end

             [szy szx c]         = size(stdFlow);
             staticMask          = zeros(szy,szx); 
             staticMask(index)   = 1;
             %%% remove single pixels and close holes
             nloops = 100;
             staticMask = bwmorph(staticMask,'majority',nloops);
             
         dataStruct.statMaskFull = staticMask;
         dataStruct.statMask     = staticMask.*dataStruct.eddyMask;
  
     %% end of: Separate static regions from noise and blood-flow
     if sum(dataStruct.statMask(:))>0     
         %% Second order correction: eddy currents + concomitant fields
            % Fit paraboloid to static regions (least squares) as:
            % w = params(1)*X + params(2)*Y + params(3)*(X.^2) + params(4)*(Y.^2) + params(5)*(X.*Y) + params(6);
            % use last phase (late diastole) for correction
            % apply correction only to non-zero pixels
         if strcmp(dataStruct.statRegStatus, 'off')
         % needed for least square fit (X-Y independent variables)
            X = double(repmat(1:szx,[szy 1]));     
            Y = double(repmat(fliplr(1:szy)',[ 1 szx]));
            for k=1:nFE
                params              = fitParaboloid(imaFLOW(:,:,numPhases,k),dataStruct.statMas,X,Y);

                tmpMask             =(not(imaFLOW(:,:,numPhases,k) == 0));

                factor              = params(1)*X + params(2)*Y + params(3) * (X.^2) + params(4) * (Y.^2) + params(5) * (X.*Y) + params(6);
                factor3D            = repmat(factor,[1 1 numPhases]);
                zeroMask            = ones(szy,szx);

                zeroMask            = zeroMask.*tmpMask;
                zeroMask3D          = repmat(zeroMask,[1 1 numPhases]);
                imaFLOW(:,:,:,k)    = imaFLOW(:,:,:,k) - factor3D.*zeroMask3D;
            end        
            dataStruct.imaFLOW      = imaFLOW ;
         end 
     end
%%%% end of : function for second order correction

        
function threslim1(hObj,event, volume_pcmra, pcmra_thresval) 
    % Called to set zlim of surface in figure axes
    % when user moves the slider control
    pcmra_thresval(1) = get(hObj,'Value');   
    subplot(1,3,1);
    cla
    p1 = patch(isosurface(volume_pcmra(:,:,:,1), pcmra_thresval(1)),...
    'FaceColor','red','EdgeColor','none');
    isonormals(volume_pcmra(:,:,:,1),p1)
    axis off
    axis tight
    camlight left 
    lighting phong;
   
%%%% end of : function for pcmra threshold1

function threslim2(hObj,event, volume_pcmra, pcmra_thresval) 
    % Called to set zlim of surface in figure axes
    % when user moves the slider control
    pcmra_thresval(2) = get(hObj,'Value');
    subplot(1,3,2);
    cla
    p2 = patch(isosurface(volume_pcmra(:,:,:,2), pcmra_thresval(2)),...
    'FaceColor','red','EdgeColor','none');
    isonormals(volume_pcmra(:,:,:,2),p2)
    axis off
    axis tight
    camlight left 
    lighting phong;
    
%%%% end of : function for pcmra threshold2

function threslim3(hObj,event, volume_pcmra, pcmra_thresval) 
    % Called to set zlim of surface in figure axes
    % when user moves the slider control
    pcmra_thresval(3) = get(hObj,'Value');
    subplot(1,3,3);
    cla
    p3 = patch(isosurface(volume_pcmra(:,:,:,3), pcmra_thresval(3)),...
    'FaceColor','red','EdgeColor','none');
    isonormals(volume_pcmra(:,:,:,3),p3)
    axis off
    axis tight
    camlight left 
    lighting phong;
   
%%%% end of : function for pcmra threshold3

function [volume_pcmra] = proc_pcmra (dataStruct)
        %%       
        volume_pcmraSumSquares = dataStruct.volume_pcmraSumSquares;
        volume_pcmraMeanAbsVel = dataStruct.volume_pcmraMeanAbsVel;
        volume_pcmraPseudoComplDiff = dataStruct.volume_pcmraPseudoComplDiff;        

        %       x=1000; % # bins
        %        vss = find(volume_pcmraSumSquares(:)>0);
        %        vms = find(volume_pcmraMeanAbsVel(:)>0);
        %        vcd = find(volume_pcmraPseudoComplDiff(:)>0);
        %        volume_pcmraSumSquaresr = volume_pcmraSumSquares(vss);
        %        volume_pcmraMeanAbsVelr = volume_pcmraMeanAbsVel(vms);
        %        volume_pcmraPseudoComplDiffr = volume_pcmraPseudoComplDiff(vcd);
        %        
        %        figure, hist(volume_pcmraSumSquaresr,x),title('Orignal PCMRA 1(mean of squares)');
        %        figure, hist(volume_pcmraMeanAbsVelr,x),title('Orignal PCMRA 2(systolic mean of squares)');
        %        figure, hist(volume_pcmraPseudoComplDiffr,x),title('Orignal PCMRA 3(pseudo complex difference)');
        % 
        if strcmp(dataStruct.pcmraMfilterq,'Yes')
           % 3D median filter
           volume_pcmraSumSquares1 = ordfilt3(volume_pcmraSumSquares,'med',3);
           volume_pcmraMeanAbsVel1 = ordfilt3(volume_pcmraMeanAbsVel,'med',3);
           volume_pcmraPseudoComplDiff1 = ordfilt3(volume_pcmraPseudoComplDiff,'med',3);

        %            vss = find(volume_pcmraSumSquares1(:)>0);
        %            vms = find(volume_pcmraMeanAbsVel1(:)>0);
        %            vcd = find(volume_pcmraPseudoComplDiff1(:)>0);
        %            volume_pcmraSumSquaresr = volume_pcmraSumSquares1(vss);
        %            volume_pcmraMeanAbsVelr = volume_pcmraMeanAbsVel1(vms);
        %            volume_pcmraPseudoComplDiffr = volume_pcmraPseudoComplDiff1(vcd);
        % 
        %            figure, hist(volume_pcmraSumSquaresr,x),title('Median Filter PCMRA 1(mean of squares)');
        %            figure, hist(volume_pcmraMeanAbsVelr,x),title('Median Filter PCMRA 2(systolic mean of squares)');
        %            figure, hist(volume_pcmraPseudoComplDiffr,x),title('Median Filter PCMRA 3(pseudo complex difference)');        
        else
            volume_pcmraSumSquares1 = volume_pcmraSumSquares;
            volume_pcmraMeanAbsVel1 = volume_pcmraMeanAbsVel;
            volume_pcmraPseudoComplDiff1 = volume_pcmraPseudoComplDiff;
        end
        %       
        if strcmp(dataStruct.pcmraEqu,'Yes')
            % histogram equalization        
            % Calculate over full data excluding zeros
            [ys, xs, zs] = size(volume_pcmraSumSquares1);                   
            ps = ys*xs*zs;
            xd = 1/ps:1/ps:1;
            p = raylpdf(xd,0.5);
            
            volume_pcmraSumSquares2r = volume_pcmraSumSquares1(:);
            volume_pcmraMeanAbsVel2r = volume_pcmraMeanAbsVel1(:);
            volume_pcmraPseudoComplDiff2r = volume_pcmraPseudoComplDiff1(:);

            vss = find(volume_pcmraSumSquares2r(:)>0);
            vms = find(volume_pcmraMeanAbsVel2r(:)>0);
            vcd = find(volume_pcmraPseudoComplDiff2r(:)>0);          
            
            volume_pcmraSumSquares2requ = histeq(volume_pcmraSumSquares1(vss));            
            volume_pcmraMeanAbsVel2requ = histeq(volume_pcmraMeanAbsVel1(vms),p);
            volume_pcmraPseudoComplDiff2requ = histeq(volume_pcmraPseudoComplDiff1(vcd),p);
            
            volume_pcmraSumSquares2r(vss) = volume_pcmraSumSquares2requ;
            volume_pcmraMeanAbsVel2r (vms) = volume_pcmraMeanAbsVel2requ;
            volume_pcmraPseudoComplDiff2r (vcd) = volume_pcmraPseudoComplDiff2requ;
            
            volume_pcmraSumSquares2 = reshape(volume_pcmraSumSquares2r, ys, xs, zs)/max(volume_pcmraSumSquares2r);
            volume_pcmraMeanAbsVel2 = reshape(volume_pcmraMeanAbsVel2r, ys, xs, zs)/max(volume_pcmraMeanAbsVel2r);
            volume_pcmraPseudoComplDiff2 = reshape(volume_pcmraPseudoComplDiff2r, ys, xs, zs)/max(volume_pcmraPseudoComplDiff2r);

        %             vss = find(volume_pcmraSumSquares2(:)>0);
        %             vms = find(volume_pcmraMeanAbsVel2(:)>0);
        %             vcd = find(volume_pcmraPseudoComplDiff2(:)>0);
        %             volume_pcmraSumSquaresr = volume_pcmraSumSquares2(vss);
        %             volume_pcmraMeanAbsVelr = volume_pcmraMeanAbsVel2(vms);
        %             volume_pcmraPseudoComplDiffr = volume_pcmraPseudoComplDiff2(vcd);
        %            
        %            figure, hist(volume_pcmraSumSquaresr(:),x),title('Equalized PCMRA 1(mean of squares)');
        %            figure, hist(volume_pcmraMeanAbsVelr(:),x),title('Equalized PCMRA 2(systolic mean of squares)');
        %            figure, hist(volume_pcmraPseudoComplDiffr(:),x),title('Equalized PCMRA 3(pseudo complex difference)');
        else
            volume_pcmraSumSquares2 = volume_pcmraSumSquares1;
            volume_pcmraMeanAbsVel2 = volume_pcmraMeanAbsVel1;
            volume_pcmraPseudoComplDiff2 = volume_pcmraPseudoComplDiff1;
        end
        
        volume_pcmra(:,:,:,1) = volume_pcmraSumSquares2;
        volume_pcmra(:,:,:,2) = volume_pcmraMeanAbsVel2;
        volume_pcmra(:,:,:,3) = volume_pcmraPseudoComplDiff2;

%%%% end of : function for process pcmra

%% MM: local functions for pressure difference calculations
%% ---------------------------------------------------------

% Copyright (c) July 24th, 2006 by University of Freiburg, 
% Dept. of Radiology, Section of Medical Physics


