%% //MM
function dataStruct = calculate_pcmra(handleStruct,dataStruct,modPCMRA,vencThreshold,tPCMRA_min,tPCMRA_max,medFilterPCMRA_Flag,histEqualPCMRA_Flag,argvarin)
    
%% calculate velocity magnitude --> sqrt(Vx^2+Vy^2+Vz^2)
    imaSPEED = sqrt(mean(dataStruct.dataFlow.^2,5)); %LiliMa: this is not actually speed. This is the sqrt of mean sum squares (a scaled version of speed) 
    
    imaMAG   = dataStruct.dataMag;
    minV = min(0.7*(imaMAG(:)));
    maxV = max(0.7*(imaMAG(:)));
    tmp  = imaMAG > maxV;    
    imaMAG =(imaMAG.*not(tmp))+tmp*maxV;   
    tmp  = imaMAG < minV;    
    imaMAG =(imaMAG.*not(tmp))+tmp*minV;    
    imaMAG = (imaMAG-minV)/(maxV-minV);
	
	% creat maks for entire 3D volume
	% allocate memory
	pcMRa_mask3D = zeros(dataStruct.imageSize(1),dataStruct.imageSize(2),dataStruct.numSlices);
	for slc = 1:dataStruct.numSlices
		pcMRa_mask3D(:,:,slc) = dataStruct.pcmraMask;
	end
		
    	   
    for k=1:size(modPCMRA,1)
        switch (modPCMRA{k,1})
            case {'modeSumSquares'}
                %% mean over all time frames  
                pcmraData = imaSPEED.*imaMAG;
                pcmraData = squeeze(mean(pcmraData.^2,4)).*pcMRa_mask3D;  % multiply with mask here
                %% scale to [0 1]
                dataStruct.volume_pcmraSumSquares = pcmraData/max(pcmraData(:));                
            case {'modeMeanAbsVel'}
                %% mean over selected frames
                pcmraData = imaSPEED.*imaMAG;
                pcmraData = squeeze(mean(pcmraData(:,:,:,tPCMRA_min:tPCMRA_max).^2,4)).*pcMRa_mask3D;  % multiply with mask here; 
                %% scale to [0 1]
                dataStruct.volume_pcmraMeanAbsVel = (pcmraData)/max(pcmraData(:)); 
            case {'modePseudoComplexDiff'}
                %% PCMRA = Magnitude for Speed >= venc threshold
                maskSpeed = (imaSPEED >= vencThreshold);             
                pcmraData = maskSpeed.*(dataStruct.dataMag)+not(maskSpeed).*(dataStruct.dataMag).*sin(pi*(imaSPEED)/(2*vencThreshold));
                pcmraData = mean(pcmraData,4).*pcMRa_mask3D;  % multiply with mask here;
                %% scale to [0 1]
                dataStruct.volume_pcmraPseudoComplDiff = (pcmraData)/max(pcmraData(:)); 
			case {'modeTimeAveMag'}
                %% time averaged 3D magnitude data  
                pcmraData = imaMAG;
                pcmraData = squeeze(mean(pcmraData,4)).*pcMRa_mask3D;  % multiply with mask here
                %% scale to [0 1]
                dataStruct.volume_pcmraTimeAveMag = pcmraData/max(pcmraData(:));   
            case {'modeSqrtSumSquares'} %% LiliMa
                %% mean over all time frames  
                pcmraData = imaSPEED.*imaMAG;
                pcmraData = sqrt(squeeze(mean(pcmraData.^2*size(dataStruct.dataFlow, 5),4)).*pcMRa_mask3D);  % multiply with mask here; LiliMa: fixed the scaling from the mean rather than sum being used for imaSPEED  Aug 2017
                %% scale to [0 1]
                dataStruct.volume_pcmraSqrtSumSquares = pcmraData/max(pcmraData(:));
            case {'modeVelCorr'} %% MBS
                [pcmraData,~] = corrPCMRA(dataStruct.dataFlow,6);
                % Scale to 0 to 1
%                 pcmraData = pcmraData - min(pcmraData(:));
%                 pcmraData = pcmraData ./ max(pcmraData(:));
                % Remove negative values, the remainder should be scaled [0
                % 1] since it is correlation coefficients
                %pcmraData(pcmraData < 0) = 0;
                % Apply the mask
                pcmraData = pcmraData .* pcMRa_mask3D;
                dataStruct.volume_pcmraVelCorr = pcmraData;
            otherwise
                %% do nothing
        end
    end
	
	%% apply median filter
	if medFilterPCMRA_Flag == 1
		dataStruct.volume_pcmraSumSquares      = medfilt3(dataStruct.volume_pcmraSumSquares);
		dataStruct.volume_pcmraMeanAbsVel      = medfilt3(dataStruct.volume_pcmraMeanAbsVel);
		dataStruct.volume_pcmraPseudoComplDiff = medfilt3(dataStruct.volume_pcmraPseudoComplDiff);
		dataStruct.volume_pcmraTimeAveMag      = medfilt3(dataStruct.volume_pcmraTimeAveMag);
        dataStruct.volume_pcmraSqrtSumSquares  = medfilt3(dataStruct.volume_pcmraSqrtSumSquares);  %% LiliMa
        dataStruct.volume_pcmraVelCorr         = medfilt3(dataStruct.volume_pcmraVelCorr);  %% MBS
	end
	
	%% apply histogram equalization
	if histEqualPCMRA_Flag == 1
	   %% currently not implemented
	end

    % end of: calculate mean PC-MRA data set