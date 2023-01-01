function maskData = createPDmask(dataStruct,noiseThreshold)
% calculate pressure difference map from velocity data

%%% load from each slice the last/ first time frame
timeFrame= 1;

tmp_flowVolumeData = zeros(dataStruct.imageSize(1),dataStruct.imageSize(2),dataStruct.numSlices,dataStruct.nFE);
magVolumeData  = zeros(dataStruct.imageSize(1),dataStruct.imageSize(2),dataStruct.numSlices);

for k=1:dataStruct.numSlices
    actFileName = squeeze(dataStruct.fileNamesFlow(k,timeFrame,1,:))';
    actFile = fullfile(dataStruct.dataDirectory,'flow',actFileName);
    tmp_flowVolumeData(:,:,k,1)=((dicom_read_singlefile(actFile,0))-2048)*(dataStruct.venc(1))/2048;
    
    actFileName = squeeze(dataStruct.fileNamesFlow(k,timeFrame,2,:))';
    actFile = fullfile(dataStruct.dataDirectory,'flow',actFileName);
    tmp_flowVolumeData(:,:,k,2)=((dicom_read_singlefile(actFile,0))-2048)*(dataStruct.venc(2))/2048;
    
    actFileName = squeeze(dataStruct.fileNamesFlow(k,timeFrame,3,:))';
    actFile = fullfile(dataStruct.dataDirectory,'flow',actFileName);
    tmp_flowVolumeData(:,:,k,3)=((dicom_read_singlefile(actFile,0))-2048)*(dataStruct.venc(3))/2048;
    
    actFileName=squeeze(dataStruct.fileNamesMag(k,timeFrame,:))';
    actFile = fullfile(dataStruct.dataDirectory,'mag',actFileName);
    magVolumeData(:,:,k)=dicom_read_singlefile(actFile,0);
end


%%% calculate velocity magnitude --> sqrt(Vx^2+Vy^2+Vz^2)
flowVolumeData = sqrt(sum(tmp_flowVolumeData.^2,4));
clear('tmp_flowVolumeData');

maskData = ((magVolumeData/max(magVolumeData(:)))>noiseThreshold);

flowVolumeData = flowVolumeData.*maskData;
magVolumeData  = magVolumeData.*maskData;

%%% calculate sum of squares
%maskData = ((flowVolumeData.^2)*magVolumeData);

%%%scale volume data

maskData=double( flowVolumeData );%> max(flowVolumeData(:)) * 0.30 );

% add border region with zeros 
[sx sy sz] = size(maskData);
maskData(1,:,:)  = 0;
maskData(:,1,:)  = 0;
maskData(:,:,1)  = 0;
maskData(sx,:,:) = 0;
maskData(:,sy,:) = 0;
maskData(:,:,sz) = 0;

   