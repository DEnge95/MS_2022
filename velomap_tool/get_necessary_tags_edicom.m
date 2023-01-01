function  [dataInfoStruct,numOfFilesMag,eDicom] = get_necessary_tags_edicom(dataInfoStruct)
    eShareInfoStruct = dataInfoStruct.SharedFunctionalGroupsSequence.Item_1;
    eFrameInfoStruct = dataInfoStruct.PerframeFunctionalGroupsSequence.Item_1;
    % add some important fields
    dataInfoStruct.RepetitionTime          = eShareInfoStruct.MRTimingandRelatedParametersSequence.Item_1.RepetitionTime;
    dataInfoStruct.PhaseEncodingDirection  = eShareInfoStruct.MRFOVGeometrySequence.Item_1.PhaseEncodingDirection;
    dataInfoStruct.EchoTime                = eFrameInfoStruct.Private_2005_140F.Item_1.EchoTime;
    dataInfoStruct.SliceThickness          = eFrameInfoStruct.Private_2005_140F.Item_1.SliceThickness;
    dataInfoStruct.SpacingBetweenSlices    = eFrameInfoStruct.Private_2005_140F.Item_1.SpacingBetweenSlices;
    dataInfoStruct.PixelSpacing            = eFrameInfoStruct.Private_2005_140F.Item_1.PixelSpacing; 
    dataInfoStruct.ImageOrientationPatient = eFrameInfoStruct.PlaneOrientationSequence.Item_1.ImageOrientationPatient;
    dataInfoStruct.ImagePositionPatient    = eFrameInfoStruct.PlanePositionSequence.Item_1.ImagePositionPatient;
    dataInfoStruct.NominalInterval         = eFrameInfoStruct.CardiacSynchronizationSequence.Item_1.RRIntervalTimeNominal;
    
    flds_perframe = fieldnames(dataInfoStruct.PerframeFunctionalGroupsSequence);
    num_max = size(flds_perframe,1);
    dims = dataInfoStruct.PerframeFunctionalGroupsSequence.(flds_perframe{num_max}).FrameContentSequence.Item_1.DimensionIndexValues;
    dataInfoStruct.CardiacNumberOfImages = dims(end);
    %if contains(dataInfoStruct.SeriesDescription,'FH')
    if max(size(dims)) == 5
        numOfFilesMag = num_max/2;
        if mod(num_max,2) ~= 0
            error('(eDicom): The total number of images in RH should be even.');
        end
         
        for idata = 1 : dataInfoStruct.CardiacNumberOfImages
            dataInfoStruct.TimeStamps(1, idata) = dataInfoStruct.PerframeFunctionalGroupsSequence.(flds_perframe{idata}).Private_2005_140F.Item_1.TriggerTime;
            dataInfoStruct.TimeStamps(2, idata) = dataInfoStruct.PerframeFunctionalGroupsSequence.(flds_perframe{numOfFilesMag+idata}).Private_2005_140F.Item_1.TriggerTime;
        end
    else
        numOfFilesMag = 0;
        for idata = 1 : dataInfoStruct.CardiacNumberOfImages
            dataInfoStruct.TimeStamps(1, idata) = dataInfoStruct.PerframeFunctionalGroupsSequence.(flds_perframe{idata}).Private_2005_140F.Item_1.TriggerTime;
        end
    end         
    eDicom = true;
