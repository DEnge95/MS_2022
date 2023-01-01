%Sample inputs
% STLname = 'E:\PhD\Danny_mixedVA\Scans\Training_data_set\Aorta_segmentations\70_2\Aorta_70_2.stl';
% mrStructname = 'E:\PhD\Danny_mixedVA\Scans\Training_data_set\Aorta_segmentations\70_2\mag_struct.mat';
% saveMRstructFlag = 'E:\PhD\Danny_mixedVA\Scans\Training_data_set\Aorta_segmentations\70_2';
STLnames = dir('F:\Results\PWV_Batch_April2022_baselinepaper');
% STLname = 'F:\Scans\MIXF_050_2_170356\3dpc\results_user006\Slicer\Aorta.stl';
% mrStructname = 'F:\Scans\MIXF_050_2_170356\3dpc\results_user006\mrstruct_subject_20190903_user006\mag_struct.mat';
% saveMRstructFlag = 'F:\Scans\MIXF_050_2_170356\3dpc\results_user006\mrstruct_subject_20190903_user006';

% STLsize = 101.6*[1 1 1]; %Size of STL volume in mm
% voxnums = 256*[1 1 1]; %Dimensions of desired mrStruct in voxels
% offsets = [49, 46, -30]; %Offset from 0 coords desired in voxels
for i = 9 %3:(length(STLnames)-1)
    loc = strcat(STLnames(i).folder,'\',STLnames(i).name);
    STLname = strcat(loc,'\Aorta.stl');
    mrStructname = strcat(loc,'\mag_struct.mat');
    saveMRstructFlag = loc;
    stl_to_mrStruct(STLname, mrStructname, saveMRstructFlag)    
end

