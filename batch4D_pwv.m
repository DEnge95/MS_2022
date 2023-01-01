% Batch code for running PWV analysis
%
% Loading:
%   -Single subject (for testing), code will find newest mask
%   -Batch, excel file with 2 columns (SubjectID, mask full path)  
% Branch removal:
%   - Simple, two routines
%       1) Erode then dilate into the original mask
%       2) Aggressive gaussian smoothing, then thresholding
%   - See code itself for more details, both are pretty crudely implemented
% Midline fitting:
%   - Currently relies on skeleton3d from the mathworks file exchange. Mike
%   may also have code using a fast-marching algorithm someplace... but it
%   isn't quite finished right now.
% Plane fitting:
%   - Option 1: manual placement. This will pop-up a GUI [defined in
%   planeSelect.fig and planeSelect.m] and allow placement of an arbitrary
%   number of planes along the aorta.
%   - Option 2: automatic placement. Set the "spacing" parameter to how
%   often you want a plane placed.
% Plane analysis:
%   - Currently calculates the velocity parallel to the midline on each
%   plane, interpolated to a 1 mm gridded plane. The code currently doesn't
%   work with other interpolated resolutions, but could probably be easily
%   modified to do so...
%   - Code is commented out for interpolating each velocity component. 
%   - Interpolated quantities are masked using an interpolated mask
%   - Could easily interpolate a mag_struct or other quantitity you've
%   calculated.
% Other analysis:
%   - Currently does some pulse-wave velocity calculations, but others are
%   possible
% Export
%   - Example of exporting to excel
%   - Several plots are generated, showing modified masks, placed planes
%   and midlines, etc.
%   - For this code an output folder is generated at the specified path,
%   but alteratively some outputs could be written to each subject folder.
%   - One might consider saving the pts structure after analysis is
%   complete in order to simplify re-calculation or further analysis
%   without requiring re-running the entire batch analysis.
%
% Kelly Jarvis, Michael Scott 2019
% Markl Lab, Northwestern

clear variables; clc;
PLOT_MASKS      = true;
PLOT_MIDLINE    = true;
PLOT_REGIONS    = true;
PLOT_BOXPLOT    = true;
SAVE_RESULTS    = true; % Only alters the saving of the data struct, figures will be saved

fprintf('============================================================\n');
fprintf('===                 Batch Processing Code                ===\n');
fprintf('============================================================\n');

% Add functions folder and all subfolders to the path
addpath(genpath([pwd filesep() 'functions']));

tic;

% Make a save folder
saveFolderParent = '';
if isempty(saveFolderParent) 
        saveFolderParent = uigetdir(pwd,'Select a folder to save data');
end
saveFolderRoot = [saveFolderParent filesep() 'Batch_Analysis' filesep() datestr(now,'yyyymmdd_HHMMSS') '_Subject_Results'];
mkdir(saveFolderRoot);
fprintf('- Save folder created:\n   %s\n',saveFolderRoot)

%% Choose subjects

% Single subject / batch entire folder, select in GUI
pts = getSubjectDataPaths();

% % Multiple subjects from excel
% xlpath='';
% if isempty(xlpath) 
%     [file,path] = uigetfile({'*.xlsx;*.xls'},'Select the Excel File');
%     xlpath = [path filesep() file];
%     clear file path
% end
% [~,~,xldata] = xlsread(xlpath);
% xldata(1,:)=[];
% 
% % Get mask path
% for i=1:size(xldata,1)
%     pts(i).Name=xldata{i,1};
%     
%     %%use the following for regular load
%     pts(i).MaskPath=xldata{i,2};
%     %%end of regular load
%     
% %     %%use the following for APS study
% %     if exist([xldata{i,2},'\Aorta_nobranch_grayvalues_mask_struct.mat'],'file')==2
% %         pts(i).MaskPath=[xldata{i,2},'\Aforta_nobranch_grayvalues_mask_struct.mat'];
% %     else
% %         pts(i).MaskPath='not found';
% %     end
% %     %%end of APS study load
% end

% Loop through the subjects
for pt = 1:size(pts,2)
    fprintf('Starting patient %i of %i: %s\n',pt,size(pts,2),pts(pt).Name);
    pts(pt).MaskPath = pts(pt).AoMaskPath;
    try
        %% Loading
        % Find the index of the last slash
        lastSlashIdx = find(pts(pt).MaskPath == filesep(), 1, 'last');
        % Find the mrStruct folder (assume mask is in this folder)
        mrStructFolder = pts(pt).MaskPath(1:lastSlashIdx-1);
        
        % Load the velStruct
        temp = load([mrStructFolder filesep() 'vel_struct.mat']);
        fields = fieldnames(temp);
        vel_struct = temp.(fields{1});
        clear temp fields
        fprintf('  - vel_struct loaded\n')
        
        % Load the mask
        temp = load(pts(pt).MaskPath);
        fields = fieldnames(temp);
        mask_struct = temp.(fields{1});
        clear temp fields
        fprintf('  - mask_struct loaded\n')
        
        % Make the patient save folder
        savefolder = [saveFolderRoot filesep() pts(pt).Name];
        mkdir(savefolder)
        
        %% Find the velocity magnitude
        % Get the velocity magnitude
        velmag = squeeze((vel_struct.dataAy(:,:,:,1,:).^2 + vel_struct.dataAy(:,:,:,2,:).^2 + vel_struct.dataAy(:,:,:,3,:).^2).^0.5);
        
        % Apply the mask to the velmag
        for ii = 1:size(velmag,4)
           velmag(:,:,:,ii) = velmag(:,:,:,ii) .* mask_struct.dataAy; 
        end
        
        %% Work on the mask
        % This section erodes the mask twice,then checks if the erosion was 
        % too aggressive (in which case it erodes only once). Then the
        % resulting mask is dilated twice, checking to make sure that the
        % mask has not exceeded its original bounds after each dilation.
        % The goal here is to eliminate branches (which tend to be
        % thinner), but this is a pretty crude morphological solution and
        % fails ~5-10% of the time. The usual solution is to erode more
        % gently or not at all (if no pre-processing is done in general,
        % the midline fitting is slower and more prone to errors). 
        % Load the mask
        mask = mask_struct.dataAy;
        vox = vel_struct.vox;
        
        % Remove any small islands
        mask = bwareaopen(mask,1000);
        
        % TEST
        N = 5;
        kernel = ones(N, N, N) / N^3;
        smoothedMask = convn(double(mask), kernel, 'same');
        newMask = smoothedMask > 0.75;
        clear smoothedMask
        
        % END TEST
        
%         % Erode the mask twice
%         se = ones(3,3,3);
%         emask = imerode(imerode(mask,se),se);
%         % Remove any small islands
%         emask = bwareaopen(emask,1000);
%         if sum(emask(:) < 1500)
%             % Erosion too aggressive, only do it once!
%             emask = imerode(mask,se);
%             % Remove any small islands
%             emask = bwareaopen(emask,1000);
%         end
%         % Grow the mask twice, but check that it does not outgrow the mask each
%         % time
%         gmask = imdilate(emask,se);
%         % Check the mask hasn't grown too large
%         gmask = gmask & mask;
%         % Grow again
%         gmask = imdilate(gmask,se);
%         % Check again
%         gmask = gmask & mask;
%         
%         % Make a velocity and WSS mask:
%         velmask = imdilate(gmask,se);
%         % Check again
%         velmask = velmask & mask;
%         
%         % Set the new mask
%         newMask = gmask;
        
        if PLOT_MASKS
            figure('units','normalized','outerposition',[0 0 1 1])
            subplot(1,2,1)
            p1iso = isosurface(mask,0.5);
            p1iso.vertices(:,1) = p1iso.vertices(:,1)*vox(2);
            p1iso.vertices(:,2) = p1iso.vertices(:,2)*vox(1);
            p1iso.vertices(:,3) = p1iso.vertices(:,3)*vox(3);
            p1 = patch(p1iso);
            p1.FaceColor = 'red';
            p1.EdgeColor = 'none';
            p1.FaceAlpha = 0.65;
            axis equal
            axis ij
            axis off
            title('Original Mask')
            subplot(1,2,2)
            p2iso = isosurface(newMask,0.5);
            p2iso.vertices(:,1) = p2iso.vertices(:,1)*vox(2);
            p2iso.vertices(:,2) = p2iso.vertices(:,2)*vox(1);
            p2iso.vertices(:,3) = p2iso.vertices(:,3)*vox(3);
            p2 = patch(p2iso);
            p2.FaceColor = 'blue';
            p2.EdgeColor = 'none';
            p2.FaceAlpha = 0.65;
            axis equal
            axis ij
            axis off
            title('Edited Mask')
            
            saveas(gcf,[savefolder filesep() 'mask_visualization.tif'])
            close(gcf)
            clear p1 p1iso p2 p2iso
        end
        
        %% Fit a midline
        [skel, interpolatedMidline, nodefinal, linkfinal] = aorta_midline(newMask,vox);
        fprintf('  - Midline fit\n')
        if PLOT_MIDLINE
            % Plot interpolated midline
            figure('units','normalized','outerposition',[0 0 1 1]);
                     
            % Put a top/bottom on the mask. This closes any surface holes that are due
            % to the lumen leaving the FOV
            top = false(size(mask,1),size(mask,2));
            temp_mask = cat(3,top,mask,top);
            maskpatch = isosurface(temp_mask,0.5);
            
            % Scale the patch based on the voxel size
            maskpatch.vertices(:,1) = maskpatch.vertices(:,1)*vox(2);
            maskpatch.vertices(:,2) = maskpatch.vertices(:,2)*vox(1);
            maskpatch.vertices(:,3) = maskpatch.vertices(:,3)*vox(3);
            
            % Plot the mask
            maskplot = patch(maskpatch);
            maskplot.FaceColor = [0.85 0.85 0.85];
            maskplot.EdgeColor = 'none';
            maskplot.FaceAlpha = 0.5;
            
            axis equal
            axis ij
            axis off
            camlight('right')
            hold on;
            plot3(interpolatedMidline(:,2),interpolatedMidline(:,1),interpolatedMidline(:,3),'LineWidth',2)
            set(gcf,'Color','white');
            title(sprintf('Subject %s midline',regexprep(pts(pt).Name,'_',' ')))
            saveas(gcf,[savefolder filesep() 'Midline.tif'])
            close(gcf)
            clear top temp_mask maskpatch maskplot
        end
        
%         %% Plane placement (manual)
%         data.midline    = interpolatedMidline;
%         data.mask       = mask;
%         data.vox        = vox;
%         
%         % Open the GUI
%         [~,planeIdx] = planeSelect('UserData',data);
%         
%         % Calculate the distance between each point
%         distances = zeros(size(interpolatedMidline,1),1);
%         for ii = 2:numel(distances)
%             distances(ii) = sum((interpolatedMidline(ii-1,:)-interpolatedMidline(ii,:)).^2)^0.5;
%         end
%         
%         % Calculate the cumulative sum
%         cumulativeDistance = cumsum(distances);
        
        %% Find the plane indices (automated placement!)
        spacing = 4; % mm
        [planeIdx,cumulativeDistance] = getMidlineIdx(interpolatedMidline,spacing);
        
        %% Save the midline and the cumulative distance
        pts(pt).interpolatedMidline = interpolatedMidline;
        pts(pt).planeIdx = planeIdx;
        pts(pt).cumulativeDistance = cumulativeDistance;
        
        %% Analyze the planes
        fprintf('  - Beginning plane analysis...\n')
        [pts(pt).subject] = planeAnalysis(mask,vel_struct,interpolatedMidline,planeIdx);
        %[pts(pt).subject] = planeAnalysis(newMask,vel_struct,interpolatedMidline,planeIdx);
        fprintf('  - Plane analysis finished!\n')
        
        %% Plot the planes
        figure('units','normalized','outerposition',[0 0 1 1]);
        % Put a top/bottom on the mask. This closes any surface holes that are due
        % to the lumen leaving the FOV
        top = false(size(mask,1),size(mask,2));
        temp_mask = cat(3,top,mask,top);
        maskpatch = isosurface(temp_mask,0.5);
        
        % Scale the patch based on the voxel size
        maskpatch.vertices(:,1) = maskpatch.vertices(:,1)*vox(2);
        maskpatch.vertices(:,2) = maskpatch.vertices(:,2)*vox(1);
        maskpatch.vertices(:,3) = maskpatch.vertices(:,3)*vox(3);
        
        % Plot the mask
        maskplot = patch(maskpatch);
        maskplot.FaceColor = [0.85 0.85 0.85];
        maskplot.EdgeColor = 'none';
        maskplot.FaceAlpha = 0.5;
        
        axis equal
        axis ij
        axis off
        camlight('right')
        hold on;
        plot3(interpolatedMidline(:,2),interpolatedMidline(:,1),interpolatedMidline(:,3),'LineWidth',2)
        % Plot the planes
        for ii = 1:size(pts(pt).subject,2)
            pc = zeros(2,2,3);
            pc(1,1,:) = pts(pt).subject(ii).PlaneCoords(1,1,:);
            pc(1,2,:) = pts(pt).subject(ii).PlaneCoords(1,end,:);
            pc(2,1,:) = pts(pt).subject(ii).PlaneCoords(end,1,:);
            pc(2,2,:) = pts(pt).subject(ii).PlaneCoords(end,end,:);
            
            plotdata.xvec = [pc(1,1,1) pc(1,2,1) pc(2,2,1) pc(2,1,1)];
            plotdata.yvec = [pc(1,1,2) pc(1,2,2) pc(2,2,2) pc(2,1,2)];
            plotdata.zvec = [pc(1,1,3) pc(1,2,3) pc(2,2,3) pc(2,1,3)];
            
            patch(plotdata.yvec,plotdata.xvec,plotdata.zvec,'b');
            % Plot the normal vector
            %quiver3(maskInfo.midline(maskInfo.plane1,2),maskInfo.midline(maskInfo.plane1,1),maskInfo.midline(maskInfo.plane1,3),scaling*AAoNhat(2),scaling*AAoNhat(1),scaling*AAoNhat(3),0,'black','LineWidth',3,'MaxHeadSize',3)
        end
        set(gcf,'Color','white');
        title(sprintf('%s fit planes',regexprep(pts(pt).Name,'_',' ')))
        saveas(gcf,[savefolder filesep() 'Planes.tif'])
        close(gcf)
        clear top temp_mask maskpatch maskplot
        
        %% Prep for PWV analysis
        fprintf('  - PWV analysis time...\n')
        % flowdata: each row is one plane, column 1 is the length in mm, 
        % columns 2-tt are the flowrate at each timepoint
        flowdata = zeros(size(pts(pt).subject,2),size(vel_struct.dataAy,5)+1);
        flowdata(:,1) = cumulativeDistance(planeIdx);
        for ii = 1:size(pts(pt).subject,2)
           flowdata(ii,2:end) = pts(pt).subject(ii).FlowPerTimepoint;
           %flowdata(ii,2:end) = pts(pt).subject(ii).mid50waveform;
        end
        
        %% Calculate the PWV
        [pwv,rmse] = ao_stiffness('',[savefolder filesep() 'results_pwv'],vel_struct.tr,flowdata);
        
        % Store: in pwv and rmse:
        %   1: TTF
        %   2: XCOR
        %   3: XCOR2
        pts(pt).pwv = pwv;
        pts(pt).rmse = rmse;
        
        close all
        
        fprintf('  - Patient complete, data saved to:\n   %s\n',savefolder)
    catch
        warning('  - Error in patient %s, check this patient later.\n',pts(pt).Name)
    end
    
    % Clear the loaded data to make sure the next one loads correctly.
    clear vel_struct mask_struct mag_struct
end

fprintf('\n\nAll subjects analyzed, attempting to write to excel.\n')
% Make a cell array to output to excel
output = cell(size(pts,2)+1,15);
output(1,:) = {'Patient','PWV TTF','RMSE TTF','PWV XCOR','RMSE XCOR','PWV XCOR2','RMSE XCOR2','PWV XCOR3','RMSE XCOR3','MIN PWV XCOR','MAX PWV XCOR','RANGE PWV XCOR','IQR PWV XCOR','STD PWV XCOR','NormDist PWV XCOR'};
for pt = 1:size(pts,2)
    if ~isempty(pts(pt).Name)
        output{1+pt,1} = pts(pt).Name;
    end
    if ~isempty(pts(pt).pwv)
        output{1+pt,2} = pts(pt).pwv(1);
        output{1+pt,4} = pts(pt).pwv(2);
        output{1+pt,6} = pts(pt).pwv(3);
        output{1+pt,8} = pts(pt).pwv(4);
        output{1+pt,10} = pts(pt).pwv(5);%min
        output{1+pt,11} = pts(pt).pwv(6);%max
        output{1+pt,12} = pts(pt).pwv(7);%range
        output{1+pt,13} = pts(pt).pwv(8);%iqr
        output{1+pt,14} = pts(pt).pwv(9);%std
        output{1+pt,15} = pts(pt).pwv(10);%normDist
    end
    if ~isempty(pts(pt).rmse)
        output{1+pt,3} = pts(pt).rmse(1);
        output{1+pt,5} = pts(pt).rmse(2);
        output{1+pt,7} = pts(pt).rmse(3);
        output{1+pt,9} = pts(pt).rmse(4);
    end
end
xlswrite([saveFolderRoot filesep() 'Results.xlsx'],output)

toc;

fprintf('\n\nExcel written to:\n   %s\n',[saveFolderRoot filesep() 'Results.xlsx'])
fprintf('Consider saving variable pts for future analysis.\n\n')

fprintf('============================================================\n');
fprintf('===                  Execution complete!                 ===\n');
fprintf('==========================================================(:\n');

% Clean up
clearvars -except pts output
close all
