function create_slicer_mrml(mrStruct,labelfnames,volfnames,outname,outpath)
% Output Slicer restart file
% 11/27/2019 Takashi Fujiwara

% parameter check
if nargin == 0 % mrstruct read
    mrStruct = mrstruct_read;
end

if nargin < 2
    [labelfnames,labelpath] = uigetfile('*.nrrd','Select label map files','MultiSelect','on');
end    

if nargin < 3
    [volfnames,volpath] = uigetfile('*.nrrd','Select volume files','MultiSelect','on');
end    

if nargin < 4
    uoutname = inputdlg({'Specify outputfile mrml name (without extension):'},'Input file name',[1 10],{'template'});
    outname = uoutname{1};
    if isempty(outname)
        outname = 'template';
    end
end

if nargin < 5
    outpath = uigetdir('Select output directory');
end

% when dim of cell array = 1, forcibly convert to 1*1 cell array
if ischar(labelfnames) == 1
    labelfnames = {labelfnames};
end

if ischar(volfnames) == 1
    volfnames = {volfnames};
end

nlabel = length(labelfnames);
nvol   = length(volfnames);
[~, lpartnames, ~] = cellfun(@fileparts,labelfnames,'UniformOutput',false); 
[~, vpartnames, ~] = cellfun(@fileparts,volfnames,'UniformOutput',false); 

% sort vectors for Slicer view
vec = [mrStruct.edges(:,2) mrStruct.edges(:,1) mrStruct.edges(:,3) mrStruct.edges(:,4)];
% voxel sizes
vox = mrStruct.vox;
% normal vectors
vec(1:3,1) = vec(1:3,1)/norm(vec(1:3,1));
vec(1:3,2) = vec(1:3,2)/norm(vec(1:3,2));
vec(1:3,3) = vec(1:3,3)/norm(vec(1:3,3));

% Open output file
fid = fopen(fullfile(outpath,[outname,'.mrml']),'w');

% Headers
fprintf(fid,'<?xml version="1.0" encoding="ISO-8859-1"?>\r\n');
fprintf(fid,'<MRML  version="Slicer4.4.0" userTags="">\r\n');

% window settings
fprintf(fid,'<Crosshair\r\n  id="vtkMRMLCrosshairNodedefault" name="Crosshair" hideFromEditors="true" selectable="true" selected="false" singletonTag="default" crosshairMode="NoCrosshair" navigation="false" crosshairBehavior="OffsetJumpSlice" crosshairThickness="Fine" crosshairRAS="0 0 0"></Crosshair>\r\n');
fprintf(fid,' <Selection\r\n  id="vtkMRMLSelectionNodeSingleton" name="Selection" hideFromEditors="true" selectable="true" selected="false" singletonTag="Singleton" frequencyUnitNodeRef="vtkMRMLUnitNodeApplicationFrequency" intensityUnitNodeRef="vtkMRMLUnitNodeApplicationIntensity" lengthUnitNodeRef="vtkMRMLUnitNodeApplicationLength" timeUnitNodeRef="vtkMRMLUnitNodeApplicationTime" velocityUnitNodeRef="vtkMRMLUnitNodeApplicationVelocity" references="unit/frequency:vtkMRMLUnitNodeApplicationFrequency;unit/intensity:vtkMRMLUnitNodeApplicationIntensity;unit/length:vtkMRMLUnitNodeApplicationLength;unit/time:vtkMRMLUnitNodeApplicationTime;unit/velocity:vtkMRMLUnitNodeApplicationVelocity;" activeVolumeID="vtkMRMLScalarVolumeNode1" secondaryVolumeID="NULL" activeLabelVolumeID="vtkMRMLLabelMapVolumeNode1" activeFiducialListID="NULL" activePlaceNodeID="NULL" activePlaceNodeClassName="NULL" activeROIListID="NULL" activeCameraID="NULL" activeTableID="NULL" activeViewID="NULL" activeLayoutID="NULL" activePlotChartID="NULL" ></Selection>\r\n');
fprintf(fid,' <Interaction\r\n  id="vtkMRMLInteractionNodeSingleton" name="Interaction" hideFromEditors="true" selectable="true" selected="false" singletonTag="Singleton" currentInteractionMode="ViewTransform" placeModePersistence="false" lastInteractionMode="ViewTransform" ></Interaction>\r\n');
fprintf(fid,' <View\r\n  id="vtkMRMLViewNode1" name="View1" hideFromEditors="false" selectable="true" selected="false" singletonTag="1" attributes="MappedInLayout:1" layoutLabel="1" layoutName="1" active="false" visibility="true" backgroundColor="0.756863 0.764706 0.909804" backgroundColor2="0.454902 0.470588 0.745098" layoutColor="0.454902 0.513725 0.913725" orientationMarkerType="none" orientationMarkerSize="medium" rulerType="none" AxisLabels="L;R;P;A;I;S" fieldOfView="200" letterSize="0.05" boxVisible="true" fiducialsVisible="true" fiducialLabelsVisible="true" axisLabelsVisible="true" axisLabelsCameraDependent="true" animationMode="Off" viewAxisMode="LookFrom" spinDegrees="2" spinMs="5" spinDirection="YawLeft" rotateDegrees="5" rockLength="200" rockCount="0" stereoType="NoStereo" renderMode="Perspective" useDepthPeeling="0" gpuMemorySize="0" expectedFPS="8" volumeRenderingQuality="Adaptive" raycastTechnique="Composite" volumeRenderingSurfaceSmoothing="0" volumeRenderingOversamplingFactor="2" linkedControl="0" ></View>\r\n');
fprintf(fid,' <Slice\r\n  id="vtkMRMLSliceNodeRed" name="Red" hideFromEditors="false" selectable="true" selected="false" singletonTag="Red" attributes="MappedInLayout:1" layoutLabel="R" layoutName="Red" active="false" visibility="true" backgroundColor="0 0 0" backgroundColor2="0 0 0" layoutColor="0.952941 0.290196 0.2" orientationMarkerType="none" orientationMarkerSize="medium" rulerType="none" AxisLabels="L;R;P;A;I;S" fieldOfView="302.787 335.699 2.5" dimensions="460 510 1" xyzOrigin="0 0 0" sliceResolutionMode="1" uvwExtents="302.787 335.699 2.5" uvwDimensions="256 256 1" uvwOrigin="0 0 0" activeSlice="0" layoutGridRows="1" layoutGridColumns="1" sliceToRAS="-1 0 0 -14.6404 0 1 0 14.2311 0 0 1 59.8643 0 0 0 1" orientationMatrixAxial="-1 0 0 0 1 0 0 0 1" orientationMatrixSagittal="0 0 1 -1 0 0 0 1 0" orientationMatrixCoronal="-1 0 0 0 0 1 0 1 0" orientation="Axial" orientationReference="Axial" jumpMode="1" sliceVisibility="false" widgetVisibility="false" widgetOutlineVisibility="true" useLabelOutline="false" sliceSpacingMode="0" prescribedSliceSpacing="1 1 1" ></Slice>\r\n');
fprintf(fid,' <Slice\r\n  id="vtkMRMLSliceNodeYellow" name="Yellow" hideFromEditors="false" selectable="true" selected="false" singletonTag="Yellow" attributes="MappedInLayout:1" layoutLabel="Y" layoutName="Yellow" active="false" visibility="true" backgroundColor="0 0 0" backgroundColor2="0 0 0" layoutColor="0.929412 0.835294 0.298039" orientationMarkerType="none" orientationMarkerSize="medium" rulerType="none" AxisLabels="L;R;P;A;I;S" fieldOfView="335.699 371.381 2.50001" dimensions="461 510 1" xyzOrigin="0 0 0" sliceResolutionMode="1" uvwExtents="335.699 371.381 2.50001" uvwDimensions="256 256 1" uvwOrigin="0 0 0" activeSlice="0" layoutGridRows="1" layoutGridColumns="1" sliceToRAS="0 0 1 -13.5923 -1 0 0 14.2311 0 1 0 60.5346 0 0 0 1" orientationMatrixAxial="-1 0 0 0 1 0 0 0 1" orientationMatrixSagittal="0 0 1 -1 0 0 0 1 0" orientationMatrixCoronal="-1 0 0 0 0 1 0 1 0" orientation="Sagittal" orientationReference="Sagittal" jumpMode="1" sliceVisibility="false" widgetVisibility="false" widgetOutlineVisibility="true" useLabelOutline="false" sliceSpacingMode="0" prescribedSliceSpacing="1 1 1" ></Slice>\r\n');
fprintf(fid,' <Slice\r\n  id="vtkMRMLSliceNodeGreen" name="Green" hideFromEditors="false" selectable="true" selected="false" singletonTag="Green" attributes="MappedInLayout:1" layoutLabel="G" layoutName="Green" active="false" visibility="true" backgroundColor="0 0 0" backgroundColor2="0 0 0" layoutColor="0.431373 0.690196 0.294118" orientationMarkerType="none" orientationMarkerSize="medium" rulerType="none" AxisLabels="L;R;P;A;I;S" fieldOfView="292.092 323.841 2.5" dimensions="460 510 1" xyzOrigin="0 0 0" sliceResolutionMode="1" uvwExtents="292.092 323.841 2.5" uvwDimensions="256 256 1" uvwOrigin="0 0 0" activeSlice="0" layoutGridRows="1" layoutGridColumns="1" sliceToRAS="-1 0 0 -14.6404 0 0 1 15.1317 0 1 0 60.5346 0 0 0 1" orientationMatrixAxial="-1 0 0 0 1 0 0 0 1" orientationMatrixSagittal="0 0 1 -1 0 0 0 1 0" orientationMatrixCoronal="-1 0 0 0 0 1 0 1 0" orientation="Coronal" orientationReference="Coronal" jumpMode="1" sliceVisibility="false" widgetVisibility="false" widgetOutlineVisibility="true" useLabelOutline="false" sliceSpacingMode="0" prescribedSliceSpacing="1 1 1" ></Slice>\r\n');
fprintf(fid,' <Layout\r\n  id="vtkMRMLLayoutNodevtkMRMLLayoutNode" name="Layout" hideFromEditors="true" selectable="true" selected="false" singletonTag="vtkMRMLLayoutNode" currentViewArrangement="2" guiPanelVisibility="1" bottomPanelVisibility ="1" guiPanelLR="0" collapseSliceControllers="0"\r\n numberOfCompareViewRows="1" numberOfCompareViewColumns="1" numberOfLightboxRows="6" numberOfLightboxColumns="6" mainPanelSize="400" secondaryPanelSize="400" ></Layout>\r\n');
fprintf(fid,' <SliceComposite\r\n  id="vtkMRMLSliceCompositeNodeRed" name="SliceComposite" hideFromEditors="true" selectable="true" selected="false" singletonTag="Red" backgroundVolumeID="vtkMRMLScalarVolumeNode1" foregroundVolumeID="vtkMRMLLabelMapVolumeNode2" labelVolumeID="vtkMRMLLabelMapVolumeNode1" compositing="0" foregroundOpacity="1" labelOpacity="1" linkedControl="0" fiducialVisibility="1" fiducialLabelVisibility="1" sliceIntersectionVisibility="0" layoutName="Red" annotationSpace="IJKAndRAS" annotationMode="All" doPropagateVolumeSelection="1" ></SliceComposite>\r\n');
fprintf(fid,' <SliceComposite\r\n  id="vtkMRMLSliceCompositeNodeYellow" name="SliceComposite_1" hideFromEditors="true" selectable="true" selected="false" singletonTag="Yellow" backgroundVolumeID="vtkMRMLScalarVolumeNode1" foregroundVolumeID="vtkMRMLLabelMapVolumeNode2" labelVolumeID="vtkMRMLLabelMapVolumeNode1" compositing="0" foregroundOpacity="1" labelOpacity="1" linkedControl="0" fiducialVisibility="1" fiducialLabelVisibility="1" sliceIntersectionVisibility="0" layoutName="Yellow" annotationSpace="IJKAndRAS" annotationMode="All" doPropagateVolumeSelection="1" ></SliceComposite>\r\n');
fprintf(fid,' <SliceComposite\r\n  id="vtkMRMLSliceCompositeNodeGreen" name="SliceComposite_2" hideFromEditors="true" selectable="true" selected="false" singletonTag="Green" backgroundVolumeID="vtkMRMLScalarVolumeNode1" foregroundVolumeID="vtkMRMLLabelMapVolumeNode2" labelVolumeID="vtkMRMLLabelMapVolumeNode1" compositing="0" foregroundOpacity="1" labelOpacity="1" linkedControl="0" fiducialVisibility="1" fiducialLabelVisibility="1" sliceIntersectionVisibility="0" layoutName="Green" annotationSpace="IJKAndRAS" annotationMode="All" doPropagateVolumeSelection="1" ></SliceComposite>\r\n');
fprintf(fid,' <Camera\r\n  id="vtkMRMLCameraNode1" name="Camera" hideFromEditors="false" selectable="true" selected="false" userTags="" position="0 500 0" focalPoint="0 0 0" viewUp="0 0 1" parallelProjection="false" parallelScale="1" viewAngle="30" activetag="vtkMRMLViewNode1" appliedTransform="1 0 0 0  0 1 0 0  0 0 1 0  0 0 0 1" ></Camera>\r\n');

% image data
fprintf(fid,' <SubjectHierarchy\r\n  id="vtkMRMLSubjectHierarchyNode1" name="SubjectHierarchy" hideFromEditors="false" selectable="true" selected="false" attributes="SubjectHierarchyVersion:2" >\r\n');
fprintf(fid,'   <SubjectHierarchyItem id="3" name="Scene" parent="0" type="" expanded="true" attributes="Level^Scene|">\r\n');
% labels (belong to id = 3)
for ilabel = 1 : nlabel
    fprintf(fid,'     <SubjectHierarchyItem id="%d" dataNode="vtkMRMLLabelMapVolumeNode%d" parent="3" type="LabelMaps" expanded="true"></SubjectHierarchyItem>\r\n',10+3*ilabel,ilabel);
end
% volumes (belong to id = 3)
for ivol = 1 : nvol
    if ivol ~= nvol
        fprintf(fid,'     <SubjectHierarchyItem id="%d" dataNode="vtkMRMLScalarVolumeNode%d" parent="3" type="Volumes" expanded="true"></SubjectHierarchyItem>\r\n',10+3*nlabel+3*ivol,ivol);
    else
        fprintf(fid,'     <SubjectHierarchyItem id="%d" dataNode="vtkMRMLScalarVolumeNode%d" parent="3" type="Volumes" expanded="true"></SubjectHierarchyItem>',10+3*nlabel+3*ivol,ivol);
    end
end
fprintf(fid,'</SubjectHierarchyItem></SubjectHierarchy>\r\n');
% end image data

% data state
fprintf(fid,' <ClipModels\r\n  id="vtkMRMLClipModelsNodevtkMRMLClipModelsNode" name="ClipModels" hideFromEditors="true" selectable="true" selected="false" singletonTag="vtkMRMLClipModelsNode" clipType="0" redSliceClipState="0" yellowSliceClipState="0" greenSliceClipState="0" ></ClipModels>\r\n');
fprintf(fid,' <ScriptedModule\r\n  id="vtkMRMLScriptedModuleNodeDataProbe" name="ScriptedModule" hideFromEditors="true" selectable="true" selected="false" singletonTag="DataProbe" ModuleName ="DataProbe" ></ScriptedModule>\r\n');

% details of the labels
for ilabel = 1 : nlabel
    fprintf(fid,' <LabelMapVolumeDisplay\r\n  id="vtkMRMLLabelMapVolumeDisplayNode%d" name="LabelMapVolumeDisplay" hideFromEditors="true" selectable="true" selected="false" color="0.5 0.5 0.5" edgeColor="0 0 0" selectedColor="1 0 0" selectedAmbient="0.4" ambient="0" diffuse="1" selectedSpecular="0.5" specular="0" power="1" opacity="1" sliceIntersectionOpacity="1" pointSize="1" lineWidth="1" representation="2" lighting="true" interpolation="1" shading="true" visibility="true" edgeVisibility="false" clipping="false" sliceIntersectionVisibility="false" sliceIntersectionThickness="3" frontfaceCulling="false" backfaceCulling="true" scalarVisibility="false" vectorVisibility="false" tensorVisibility="false" interpolateTexture="false" scalarRangeFlag="UseData" scalarRange="0 100" colorNodeID="vtkMRMLColorTableNodeFileGenericColors.txt"  ></LabelMapVolumeDisplay>\r\n',ilabel);
    fprintf(fid,' <VolumeArchetypeStorage\r\n  id="vtkMRMLVolumeArchetypeStorageNode%d" name="VolumeArchetypeStorage_2" hideFromEditors="true" selectable="true" selected="false" fileName="%s" useCompression="1" defaultWriteFileExtension="nrrd" readState="0" writeState="0" centerImage="0" UseOrientationFromFile="1" ></VolumeArchetypeStorage>\r\n',ilabel+2,strrep(labelfnames{ilabel},' ','%20'));
    fprintf(fid,' <LabelMapVolume\r\n  id="vtkMRMLLabelMapVolumeNode%d" name="%s" hideFromEditors="false" selectable="true" selected="false" displayNodeRef="vtkMRMLLabelMapVolumeDisplayNode%d" storageNodeRef="vtkMRMLVolumeArchetypeStorageNode%d" references="display:vtkMRMLLabelMapVolumeDisplayNode%d;storage:vtkMRMLVolumeArchetypeStorageNode%d;" userTags="" ijkToRASDirections="%.6f   %.6f   %.6f   %.6f   %.6f   %.6f   %.6f   %.6f   %.6f " spacing="%f %f %f" origin="%.3f %.3f %.3f" ></LabelMapVolume>\r\n',...
              ilabel,lpartnames{ilabel},ilabel,ilabel+2,ilabel,ilabel+2,-vec(1,1:3),-vec(2,1:3),vec(3,1:3),vox(1:3),-vec(1,4),-vec(2,4),vec(3,4));
end
for ivol = 1 : nvol
    fprintf(fid,' <VolumeDisplay\r\n  id="vtkMRMLScalarVolumeDisplayNode%d" name="VolumeDisplay" hideFromEditors="true" selectable="true" selected="false" color="0.5 0.5 0.5" edgeColor="0 0 0" selectedColor="1 0 0" selectedAmbient="0.4" ambient="0" diffuse="1" selectedSpecular="0.5" specular="0" power="1" opacity="1" sliceIntersectionOpacity="1" pointSize="1" lineWidth="1" representation="2" lighting="true" interpolation="1" shading="true" visibility="true" edgeVisibility="false" clipping="false" sliceIntersectionVisibility="false" sliceIntersectionThickness="1" frontfaceCulling="false" backfaceCulling="true" scalarVisibility="false" vectorVisibility="false" tensorVisibility="false" interpolateTexture="false" scalarRangeFlag="UseData" scalarRange="0 100" colorNodeID="vtkMRMLColorTableNodeGrey"  window="5396" level="2698" upperThreshold="32767" lowerThreshold="-32768" interpolate="1" windowLevelLocked="false" autoWindowLevel="1" applyThreshold="0" autoThreshold="0" ></VolumeDisplay>\r\n',ivol);
    fprintf(fid,' <VolumeArchetypeStorage\r\n  id="vtkMRMLVolumeArchetypeStorageNode%d" name="VolumeArchetypeStorage_3" hideFromEditors="true" selectable="true" selected="false" fileName="%s" useCompression="1" defaultWriteFileExtension="nrrd" readState="0" writeState="0" centerImage="0" UseOrientationFromFile="1" ></VolumeArchetypeStorage>\r\n',nlabel+2+ivol,strrep(volfnames{ivol},' ','%20'));
    fprintf(fid,' <Volume\r\n  id="vtkMRMLScalarVolumeNode%d" name="%s" hideFromEditors="false" selectable="true" selected="false" displayNodeRef="vtkMRMLScalarVolumeDisplayNode%d" storageNodeRef="vtkMRMLVolumeArchetypeStorageNode%d" references="display:vtkMRMLScalarVolumeDisplayNode%d;storage:vtkMRMLVolumeArchetypeStorageNode%d;" userTags="" ijkToRASDirections="%.6f   %.6f   %.6f   %.6f   %.6f   %.6f   %.6f   %.6f   %.6f " spacing="%f %f %f" origin="%.3f %.3f %.3f" ></Volume>\r\n',...
              ivol,vpartnames{ivol},ivol,nlabel+2+ivol,ivol,nlabel+2+ivol,-vec(1,1:3),-vec(2,1:3),vec(3,1:3),vox(1:3),-vec(1,4),-vec(2,4),vec(3,4));
end
fprintf(fid,'</MRML>\r\n');

fclose(fid);
end