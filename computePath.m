function pathLength = computePath

[file2D,path2D]   = uigetfile('*.dcm');
[fileCen,pathCen] = uigetfile('*.dat');
info = dicominfo(fullfile(path2D,file2D));
cenAlltab = readtable(fullfile(pathCen,fileCen));

% sort due to the data structure of slicer output
cenAlltab = cenAlltab(end:-1:1,:);
% conversion to a matrix
cenAll = table2array(cenAlltab);
cenAll(:,1:2) = -cenAll(:,1:2);
% conversion to a cell array
cenAllCell = table2cell(cenAlltab);
% conversion of each cell components to string to search branch centerlines
tmpmat_str = cellfun(@num2str,cenAllCell,'UniformOutput',false);

tmpmat_branchstartInd = cellfun(@(x) strcmp(x,tmpmat_str{1,3}),tmpmat_str,'UniformOutput', false );
[row_startInd,~] = find(cell2mat(tmpmat_branchstartInd));

% centerline for the main aorta.
numPoints = diff([row_startInd;size(cenAll,1)]);
[~,maxPointsInd] = max(numPoints);
cen = cenAll(row_startInd(maxPointsInd):row_startInd(maxPointsInd)+numPoints(maxPointsInd)-1,:);
        
fnplt(cscvn(cen'),'g',2);
hold on;

% compute determinant of three vectors (imageOrient[1],imageOrient[2],
% vector from "imagePosition" to Centerline)) to identify closest points
% from the 2D-PC.
orgn(:,1) = info.ImagePositionPatient;
vecs    = info.ImageOrientationPatient;
% Check thickness
thickness = info.SliceThickness;
disp(['Slice thickness = ',num2str(thickness),'mm.']);
n = cross(vecs(1:3),vecs(4:6))/norm(cross(vecs(1:3),vecs(4:6)));
% Taking into account the slice thickness (origin(1): baseline, origin(2):upper limit, origin(3): upper limit)
orgn(:,2) = orgn(:,1) + n*thickness/2;
orgn(:,3) = orgn(:,1) - n*thickness/2;

vecs = reshape(vecs,[3,2]);

for in = 1 : 3
    dets = zeros(size(cen,1),1);
    for i = 1 : size(cen,1)
        dets(i,1) = det([cen(i,:)'-orgn(:,in) vecs]);
    end
    
    %Determinant is close to zero as approaching to the plane. Check the sign of determinants.
    %The two points where the sign changes should be closest points to the pldane.
    detsCheck = dets(1:end-1).*dets(2:end);
    detsCheck(detsCheck>0) = 0;
    interPoint = find(detsCheck);
    cen_path = cen(interPoint(1):interPoint(2)+1,1:3);
    
    %Compute the intersection points (Linear interpolation)
    t = ((orgn(1,in)-cen_path(1,1))*n(1)+(orgn(2,in)-cen_path(1,2))*n(2)+(orgn(3,in)-cen_path(1,3))*n(3))/...
        ((cen_path(2,1)-cen_path(1,1))*n(1)+(cen_path(2,2)-cen_path(1,2))*n(2)+(cen_path(2,3)-cen_path(1,3))*n(3));
    
    p1 = cen_path(1,:)+t*(cen_path(2,:)-cen_path(1,:));
    
    t = ((orgn(1,in)-cen_path(end-1,1))*n(1)+(orgn(2,in)-cen_path(end-1,2))*n(2)+(orgn(3,in)-cen_path(end-1,3))*n(3))/...
        ((cen_path(end,1)-cen_path(end-1,1))*n(1)+(cen_path(end,2)-cen_path(end-1,2))*n(2)+(cen_path(end,3)-cen_path(end-1,3))*n(3));
    
    p2 = cen_path(end-1,:)+t*(cen_path(end,:)-cen_path(end-1,:));
    
    % replace the first and last point of the path
    cen_path(1  ,:) = p1;
    cen_path(end,:) = p2;
    
    if in == 1
        fnplt(cscvn(cen_path'),'r',2);
        axis equal;
    end
    % Calculate length
    %Diff = diff(cen_path,1);
    %pathLength = sum(sum(Diff.*Diff,2),1);
    fn = cscvn(cen_path');
    fnprime = fnder(fn);
    Lfun = @(s) sqrt(sum(fnval(fnprime,s).^2,1));
    pathLength(in) = integral(Lfun,fn.breaks(1),fn.breaks(end));
    
    disp(['Intersection1:',num2str(p1(1),'%10.4f'),',',num2str(p1(2),'%10.4f'),',',num2str(p1(3),'%10.4f')]);
    disp(['Intersection2:',num2str(p2(1),'%10.4f'),',',num2str(p2(2),'%10.4f'),',',num2str(p2(3),'%10.4f')]);
end
end

