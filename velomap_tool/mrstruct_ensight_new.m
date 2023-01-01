function [res, errStr, oArg]= mrstruct_ensight_new(varargin)
%function [res, errStr, oArg]= mrstruct_ensight(mrStruct, comStr, [, op1[, op2[... [, opN]]]])
%
%   command:    {'geoFile' | 'dataVolume' | 'dataVector' | '' | 
%
%   'geoFile': fName, partNo, transf_My, caFlag
%   'dataVolume': fName, partNo, commentStr
%   'dataVector': fName, partNo, commentStr, transMy, NormFlag
%
%           transf_My: 4x4 Matrix; M*voxIndex -> world system; Default: edge entry of the mrStruct
%           caFlag: 'create' | 'append'; On 'create' a new geo file will be created
%                                       On 'append' new geofile will be appended to existing geofile
%           partNo: Id, to link vector or volume datasets to a geometriefile entry
%           normFlag: 'none' | 'L1' | 'L2'; Norm methode for vector data
%
% Bjoern W. Kreher
% 08/05
%
% UNIX


res= []; errStr= ''; oArg= [];
maxArg= 10;

% Hauptparameter kontrolle
if (nargin >= 1) & mrstruct_istype(varargin{1})
    mrStruct= varargin{1};
else
    errStr= strcat(mfilename, ' (error): First argument have to be of the type mrStruct');
    return
end

if (nargin >= 2) & isstr(varargin{2})
    comStr= varargin{2};
else
    errStr= strcat(mfilename, ' (error): Second argument have to be of the type string');
    return
end

% Nebenparameter übergabe für beseres handling
argCell= cell(1, maxArg);
for i= 1:(nargin - 2)
    argCell{i}= varargin{2 + i};
end

% Parsen des commando Strings
if strcmp(comStr, 'geoFile')            %'geoFile'
    if isempty(argCell{1})  % fName
        [fNameStr, pathStr]= uiputfile;
        if isstr(fNameStr)
            fNameStr= fullfile(pathStr, fNamestr);
        else
            errStr= strcat(mfilename, ' (message): Aborted by user');
            return
        end
    else
        fNameStr= argCell{1};
    end
    if isempty(argCell{2})      % partNo
        partNo= 1;
    elseif isnumeric(argCell{2}) & (prod(size(argCell{2})) == 1) & (round(argCell{2}) == argCell{2})
        partNo= argCell{2};
    else
        errStr= strcat(mfilename, ' (error): partNo should be empty or a natural number');
        return
    end
    if isequal(size(argCell{3}), [4 4])  % transMy
        transMy= argCell{3};
    elseif isstr(argCell{3})
        [transMy, errStr]= private_getTransMy(mrStruct, argCell{3});
    elseif isempty(argCell{3})
        [transMy, errStr]= private_getTransMy(mrStruct);
    else
        errStr= strcat(mfilename, ' (error): transMy should be empty or a 4x4 matrix');
    end
    if isempty(transMy)
        return
    end
    if isempty(argCell{4})      % appFlag: {'create' | 'append'}
        appFlag= 'create';
    elseif isstr(argCell{4}) & (strcmp(argCell{4}, 'create') | strcmp(argCell{4}, 'append'))
        appFlag= argCell{4};
    else
        errStr= strcat(mfilename, ' (error): appFlag should be either ''create'' or ''append''');
        return
    end
    
    [res, errStr]= local_createGeoFile(mrStruct, fNameStr, partNo, transMy, appFlag);
    
elseif strcmp(comStr, 'dataVolume') %     'dataFile'
    if isempty(argCell{1})          % fName
        [fNameStr, pathStr]= uiputfile;
        if isstr(fNameStr)
            fNameStr= fullfile(pathStr, fNamestr);
        else
            errStr= strcat(mfilename, ' (message): Aborted by user');
            return
        end
    else
        fNameStr= argCell{1};
    end
    if isempty(argCell{2})          % partNo
        partNo= 1;
    elseif isnumeric(argCell{2}) & (prod(size(argCell{2})) == 1) & (round(argCell{2}) == argCell{2})
        partNo= argCell{2};
    else
        errStr= strcat(mfilename, ' (error): partNo should be empty or a natural number');
        return
    end
    if isempty(argCell{3}) | isstr(argCell{3}) %commentStr
        commentStr= argCell{3};
    else
        errStr= strcat(mfilename, ' (error): commentStr should be empty or a string');
        return        
    end
    [res, errStr]= local_createMagDataFile(mrStruct, fNameStr, partNo, commentStr);   
elseif strcmp(comStr, 'dataVector') %     'dataFile'
    if isempty(argCell{1})          % fName
        [fNameStr, pathStr]= uiputfile;
        if isstr(fNameStr)
            fNameStr= fullfile(pathStr, fNamestr);
        else
            errStr= strcat(mfilename, ' (message): Aborted by user');
            return
        end
    else
        fNameStr= argCell{1};
    end
    if isempty(argCell{2})          % partNo
        partNo= 1;
    elseif isnumeric(argCell{2}) & (prod(size(argCell{2})) == 1) & (round(argCell{2}) == argCell{2})
        partNo= argCell{2};
    else
        errStr= strcat(mfilename, ' (error): partNo should be empty or a natural number');
        return
    end
    if isempty(argCell{3}) | isstr(argCell{3}) % CommentStr
        commentStr= argCell{3};
    else
        errStr= strcat(mfilename, ' (error): commentStr should be empty or a string');
        return        
    end
    if isequal(size(argCell{4}), [4 4])  % transMy
        transMy= argCell{4};
    elseif isstr(argCell{4})
        [transMy, errStr]= private_getTransMy(mrStruct, argCell{4});
    elseif isempty(argCell{4})
        [transMy, errStr]= private_getTransMy(mrStruct);
    else
        errStr= strcat(mfilename, ' (error): transMy should be empty or a 4x4 matrix');
    end
    if isempty(transMy)
        return
    end
    if isempty(argCell{5}) | isstr(argCell{5}) % normFlag
        normFlag= argCell{5};
    else
        errStr= strcat(mfilename, ' (error): normFlag should be a string like ''none'', ''L0'', ''L1'', ''L2'', etc');
        return        
    end
    [res, errStr]= local_createVectorDataFile(mrStruct, fNameStr, partNo, transMy, normFlag, commentStr);   
    
else
    errStr= strcat(mfilename, ' (error): Command ''', comStr, ''' is not supported yet');
    return    
end






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%
%
%  START:
%       [res, errStr]= local_createGeoFile(mrStruct, transMy, fullNameStr)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [res, errStr]= local_createGeoFile(mrStruct, fullNameStr, partNo, transMy, appFlag)

res= [];    errStr= '';

[pathStr, fNameStr, extStr]= fileparts(fullNameStr);
geoPathName= fullfile(pathStr, strcat(fNameStr, '.geo'));

sizeAy= mrstruct_query(mrStruct, 'sizeAy');
if prod(sizeAy) == 0
    errStr= strcat(mfilename, '::local_createGeoFile (error): mrStruct should not be empty');
    return;
end
if length(sizeAy) < 3
    sizeAy(end:3)= 1;
end

%[pX, pY, pZ]= meshgrid(1:sizeAy(1), 1:sizeAy(2), 1:sizeAy(3));
%posVc= transMy*[reshape(pX, [1, voxNo]); reshape(pY, [1, voxNo]); reshape(pZ, [1, voxNo]); ones(1, voxNo)];


%posVc= transMy*(double([pX; pY; pZ; int16(ones(1, voxNo))]));

%%% write geoFile
if strcmp(appFlag, 'append') & ~exist(geoPathName, 'file')
    errStr= strcat(mfilename, '::local_createGeoFile (error): Cannot append geo data; file does not exist');
    return
elseif strcmp(appFlag, 'create')
    fidGeo  = fopen(geoPathName, 'w+', 'ieee-be');
    %matSize = szx * szy;
    
    % generate and write header
    geoheaderStr(1:8) = 'C Binary';
    geoheaderStr(9:80) = ' ';
    fwrite(fidGeo,geoheaderStr,'char');
    geoheaderStr(1:21) = 'Ensight Geometry File';
    geoheaderStr(22:80) = ' ';
    fwrite(fidGeo,geoheaderStr,'char');
    geoheaderStr(1:44) = 'Created by MATLAB routine, (c) M. Markl 2005';
    geoheaderStr(45:80) = ' ';
    fwrite(fidGeo,geoheaderStr,'char');
    geoheaderStr(1:14) = 'node id assign';
    geoheaderStr(15:80) = ' ';
    fwrite(fidGeo,geoheaderStr,'char');
    geoheaderStr(1:17) = 'element id assign';
    geoheaderStr(18:80) = ' ';
    fwrite(fidGeo,geoheaderStr,'char');    
else
    fidGeo  = fopen(geoPathName, 'a+', 'ieee-be');
    fseek(fidGeo, 0, 'eof');
end

geoheaderStr(1:4) = 'part';
geoheaderStr(5:80) = ' ';
fwrite(fidGeo,geoheaderStr,'char');

fwrite(fidGeo, partNo, 'int');              %

geoheaderStr(1:12) = 'Total Volume';
geoheaderStr(13:80) = ' ';
fwrite(fidGeo,geoheaderStr,'char');
geoheaderStr(1:5) = 'block';
geoheaderStr(6:80) = ' ';
fwrite(fidGeo,geoheaderStr,'char');
%clear('mrStruct')
%%%
%%% create and transform datapoints
voxNo= prod(sizeAy(1:3));


% % % % % pX= reshape((1:sizeAy(1))'*ones(1, voxNo/sizeAy(1)), [1 voxNo]);
% % % % % pY= reshape(reshape(ones(sizeAy(1), 1)*(1:sizeAy(2)), [voxNo/sizeAy(3) 1])*ones(1, sizeAy(3)), [1, voxNo]);
% % % % % pZ= reshape(ones(voxNo/sizeAy(3), 1)*(1:sizeAy(3)), [1 voxNo]);

%%%try to allocate memory
try
    posVc = ones(1,voxNo);
catch
    errordlg('An error is occured while trying to allocate memory!')
    return
end


%%% write data to geo file
fwrite(fidGeo,sizeAy(1),'int');         %
fwrite(fidGeo,sizeAy(2),'int');         %
fwrite(fidGeo,sizeAy(3),'int');         %

for n=1:3
    for m=1:sizeAy(3)
        %% create matrix
        pX= reshape((1:sizeAy(1))'*ones(1, sizeAy(2)),[1 voxNo/sizeAy(3)]);
        pY= reshape(reshape(ones(sizeAy(1), 1)*(1:sizeAy(2)), [voxNo/sizeAy(3) 1]), [1, voxNo/sizeAy(3)]);
        pZ= reshape(ones(voxNo/sizeAy(3), 1)*(m), [1 voxNo/sizeAy(3)]);
        startPos=(m-1)*voxNo/sizeAy(3)+1;
        endPos  = m*voxNo/sizeAy(3);
        posVc(:,startPos:endPos)= squeeze(transMy(n,:))*([pX; pY; pZ; (ones(1, voxNo/sizeAy(3)))]);
    end
    fwrite(fidGeo, posVc, 'float');
end
% posVc= squeeze(transMy(1,:))*([pX; pY; pZ; (ones(1, voxNo))]);
%  % x-coordinate
% clear('posVc');
% posVc= squeeze(transMy(2,:))*([pX; pY; pZ; (ones(1, voxNo))]);
% fwrite(fidGeo, posVc, 'float'); % y-coordinate
% clear('posVc');
% posVc= squeeze(transMy(3,:))*([pX; pY; pZ; (ones(1, voxNo))]);
% fwrite(fidGeo, posVc, 'float'); % z-coordinate
% clear('posVc');

% % original version von Björn
% fwrite(fidGeo, posVc(1, :), 'float'); % x-coordinate
% fwrite(fidGeo, posVc(2, :), 'float'); % y-coordinate
% fwrite(fidGeo, posVc(3, :), 'float'); % z-coordinate
% fclose(fidGeo); % close file

res= geoPathName;

% end generate geo file               


%
%
%  START:
%       [res, errStr]= local_createDataFile(mrStruct, dataNameStr, transMy, fullNameStr)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [res, errStr]= local_createMagDataFile(mrStruct, fullNameStr, partNo, commentStr)
res= [];    errStr= '';

sizeAy= size(squeeze(mrStruct.dataAy));
if prod(sizeAy) == 0
    errStr= strcat(mfilename, '::local_createMagDataFile (error): mrstruct doesen''t contain data');
    return
end
if (length(sizeAy) > 3)
    errStr= strcat(mfilename, '::local_createMagDataFile (error): mrstruct doesen''t describe a 3D magnitude data');
    return
end
sizeVol= ones(1, 3);
sizeVol(1:length(sizeAy))= sizeAy(1:end);

% Tranform data
dataVc= mrStruct.dataAy(:);
if isempty(commentStr)
    commentStr= sprintf('Magnitudedata [%dx%dx%d]', sizeVol(1), sizeVol(2), sizeVol(3));
end

[res, errStr]= private_write_data(fullNameStr, partNo, dataVc, commentStr);


%
%
%  START:
%       [res, errStr]= local_createVectorDataFile(mrStruct, fullNameStr, partNo, transMy, normFlag, commentStr);   
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [res, errStr]= local_createVectorDataFile(mrStruct, fullNameStr, partNo, transMy, normFlag, commentStr)
res= [];    errStr= '';

sizeAy= size(squeeze(mrStruct.dataAy));
if prod(sizeAy) == 0
    errStr= strcat(mfilename, '::local_createVectorDataFile (error): mrstruct doesen''t contain data');
    return
end
if (length(sizeAy) > 4) | (sizeAy(end) ~= 3)
    errStr= strcat(mfilename, '::local_createVectorDataFile (error): mrstruct doesen''t describe a 3D vectorfield');
    return
end
sizeVol= ones(1, 3);
idx= find(sizeAy == 3); 
sizeVol(1:(idx-1))= sizeAy(1:(idx-1));
voxNo= prod(sizeVol);

% Tranform data
dataVc= reshape(mrStruct.dataAy, [voxNo 3]);
% Only rotation and scaling)
P= private_spm_imatrix(transMy);
if ~isequal([P(4:6) sign(P(7:9))], [0 0 0 1 1 1])
    rotMy= private_spm_matrix([0 0 0 P(4:6) sign(P(7:9)) 0 0 0]); %%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%SCALING?????????
    dataVc= (rotMy(1:3, 1:3)*(dataVc'))';
end
if ~isempty(normFlag) & ~strcmp(normFlag, 'none')
    idx= find(sum(abs(dataVc), 2) > 0);
    if strcmp(normFlag, 'L2')
        dataVc(idx, :)= dataVc(idx, :)./(sqrt(sum(dataVc(idx, :).^2, 2))*ones(1, 3));
    elseif strcmp(normFlag, 'L1')
        dataVc(idx, :)= dataVc(idx, :)./(sum(abs(dataVc(idx, :)), 2)*ones(1, 3));
    else
        errStr= strcat(mfilename, '::local_createDataFile (error): unsupported norm');
        return
    end
else
    normFlag= 'none';
end

if isempty(commentStr)
    commentStr= sprintf('Vectordata [%dx%dx%d] normalized by %s', sizeVol(1), sizeVol(2), sizeVol(3), normFlag);
end
[res, errStr]= private_write_data(fullNameStr, partNo, dataVc, commentStr);


%
%
%  START:
%       [res, errStr]= private_write_data(fName, partNo, timeId, dataVc)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [res, errStr]= private_write_data(fullFName, partNo, dataVc, commentStr)

res= []; errStr= '';

% write header and data
fid= fopen(fullFName, 'w+', 'ieee-be');
if fid < 0
    errStr= strcat(mfilename, '::private_write_data (error): cant''t openfile (', lasterr, ')');
    return
end

% write data header
if length(commentStr) >= 80
    commentStr= commentStr(1:79);
end
dataheaderStr(1:length(commentStr)) = commentStr;
dataheaderStr(length(commentStr)+1:80) = ' ';
fwrite(fid,dataheaderStr, 'char');
dataheaderStr(1:4) = 'part';
dataheaderStr(5:80) = ' ';
fwrite(fid, dataheaderStr, 'char');
fwrite(fid, partNo, 'int');     

dataheaderStr(1:5) = 'block';
dataheaderStr(6:80) = ' ';
fwrite(fid, dataheaderStr, 'char');
% write data body
fwrite(fid, dataVc, 'float');
% close file and return
fclose(fid);
res= fullFName;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%
%
%  START:
%       [transMy, errStr]= private_getTransMy(mrStruct)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [transMy, errStr]= private_getTransMy(mrStruct, scaleOnly)

transMy= [];    errStr= '';

if ~exist('scaleOnly') | isempty(scaleOnly)
    scaleOnly= 'edge';
end

if ~mrstruct_istype(mrStruct)
    errStr= strcat(mfilename, '::private_getTransMy (error): mrStruct should be a mrStruct');
    return;
end

if isequal(size(mrStruct.edges), [4 4]) & ~strcmp(scaleOnly, 'vox')% wenn edge eintrag vorhanden ist alles ok....
    transMy= mrStruct.edges;
    return
else
    sizeAy= mrstruct_query(mrStruct, 'sizeAy');
    if length(sizeAy) < 3
        errStr= strcat(mfilename, '::private_getTransMy (error): mrStruct have to be 3dim at least');
        return;
    end
    
    vox= mrstruct_query(mrStruct, 'vox');
    if length(vox) == 4
        vox= [vox(1) vox(2) vox(3) + vox(4)];
    elseif length(vox) ~= 3
        errStr= strcat(mfilename, '::private_getTransMy (error): mrStruct contain a vox size');
        return;
    end
    
    transMy= diag([vox 1]);
    
    if strcmp(mrStruct.orient, 'axial')
        %do nothing
    elseif strcmp(mrStruct.orient, 'saggital')
        %do nothing
    elseif strcmp(mrStruct.orient, 'coronar')
        transMy(3, 4)= transMy(3, 3)*(sizeAy(3) - 1);
        transMy(3, 3)= -transMy(3, 3);
    else        
        %do nothing
    end
end


function P = private_spm_imatrix(M)
% returns the parameters for creating an affine transformation
% FORMAT P = spm_imatrix(M)
% M      - Affine transformation matrix
% P      - Parameters (see spm_matrix for definitions)
%___________________________________________________________________________
% @(#)spm_imatrix.m	2.1 John Ashburner & Stefan Kiebel 98/12/18

% Translations and zooms
%-----------------------------------------------------------------------
R         = M(1:3,1:3);
C         = chol(R'*R);
P         = [M(1:3,4)' 0 0 0  diag(C)'  0 0 0];
if det(R)<0, P(7)=-P(7);end % Fix for -ve determinants

% Shears
%-----------------------------------------------------------------------
C         = diag(diag(C))\C;
P(10:12)  = C([4 7 8]);
R0        = private_spm_matrix([0 0 0  0 0 0 P(7:12)]);
R0        = R0(1:3,1:3);
R1        = R/R0;

% This just leaves rotations in matrix R1
%-----------------------------------------------------------------------
%[          c5*c6,           c5*s6, s5]
%[-s4*s5*c6-c4*s6, -s4*s5*s6+c4*c6, s4*c5]
%[-c4*s5*c6+s4*s6, -c4*s5*s6-s4*c6, c4*c5]

P(5) = asin(rang(R1(1,3)));
if (abs(P(5))-pi/2).^2 < 1e-9,
	P(4) = 0;
	P(6) = atan2(-rang(R1(2,1)), rang(-R1(3,1)/R1(1,3)));
else,
	c    = cos(P(5));
	P(4) = atan2(rang(R1(2,3)/c), rang(R1(3,3)/c));
	P(6) = atan2(rang(R1(1,2)/c), rang(R1(1,1)/c));
end;
return;

% There may be slight rounding errors making b>1 or b<-1.
function a = rang(b)
a = min(max(b, -1), 1);
return;


function [A] = private_spm_matrix(P)
% returns an affine transformation matrix
% FORMAT [A] = spm_matrix(P)
% P(1)  - x translation
% P(2)  - y translation
% P(3)  - z translation
% P(4)  - x rotation about - {pitch} (radians)
% P(5)  - y rotation about - {roll}  (radians)
% P(6)  - z rotation about - {yaw}   (radians)
% P(7)  - x scaling
% P(8)  - y scaling
% P(9)  - z scaling
% P(10) - x affine
% P(11) - y affine
% P(12) - z affine
%
% A     - affine transformation matrix
%___________________________________________________________________________
%
% spm_matrix returns a matrix defining an orthogonal linear (translation,
% rotation, scaling or affine) transformation given a vector of
% parameters (P).  The transformations are applied in the following order
% (i.e., the opposite to which they are specified):
%
% 1) shear
% 2) scale
% 3) rotation - yaw, roll & pitch
% 4) translation
%
% SPM uses a PRE-multiplication format i.e. Y = A*X where X and Y are 4 x n
% matrices of n coordinates.
%
%__________________________________________________________________________
% @(#)spm_matrix.m	2.1 02/08/14

% pad P with 'null' parameters
%---------------------------------------------------------------------------
q  = [0 0 0 0 0 0 1 1 1 0 0 0];
P  = [P q((length(P) + 1):12)];

T  =   [1 	0 	0 	P(1);
        0 	1 	0 	P(2);
        0 	0 	1 	P(3);
        0 	0 	0 	1];

R1  =  [1    0   	0   	   0;
        0    cos(P(4))  sin(P(4))  0;
        0   -sin(P(4))  cos(P(4))  0;
        0    0    	0   	   1];

R2  =  [cos(P(5))  0   	sin(P(5))  0;
        0    	   1    0  	   0;
       -sin(P(5))  0  	cos(P(5))  0;
        0          0    0   	   1];

R3  =  [cos(P(6))   sin(P(6))   0  0;
       -sin(P(6))   cos(P(6))   0  0;
        0           0           1  0;
        0     	    0    	0  1];

Z   =  [P(7) 	0   	0    	0;
        0    	P(8) 	0    	0;
        0    	0    	P(9) 	0;
        0    	0    	0    	1];

S   =  [1   	P(10)   P(11)   0;
        0   	1 	P(12)   0;
        0   	0   	1	0;
        0    	0    	0    	1];

A = T*R1*R2*R3*Z*S;



% %% -------------------------------------------------
%  %% Init EnSight file conversion
%  %% -------------------------------------------------    
%  if enSightConvFlag
%       
%      if nFE == 3
%          
%         % update text output
%         subplot('position',[0.15 0.05 0.8 0.05])
%         outStr = ['Dicom --> EnSight, generating case & geo files ... '];
%         cla;
%         text(0.02,0.5,outStr,'FontSize',8,'Interpreter','none');
%         axis off;
%         colormap(gray);
%         drawnow;
%         
%         % create new directory for EnSight & matlab files
%         ensightDirPath  = sprintf('%s%s%s%s',dirStrData,'\','EnSight_',patientName);
%         [s,mess,messid] = mkdir(ensightDirPath);
%                  
%         % generate case file 
%         % ------------------------------------------------
%         casePathName  = sprintf('%s%s%s%s%s',ensightDirPath,'\','EnSight_',patientName,'.case');
%         geoPathName   = sprintf('%s%s%s%s%s',ensightDirPath,'\','EnSight_',patientName,'.geo');
%         dataPathName  = sprintf('%s%s%s%s%s',ensightDirPath,'\','EnSight_',patientName,'_');
%         
%         geoFileName   = sprintf('EnSight_%s%s%',patientName,'.geo');
%         dataFileName  = sprintf('EnSight_%s%s%',patientName,'_');
% 
%          RRinterval= 1; numPhases= 120;
%          timeStamps   = [1:numPhases]*RRinterval/numPhases - RRinterval/numPhases/2;
% %         
%         % open text file an write data
%         fidCase  = fopen(casePathName, 'wt');
%         fprintf(fidCase,'FORMAT\n');
%         fprintf(fidCase,'type:	 ensight gold\n');
%         fprintf(fidCase,'GEOMETRY\n');
%         fprintf(fidCase,'model:	 %s\n',geoFileName);
%         fprintf(fidCase,'VARIABLE\n');
%         fprintf(fidCase,'scalar per node:	 Magnitude	 %s**.mag\n',dataFileName );
%         if pcmraFlag % add speed data if selected by user
%             fprintf(fidCase,'scalar per node:	 Speed	     %s**.spd\n',dataFileName );
%         end
%         fprintf(fidCase,'vector per node:	 Velocity	 %s**.vel\n',dataFileName );
%         fprintf(fidCase,'TIME\n');
%         fprintf(fidCase,'time set:		 1\n');
%         fprintf(fidCase,'number of steps:	 %s\n',num2str(numPhases));
%         fprintf(fidCase,'filename start number:	 0\n');
%         fprintf(fidCase,'filename increment:	 1\n');
%         fprintf(fidCase,'time values:\n');
%         fprintf(fidCase,'%.3f\n',timeStamps);
%         fclose(fidCase); % close file
%         % end generate case file
%                 
%         
%         % generate geo file
%         % ------------------------------------------------
%         % open binary file an write data
%         fidGeo  = fopen(geoPathName, 'w', 'ieee-be');
%         matSize = szx * szy;
%         part    = 1;
%         
%         % generate and write header
%         geoheaderStr(1:8) = 'C Binary';
%         geoheaderStr(9:80) = ' ';
%         fwrite(fidGeo,geoheaderStr,'char');
%         geoheaderStr(1:21) = 'Ensight Geometry File';
%         geoheaderStr(22:80) = ' ';
%         fwrite(fidGeo,geoheaderStr,'char');
%         geoheaderStr(1:44) = 'Created by MATLAB routine, (c) M. Markl 2005';
%         geoheaderStr(45:80) = ' ';
%         fwrite(fidGeo,geoheaderStr,'char');
%         geoheaderStr(1:14) = 'node id assign';
%         geoheaderStr(15:80) = ' ';
%         fwrite(fidGeo,geoheaderStr,'char');
%         geoheaderStr(1:17) = 'element id assign';
%         geoheaderStr(18:80) = ' ';
%         fwrite(fidGeo,geoheaderStr,'char');
%         geoheaderStr(1:4) = 'part';
%         geoheaderStr(5:80) = ' ';
%         fwrite(fidGeo,geoheaderStr,'char');
%         
%         fwrite(fidGeo,part,'int');         
%         
%         geoheaderStr(1:12) = 'Total Volume';
%         geoheaderStr(13:80) = ' ';
%         fwrite(fidGeo,geoheaderStr,'char');
%         geoheaderStr(1:5) = 'block';
%         geoheaderStr(6:80) = ' ';
%         fwrite(fidGeo,geoheaderStr,'char');
%         
%         fwrite(fidGeo,szx,'int');
%         fwrite(fidGeo,szy,'int');
%         fwrite(fidGeo,numSlices,'int');
%         
%         % needed for coordiante definition
%         for l=1:szy
%             %Xtmp(1:szx,l) = [1:szx];
%             Xtmp(l,1:szx) = [1:szx];
%         end; 
%         
%         for l=1:szx
%             %Ytmp(l,1:szy) = fliplr(1:szy);
%             Ytmp(1:szx,l) = fliplr(1:szy);
%         end; 
%          
%         for k = 1:numSlices
%             Xcor((k-1)*matSize+1:k*matSize) = Xtmp(:) * dx;
%             Ycor((k-1)*matSize+1:k*matSize) = Ytmp(:) * dy;
%             if orientation == 'cor'
%                 Zcor((k-1)*matSize+1:k*matSize) = (numSlices - k + 1) * dz;
%             elseif orientation == 'sag'
%                 Zcor((k-1)*matSize+1:k*matSize) = k * dz;
%             elseif orientation == 'tra'
%                 Zcor((k-1)*matSize+1:k*matSize) = k * dz;                        
%             end
%         end
%         
%         fwrite(fidGeo,Xcor,'float'); % x-coordinate
%         fwrite(fidGeo,Ycor,'float'); % y-coordinate
%         fwrite(fidGeo,Zcor,'float'); % z-coordinate
%         fclose(fidGeo); % close file
%         % end generate geo file               
%         
%      end
%  end
% 
%  
%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%  
%  
%      %% -------------------------------------------------
%     %% EnSight file conversion, save intermediate data
%     %% -------------------------------------------------    
%     if enSightConvFlag
%       
%         if nFE == 3
%         
%             % update text output
%             subplot('position',[0.15 0.05 0.8 0.05])
%             outStr = ['Dicom --> EnSight, saving data to disk , slice ',num2str(k),' ...'];
%             cla;
%             text(0.02,0.5,outStr,'FontSize',8,'Interpreter','none');
%             axis off;
%             colormap(gray);
%             drawnow;
% 
%             % for each slice/partition 
%             % write corrected and filtered data to disk for later use
%             fileStrMag(k,:)  = sprintf('%s%s%s%s%s%s%s',ensightDirPath,'\',patientName,'_sl',num2str(k,'%.3d'),'_mag','.mat');
%             fileStrFlow(k,:) = sprintf('%s%s%s%s%s%s%s',ensightDirPath,'\',patientName,'_sl',num2str(k,'%.3d'),'_flow','.mat');
%             save(fileStrMag(k,:),'imaMAG','k','numPhases','numSlices');
%             save(fileStrFlow(k,:),'imaFLOW');  
%             % save speed images if selected by user
%             if pcmraFlag
%                 fileStrSpeed(k,:)  = sprintf('%s%s%s%s%s%s%s',ensightDirPath,'\',patientName,'_sl',num2str(k,'%.3d'),'_spd','.mat');
%                 save(fileStrSpeed(k,:),'imaSPEED');  
%             end
%         
%         end  
%     end
% 
%     
%     
%     
%     
%     
%     
%     
%     
%     
%     
%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%     
%     
%     
%     
%  %% ------------------------------------------------------
%  %% EnSight file conversion, generate EnSight data files
%  %% ------------------------------------------------------    
%  if enSightConvFlag
%      
%      if (nFE == 3)
%         
%         % arrange / mirror velocity data according to the main
%         % orientation and in-plane phase encoding direction
%         if orientation == 'cor'  
%             signVx = 1;
%             signVy = -1;
%             signVz = -1;
%         elseif orientation == 'sag'
%             signVx = 1;
%             signVy = 1;
%             signVz = 1;
%         elseif orientation == 'tra'
%             signVx = -1;
%             signVy = -1;
%             signVz = -1;
%         else
%             subplot('position',[0.15 0.05 0.8 0.05])
%             outStr = ['Dicom --> EnSight Error: Main orientation not recognized - aborting !! '];
%             cla;
%             text(0.02,0.5,outStr,'FontSize',8,'Interpreter','none');
%             axis off;
%             colormap(gray);
%             drawnow;
%             return;
%         end
%                 
%         for m = 1:numPhases
%  
%             % update text output
%             subplot('position',[0.15 0.05 0.8 0.05])
%             outStr = ['Dicom --> EnSight, generating EnSight data files , phase ',num2str(m),' ...'];
%             cla;
%             text(0.02,0.5,outStr,'FontSize',8,'Interpreter','none');
%             axis off;
%             colormap(gray);
%             drawnow;
%                
%             % generate data files
%             % ------------------------------------------------
%             % open binary file an write data
%             dataPathMag  = sprintf('%s%s%s',dataPathName,num2str(m-1,'%02d'),'.mag');
%             dataPathFlow = sprintf('%s%s%s',dataPathName,num2str(m-1,'%02d'),'.vel');
%             fidMag       = fopen(dataPathMag, 'w', 'ieee-be');
%             fidFlow      = fopen(dataPathFlow, 'w', 'ieee-be');
%             part         = 1;
%             
%             % write data header
%             [dummy szStr]= size(sprintf('%s%s','Magnitude, time ',num2str(m-1)));
%             dataheaderStr(1:szStr) = sprintf('%s%s','Magnitude, time ',num2str(m-1));
%             dataheaderStr(szStr+1:80) = ' ';
%             fwrite(fidMag,dataheaderStr,'char');
%             [dummy szStr]= size(sprintf('%s%s','Velocity, time ',num2str(m-1)));
%             dataheaderStr(1:szStr) = sprintf('%s%s','Velocity, time ',num2str(m-1));
%             dataheaderStr(szStr+1:80) = ' ';
%             fwrite(fidFlow,dataheaderStr,'char');
%             dataheaderStr(1:4) = 'part';
%             dataheaderStr(5:80) = ' ';
%             fwrite(fidMag,dataheaderStr,'char');
%             fwrite(fidFlow,dataheaderStr,'char');
% 
%             fwrite(fidMag,part,'int');     
%             fwrite(fidFlow,part,'int');      
% 
%             geoheaderStr(1:5) = 'block';
%             geoheaderStr(6:80) = ' ';
%             fwrite(fidMag,geoheaderStr,'char');
%             fwrite(fidFlow,geoheaderStr,'char');
%             
%             
%                        
%             % load 2D/3D data volume for each time frame (phase)
%             for k = 1:numSlices
%         
%                 % for each phase load previously saved matlab files
%                 fileStrMag(k,:)  = sprintf('%s%s%s%s%s%s%s',ensightDirPath,'\',patientName,'_sl',num2str(k,'%.3d'),'_mag','.mat');
%                 fileStrFlow(k,:) = sprintf('%s%s%s%s%s%s%s',ensightDirPath,'\',patientName,'_sl',num2str(k,'%.3d'),'_flow','.mat');
%                 load(fileStrMag(k,:));
%                 load(fileStrFlow(k,:)); 
%                 if pcmraFlag
%                     fileStrSpeed(k,:) = sprintf('%s%s%s%s%s%s%s',ensightDirPath,'\',patientName,'_sl',num2str(k,'%.3d'),'_spd','.mat');
%                     load(fileStrSpeed(k,:));
%                     % write speed data
%                     dataTmpSpeed(:,:) = squeeze(imaSPEED(:,:,m));
%                     fwrite(fidSpeed,dataTmpSpeed(:),'float');
%                 end
%                 
%                 % write magnitude data
%                 dataTmpMag(:,:) = squeeze(imaMAG(:,:,m));
%                 fwrite(fidMag,dataTmpMag(:),'float');
%                 
%                 % construct velocity data, convert to m/s and adjust sign & orientation
%                 if peDir == 'x'
%                     dataTmpFlow(:,:,k,1) = squeeze(imaFLOW(:,:,m,2)) * 0.01 * signVx;
%                     dataTmpFlow(:,:,k,2) = squeeze(imaFLOW(:,:,m,1)) * 0.01 * signVy;
%                 elseif peDir == 'y'
%                     dataTmpFlow(:,:,k,1) = squeeze(imaFLOW(:,:,m,1)) * 0.01 * signVx;
%                     dataTmpFlow(:,:,k,2) = squeeze(imaFLOW(:,:,m,2)) * 0.01 * signVy;
%                 else
%                     subplot('position',[0.15 0.05 0.8 0.05])
%                     outStr = ['Dicom --> EnSight Error: PE-direction not recognized - aborting !! '];
%                     cla;
%                     text(0.02,0.5,outStr,'FontSize',8,'Interpreter','none');
%                     axis off;
%                     colormap(gray);
%                     drawnow;
%                     return;  
%                 end
%                 dataTmpFlow(:,:,k,3) = squeeze(imaFLOW(:,:,m,3)) * 0.01 * signVz;
%                                 
%                 clear imaMAG imaFLOW imaSPEED
%               
%             end % numslices
%             
%             % write velocity data (m/s)
%             fwrite(fidFlow,dataTmpFlow(:),'float');
%             
%             fclose(fidMag);
%             fclose(fidFlow);
%             if pcmraFlag
%                 fclose(fidSpeed);
%             end
%             
%         end; % end phase loop
%             
%         % delete intermediate matlab files
%         fileStrDelete = sprintf('%s%s%s',ensightDirPath,'\','*.mat');
%         delete(fileStrDelete);
%             
%     end % if nFE
% end % end EnSight conversion
