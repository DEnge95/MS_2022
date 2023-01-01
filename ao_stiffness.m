function [pwv,rmse]=ao_stiffness(strMrstructMag,savePath,tr,flowdata)

%AO_STIFFNESS - Evaluate aortic stiffness by pulse wave velocity
%
%   Syntax:  [pwv,rmse]=ao_stiffness(strMrstructMag)
%
%   Example: 
%       [pwv,rmse]=ao_stiffness('') for selection of mag file
%       [pwv,rmse]=ao_stiffness('',savePath,tr,flowdata) for batch
%
%   Other m-files required:
%   Subfunctions:
%   MAT-files required:
%   
% Kelly Jarvis
% 2019

%% error check
pwv = 0;
rmse = 0;

if nargin < 1
   strMrstructMag = '';
end

%% load data

if nargin < 2
    if ~exist(strMrstructMag) == 2 || isempty(strMrstructMag) 
        [magName, magPath] = uigetfile('*.mat','Load the mag mrstruct file','Multiselect','Off');
    else
        [magPath, name, ext] = fileparts(strMrstructMag);
        magName = [name,ext];
    end
    X = xlsread(fullfile(magPath,'transient_flowdata.xlsx'));
    [magStruct,path] = mrstruct_read(fullfile(magPath,magName));
    pStruct = struct('length',X(2:end-1,1),'times','','flow',X(2:end-1,2:1+(end-1)/4)); %remove first plane and last
    clear X;
    
    newSubFolder='results_pwv';
    %savePath=fullfile(magPath,newSubFolder); %save in mrstruct folder
    idcs   = strfind(magPath,'\'); magPath1 = magPath(1:idcs(end)-1);
    savePath=fullfile(magPath1,newSubFolder); %save one above mrstruct folder
    
    tr=magStruct.tr;
else
    pStruct = struct('length',flowdata(2:end-1,1),'times','','flow',flowdata(2:end-1,2:end)); %remove first plane and last
end

if ~exist(savePath, 'dir')
  mkdir(savePath);
end

%% interpolate data

x1 = 1:size(pStruct.flow,1);
x2 = tr/2:tr:tr*size(pStruct.flow,2)-tr/2;
x1q = 1:size(pStruct.flow,1); %no interp
x2q = tr/2:tr*size(pStruct.flow,2)-tr/2; % interp to 1 ms intervals
F = griddedInterpolant({x1,x2},pStruct.flow,'linear');
pStruct.flow = F({x1q,x2q});
pStruct.times = x2q;

%% smooth data
%pStruct.flow=smoothdata(pStruct.flow,1);

%% plot all flow waveforms
figure
for p=1:size(pStruct.flow,1) 
    plot(pStruct.times,pStruct.flow(p,:),'-b')
    hold on
end
hold off
title('Flow Waveforms')
xlabel('time [ms]')
ylabel('flow [ml/s]')
%ylabel('normalized flow [(ml/s)/(ml/s)]')
savefig(fullfile(savePath,'Flow Waveforms'))

%% pulse wave velocity

% find peakflow: flow idx time
[pStruct.peakflow(:,1),pStruct.peakflow(:,2)] = max(pStruct.flow,[],2);
pStruct.peakflow(:,3)=pStruct.times(pStruct.peakflow(:,2));

% find upstroke points
tempFlow=pStruct.flow;
for i=1:size(tempFlow,1)
    tempFlow(i,pStruct.peakflow(i,2):end)=nan; % remove timepts at peak and after
    tempFlow(i,tempFlow(i,:)>0.8*pStruct.peakflow(i,1))=nan;
    for j=pStruct.peakflow(i,2)-1:-1:1
        if tempFlow(i,j)<0.2*pStruct.peakflow(i,1) || tempFlow(i,j)>tempFlow(i,j+1)
            tempFlow(i,1:j)=nan;
            break
        end
    end
end
pStruct.upstroke = tempFlow; clear tempFlow;

% find TTF
pStruct.TTF=zeros(size(pStruct.flow,1),1);
for i=1:size(pStruct.upstroke,1)
    keep = ~(isnan(pStruct.upstroke(i,:))) ;
    pStruct.fit(i,:) = polyfit(pStruct.times(keep), pStruct.upstroke(i,keep),1);
    pStruct.TTF(i) = roots(pStruct.fit(i,:));
end
%
figure
for i=1:round(size(pStruct.upstroke,1)/5):size(pStruct.upstroke,1) %plot some of the planes
    plot(pStruct.times,pStruct.flow(i,:),'-b');
    keep = ~(isnan(pStruct.upstroke(i,:))) ;
    hold on
    plot(pStruct.times(keep),pStruct.flow(i,keep),'Color','b','Marker','o')
    f = polyval(pStruct.fit(i,:),pStruct.times(keep));
    plot(pStruct.times(keep),f,'-r')
    plot(pStruct.TTF(i),0,'*r','MarkerSize',7) %default size is 6
end
clear f
title('Time to Foot (TTF)')
xlabel('time [ms]')
ylabel('flow [ml/s]')
hold off
savefig(fullfile(savePath,'TTF'))

% get pulse wave velocity TTF
fit = polyfit(pStruct.length,pStruct.TTF,1);
pStruct.pwvTTF = 1/fit(1);
lm = fitlm(pStruct.length,pStruct.TTF,'linear');
pStruct.rmseTTF=lm.RMSE;
%
figure
scatter(pStruct.length,pStruct.TTF,20,'k','o')
hold on
f = polyval(fit,pStruct.length);
plot(pStruct.length,f,'k')
title('Pulse Wave Velocity TTF')
xlabel('length [mm]')
ylabel('time [ms]')
text1=xlim; text2=ylim;
text(0.05*(text1(2)-text1(1)),.95*text2(2),strjoin({'PWV = ',num2str(round(pStruct.pwvTTF,1),'%.1f'),'m/s'}))
hold off
clear fit f lm
savefig(fullfile(savePath,'PWV_TTF'))

% %upstroke only
% pStruct.flow=pStruct.upstroke;
% pStruct.flow(isnan(pStruct.flow))=0;
% for j=1:size(pStruct.flow,1)
%     pStruct.flow(j,:)=pStruct.flow(j,:)./max(pStruct.flow(j,:));
% end
% %end upstroke only

%get pulse wave velocity XCOR
Ts=pStruct.times(2)-pStruct.times(1); %sampling period
%figure
for p=1:size(pStruct.flow,1) 
    [acor,lag] = xcorr(pStruct.flow(1,:),pStruct.flow(p,:)); %compare each plane to first one
    %[~,Id] = max(abs(acor)); %xcor is maximum at lag=delay
    [~,Id] = max(acor); % raw correlation to avoid inverted flow match
    lagDiff = lag(Id);
    pStruct.XCORdelay(p,1) = lagDiff*Ts*-1; %timeDiff(1)=0
    pStruct.acor(p,1)=acor(Id); %report cross-correlation for each plane
    %plot(lag,acor)
end
fit = polyfit(pStruct.length,pStruct.XCORdelay,1);
pStruct.pwvXCOR = 1/fit(1); 
lm = fitlm(pStruct.length,pStruct.XCORdelay,'linear');
pStruct.rmseXCOR=lm.RMSE;
%
figure
scatter(pStruct.length,pStruct.XCORdelay,20,'k','o')
hold on
f = polyval(fit,pStruct.length);
plot(pStruct.length,f,'k')
title('Pulse Wave Velocity XCOR')
xlabel('length [mm]')
ylabel('time [ms]')
text1=xlim; text2=ylim;
text(0.05*(text1(2)-text1(1)),.95*text2(2),strjoin({'PWV = ',num2str(round(pStruct.pwvXCOR,1),'%.1f'),'m/s'}))
hold off
clear fit f lm acor lag Id lagDiff
savefig(fullfile(savePath,'PWV_XCOR'))

%get pulse wave velocity XCOR2
%figure
for p=1:size(pStruct.flow,1) 
    [acor,lag] = xcorr(pStruct.flow(round(size(pStruct.flow,1)/2),:),pStruct.flow(p,:)); %compare each plane to middle one
    %[~,Id] = max(abs(acor)); %xcor is maximum at lag=delay
    [~,Id] = max(acor); % raw correlation to avoid inverted flow match
    lagDiff = lag(Id);
    pStruct.XCORdelay2(p,1) = lagDiff*Ts*-1; %timeDiff(middle)=0
    pStruct.acor2(p,1)=acor(Id); %report cross-correlation for each plane
    %plot(lag,acor)
end
fit = polyfit(pStruct.length,pStruct.XCORdelay2,1);
pStruct.pwvXCOR2 = 1/fit(1); 
lm = fitlm(pStruct.length,pStruct.XCORdelay2,'linear');
pStruct.rmseXCOR2=lm.RMSE;
%
figure
scatter(pStruct.length,pStruct.XCORdelay2,20,'k','o')
hold on
f = polyval(fit,pStruct.length);
plot(pStruct.length,f,'k')
title('Pulse Wave Velocity XCOR2')
xlabel('length [mm]')
ylabel('time [ms]')
text1=xlim; text2=ylim;
text(0.05*(text1(2)-text1(1)),.95*text2(2),strjoin({'PWV = ',num2str(round(pStruct.pwvXCOR2,1),'%.1f'),'m/s'}))
hold off
clear fit f
savefig(fullfile(savePath,'PWV_XCOR2'))

%get pulse wave velocity XCOR3
for pr=1:size(pStruct.flow,1)
    for p=1:size(pStruct.flow,1) 
        [acor,lag] = xcorr(pStruct.flow(pr,:),pStruct.flow(p,:)); %compare each plane
        %[~,Id] = max(abs(acor)); %xcor is maximum at lag=delay
        [~,Id] = max(acor); % raw correlation to avoid inverted flow match
        lagDiff = lag(Id);
        pStruct.XCORdelay2(p,1) = lagDiff*Ts*-1; %timeDiff(middle)=0
        pStruct.acor2(p,1)=acor(Id); %report cross-correlation for each plane
        %plot(lag,acor)
    end
    fit = polyfit(pStruct.length,pStruct.XCORdelay2,1);
    pStruct.pwvXCOR3all(pr) = 1/fit(1); 
    lm = fitlm(pStruct.length,pStruct.XCORdelay2,'linear');
    pStruct.rmseXCOR3all(pr)=lm.RMSE;
end

normDist=1-lillietest(pStruct.pwvXCOR3all);
if normDist==1
    pStruct.pwvXCOR3=mean(pStruct.pwvXCOR3all);
    pStruct.rmseXCOR3=mean(pStruct.rmseXCOR3all);
else
    pStruct.pwvXCOR3=median(pStruct.pwvXCOR3all);
    pStruct.rmseXCOR3=median(pStruct.rmseXCOR3all);
    %[y, medianIDX] = min(abs(pStruct.pwvXCOR3all-pStruct.pwvXCOR3));%if odd, finds idx of first occurance of median; if even, finds idx of value closest to median 
    %pStruct.rmseXCOR3=pStruct.rmseXCOR3all(medianIDX);
end

%get descriptive info on PWV
pStruct.minpwvXCOR=min(pStruct.pwvXCOR3all);
pStruct.maxpwvXCOR=max(pStruct.pwvXCOR3all);
pStruct.rangepwvXCOR=pStruct.maxpwvXCOR - pStruct.minpwvXCOR;
pStruct.iqrpwvXCOR=iqr(pStruct.pwvXCOR3all);
pStruct.stdpwvXCOR=std(pStruct.pwvXCOR3all);
pStruct.normDistpwvXCOR=normDist;

figure
%
subplot(2,1,1)
plot(1:size(pStruct.pwvXCOR3all,2),pStruct.pwvXCOR3all,'-k.');
hold on;
plot([1,size(pStruct.pwvXCOR3all,2)],[pStruct.pwvXCOR3,pStruct.pwvXCOR3],'-b');
%plot(medianIDX,pStruct.pwvXCOR3,'bo') 
hold off;
title('PWV Analysis by XCOR')
ylabel('PWV [m/s]')
%
subplot(2,1,2)
plot(1:size(pStruct.rmseXCOR3all,2),pStruct.rmseXCOR3all,'-k.'); 
hold on;
plot([1,size(pStruct.rmseXCOR3all,2)],[pStruct.rmseXCOR3,pStruct.rmseXCOR3],'-b');
hold off;
ylabel('RMSE [m/s]')
xlabel('Reference Plane (Ao: Prox. to Dist.)')
savefig(fullfile(savePath,'PWV_XCORall'))

% %% TEMP
% %get pulse wave velocity XCOR
% Ts=pStruct.times(2)-pStruct.times(1); %sampling period
% %figure
% for p=1:size(pStruct.flow,1)
%     for p2 = 1:p
%         [acor,lag] = xcorr(pStruct.flow(p2,:),pStruct.flow(p,:)); %compare each plane to first one
%         %[~,Id] = max(abs(acor)); %xcor is maximum at lag=delay
%         [~,Id] = max(acor); % raw correlation to avoid inverted flow match
%         lagDiff = lag(Id);
%         pStruct.XCORdelayImage(p,p2) = lagDiff*Ts*-1; %timeDiff(1)=0
%         pStruct.acorImage(p,p2)=acor(Id); %report cross-correlation for each plane
%         %plot(lag,acor)
%     end
% end
% for ii = 1:size(pStruct.flow,1)-3 %3 arbitrary
%     fit = polyfit(pStruct.length(ii:end),pStruct.XCORdelayImage(ii:end,ii),1);
%     pStruct.pwvXCORvec(ii) = 1/fit(1);
%     lm = fitlm(pStruct.length,pStruct.XCORdelay,'linear');
%     pStruct.rmseXCOR=lm.RMSE;
% end
% figure
% pStruct.XCORdelayImage(pStruct.XCORdelayImage==0)=nan;
% surf(pStruct.XCORdelayImage);
% figure
% pStruct.acorImage(pStruct.acorImage==0)=nan;
% surf(pStruct.acorImage);
% 
% %% TEMP

pwv(1)=pStruct.pwvTTF; rmse(1)=pStruct.rmseTTF;
pwv(2)=pStruct.pwvXCOR; rmse(2)=pStruct.rmseXCOR;
pwv(3)=pStruct.pwvXCOR2; rmse(3)=pStruct.rmseXCOR2;
pwv(4)=pStruct.pwvXCOR3; rmse(4)=pStruct.rmseXCOR3;
pwv(5)=pStruct.minpwvXCOR;
pwv(6)=pStruct.maxpwvXCOR;
pwv(7)=pStruct.rangepwvXCOR;
pwv(8)=pStruct.iqrpwvXCOR;
pwv(9)=pStruct.stdpwvXCOR;
pwv(10)=pStruct.normDistpwvXCOR;

save (fullfile(savePath,'pStruct.mat'),'pStruct')
tic;

%% misc
   
% for plotting centerline and iso and mag slices
% a = magStruct.dataAy;
% size(a) ans = 160   120    26    40
% figure; slice(a(:,:,:,10),80,60,13)
% b = centerline1Aorta;
% hold on; plot3(b(:,1),b(:,2),b(:,3));
% hold on; scatter3(b(:,1),b(:,2),b(:,3));
% load('Aorta_nobranch_grayvalues_mask_struct.mat')
% c = SegStruct;
% c = SegStruct.dataAy;
% hold on; patch(isosurface(c,0.5));
% axes equal

