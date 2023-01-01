%% Arch PWV Analysis Script
%% Analysis
close all,clear all,clc
Files = dir ("*.xlsx");
%Get TR from dicom headers. Can differ between patients.

% n=0;

for ii = 23:length(Files)-1
%     n = n+1;
DescIndexRange = 1; DescRange = 1; AscRange = 1; AscIndexRange = 1; Asc_rsq = 0; Desc_rsq = 0;
num = importfile(Files(ii).name);

time = getfield(num,'Times');
flow = getfield(num,'FlowRatemls');
area = getfield(num,'Areamm');
label(ii,1) = convertCharsToStrings(Files(ii).name);
% Ascending aorta is first half of data, descending is second

% time = time(1:end/2); %Ascending aorta
ascending = flow(1:end/2);
% area = area(1:end/2);

time = time(((end/2)+1):end); %Descending aorta
descending = flow(((end/2)+1):end);
area = area(((end/2)+1):end);



% ascending = input('Please copy column vector of ascending aorta flow waveform: (use brackets)');
% descending = input('Please copy column vector of descending aorta flow waveform: (use brackets)');
% Projections = readmda('projections.mda');
% Projections = Projections';
% References = readmda('references.mda');
% hr=40; %beats per minute
% tr = 4.1752;
% % tr = 4.074699878;
% % **Important note: If not all (or too many) of the peaks are detected,
% % change the hr
% %%
% interval = round((60/hr)*(1000/tr));
% close all
% figure('units','normalized','outerposition',[0 0 0.3 1]), imagesc(References), colormap gray, axis image
% figure('units','normalized','outerposition',[0.3 0 0.7 1]), imagesc(Projections(800:1600,:)'), colormap gray
% %% Selection of proximal and distal aorta locations in projections
% figure(2), title('select proximal axis')
% [aMask,xa,ya] = roipoly;
% 
% ya1 = zeros(1,length(Projections(:,1)));
% ya1(:) = round(ya(1));
% ya2 = zeros(1,length(Projections(:,1)));
% ya2(:) = round(ya(2));
% %
% figure(2)
% hold on, plot(1:length(ya1),ya1,'r')
% hold on, plot(1:length(ya1),ya2,'r')
% figure(1), hold on,
% plot(1:size(References,2),ya1(1:size(References,2)),'r')
% plot(1:size(References,2),ya2(1:size(References,2)),'r')
% 
% figure(2), title('select distal axis')
% [dMask,xd,yd] = roipoly;
% 
% yd1 = zeros(1,length(Projections(:,1)));
% yd1(:) = round(yd(1));
% yd2 = zeros(1,length(Projections(:,1)));
% yd2(:) = round(yd(2));
% 
% hold on, plot(1:length(ya1),yd1,'b')
% hold on, plot(1:length(ya1),yd2,'b')
% 
% figure(1), hold on,
% plot(1:size(References,2),yd1(1:size(References,2)),'r')
% plot(1:size(References,2),yd2(1:size(References,2)),'r')
% %%
% 
% close all
% figure('units','normalized','outerposition',[0 0 1 1]), subplot(1,2,1), imagesc(References), colormap gray
% hold on
% subplot(1,2,1), plot(1:length(ya1),ya1,'r')
% subplot(1,2,1), plot(1:length(ya1),ya2,'r')
% subplot(1,2,1), plot(1:length(ya1),yd1,'b')
% subplot(1,2,1), plot(1:length(ya1),yd2,'b')
% 
% figure(1), subplot(1,2,2), imagesc(Projections(100:750,:)'), colormap gray
% hold on
% subplot(1,2,2), plot(1:length(ya1),ya1,'r')
% subplot(1,2,2), plot(1:length(ya1),ya2,'r')
% subplot(1,2,2), plot(1:length(ya1),yd1,'b')
% subplot(1,2,2), plot(1:length(ya1),yd2,'b')
% 
% for k = 1:length(Projections(:,1))
%     ascending(k) = mean(Projections(k,ya1(1):ya2(1)));
%     descending(k) = mean(Projections(k,yd1(1):yd2(1)));
% end
%%
% figure(4), plot(ascending,'r'), hold on, plot(descending,'b')
% asc_avg = zeros(1,length(ascending)-1);
% desc_avg = zeros(1,length(descending)-1);
%
% for k = 1:length(ascending)-1
%     asc_avg(k) = (ascending(k)+ascending(k+1))/2;
%     desc_avg(k) = (descending(k)+descending(k+1))/2;
% end
% interval = 
%% READ ASCENDING AND DESCENDING

ascending2(1:20) = ascending(21:40);
ascending2(21:40) = ascending(1:20);
descending2(1:20) = -descending(21:40);
descending2(21:40) = -descending(1:20);

MPD = 250;%round(0.5*interval*10); % Set minimum peak distance to be at least 80% of RR interval

% AscAvgInterp = interp(asc_avg,10); % Interpolate ascending and descending data
% DescAvgInterp = interp(desc_avg,10);
x = 1:length(ascending2);
xi = 1:length(ascending2)*10;
xi = xi/10;

AscAvgInterp1 = interp1(x,ascending2,xi);
DescAvgInterp1 = interp1(x,descending2,xi); % New interpolation = resampling

for k = 1:length(AscAvgInterp1)-10
    AscAvgInterp(k) = mean(AscAvgInterp1(k:k+10));
    DescAvgInterp(k) = mean(DescAvgInterp1(k:k+10));
end
%
% figure, hold on
% plot(xi,AscAvgInterp1,'r')
% plot(xi(13:end-13),AscAvgInterp)
% plot(x,asc_avg,'k.')
% ginput(1);
%
close all, figure(1)
[AscPeaks,AscLocs] = findpeaks(AscAvgInterp,'minpeakdistance',MPD);
[DescPeaks,DescLocs] = findpeaks(DescAvgInterp,'minpeakdistance',MPD);

figure('units','normalized','outerposition',[0 0.75 0.25 0.25])
plot(AscAvgInterp,'r'), hold on, plot(DescAvgInterp)
plot(AscLocs,AscPeaks,'k.','MarkerSize',12)
plot(DescLocs,DescPeaks,'k.','MarkerSize',12)


mean_Asc = AscAvgInterp;
mean_Desc = DescAvgInterp;

%
close all
LowerLimits = [0.1:0.05:0.4, 0.4, 0.2, 0.2, 0.2, 0.25, 0.25, 0.25, 0.25]; %Range for selection
UpperLimits = [0.9:-0.05:0.55, 0.6, 0.55, 0.5, 0.65, 0.6, 0.55, 0.5];
for irange = 1:length(LowerLimits)
    LowerLimit = LowerLimits(irange);
    UpperLimit = UpperLimits(irange);
    if irange == 1
        figure
        figure('units','normalized','outerposition',[0 0.75 0.25 0.25])
        plot(AscAvgInterp,'r.'), hold on, plot(DescAvgInterp,'b.')
        plot(AscLocs,AscPeaks,'k.','MarkerSize',12)
        plot(DescLocs,DescPeaks,'k.','MarkerSize',12)

        figure('units','normalized','outerposition',[0.75 0.75 0.25 0.25])
        figure('units','normalized','outerposition',[0.3 0.5 0.4 0.5])
        figure('units','normalized','outerposition',[0.3 0 0.4 0.5])
    end
    [~,idx] = max(mean_Desc);
    loc1 = find(isnan(mean_Desc)-1,1);
    
    LocalAscBase = mean_Asc(loc1+10:loc1+100);
    LocalDescBase = mean_Desc(loc1+10:loc1+100);
    
    figure(3), plot(LocalAscBase,'r'), hold on, plot(LocalDescBase,'b')
    LocalAscBaseAvg = mean(LocalAscBase);
    LocalDescBaseAvg = mean(LocalDescBase);
    LocalAscNorm = mean_Asc;%-LocalAscBaseAvg;
    LocalDescNorm = mean_Desc;%-LocalDescBaseAvg;
    LocalAscPeak = max(LocalAscNorm);
    LocalAscPeakLoc = find(LocalAscNorm == max(LocalAscNorm));
    LocalDescPeak = max(LocalDescNorm);
    LocalDescPeakLoc = find(LocalAscNorm == max(LocalDescNorm));
    LocalAscNormScaled = LocalAscNorm*(1/LocalAscPeak);
    LocalDescNormScaled = LocalDescNorm*(1/LocalDescPeak);
    figure(4),hold off
    plot(LocalAscNormScaled,'r.'), hold on, plot(LocalDescNormScaled,'b.')
    hold on,plot([0 length(LocalDescNormScaled)], [LowerLimit,LowerLimit],'k')
    hold on, plot([0 length(LocalDescNormScaled)], [UpperLimit,UpperLimit],'k')
    AscLL = LocalAscNormScaled(find(LocalAscNormScaled>LowerLimit));
    AscIndexLL = find(LocalAscNormScaled>LowerLimit);
    DescLL = LocalDescNormScaled(find(LocalDescNormScaled>LowerLimit));
    DescIndexLL = find(LocalDescNormScaled>LowerLimit);
    AscRange2 = AscLL(find(AscLL<UpperLimit));
    AscIndexRange2 = AscIndexLL(find(AscLL<UpperLimit));
    DescRange2 = DescLL(find(DescLL<UpperLimit));
    DescIndexRange2 = DescIndexLL(find(DescLL<UpperLimit));
    clear AscRange AscIndexRange DescIndexRange DescRange
    DescIndexRange = 1; DescRange = 1; AscRange = 1; AscIndexRange = 1;
    for j = 1:length(AscIndexRange2)
        if j<length(AscIndexRange2)
            if AscIndexRange2(j)+1==AscIndexRange2(j+1)
                AscIndexRange(j) = AscIndexRange2(j);
                AscRange(j) = AscRange2(j);
            else
                break
            end
        elseif j==length(AscIndexRange2)
            AscIndexRange(j) = AscIndexRange2(j);
            AscRange(j) = AscRange2(j);
        end
    end
    
    for j = 1:length(DescIndexRange2)
        if j<length(DescIndexRange2)
            if DescIndexRange2(j)+1==DescIndexRange2(j+1)
                DescIndexRange(j) = DescIndexRange2(j);
                DescRange(j) = DescRange2(j);
            else
                break
            end
        elseif j==length(DescIndexRange2)
            DescIndexRange(j) = DescIndexRange2(j);
            DescRange(j) = DescRange2(j);
        end
    end
    
    tAsc = min(AscIndexRange):max(AscIndexRange);
    tDesc = min(DescIndexRange):max(DescIndexRange);
    AscPolyFit = polyfit(AscIndexRange,AscRange,1);
    DescPolyFit = polyfit(DescIndexRange,DescRange,1);
    
    if length(AscRange)>3       % Added by Danny to bypass ranges with too few points
    Asc_s = regstats(AscRange,AscIndexRange,'linear');
    %                 Asc_yhat(:,k) = Asc_s(k).yhat;
    Asc_rsq = Asc_s.rsquare;
    end
    if length(DescRange)>3
    Desc_s = regstats(DescRange,DescIndexRange,'linear');
    %                 Desc_yhat(:,k) = Desc_s(k).yhat;
    Desc_rsq = Desc_s.rsquare;
    end
    
    rsq = min(Asc_rsq,Desc_rsq);
    if rsq>0.995
        
        Asc_yx(1) = 1/AscPolyFit(1);
        Asc_yx(2) = -AscPolyFit(2)/AscPolyFit(1);
        Desc_yx(1) = 1/DescPolyFit(1);
        Desc_yx(2) = -DescPolyFit(2)/DescPolyFit(1);
        
        
        AscFitL = Asc_yx(1)*LowerLimit+Asc_yx(2);
        AscFitU = Asc_yx(1)*UpperLimit+Asc_yx(2);
        DescFitL = Desc_yx(1)*LowerLimit+Desc_yx(2);
        DescFitU = Desc_yx(1)*UpperLimit+Desc_yx(2);
        hold off
        AscPolyFit(1);
        DescPolyFit(1);
        figure(5), plot(AscIndexRange,AscRange,'r.')
        hold on, plot(DescIndexRange,DescRange,'b.')
        plot([AscFitL,AscFitU],[LowerLimit,UpperLimit],'k')
        plot([DescFitL,DescFitU],[LowerLimit,UpperLimit],'k')
        drawnow, pause(0.1), hold off
        AscFits(1) = Asc_yx(2);
        AscFits(2) = AscFitL;
        AscFits(3) = AscFitU;
        
        DescFits(1) = Desc_yx(2);
        DescFits(2) = DescFitL;
        DescFits(3) = DescFitU;
        
        deltaT_Intercept = (Desc_yx(2)-Asc_yx(2))/10;
        deltaT = ((DescFitU-AscFitU)+(DescFitL-AscFitL))/20;
        %                     [tempx,tempy,button]=ginput(1);
        %                     button_list(k) = button;
        button_list = 1;
        %                     button
        button = 1;
    else
        deltaT_Intercept = 0;
        deltaT = 0;
        button = 32;
        figure(5), title('Rsq too low')
    end
    

rsq_range(irange) = rsq;
deltaT_Intercept_range(irange) = deltaT_Intercept;
deltaT_range(irange) = deltaT;
pause(1)
% ginput(1);
end
ranges = [LowerLimits', UpperLimits'];
output = [ranges,rsq_range',deltaT_Intercept_range',deltaT_range']
save('20211215_2DPC_Results.mat') %Make sure a good output range is saved
output(13,:)
final(ii,:) = [label(ii),output(13,:)];
end