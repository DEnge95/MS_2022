close all;
clc;
clearvars -except n j output Files;
% rootstring = '/Users/ailguestuser/Desktop/Danny_mixedVA';
% Please provide path to the data folder
% rootpath = rootstring;
% cd(rootpath);
% n=0; j=1;

Files = dir ("*.xlsx");
group = [8, 14, 16, 19, 3, 4, 35, 36, 45, 47, 48];
% n=0;

% for i = [34,38,39,46,48,50,52,53,55,60]
for i = 49 % group(j) %
   % length(Files)
num = importfile(Files(i).name);

time = getfield(num,'Times');
flow = getfield(num,'FlowRatemls');
area = getfield(num,'Areamm');

% Ascending aorta is first half of data, descending is second

% time = time(1:end/2); %Ascending aorta
% flow = flow(1:end/2);
% area = area(1:end/2);

time = time(((end/2)+1):end); %Descending aorta
flow = flow(((end/2)+1):end);
area = area(((end/2)+1):end);



area = smooth(area);
flow = smooth(flow);

% if flow(5) < 0
%     flow = -flow
% end 

T = max(time);
time(1)=1;
% Creates time line with TR of 11ms up to given period T
time_i = 0:10:T;
flow_i = spline(time,flow,time_i); % newly interpolated flow waveform [U(t)]
area_i = spline(time,area,time_i); % new;y interpolated area wavefrom [P(t)]

max_A = max(area);
min_A = min(area);
RAC = (max_A-min_A)/max_A*100;

% The real TR
TR = T/length(time);

dQdt = gradient(flow_i,10); %derivative of flow curve
dAdt = gradient(area_i,10); %derivative of area curve 

figure
subplot(3,3,5)
plot(area_i,flow_i,'o')
xlabel('Area (cm^2)')
ylabel('Flow (L/min)')

subplot(3,3,6)
plot(area(1:10),flow(1:10),'o')
xlabel('Area (cm^2)')
ylabel('Flow (L/min)')

view = [area(1:10),flow(1:10)] %#ok<NOPTS> %Showing points for PWV selection
    PWV = LinearModel.fit(area(20:25),flow(20:25)); %Breakpoint for selecting 4-6 points for PWV (c)
c = -PWV.Coefficients{2,1} %#ok<NOPTS>

dQ_plus = 0.5*(dQdt + c*dAdt);
dQ_minus = 0.5*(dQdt - c*dAdt);


% for i = 1:length(dQdt)
%     dA_plus(i) = 0.5*(area_i(i) + (1/c)*dQdt(i));
%     dA_minus(i) = 0.5*(area_i(i) - (1/c)*dQdt(i));
% end 


dA_plus = 0.5*(dAdt+(1/c)*dQdt);
dA_minus = 0.5*(dAdt - (1/c)*dQdt);


dI_net = dAdt.*dQdt;
dI_plus = (c/4)*((dAdt+(dQdt/c)).^2)*1000; %To fix FCW/BCW try moving minus sign here
dI_minus = -(c/4)*((dAdt-(dQdt/c)).^2)*1000;

WRI = -trapz(time_i,dI_minus)/trapz(time_i,dI_plus);
WRI2 = -trapz(time_i,dI_plus)/trapz(time_i,dI_minus);

FCW_m = max(dI_plus);
[a,b] = max(dI_plus);
FCW_time = time_i(b);

BCW_m = min(dI_minus);
[e,f] = min(dI_minus);
BCW_time = time_i(f);

j = j+1; %advancing through list of scans
n = n+1; 

output(n,:) = [c,RAC,FCW_m,-BCW_m,WRI,WRI2];

subplot(3,3,1)
plot(time,area,'o-')
xlabel('Time (ms)')
ylabel('Area (cm^2)')

subplot(3,3,2)
plot(time,flow,'o-')
xlabel('Time (ms)')
ylabel('Flow (cm^2)')

subplot(3,3,3)
plot(time_i,area_i,'o-')
xlabel('Time (ms)')
ylabel('Area-interpolated (cm^2)')

subplot(3,3,4)
plot(time_i,flow_i,'o-')
xlabel('Time (ms)')
ylabel('Flow-interpolated (L/min)')




subplot(3,3,7)
plot(time_i,dI_plus,'--')
hold on
plot(time_i,dI_minus,'--')
xlabel('Time (ms)')
ylabel('WI (cm^5/s)')

subplot(3,3,8)
plot(time_i,dA_plus,'-')
hold on
plot(time_i(1:length(dQdt)),dA_minus,'-')
xlabel('Time (ms)')
ylabel('Delta Area (cm^2)')

subplot(3,3,9)
plot(time_i,dI_plus,'--')
hold on
plot(time_i,dI_minus,'--')
hold on 
plot(time_i,dA_plus*1000,'-')
hold on
plot(time_i,dA_minus*1000,'-')

xlabel('Time (ms)')
ylabel('WI (cm^5/s) and dA data')


end








