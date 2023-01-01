clc;
clearvars;
close all;



Path = 'F:\Aorta_results.xlsx';
[Dat,~,~] = xlsread(Path,2);
Data = Dat(:,[2:3,5:6,10:11,24:30]); %for correlation
Pre = Data([1:11,23:28],:);
Post = Data([12:22,29:34],:);

figure
[Corrpre, Ppre] = corrcoef(Pre,'Rows','complete');
corrplot(Pre,'Rows','complete','varNames',{'4DPWV','SphPWV','VO2','GIR','EDVn','ESVn','SV','HR','CO','SBP','DBP','EF','PP'}); % aortic measures

figure
[Corrpost, Ppost] = corrcoef(Post,'Rows','complete');
corrplot(Post,'Rows','complete','varNames',{'4DPWV','SphPWV','VO2','GIR','EDVn','ESVn','SV','HR','CO','SBP','DBP','EF','PP'}); % aortic measures

% [Corr, P] = corrcoef(Data,'Rows','complete');
% corrplot(Data,'Rows','complete','varNames',{'4DPWV','SphPWV','VO2','GIR','EDVn','ESVn','SV','HR','CO','SBP','DBP','EF','PP'}); % aortic measures
Conpre = Dat(1:11,[2:6,10:11]);    
Conpost = Dat(12:22,[2:6,10:11]);

Diapre = Dat(23:28,[2:6,10:11]);
Diapost = Dat(29:34,[2:6,10:11]);

NWC = [6.4581,0,0,0,0,0];
NWCe = [1.1552,0,0,0,0,0];

% % clinical measures
% Conpre = Dat(1:9,12:23);    
% Conpost = Dat(10:18,12:23);
% 
% Diapre = Dat(19:23,12:23);
% Diapost = Dat(24:28,12:23);

Cpr = nanmean(Conpre);     Cpre = nanstd(Conpre);
Cpo = nanmean(Conpost);    Cpoe = nanstd(Conpost);
Dpr = nanmean(Diapre);     Dpre = nanstd(Diapre);
Dpo = nanmean(Diapost);    Dpoe = nanstd(Diapost);
n = [2,3;4,5;5,6];
% Titles = ["Normalized EDV (mL/m^2)","Normalized ESV (mL/m^2)",...
%     "Septal wall thickness (mm)", "Free wall thickness (mm)"];
% dim = [.2 .5 .3 .3];
% str = 'Ex trained:';


%% Pre-only analysis, by disease
close all

Con = Dat(1:14,6:8);
Dia = Dat(22:29,6:8);

Cpr = mean(Con);     Cpre = std(Con)/sqrt(length(Con(:,1)));
Dpr = mean(Dia);     Dpre = std(Dia)/sqrt(length(Dia(:,1)));
n = [1,2;3,4;5,6];
Titles = ["4D PWV (m/s)","Sphygmacor PWV (m/s)","2D PWV (m/s)","RAC (%)","Distensibilty (10-3 mm Hg-1)","FCW_m (mL/s)","BCW_m (mL/s)","WRI (Ratio, 0-1)","Vo2 Peak (mL/kg/min)","GIR (mL/kg-min)","EDV (mL)","ESV (mL)","Normalized EDV (mL/m2)","Normalized ESV (mL/m2)","Septal wall thickness (mm)","Free Wall Thickness (mm)","Rad. Peak strain (%)","Circ. Peak Strain (%)","Long. Peak Strain (%)","Radial, systolic peak SR","Circumferential, systolic peak SR","Longitudinal, systolic peak SR","Radial, diastolic peak SR","Circumferential, diastolic peak SR","Longitudinal, diastolic peak SR","Age (years)","Lean body mass (kg)","Body mass index (kg/m2)","Body Fat %","HbA-1c (%)","Insulin (µU/mL) in serum","Glucose at Pre-scan (mg/dL)","Adiponectin (µg/mL)","LDL Cholesterol (mg/dL)","HDL Cholesterol (mg/dL)","Triglycerides (mg/dL)","Stroke Volume (mL)","Heart Rate (1/min)","Cardiac Output (mL/min)","Systolic BP (mmHg)","Diastolic BP (mmHg)","Ejection Fraction (Ratio, 0-1)","Pulse Pressure (mmHg)","BSA (m2)"];
dim = [.2 .5 .3 .3];

for j = 1:3
    figure
    bar([1,2],diag([Cpr(j);Dpr(j)]),'stacked');
    hold on
    annotation('textbox',dim,'String',str,'FitBoxToText','on','LineStyle','none','FontWeight','bold');
    xticklabels({'-','-'})
    ylabel(Titles(j),'FontWeight','bold')
    xlabel('Control                    Diabetes','FontWeight','bold')
    for i = 1:length(Con(:,1))
        plot(1,Con(i,j),'ko')
    end
    for k = 1:length(Dia(:,1))
        plot(2,Dia(k,j),'ko')
    end
    errorbar([1;2],[Cpr(j);Dpr(j)],[Cpre(j);Dpre(j)],'k.','LineWidth',1.5,'Capsize',10)
    [h(j),p(j)] = ttest2(Con(:,j),Dia(:,j));
end

%% Pre/post analysis, by disease, new format
Titles = ["4D Flow PWV (m/s)","SphygmoCor PWV (m/s)","2D PWV (m/s)","Vo2 Peak (mL/kg/min)","GIR (mL/kg-min)","EDV_n (mL/m^2)","ESV_n (mL/m^2)"];
dim = [.2 .5 .3 .3];
str = 'Ex trained:';
n = [2,3;4,5;5,6];

% Colors = ['b','b','r','r'];
for j = 1
    figure
    bar([0.5;2;3;4;5],[0;Cpr(j);Cpo(j);0;0],'b');
    hold on
    xticks([0.5 2 3 4 5]);
    xticklabels({'-','-','+','-','+'});
    xline(1.27)
    bar([4;5],[Dpr(j);Dpo(j)],'r')
    bar(0.5,NWC(j),'g')
    annotation('textbox',dim,'String',str,'FitBoxToText','on','LineStyle','none','FontWeight','bold','FontSize',12);
    ylabel(Titles(j),'FontWeight','bold')
    xlabel('Lean-Control (n=20)   BMI-Control (n=11)     T2D (n=6)','FontWeight','bold');
    ax = gca;
    ax.FontSize = 10;
    ylim([0 15])
    for i = 1:length(Conpre(:,1))
        plot(n(1,:),[Conpre(i,j),Conpost(i,j)],'k-o')
    end
    for k = 1:length(Diapre(:,1))
        plot(n(2,:),[Diapre(k,j),Diapost(k,j)],'k-o')
    end
    errorbar([2;3;4;5;0.5],[Cpr(j);Cpo(j);Dpr(j);Dpo(j);NWC(j)],[Cpre(j);Cpoe(j);Dpre(j);Dpoe(j);NWCe(j)],'k.','LineWidth',1.5,'Capsize',10)
    [~,pc(j)] = ttest(Conpre(:,j),Conpost(:,j));
    [~,pd(j)] = ttest(Diapre(:,j),Diapost(:,j));
    [~,pa(j)] = ttest([Conpre(:,j);Diapre(:,j)],[Conpost(:,j);Diapost(:,j)]);
    [~,pp(j)] = ttest2(Conpre(:,j),Diapre(:,j));
end

%% Pre/post analysis, by disease, old format
Titles = ["4D Flow PWV (m/s)","SphygmoCor PWV (m/s)","2D PWV (m/s)","Vo2 Peak (mL/kg/min)","GIR (mL/kg-min)","EDV_n (mL/m^2)","ESV_n (mL/m^2)"];
dim = [.2 .5 .3 .3];
str = 'Ex trained:';
n = [1,2;3,4;5,6];
% Colors = ['b','b','r','r'];
for j = 2
    figure
    bar([1;2;3;4],[Cpr(j);Cpo(j);0;0],'b');
    hold on
    xticklabels({'-','+','-','+'})
    bar([3;4],[Dpr(j);Dpo(j)],'r')
%     bar(5,NWC(j),'g')
    annotation('textbox',dim,'String',str,'FitBoxToText','on','LineStyle','none','FontWeight','bold','FontSize',12);
    ylabel(Titles(j),'FontWeight','bold')
    xlabel('BMI-Control (n=10)               T2D (n=5)','FontWeight','bold')
    ax = gca;
    ax.FontSize = 10;
    for i = 1:8
        plot(n(1,:),[Conpre(i,j),Conpost(i,j)],'k-o')
    end
    for i = 10:length(Conpre(:,1))
        plot(n(1,:),[Conpre(i,j),Conpost(i,j)],'k-o')
    end
    for k = 1:length(Diapre(:,1))
        plot(n(2,:),[Diapre(k,j),Diapost(k,j)],'k-o')
    end
    errorbar([1;2;3;4],[Cpr(j);Cpo(j);Dpr(j);Dpo(j)],[Cpre(j);Cpoe(j);Dpre(j);Dpoe(j)],'k.','LineWidth',1.5,'Capsize',10)
    [~,pc(j)] = ttest(Conpre(:,j),Conpost(:,j));
    [~,pd(j)] = ttest(Diapre(:,j),Diapost(:,j));
    [~,pa(j)] = ttest([Conpre(:,j);Diapre(:,j)],[Conpost(:,j);Diapost(:,j)]);
    [~,pp(j)] = ttest2(Conpre(:,j),Diapre(:,j));
end

%% Vo2 peak & GIR - whole study
clc;
clearvars;
close all;



Path = 'E:\PhD\Danny_mixedVA\Demographics\Vo2andGIR_preandpost.csv';
[Dat,~,~] = xlsread(Path);
Con = Dat(1:56,3:6);
Conpre = Dat(:,[3,5]);    
Conpost = Dat(:,[4,6]);

Diapre = Dat(57:end,[3,5]);
Diapost = Dat(57:end,[4,6]);

Cpr = nanmean(Conpre);     Cpre = nanstd(Conpre);
Cpo = nanmean(Conpost);    Cpoe = nanstd(Conpost);
Dpr = nanmean(Diapre);     Dpre = nanstd(Diapre);
Dpo = nanmean(Diapost);    Dpoe = nanstd(Diapost);
n = [1,2;3,4;5,6];

Titles = ["Vo2 Peak (mL/kg/min)","GIR (mL/kg-min)"];
dim = [.2 .5 .3 .3];
str = 'Ex trained:';
% Colors = ['b','b','r','r'];
for j = 1:2
    figure
    bar([1;2;3;4],[Cpr(j);Cpo(j);0;0],'b');
    hold on
    xticklabels({'-','+','-','+'})
    bar([3;4],[Dpr(j);Dpo(j)],'r')
    annotation('textbox',dim,'String',str,'FitBoxToText','on','LineStyle','none','FontWeight','bold');
    ylabel(Titles(j),'FontWeight','bold')
    xlabel('Control                    Diabetes','FontWeight','bold')
    for i = 1:length(Conpre(:,1))
        plot(n(1,:),[Conpre(i,j),Conpost(i,j)],'k-o')
    end
    for k = 1:length(Diapre(:,1))
        plot(n(2,:),[Diapre(k,j),Diapost(k,j)],'k-o')
    end
    errorbar([1;2;3;4],[Cpr(j);Cpo(j);Dpr(j);Dpo(j)],[Cpre(j);Cpoe(j);Dpre(j);Dpoe(j)],'k.','LineWidth',1.5,'Capsize',10)
    [~,pc(j)] = ttest(Conpre(:,j),Conpost(:,j));
    [~,pd(j)] = ttest(Diapre(:,j),Diapost(:,j));
end

%% Sphygmocor PWV - whole study
clc;
clearvars;
close all;



Path = 'E:\PhD\Danny_mixedVA\Demographics\VAMIXED-PWVComparison_DATA_2020-09-08_1027.csv';
[Dat,~,~] = xlsread(Path);

Conpre = Dat(1:56,3);    
Conpost = Dat(1:56,4);

Diapre = Dat(57:end,3);
Diapost = Dat(57:end,4);

Cpr = nanmean(Conpre);     Cpre = nanstd(Conpre);
Cpo = nanmean(Conpost);    Cpoe = nanstd(Conpost);
Dpr = nanmean(Diapre);     Dpre = nanstd(Diapre);
Dpo = nanmean(Diapost);    Dpoe = nanstd(Diapost);
n = [1,2;3,4;5,6];

Titles = ["Sphygmocor PWV (m/s)"];
dim = [.2 .5 .3 .3];
str = 'Ex trained:';
% Colors = ['b','b','r','r'];
for j = 1
    figure
    bar([1;2;3;4],[Cpr(j);Cpo(j);0;0],'b');
    hold on
    xticklabels({'-','+','-','+'})
    bar([3;4],[Dpr(j);Dpo(j)],'r')
    annotation('textbox',dim,'String',str,'FitBoxToText','on','LineStyle','none','FontWeight','bold');
    ylabel(Titles(j),'FontWeight','bold')
    xlabel('Control                    Diabetes','FontWeight','bold')
    for i = 1:length(Conpre(:,1))
        plot(n(1,:),[Conpre(i,j),Conpost(i,j)],'k-o')
    end
    for k = 1:length(Diapre(:,1))
        plot(n(2,:),[Diapre(k,j),Diapost(k,j)],'k-o')
    end
    errorbar([1;2;3;4],[Cpr(j);Cpo(j);Dpr(j);Dpo(j)],[Cpre(j);Cpoe(j);Dpre(j);Dpoe(j)],'k.','LineWidth',1.5,'Capsize',10)
    [~,pc(j)] = ttest(Conpre(:,j),Conpost(:,j));
    [~,pd(j)] = ttest(Diapre(:,j),Diapost(:,j));
end