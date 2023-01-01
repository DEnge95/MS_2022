 close all,clear all,clc

for i = 1:length(dir("F:\Scans"))
    file = uigetdir('F:\Scans\'); %Select 2D PC-MRI scan folder
    dire = dir(file);
    label = dire(1).folder; %#ok<*SAGROW>
    n = 1;
    for j = 4:3:13
        info = dicominfo(strcat(file,'\',dire(4).name));
        TRa(n) = info.TriggerTime/(j-3); 
        n = n+1;
        TR = mean(TRa);
    end
    output(i,:) = [string(info.Filename), TR]
end

info = dicominfo