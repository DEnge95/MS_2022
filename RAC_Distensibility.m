clc;
clearvars;
close all;


Files = dir ("*.txt");
for i = 1:length(Files)

    A = importfile1(Files(i).name,84,134);
    Asc(:,i) = A.VarName7(2:end);
    D = importfile(Files(i).name,178,228);
    Desc(:,i) = D.VarName7(2:end);
    %RACa(i) = (max(Asc(2:end)))-min(Asc(2:end))/max(Asc(2:end));
end

