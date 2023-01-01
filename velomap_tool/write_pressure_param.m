function write_pressure_param(dataStruct)

%%% build cell array to write pressure param in existing excel document

excelFileStr = sprintf('%s%s%s',dataStruct.dataDirectory,filesep,'functions_settings.xls');
outExcel = cell(40,6);
numTF=dataStruct.numPhases;

outExcel(1,:)   ={'Static Pressure Variables','','','time frame','num of iterations steps','error'}; 
%%% write static pressure variables
outExcel(2,1:2) ={'Viscosity',dataStruct.pdStruct.visc}; 
outExcel(3,1:2) ={'Density',dataStruct.pdStruct.dens}; 
outExcel(4,1:2) ={'Max Iter',dataStruct.pdStruct.max_iter}; 
outExcel(5,1:2) ={'Alpha',dataStruct.pdStruct.alpha}; 
outExcel(6,1:2) ={'Poly Num',dataStruct.pdStruct.poly_num}; 
outExcel(7,1:2) ={'Max error',dataStruct.pdStruct.max_error};     
outExcel(8,1:2) ={'DelX',dataStruct.pdStruct.delX};     
outExcel(9,1:2) ={'DelY',dataStruct.pdStruct.delY};     
outExcel(10,1:2)={'DelZ',dataStruct.pdStruct.delZ};  
outExcel(11,1:2)={'Tres',dataStruct.pdStruct.tres}; 
outExcel(12,1:2)={'Num Pts',dataStruct.pdStruct.npts}; 


for tf=1:numTF
    outExcel(tf+1,4) = {tf};
    outExcel(tf+1,5) = {dataStruct.pdStruct.errorAy(tf,1)};
    outExcel(tf+1,6) = {dataStruct.pdStruct.errorAy(tf,2)};
end

warning off MATLAB:xlswrite:AddSheet;
status = xlswrite(excelFileStr, outExcel, 'functions settings', 'E1');