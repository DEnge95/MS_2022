function var =load_into_variable(fileName)
% load from path directly into the desired variable
tmp     = load(fileName);            
names   = fieldnames(tmp);
var     = getfield(tmp,names{1}); %#ok<GFLD>