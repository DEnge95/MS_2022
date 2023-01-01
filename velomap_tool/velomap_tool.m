%		velomap_tool
% 
%       Import and visualization of dicom data(PC-MRI-data)
%  
%       Modification of data (noise filtering,anti-aliasing,eddy current correction) 
%  
%       Conversion of modified data to EnSight, avi-movie or mrstruct.
%
%       Originally developed by Jelena Bock, University of Freiburg, Germany
%       Modified and updated by: Alex Barker, Michael Markl, Northwestern University, Chicago, USA
% 



function [mag_struct,flow_struct] = velomap_tool(config)



% % %%% return state structure if demanded
% % if nargin==1 & ischar(config) & strcmp(config,'struct')   
% %    ret = velomap_model('struct');
% %    return;      
% % end
% % %%% End of: return state structure if demanded



%%% create gui and init structures

fig = findobj('tag','velomap_figure');
if (nargin ~= 1) || ~isstr(config) || ~strcmp(config, 'block')
    if ~isempty(fig)
        delete(fig);
    end
    clear velomap_model;
    clear velomap_gui;
    fig= [];
end

if isempty(fig)
    warning off;
    fig=openfig('velomap_gui'); 
    velomap_model('init');
    %warning on;
end
%%% End of: create gui and init structures



%%% read config m-file or structure
if (nargin==1)
   
   if ischar(config)
      if strcmp(config, 'block')
      elseif exist(config)~=2
         %warning('Config file does not exist');
      else
         eval(config);
      end
   elseif isstruct(config)
      velomap_model('set','all',config);      

   else
      warning('Input argument is not a string nor a structure. Will be ignored');
   end
   
end

%%% End of: read config m-file or structure




%%% read data structure and return value of interesst if expected
if nargout~=0 || ((nargin == 1) && isstr(config) && strcmp(config, 'returnstruct'))
   velomap_model('set','selStruct');
   waitfor(fig);
   mag_struct  = [];
   flow_struct = [];
   ret  = velomap_model('return');
   if ~isempty(ret)
    mag_struct   = mrstruct_read(ret(1,:));
    flow_struct  = mrstruct_read(ret(2,:));
   end
   
end

%%% End of: read data structure and return value of interesst if expected
