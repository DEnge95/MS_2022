
% see preprocessing_tool for a full description of the tool.
%
%		Jelena Bock
%		2/06

function ret = velomap_model(commandStr,entryStr,dataValue)


%%%%% declare persistent structures and variables
persistent handleStruct;
persistent dataStruct;
persistent velomap_semaphore;
%%%%% End of: declare persistent structures and variables


%%%%% init unused arguments
if nargin<3
   dataValue = [];
end
if nargin<2
   entryStr = '';
end
if nargin==0
   commandStr = 'init';
end
%%%%% End of: init unused arguments


%%%%% is application still busy ?
if strcmp(commandStr,'gui') && strcmp(entryStr,'gui_destroy') 
   velomap_semaphore=0;
end
if ~isempty(velomap_semaphore) && velomap_semaphore==1 
   msgbox('Application is still busy. Please retry later ....');
   return;
else
   velomap_semaphore=1;
end
%%%%% End of: is application still busy ?


%%%%% init structure of figure handles and data items
if strcmp(commandStr,'init')
   try
      handleStruct = local_init_handleStruct;
      %% initializing tab groups
      TabGroup    = 'TabGroupPreprocessing';
      PanelNames  = {'noise_masking','pcmra','eddy_current_correction','anti_aliasing','pressure_diff','nonuniformity'};
      handleStruct = tabpanelfcn('make_groups',TabGroup, PanelNames, handleStruct, 1);
      
      TabGroup    = 'TabGroupConversion';
      PanelNames  = {'ensight','mrstruct'};
      handleStruct = tabpanelfcn('make_groups',TabGroup, PanelNames, handleStruct, 1);
      %% end of: initializing tab groups
      dataStruct   = local_init_dataStruct;
   catch
      warndlg(lasterr,'An error occured while initializing velomap_model');
   end
%%%%% End of: init structure of figure handles and data items


%%%%% return dataStruct
elseif strcmp(commandStr,'struct')
   try
      if isempty(dataStruct)
         ret = local_init_dataStruct;
      else
         ret = dataStruct;
      end
   catch
      warndlg(lasterr,'An error occured while reading structure of velomap_model');
   end
%%%%% End of: return dataStruct


%%%%% read internal data
elseif strcmp(commandStr,'get')
   try
      ret = local_get_dataStruct_entry(dataStruct,entryStr);
   catch
      warndlg(lasterr,'An error occured while reading an entry in velomap_model');
   end
%%%%% End of: read internal data


%%%%% set interface data   
elseif strcmp(commandStr,'set')
   try
      dataStruct = local_set_interface_entry(handleStruct,dataStruct,entryStr,dataValue);
   catch
      warndlg(lasterr,'An error occured while setting an entry in velomap_model');
   end
%%%%% End of: set interface data


%%%%% calls from gui
elseif strcmp(commandStr,'gui')
%   try
      [dataStruct,handleStruct] = local_handle_gui_event(handleStruct,dataStruct,entryStr);
%    catch
%       warndlg(lasterr,'An error occured while processing gui events for velomap_model');
%    end
%%%%% End of: calls from gui


%%%%% return internal data to main tool, empty persistent variables
elseif strcmp(commandStr,'return')
   try       
      ret = dataStruct.returnPath;
      dataStruct   =[];
      handleStruct =[];
   catch
      warndlg(lasterr,'An error occured while returning results from velomap_model');
   end
%%%%% End of: return internal data to main tool, empty persistent variables


%%%%% detect invalid commandStr
else 
   warndlg(lasterr,'Sorry, command switch was not recognized');
end
%%%%% End of: detect invalid commandStr


%%%%% reset semaphore
 velomap_semaphore=0;
%%%%% End of: reset semaphore


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  local functions  (main level)                                            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 
%%%%% retrieve and save handles of gui controls
function  handleStruct = local_init_handleStruct


handleStruct.figure                 = findobj('tag','velomap_figure'); 
handleStruct.axesMag                = findobj('tag','velomap_axes_mag');
handleStruct.axesFlowUp             = findobj('tag','velomap_axes_vx');
handleStruct.axesFlowDown           = findobj('tag','velomap_axes_vy');
handleStruct.axesFlowBig            = findobj('tag','velomap_axes_vz');
handleStruct.axesFlowUpMod          = findobj('tag','velomap_axes_vx_new');
handleStruct.axesFlowDownMod        = findobj('tag','velomap_axes_vy_new');
handleStruct.axesFlowBigMod         = findobj('tag','velomap_axes_vz_new');
handleStruct.changeViewPop          = findobj('tag','velomap_changeViewPop');
handleStruct.patientTxt             = findobj('tag','velomap_patientTxt'); 
handleStruct.fileTypePop            = findobj('tag','velomap_fileTypePop'); 
handleStruct.encodingTypePop        = findobj('tag','velomap_encodingTypePop');
handleStruct.userID                 = findobj('tag','velomap_userID');
handleStruct.encodingTypeString     = findobj('tag','velomap_encodingTypeString');
handleStruct.openMRDataBut          = findobj('tag','velomap_openMRDataBut'); 
handleStruct.closeBut               = findobj('tag','velomap_closeBut'); 
handleStruct.saveImageBut           = findobj('tag','velomap_saveImageBut'); 
handleStruct.statusEdt              = findobj('tag','velomap_statusEdt');
handleStruct.SlicesNumSlide         = findobj('tag','velomap_SlicesNumSlide');
handleStruct.SlicesNumEdit          = findobj('tag','velomap_SlicesNumEdit');
handleStruct.PhasesNumSlide         = findobj('tag','velomap_PhasesNumSlide');
handleStruct.PhasesNumEdit          = findobj('tag','velomap_PhasesNumEdit');
handleStruct.VorigBigTxt            = findobj('tag','velomap_V1origTxt');
handleStruct.VorigUpTxt             = findobj('tag','velomap_V2origTxt');
handleStruct.VorigDownTxt           = findobj('tag','velomap_V3origTxt');
handleStruct.VmodBigTxt             = findobj('tag','velomap_V1modTxt');
handleStruct.VmodUpTxt              = findobj('tag','velomap_V2modTxt');
handleStruct.VmodDownTxt            = findobj('tag','velomap_V3modTxt');
handleStruct.shutdownCheck          = findobj('tag','velomap_shutdownCheck');

%%%  -------------------------------------------------------------------------
%%% handles for preprocesing functions
handleStruct.TabGroupPreprocessing  = findobj('tag','velomap_TabGroupPreprocessing');
handleStruct.last_tab   = [];
%%% handles for non-uniformity correction
handleStruct.nonuniformityCheck    = findobj('tag','velomap_nonuniformityCheck');
handleStruct.backgroundSlide       = findobj('tag','velomap_backgroundSlide');
handleStruct.backgroundEdt         = findobj('tag','velomap_backgroundEdt');
handleStruct.previewBackgroundBut  = findobj('tag','velomap_previewBackgroundBut');
handleStruct.nonuniformityTab       = findobj('String','NON-UNIFORMITY');
handleStruct.backgroundThresholdCaption = findobj('tag','velomap_backgroundThresholdCaption');
%%% set user data for noise masking  object
%set(handleStruct.nonuniformityCheck,    'UserData', 'nonuniformity');
set(handleStruct.backgroundSlide,       'UserData', 'nonuniformity');
set(handleStruct.backgroundEdt,         'UserData', 'nonuniformity');
set(handleStruct.previewBackgroundBut,  'UserData', 'nonuniformity');
set(handleStruct.nonuniformityTab,       'UserData', 'nonuniformity');
set(handleStruct.backgroundThresholdCaption,'UserData', 'nonuniformity');

%%% find noise masking objects 
handleStruct.filterCheck            = findobj('tag','velomap_filterCheck');
handleStruct.radioANDBut            = findobj('tag','velomap_radioANDBut');
handleStruct.radioORBut             = findobj('tag','velomap_radioORBut');
handleStruct.radioButGroup          = findobj('tag','velomap_radioButGroup');
handleStruct.noiseValueEdt          = findobj('tag','velomap_noiseValueEdt');
handleStruct.noiseValueSlide        = findobj('tag','velomap_noiseValueSlide');
handleStruct.stdevNoiseValueEdt     = findobj('tag','velomap_stdevNoiseValueEdt');
handleStruct.stdevNoiseValueSlide   = findobj('tag','velomap_stdevNoiseValueSlide');
handleStruct.derivatNoiseValueEdt   = findobj('tag','velomap_derivatNoiseValueEdt');
handleStruct.derivatNoiseValueSlide = findobj('tag','velomap_derivatNoiseValueSlide');
handleStruct.chooseFilterList       = findobj('tag','velomap_chooseFilterList');
handleStruct.noiseMaskingTab        = findobj('String','NOISE MASKING');
handleStruct.staticTissueCheck      = findobj('tag','velomap_staticTissueCheck');
handleStruct.derivFilterThresholdCaption   = findobj('tag','velomap_derivFilterThresholdCaption');
handleStruct.noiseFilterThresholdCaption   = findobj('tag','velomap_noiseFilterThresholdCaption');
handleStruct.stdevFilterThresholdCaption   = findobj('tag','velomap_stdevFilterThresholdCaption');

%%% set user data for noise masking  objects
set(handleStruct.radioANDBut,                   'UserData', 'noise_masking');
set(handleStruct.radioORBut,                    'UserData', 'noise_masking');
set(handleStruct.radioButGroup,                 'UserData', 'noise_masking');
set(handleStruct.noiseValueEdt,                 'UserData', 'noise_masking');
set(handleStruct.noiseValueSlide,               'UserData', 'noise_masking');
set(handleStruct.stdevNoiseValueEdt,            'UserData', 'noise_masking');
set(handleStruct.stdevNoiseValueSlide,          'UserData', 'noise_masking');
set(handleStruct.derivatNoiseValueEdt,          'UserData', 'noise_masking');
set(handleStruct.derivatNoiseValueSlide,        'UserData', 'noise_masking');
set(handleStruct.chooseFilterList,              'UserData', 'noise_masking');
set(handleStruct.noiseMaskingTab,               'UserData', 'noise_masking');
set(handleStruct.staticTissueCheck,             'UserData', 'noise_masking');
set(handleStruct.derivFilterThresholdCaption,   'UserData', 'noise_masking');
set(handleStruct.noiseFilterThresholdCaption,   'UserData', 'noise_masking');
set(handleStruct.stdevFilterThresholdCaption,   'UserData', 'noise_masking');

%%% find objects for eddy current correction
handleStruct.eddyCurrentCheck         = findobj('tag','velomap_eddyCurrentCheck');
handleStruct.secOrderCorrCheck        = findobj('tag','velomap_secOrderCorrCheck');
handleStruct.staticValueEdt           = findobj('tag','velomap_staticValueEdt');
handleStruct.staticValueSlide         = findobj('tag','velomap_staticValueSlide');
handleStruct.eddyCurrentCorrectionTab = findobj('String','EDDY CURRENT CORR.');
handleStruct.previewStatRegBut        = findobj('tag','velomap_previewStatRegBut');
%%% set user data for eddy current correction objects
set(handleStruct.staticValueEdt,            'UserData', 'eddy_current_correction');
set(handleStruct.staticValueSlide,          'UserData', 'eddy_current_correction');
set(handleStruct.previewStatRegBut,         'UserData', 'eddy_current_correction');
set(handleStruct.eddyCurrentCorrectionTab,  'UserData', 'eddy_current_correction');

%%% handles for calculating & previewing of pc-mra
handleStruct.pcmraCheck             = findobj('tag','velomap_pcmraCheck');
handleStruct.pcmraSumSquaresCheck   = findobj('tag','velomap_pcmraSumSquaresCheck');
handleStruct.pcmraMeanAbsVelCheck   = findobj('tag','velomap_pcmraMeanAbsVelCheck');
handleStruct.pcmraPseudoComplDiffCheck = findobj('tag','velomap_pcmraPseudoComplDiffCheck');
handleStruct.pcmraTimeAveMagCheck   = findobj('tag','velomap_pcmraTimeAveMagCheck');
handleStruct.pcmraSqrtSumSquaresCheck   = findobj('tag','velomap_pcmraSqrtSumSquaresCheck'); %%//LiliMa: August 2017
handleStruct.pcmraVelCorrCheck      = findobj('tag','velomap_pcmraVelCorrCheck'); %%//MBS
handleStruct.pcmraMedianFilterCheck = findobj('tag','velomap_pcmraMedianFilterCheck');  %% //MM
handleStruct.pcmraHistEqualCheck    = findobj('tag','velomap_pcmraHistEqualCheck');     %% //MM
handleStruct.previewPCmraBut        = findobj('tag','velomap_previewPCmraBut');
handleStruct.vencThresholdEdit      = findobj('tag','velomap_vencThresholdEdit');
handleStruct.timeFramesMinEdit = findobj('tag','velomap_timeFramesMinEdit');  %% //MM
handleStruct.timeFramesMaxEdit = findobj('tag','velomap_timeFramesMaxEdit');  %% //MM
handleStruct.vencThresholdCaption   = findobj('tag','velomap_vencThresholdCaption');
handleStruct.pcmraTab               = findobj('String','PCMRA');
handleStruct.timeFramesPCMRATxt     = findobj('tag','velomap_timeFramesPCMRA');
handleStruct.timeFramesMinTxt       = findobj('tag','velomap_timeFramesMin');
handleStruct.timeFramesMaxTxt       = findobj('tag','velomap_timeFramesMax');
%%% set user data for calculating & previewing of pc-mra
set(handleStruct.pcmraSumSquaresCheck,      'UserData', 'pcmra', 'Value', 1);
set(handleStruct.pcmraMeanAbsVelCheck,      'UserData', 'pcmra', 'Value', 0);
set(handleStruct.pcmraPseudoComplDiffCheck, 'UserData', 'pcmra', 'Value', 0);
set(handleStruct.pcmraTimeAveMagCheck,      'UserData', 'pcmra', 'Value', 0);
set(handleStruct.pcmraSqrtSumSquaresCheck,  'UserData', 'pcmra', 'Value', 0); %% LiliMa
set(handleStruct.pcmraVelCorrCheck,         'UserData', 'pcmra', 'Value', 0); %% MBS
set(handleStruct.pcmraMedianFilterCheck,    'UserData', 'pcmra'); %% //MM
set(handleStruct.pcmraHistEqualCheck,       'UserData', 'pcmra'); %% //MM
set(handleStruct.previewPCmraBut,           'UserData', 'pcmra');
set(handleStruct.vencThresholdEdit,         'UserData', 'pcmra');
set(handleStruct.timeFramesMinEdit,         'UserData', 'pcmra'); %% //MM
set(handleStruct.timeFramesMaxEdit,         'UserData', 'pcmra'); %% //MM 
set(handleStruct.vencThresholdCaption,      'UserData', 'pcmra');
set(handleStruct.pcmraTab,                  'UserData', 'pcmra');
set(handleStruct.timeFramesPCMRATxt,        'UserData', 'pcmra');
set(handleStruct.timeFramesMinTxt,          'UserData', 'pcmra');
set(handleStruct.timeFramesMaxTxt,          'UserData', 'pcmra');



%%% handles for performing of anti-aliasing
handleStruct.antiAliasCheck         = findobj('tag','velomap_antiAliasCheck');
handleStruct.unwrapManualBut        = findobj('tag','velomap_unwrapManualBut');
handleStruct.loadUnwrapDataBut      = findobj('tag','velomap_loadUnwrapDataBut');
handleStruct.numberOfIterPop        = findobj('tag','velomap_numIterPop');
handleStruct.numIterationsCaption   = findobj('tag','velomap_numIterationsCaption');
handleStruct.antiAliasingTab        = findobj('String','ANTI-ALIASING');
%%% set user data for performing of anti-aliasing
set(handleStruct.unwrapManualBut,      'UserData', 'anti_aliasing');
set(handleStruct.loadUnwrapDataBut,    'UserData', 'anti_aliasing');
set(handleStruct.numIterationsCaption, 'UserData', 'anti_aliasing');
set(handleStruct.numberOfIterPop,      'UserData', 'anti_aliasing');
set(handleStruct.antiAliasingTab,      'UserData', 'anti_aliasing');

%%% handles for calculation of pressure difference
handleStruct.pdCheck                = findobj('tag','velomap_pdCheck');
handleStruct.previewPDMaskBut       = findobj('tag','velomap_previewPDMaskBut');
handleStruct.thresholdPDmaskEdit    = findobj('tag','velomap_thresholdPDmaskEdit');
handleStruct.pressureDiffTab        = findobj('String','PRESSURE DIFF.');
%%% set user data for performing of  pressure difference
set(handleStruct.previewPDMaskBut,     'UserData', 'pressure_diff');
set(handleStruct.thresholdPDmaskEdit,  'UserData', 'pressure_diff');
set(handleStruct.pressureDiffTab,      'UserData', 'pressure_diff');


%%%------------------------------------------------------------------------

%%%------------------------------------------------------------------------
%%% conversion tab group
handleStruct.TabGroupConversion     = findobj('tag','velomap_TabGroupConversion');

%%% handles for ensight conversion
handleStruct.enSightCheck           = findobj('tag','velomap_enSightCheck');
handleStruct.appAnatomDataBut       = findobj('tag','velomap_appAnatomDataBut');
handleStruct.deleteAnatomDataBut    = findobj('tag','velomap_deleteAnatomDataBut');
handleStruct.AnatomDataList         = findobj('tag','velomap_AnatomDataList');
handleStruct.ensightTab             = findobj('String','ENSIGHT');
%%% set user data for ensight conversion
set(handleStruct.appAnatomDataBut,      'UserData', 'ensight');
set(handleStruct.deleteAnatomDataBut,   'UserData', 'ensight');
set(handleStruct.AnatomDataList,        'UserData', 'ensight');
set(handleStruct.ensightTab,            'UserData', 'ensight');

%%% handles for mrstruct conversion
handleStruct.selMRstructCheck       = findobj('tag','velomap_selMRstructCheck');
handleStruct.lowlimSlicesEdt        = findobj('tag','velomap_lowlimSlicesEdt');
handleStruct.uplimSlicesEdt         = findobj('tag','velomap_uplimSlicesEdt');
handleStruct.lowlimPhasesEdt        = findobj('tag','velomap_lowlimPhasesEdt');
handleStruct.uplimPhasesEdt         = findobj('tag','velomap_uplimPhasesEdt');
handleStruct.uplimSlicesTxt         = findobj('tag','velomap_uplimSlicesTxt');
handleStruct.uplimPhasesTxt         = findobj('tag','velomap_uplimPhasesTxt');
handleStruct.mrstructTxt            = findobj('tag','velomap_mrstructTxt');
handleStruct.mrstructTab            = findobj('String','MRSTRUCT');
%%% set user data for mrstruct conversion
set(handleStruct.lowlimSlicesEdt,   'UserData', 'mrstruct');
set(handleStruct.uplimSlicesEdt,    'UserData', 'mrstruct');
set(handleStruct.lowlimPhasesEdt,   'UserData', 'mrstruct');
set(handleStruct.uplimPhasesEdt,    'UserData', 'mrstruct');
set(handleStruct.uplimPhasesTxt,    'UserData', 'mrstruct');
set(handleStruct.uplimSlicesTxt,    'UserData', 'mrstruct');
set(handleStruct.mrstructTxt,       'UserData', 'mrstruct');
set(handleStruct.mrstructTab,       'UserData', 'mrstruct');
%%%------------------------------------------------------------------------

handleStruct.aviMovieCheck          = findobj('tag','velomap_aviMovieCheck');

handleStruct.pcmraToDICOMCheck      = findobj('tag','velomap_pcmraToDICOMCheck');
% 20190729 Takashi added
handleStruct.autoseg                = findobj('tag','velomap_autoseg');
handleStruct.site                   = findobj('tag','velomap_site');
handleStruct.autosegvers            = findobj('tag','velomap_autosegvers');

handleStruct.orderList              = findobj('tag','velomap_orderList');
handleStruct.orderUpBut             = findobj('tag','velomap_orderUpBut');
handleStruct.orderDownBut           = findobj('tag','velomap_orderDownBut');
handleStruct.vencThPlaneEdt         = findobj('tag','velomap_vencThPlaneEdt');
handleStruct.vencInPlaneXEdt        = findobj('tag','velomap_vencInPlaneXEdt');
handleStruct.vencInPlaneYEdt        = findobj('tag','velomap_vencInPlaneYEdt');
handleStruct.dataInfoEdt            = findobj('tag','velomap_dataInfoEdt');

handleStruct.queueList              = findobj('tag','velomap_queueList');

handleStruct.loadDataBut            = findobj('tag','velomap_loadDataBut');
handleStruct.previewBut             = findobj('tag','velomap_previewBut');
handleStruct.applyBut               = findobj('tag','velomap_applyBut');

handleStruct.outToExcelBut          = findobj('tag','velomap_outToExcelBut');
handleStruct.saveFuncSetBut         = findobj('tag','velomap_saveFuncSetBut');
handleStruct.takeQueueBut           = findobj('tag','velomap_takeQueueBut');
handleStruct.loadInQueueBut         = findobj('tag','velomap_loadInQueueBut');
handleStruct.startQueueBut          = findobj('tag','velomap_startQueueBut');

% freiburgLogo = imread('cvmri_logo_blacklines_lowres','jpg'); 
%ima4DFlow = imread('4dflow_example.jpg');
%CvCTIlogo = imread('Cardiovascular-Color.jpg');
handleStruct.imageMag               = image([0],'parent', handleStruct.axesMag);
handleStruct.imageFlowBig           = image([0],'parent', handleStruct.axesFlowBig);
handleStruct.imageFlowUp            = image([0],'parent', handleStruct.axesFlowUp);
handleStruct.imageFlowDown          = image([0],'parent', handleStruct.axesFlowDown);
handleStruct.imageFlowBigMod        = image([0],'parent', handleStruct.axesFlowBigMod);
handleStruct.imageFlowUpMod         = image([0],'parent', handleStruct.axesFlowUpMod);
handleStruct.imageFlowDownMod       = image([0],'parent', handleStruct.axesFlowDownMod);

%settings
axis(handleStruct.axesMag, 'image', 'off')
colormap(handleStruct.axesMag, gray(256));
set(handleStruct.axesMag, 'XLimMode', 'auto', 'YLimMode', 'auto', 'ydir', 'normal') % here is where loading logos will act funny (AJB)
set(handleStruct.imageMag, 'CDataMapping', 'scaled')

axis(handleStruct.axesFlowUp, 'image', 'off')
colormap(handleStruct.axesFlowUp, gray(256));
set(handleStruct.axesFlowUp, 'XLimMode', 'auto', 'YLimMode', 'auto', 'ydir', 'normal')
set(handleStruct.imageFlowUp, 'CDataMapping', 'scaled')

axis(handleStruct.axesFlowDown, 'image', 'off')
colormap(handleStruct.axesFlowDown, gray(256));
set(handleStruct.axesFlowDown, 'XLimMode', 'auto', 'YLimMode', 'auto', 'ydir', 'normal')
set(handleStruct.imageFlowDown, 'CDataMapping', 'scaled')

axis(handleStruct.axesFlowBig, 'image', 'off')
colormap(handleStruct.axesFlowBig, gray(256));
set(handleStruct.axesFlowBig, 'XLimMode', 'auto', 'YLimMode', 'auto', 'ydir', 'normal')
set(handleStruct.imageFlowBig, 'CDataMapping', 'scaled')

axis(handleStruct.axesFlowUpMod, 'image', 'off')
colormap(handleStruct.axesFlowUpMod, gray(256));
set(handleStruct.axesFlowUpMod, 'XLimMode', 'auto', 'YLimMode', 'auto', 'ydir', 'normal')
set(handleStruct.imageFlowUpMod, 'CDataMapping', 'scaled')

axis(handleStruct.axesFlowDownMod, 'image', 'off')
colormap(handleStruct.axesFlowDownMod, gray(256));
set(handleStruct.axesFlowDownMod, 'XLimMode', 'auto', 'YLimMode', 'auto', 'ydir', 'normal')
set(handleStruct.imageFlowDownMod, 'CDataMapping', 'scaled')

axis(handleStruct.axesFlowBigMod, 'image', 'off')
colormap(handleStruct.axesFlowBigMod, gray(256));
set(handleStruct.axesFlowBigMod, 'XLimMode', 'auto', 'YLimMode', 'auto', 'ydir', 'normal')
set(handleStruct.imageFlowBigMod, 'CDataMapping', 'scaled');

% Have to handle colorbar different in the new 2014b(8.4) version of matlab (no longer an axes handle object)
v        = ver('matlab');              % get matlab version
isnewver = str2double(v.Version)>=8.4; % test if pesky version of matlab
if isnewver % treat colorbar as its own entity
    handleStruct.colorbarAx = colorbar;
    handleStruct.colorbarAx = colorbar('peer', handleStruct.axesFlowBigMod); % here is where colorbar was getting jacked
    set(handleStruct.colorbarAx,'position',[0.739344262295082 0.49398907103825135 0.00901639344262295 0.30273224043715846]);
else % treat colorbar as a valid axes handle
%     handleStruct.colorbarAx = colorbar('peer', handleStruct.axesFlowBigMod);
%     set(handleStruct.colorbarAx,'position',[0.739344262295082 0.49398907103825135 0.00901639344262295 0.30273224043715846])
   
end

set(handleStruct.figure, 'Name', 'Tool for PreProcessing & Converting of 4D Flow MRI Data - (c) Freiburg University, Germany & Northwestern University, Chicago, USA')

%% setup the callbacks according to MR physics structure 
set(handleStruct.encodingTypePop,'Callback','velomap_model(''gui'',''encodingTypePop'')');
set(handleStruct.orderUpBut,'Callback','velomap_model(''gui'',''orderUpBut'')');
set(handleStruct.orderDownBut,'Callback','velomap_model(''gui'',''orderDownBut'')');
set(handleStruct.loadUnwrapDataBut,'Callback','velomap_model(''gui'',''loadUnwrapDataBut'')');
set(handleStruct.orderList,'KeyPressFcn','velomap_model(''gui'',''orderList'')');

%%%%% End of: retrieve and save handles of gui controls
%%%%%

%% 
%%%%% init dataStruct
function  dataStruct = local_init_dataStruct
dataStruct.SiemensFlag      = []; %% LiliMa
dataStruct.TimeStamps      = []; %% LiliMa

dataStruct.mrStruct         = [];
dataStruct.status           = 'Tool is initialized!';
dataStruct.colormap         = 'gray';

dataStruct.dataMag          = []; %% //MM
dataStruct.dataFlow         = []; %% //MM
dataStruct.fileNamesMag     = [];
dataStruct.fileNamesFlow    = [];
dataStruct.imageSize        = [];
dataStruct.signVijk         = [];
dataStruct.numSlices        = [];
dataStruct.numPhases        = [];
dataStruct.nFE              = [];
dataStruct.dataDirectory    = [];
dataStruct.subfoldersFlow   = [];%LiliMa
dataStruct.checkedArray     = [];
dataStruct.stdFlow          = [];
dataStruct.stdFlowMax       = [];
dataStruct.imaFlowOrig      = [];
dataStruct.imaFLOW          = [];
dataStruct.imaMAG           = [];
dataStruct.backgroundMask   = [];
dataStruct.previewBackground =0;
%%% pcmra data 
dataStruct.pcmraSumSquares  = [];
dataStruct.pcmraMeanAbsVel  = [];
dataStruct.pcmraPseudoComplDiff =[];
dataStruct.pcmraSqrtSumSquares  = []; %% LiliMa
dataStruct.pcmraVelCorr  = []; %% MBS
dataStruct.pcmraMask        =[];
dataStruct.volume_pcmraSumSquares  = [];
dataStruct.volume_pcmraMeanAbsVel  = [];
dataStruct.volume_pcmraPseudoComplDiff =[];
dataStruct.volume_pcmraTimeAveMag =[];
dataStruct.volume_pcmraSqrtSumSquares  = []; %% LiliMa
dataStruct.volume_pcmraVelCorr  = []; %% MBS
%%% end of pcmra data
dataStruct.meanMag          = [];
dataStruct.filterMask       = [];
dataStruct.statMask         = [];
dataStruct.statRegStatus    = 'off';
dataStruct.pcmraStatus      = 'off';
dataStruct.previewStatus    = 'off';
dataStruct.loadStatus       = 'cancelled';
dataStruct.venc             = [];
dataStruct.phaseRange       = 4096;
dataStruct.TR               = [];
dataStruct.ensightDirPath   = [];
dataStruct.mrstructDirPath  = [];
dataStruct.dataPathName     = [];
dataStruct.propMRstruct     = [];
dataStruct.actMRstruct      = [];
dataStruct.patientName      = '';
dataStruct.orientation      = '';
dataStruct.peDir            = '';
dataStruct.dataInfo         = [];
dataStruct.returnPath       = [];
dataStruct.inifilePath      = [];
dataStruct.unwrapData       = [];
dataStruct.swVersion        = '';
dataStruct.encodingType     = 'velocity';
dataStruct.frameRange       = [];
dataStruct.pcmraMfilterq    = [];
dataStruct.pcmraEqu         = [];
dataStruct.userIDStr    	= '';
dataStruct.userSubjectIDStr = '';
dataStruct.serverPath       = "//childrens/clinical/AIL/ail_group/software_and_documentation/autoseg"; 

%% MM: pressure difference calculations
%% init data structure
dataStruct.resX       = [];
dataStruct.resY       = [];
dataStruct.resZ       = [];

dataStruct.pdStruct.visc =  [];         %% Viscosity in Cp divided by 1000 to get SI units
dataStruct.pdStruct.dens = [];          %% Fluid density in kg/m^3
dataStruct.pdStruct.max_iter = [];      %% Max pressure iterations
dataStruct.pdStruct.alpha = [];         %% Under relaxation factor for iteration
dataStruct.pdStruct.poly_num = [];      %% Order of polynomial fitting
dataStruct.pdStruct.max_error = [];     %% Stopping criteria for Iterations
dataStruct.pdStruct.delX=[];            %% Grid spacing in x (mm)
dataStruct.pdStruct.delY=[];            %% Grid spacing in y (mm)
dataStruct.pdStruct.delZ=[];            %% Grid spacing in z (mm)
dataStruct.pdStruct.tres=[];            %% Time spacing in (ms)
dataStruct.pdStruct.plist = [];         %% Points to include in calcs
dataStruct.pdStruct.npts  = [];         %% Number of points
dataStruct.pdStruct.VELx = [];          %% 4D Velocity in first index
dataStruct.pdStruct.VELy = [];          %% 4D Velocity in second index
dataStruct.pdStruct.VELz = [];          %% 4D Velocity in third index
dataStruct.pdStruct.MASK = [];          %% 3D Mask Structure
dataStruct.pdStruct.GRADx = [];         %% Vectorized Gradient in first index
dataStruct.pdStruct.GRADy = [];         %% Vectorized Gradient in second index
dataStruct.pdStruct.GRADz = [];         %% Vectorized Gradient in third index
dataStruct.pdStruct.DIM = [];           %% Size of data matrix
dataStruct.pdStruct.verbose=[];         %% 1 implies verbose
dataStruct.pdStruct.tvals=[];           %% Number of time points to be use for calc
dataStruct.pdStruct.PRESSURE=[];        %% Relative pressure
dataStruct.pdStruct.xlist=[];           %% Point List in X of indexed values
dataStruct.pdStruct.ylist=[];           %% Point List in Y of indexed values
dataStruct.pdStruct.zlist=[];           %% Point List in Z of indexed values
dataStruct.pdStruct.nbs=[];             %% Details existance of neighbors 
dataStruct.pdStruct.PRESSURE_MAT=[];    %% Actual Matrix
dataStruct.pdStruct.tvals=[];           %% Time points used
dataStruct.pdStruct.GRADx_MAT=[];       %% Matrix Storage for Gradient in X
dataStruct.pdStruct.GRADy_MAT=[];       %% Matrix Storage for Gradient in Y
dataStruct.pdStruct.GRADz_MAT=[];       %% Matrix Storage for Gradient in Z
dataStruct.pdStruct.GRAD_OUT=[];        %% Whether or not to assign gradient

dataStruct.pdStruct.flood_points=[];
dataStruct.pdStruct.errorAy=[];         %% Error calculation for certain time frame stopped
dataStruct.pdStruct.tmp=[];
dataStruct.externflag = 1; %1== for extern users; 0 for intern users
%%%%% End of: init dataStruct
%%%%% 


%% read entry from dataStruct
function  dataEntry = local_get_dataStruct_entry(dataStruct,entryStr)

if strcmp(entryStr,'mrStruct')
   dataEntry = dataStruct.mrStruct;
elseif strcmp(entryStr,'pos')
   dataEntry = dataStruct.currentPos;
elseif strcmp(entryStr,'colormap')
   dataEntry = dataStruct.colormap; 
else
   warndlg(lasterr,'Sorry, entry switch was not recognized, no action performed');
end
%%%%% End of: read entry from dataStruct
%%%%% 


%% 
%%%%% set public interface value
function  dataStruct = local_set_interface_entry(handleStruct,dataStruct,entryStr,dataValue)

if strcmp(entryStr,'all')
   dataStruct = velomap_internal(handleStruct, dataStruct, 'database', 'all', dataValue);
   velomap_internal(handleStruct, dataStruct, 'gui', 'init',  []);
   velomap_internal(handleStruct, dataStruct, 'gui', 'image',  []);
elseif strcmp(entryStr,'mrStruct')
   dataStruct = velomap_internal(handleStruct, dataStruct, 'database', 'mrStruct', dataValue);
   velomap_internal(handleStruct, dataStruct, 'gui', 'init',  []);
   velomap_internal(handleStruct, dataStruct, 'gui', 'image',  []);
elseif strcmp(entryStr, 'chColormap')
    dataStruct = velomap_internal(handleStruct, dataStruct, 'database', 'chColormap', dataValue);
    velomap_internal(handleStruct, dataStruct, 'gui', 'image',  []); 
elseif strcmp(entryStr,'selStruct')
    %set(handleStruct.selMRstructCheck,'Value',1);   
else
   warndlg(lasterr,'Sorry, entry switch was not recognized, no action performed');
end
%%%%% End of: set public interface value
%%%%% 


%% 
%%%%% catch gui events
function  [dataStruct,handleStruct] = local_handle_gui_event(handleStruct,dataStruct,entryStr)

%%% load data button or a key was pressed
if strcmp(entryStr,'loadDataBut')%||strcmp(entryStr,'loadDataKeyPress')
    %% get id of the pressed button
    if strcmp(entryStr,'loadDataKeyPress')        
        keyID = double(get(handleStruct.figure,'CurrentCharacter'));        
    else
        keyID = 13;
    end
    %% start loading data
    if (keyID==13||keyID==32)
        dataStruct = velomap_internal(handleStruct, dataStruct, 'database', 'loadData',[]);
        dataStruct = velomap_internal(handleStruct, dataStruct, 'gui', 'loadData',[]);
        id = 0;
        dataStruct = velomap_internal(handleStruct, dataStruct, 'gui', 'antiAliasCheck',id);
        dataStruct = velomap_internal(handleStruct, dataStruct, 'gui', 'filterCheck',id);
        dataStruct = velomap_internal(handleStruct, dataStruct, 'gui', 'pcmraCheck',id);
        dataStruct = velomap_internal(handleStruct, dataStruct, 'gui', 'eddyCurrentCheck',id);
        dataStruct = velomap_internal(handleStruct, dataStruct, 'gui', 'secOrderCorrCheck',id);
        dataStruct = velomap_internal(handleStruct, dataStruct, 'gui', 'enSightCheck',id);
        dataStruct = velomap_internal(handleStruct, dataStruct, 'gui', 'selMRstruct',id);
        dataStruct = velomap_internal(handleStruct, dataStruct, 'gui', 'pdCheck','off');
        dataStruct = velomap_internal(handleStruct, dataStruct, 'gui', 'applyEnable',id);
%         list   = get(handleStruct.queueList,'String');
%         if isempty(list)
%             set (handleStruct.startQueueBut, 'Enable', 'off')
%         end

		% get user ID
        
% 		userIDlistStr = sprintf('%03d,',1:100);
% 		userIDlistStr = num2str(reshape(sprintf('%03d',1:100),3,100)');
        userIDlist     = importdata('txt_users.txt');
        userIDlistStr  = userIDlist;
		dataStruct.userIDStr = get(handleStruct.userID, 'String');
		if strcmp(dataStruct.userIDStr, '000')
			[resNum,v] = listdlg('PromptString','Select a User ID:','SelectionMode','single','ListString',userIDlistStr);
			resStr = num2str(resNum,'%03d');
			set(handleStruct.userID, 'String', resStr);
			dataStruct.userIDStr = resStr;
		end
		dataStruct.userSubjectIDStr = [dataStruct.patientName,'_user',dataStruct.userIDStr];
	     
        if ~isempty(dataStruct.fileNamesFlow)&& strcmp(dataStruct.loadStatus,'ok');       
           venc       = dataStruct.venc;
           set(handleStruct.orderList,'String','');
           velomap_internal(handleStruct, dataStruct, 'gui', 'image',  []);
           
           yTick = zeros(1,round(max(venc)/25));
           yTickLabel = zeros(1,round(max(venc)/25));           
           for i=0:round((max(venc)/25))
               yTick(i+1)       = dataStruct.phaseRange*i/round((max(venc)/25));
               yTickLabel(i+1)  = -max(venc)+i*50;
           end
           
           v        = ver('matlab');              % get matlab version
           isnewver = str2double(v.Version)>=8.4; % test if pesky version of matlab
           if isnewver % treat colorbar as its own entity
               delete(get(handleStruct.colorbarAx, 'children'));
               set(handleStruct.colorbarAx, 'Tag', '', 'UserData', []);
               set(handleStruct.colorbarAx,'YLimMode','manual','YLim',[0 4096],'YTickMode','manual','YTick',yTick',...
                   'YTickLabelMode','manual','YTickLabel',yTickLabel');
           end
           
        else
             velomap_internal(handleStruct, dataStruct, 'gui', 'dummy',  []); 
        end
        
    end
%% update tab panel for preprocessing functions
elseif (strcmp(entryStr,'noiseMaskingTab')||strcmp(entryStr,'pcmraTab')||strcmp(entryStr,'eddyCurrentCorrectionTab')...
        ||strcmp(entryStr,'antiAliasingTab')||strcmp(entryStr,'pressureDiffTab')...
        || strcmp(entryStr,'nonuniformityTab'))
    
   [dataStruct,handleStruct] = velomap_internal(handleStruct, dataStruct, 'gui', 'updateTabGroupPreprocessing',  entryStr);
    

%% update tab panel for conversion functions
elseif (strcmp(entryStr,'ensightTab')||strcmp(entryStr,'mrstructTab'))
    
   [dataStruct,handleStruct] = velomap_internal(handleStruct, dataStruct, 'gui', 'updateTabGroupConversion',  entryStr);
   
%% number of current slice was changed   
elseif strcmp(entryStr,'SlicesNumSlide')       
        set(handleStruct.SlicesNumEdit,'String',int2str(ceil(get(handleStruct.SlicesNumSlide, 'Value'))));
        
        if ~isempty(get(handleStruct.orderList,'String'))&& strcmp(dataStruct.previewStatus,'on')        
            dataStruct = velomap_internal(handleStruct, dataStruct, 'database', 'preProcess',[]);    
            if strcmp(dataStruct.statRegStatus, 'on')
                dataStruct = velomap_internal(handleStruct, dataStruct, 'database', 'previewStatReg',[]);
            end
            imaOP = 1;
        else
            imaOP = [];
        end
        velomap_internal(handleStruct, dataStruct, 'gui', 'image',  imaOP);
        
elseif strcmp(entryStr,'SlicesNumEdit')
    sliceNum = str2double(get(handleStruct.SlicesNumEdit, 'String'));
    if (isnan(sliceNum)|| sliceNum<1 || sliceNum >(dataStruct.numSlices))        
        sliceNum = 1;
        set(handleStruct.SlicesNumEdit,'String','1');
    end
    set(handleStruct.SlicesNumSlide,'Value',sliceNum);
    
    if ~isempty(get(handleStruct.orderList,'String'))&& strcmp(dataStruct.previewStatus,'on')        
        dataStruct = velomap_internal(handleStruct, dataStruct, 'database', 'preProcess',[]);    
        if strcmp(dataStruct.statRegStatus, 'on')
            dataStruct = velomap_internal(handleStruct, dataStruct, 'database', 'previewStatReg',[]);
        end
        imaOP = 1;
    else
        imaOP = [];
    end
    velomap_internal(handleStruct, dataStruct, 'gui', 'image',  imaOP);
    

%% number of current phase was changed       
elseif strcmp(entryStr,'PhasesNumSlide')
    set(handleStruct.PhasesNumEdit,'String',int2str(ceil(get(handleStruct.PhasesNumSlide, 'Value'))));    
    if ~isempty(get(handleStruct.orderList,'String'))&& strcmp(dataStruct.previewStatus,'on') 
        imaOP = 1;
    else 
        imaOP = [];
    end
    velomap_internal(handleStruct, dataStruct, 'gui', 'image',  imaOP);

elseif strcmp(entryStr,'PhasesNumEdit')
    phaseNum = str2double(get(handleStruct.PhasesNumEdit, 'String'));
    if (isnan(phaseNum)|| phaseNum<1 || phaseNum >(dataStruct.numPhases))        
        phaseNum = 1;
        set(handleStruct.PhasesNumEdit,'String','1');
    end
    set(handleStruct.PhasesNumSlide,'Value',phaseNum);
       
    if ~isempty(get(handleStruct.orderList,'String'))&& strcmp(dataStruct.previewStatus,'on') 
        imaOP = 1;
    else 
        imaOP = [];
    end
    velomap_internal(handleStruct, dataStruct, 'gui', 'image',  imaOP);
%%%
elseif strcmp(entryStr,'loadUnwrapDataBut')
    % search the directory for an existing unwrap data file
    fname = [dataStruct.dataDirectory, filesep,'unwrapData.mat'];
    status = exist(fname,'file');
    % if uwrap data file exists
    if status ==2
        % load file and apply changes accordingly
        dataStruct = velomap_internal(handleStruct, dataStruct, 'gui', 'loadUnwrapData',[]);
        dataStruct.unwrapData = importdata(fname);
        imaOP=1;
        dataStruct.previewStatus ='on';
        dataStruct = velomap_internal(handleStruct, dataStruct, 'database', 'preProcess',  imaOP);        
        dataStruct = velomap_internal(handleStruct, dataStruct, 'gui', 'image',  imaOP);
    else
        warndlg('There is no uwrap data file existing')
    end
%%%
elseif strcmp(entryStr,'unwrapManualBut')
    % search the directory for an existing unwrap data file 
    fname = [dataStruct.dataDirectory, filesep,'unwrapData.mat'];
    status = exist(fname,'file');
    orderList = get(handleStruct.orderList,'String');
    strFlag   = any(strcmp(orderList,'manual anti-aliasing'));
    if status ==2 && ~strFlag
        button = questdlg('An unwrap data file already exist. What do you want to do with it?','',...
                          'Overwrite existing','Load existing','Overwrite existing');
        if strcmp(button,'Load existing') % add to existing
            % load file and apply changes accordingly
            dataStruct = velomap_internal(handleStruct, dataStruct, 'database', 'loadUnwrapData',[]);
            imaOP=1;
            dataStruct = velomap_internal(handleStruct, dataStruct, 'database', 'preProcess',  imaOP);        
            dataStruct = velomap_internal(handleStruct, dataStruct, 'gui', 'image',  imaOP);
            return;
        elseif strcmp(button,'Overwrite existing')
            status = 0;
            delete (fname)
        end
    else
        status = 0;
    end
    
    if status == 0
        if get(handleStruct.unwrapManualBut,'Value')
            dataStruct.previewStatus ='on';
            dataStruct = velomap_internal(handleStruct, dataStruct, 'database', 'manualPhaseUnwrap',[]);        
            if isempty(dataStruct.unwrapData)
                %reload images
                if ~isempty(get(handleStruct.orderList,'String'))
                    imaOP = 1;
                else
                    imaOP = [];
                end
                dataStruct = velomap_internal(handleStruct, dataStruct, 'database', 'preview',imaOP);
                if ~isempty(imaOP)
                    dataStruct = velomap_internal(handleStruct, dataStruct, 'database', 'preProcess',  imaOP);
                end
            else
                imaOP = local_get_image_operations(handleStruct,dataStruct);        
            end        
            dataStruct = velomap_internal(handleStruct, dataStruct, 'gui', 'image',  imaOP);
        end
    end
   
%%% value of venc was changed
elseif strcmp(entryStr,'vencValueEdt')
    if ~isempty (dataStruct.nFE)
        id = [];
        if (dataStruct.nFE == 1)
            id(1) = str2double(get(handleStruct.vencThPlaneEdt, 'String')); 
        elseif (dataStruct.nFE == 3)
            id(1) = str2double(get(handleStruct.vencInPlaneXEdt, 'String'));
            id(2) = str2double(get(handleStruct.vencInPlaneYEdt, 'String'));
            id(3) = str2double(get(handleStruct.vencThPlaneEdt, 'String'));
        end
        dataStruct = velomap_internal(handleStruct, dataStruct, 'database', 'vencValueEdt',id);        
        dataStruct = velomap_internal(handleStruct, dataStruct, 'gui', 'image',[]);
    end
%%% change the actual bigger view
elseif strcmp(entryStr,'changeViewPop')
    imaOP = local_get_image_operations(handleStruct,dataStruct);    
    dataStruct = velomap_internal(handleStruct, dataStruct, 'gui', 'image',  imaOP);
    
%% change the current bigger view
elseif strcmp(entryStr,'encodingTypePop')
    viewID = get(handleStruct.encodingTypePop, 'Value');   
    if viewID==1 %velocity encoding
        dataStruct.encodingType ='velocity';
        set(handleStruct.encodingTypeString,'String','venc in m/s');
    else
        dataStruct.encodingType ='acceleration';
        set(handleStruct.encodingTypeString,'String','acceleration in m/s^2');
    end
    dataStruct = velomap_internal(handleStruct, dataStruct, 'gui', 'image',  []);
    
%% set user ID
elseif strcmp(entryStr,'pullDown_userID')
    dataStruct.userIDStr = get(handleStruct.pullDown_userID, 'Value');   
    
%% preview button or specified key was pressed
elseif (strcmp(entryStr,'previewBut') || strcmp(entryStr,'previewKeyPress'))    
    %% get id of the pressed button
    if strcmp(entryStr,'previewKeyPress')        
        keyID = double(get(handleStruct.figure,'CurrentCharacter'));        
    else
        keyID = 13;
    end
    
    %% start preview modus
    %RPM added the isnumeric and ~isempty portion in case of shift and such. 
    if ~isempty(keyID) && isnumeric(keyID) && (keyID==13||keyID==32)    
        if ~isempty(get(handleStruct.orderList,'String'))
            imaOP = 1;
        else
            imaOP = [];
        end
        dataStruct = velomap_internal(handleStruct, dataStruct, 'database', 'preview',imaOP);
        if ~isempty(imaOP)
            dataStruct = velomap_internal(handleStruct, dataStruct, 'database', 'preProcess',  imaOP);
        end
        dataStruct = velomap_internal(handleStruct, dataStruct, 'gui', 'image',  imaOP);
    end
%
elseif strcmp(entryStr,'previewStatRegBut')    
    dataStruct = velomap_internal(handleStruct, dataStruct, 'database', 'previewStatReg',[]);
    %imaOP      = local_get_image_operations(handleStruct,dataStruct);
    imaOP = get(handleStruct.previewStatRegBut,'Value');
    if imaOP==1% set big view to Vi
        set(handleStruct.changeViewPop,'Value',1)
    end
    dataStruct = velomap_internal(handleStruct, dataStruct, 'gui', 'image',  imaOP);
elseif strcmp(entryStr,'previewBackgroundBut') 
    dataStruct = velomap_internal(handleStruct, dataStruct, 'database', 'previewBackground',[]);
    imaOP = get(handleStruct.previewBackgroundBut,'Value');
    dataStruct.previewBackground = imaOP;
    if imaOP==1% set big view to Vi
        set(handleStruct.changeViewPop,'Value',1)
    end
    dataStruct = velomap_internal(handleStruct, dataStruct, 'gui', 'image',  imaOP);
%%
elseif strcmp(entryStr,'previewPCmraBut')    
    dataStruct = velomap_internal(handleStruct, dataStruct, 'database', 'previewPCmra',[]);    
    %% imaOP      = local_get_image_operations(handleStruct,dataStruct);
    %% dataStruct = velomap_internal(handleStruct, dataStruct, 'gui', 'image',  imaOP);
%%
elseif strcmp(entryStr,'saveFuncSetBut')
    dataStruct = velomap_internal(handleStruct, dataStruct, 'database', 'saveFuncSetting',[]);
    velomap_internal(handleStruct, dataStruct, 'gui', 'dummy',  []);
%% take data set in queue
elseif strcmp(entryStr,'takeQueueBut')
    dataStruct = velomap_internal(handleStruct, dataStruct, 'database', 'takeInQueue',[]);
    dataStruct = velomap_internal(handleStruct, dataStruct, 'database', 'saveFuncSetting',[]);
    dataStruct = velomap_internal(handleStruct, dataStruct, 'database', 'outToExcel',[]);
    velomap_internal(handleStruct, dataStruct, 'gui', 'dummy',  []);
%% load further settings to the actual queue
elseif strcmp(entryStr,'loadInQueueBut')
    set(handleStruct.queueList,'Enable','on')
    dataStruct = velomap_internal(handleStruct, dataStruct, 'database', 'loadInQueue',[]);    
    velomap_internal(handleStruct, dataStruct, 'gui', 'dummy',  []);

%% start dequeueing
elseif strcmp(entryStr, 'startQueueBut')
    %% ask for confirmation
    button = questdlg('DO YOU REALLY WANT TO START DEQUEUEING?','','YES','CANCEL','YES');
    
    if strcmp(button, 'YES')       
        %% get numbers of data sets to be preprocessed
        pathList = get(handleStruct.queueList,'String');
        numSet   = size(pathList);
        if ~isempty(pathList)            
            for i=1:numSet(1)
                %% counter                
                set(handleStruct.statusEdt, 'String',sprintf('data set # %d of %d is in process', i, numSet(1)));
               
                %% set path for status report
                if i==1
                    if ~isempty(dataStruct.dataDirectory)
                        directory_name  = uigetdir(dataStruct.dataDirectory,'select directory for status report');
                    else
                        directory_name ='C:\';
                    end
                    if directory_name==0
                       directory_name  = 'C:\'; 
                    end
                    c = clock;
                    xls_name = sprintf('data_conversion_%04d%02d%02d_%02d%02d',c(1),c(2),c(3),c(4),c(5)); 
                    excelFileStr    = sprintf('%s%s%s%s',directory_name,filesep, xls_name,'.xls');                   
                end
                %% load data                
                dataStruct.inifilePath = char(pathList{i}); 
                readKeys    = {'data path','','dirPath','','none'};
                dirPath     = inifile(dataStruct.inifilePath,'read', readKeys);                
                if strcmp(dirPath,'none')
                    %status ='invalid or nonexistent data path';
                else %% go on
                    dataStruct.dataDirectory = char(dirPath);
                    %% load data
                    dataStruct = velomap_internal(handleStruct, dataStruct, 'database', 'loadData','queue');
                    dataStruct = velomap_internal(handleStruct, dataStruct, 'database', 'outToExcel',[]);
                    %% 
                    dataStruct = velomap_internal(handleStruct, dataStruct, 'database', 'applyTodata','dequeue');
                    %% 
                    dataStruct = velomap_internal(handleStruct, dataStruct, 'database', 'dequeueing',  []);
                    
                    %% calculate pressure difference maps
                    %% if ensight conversion desired
                    readKeys     = {'conversion','','pdmaps','','0'};
                    status = char(inifile(dataStruct.inifilePath,'read', readKeys));
                    if strcmp(status,'1')
                        dataStruct = velomap_internal(handleStruct, dataStruct,'database', 'createPDmaps','dequeue');   
                    end
                    
                    %% if ensight conversion desired
                    readKeys     = {'conversion','','ensight','','0'};
                    status = char(inifile(dataStruct.inifilePath,'read', readKeys));
                    if strcmp(status,'1')
                        dataStruct = velomap_internal(handleStruct, dataStruct, 'database', 'enSight','dequeue');
                        %% read if additional angiography should be
                        %% converted to ensight
                        readKeys   = {'conversion','','angio','','0'};
                        angioCheck = char(inifile(dataStruct.inifilePath,'read', readKeys));
                        if strcmp(angioCheck,'1')
                            dataStruct = velomap_internal(handleStruct, dataStruct, 'database', 'appendToEnsight','dequeue');
                        end
                    end
                    %% read whether avi conversion is desired
                    readKeys   = {'conversion','','avi_movie','','0'};
                    status   = char(inifile(dataStruct.inifilePath,'read', readKeys));
                    if strcmp(status,'1')
                        dataStruct = velomap_internal(handleStruct, dataStruct,'database', 'avimovie',[]);
                    end
                    
                    %% read whether mrstruct should be returned
                    readKeys = {'conversion','','mrstruct','','0'};
                    status   = char(inifile(dataStruct.inifilePath,'read', readKeys));
                    if strcmp(status,'1')
                        dataStruct = velomap_internal(handleStruct, dataStruct,'database', 'selStruct','dequeue');
                    end 
                    %% read whether mrstruct should be returned
                    readKeys = {'conversion','','mra2DICOM','','0'};
                    status   = char(inifile(dataStruct.inifilePath,'read', readKeys));
                    if strcmp(status,'1')
                        dataStruct = velomap_internal(handleStruct, dataStruct,'database', 'pcmraToDICOM','dequeue');
                    end
                    
                end
                
                %% write status to an excel sheet                
                keysExcel   = inifile(dataStruct.inifilePath,'readall');
                warning off MATLAB:xlswrite:AddSheet;
                sheetName = sprintf('%d%s%s',i,'_',dataStruct.patientName);
                xlswrite(excelFileStr, keysExcel, sheetName, 'A1');               
                
            end
            set(handleStruct.statusEdt, 'String',sprintf('conversion of %d data sets is complete',numSet(1)));
        else
            msg = sprintf('%s\n%s','There is no ini files in the list.',...
                                   'No action will be performed'); 
            uiwait(warndlg(msg,'','modal'))            
        end
        %% shutdown pc if choosen;
        if get(handleStruct.shutdownCheck,'Value')
            eval(['!shutdown -s -f -t ' num2str(60)]);
        end
    end
%%%%%%%%%%%%%%%%%%%%

%% apply chosen preprocessing functions to the data
elseif strcmp(entryStr,'applyBut')
    %% check parameters for pcmra calculation
    if get(handleStruct.pcmraCheck,'Value')
        order_list = char(get(handleStruct.orderList,'String'));
        if isempty(strmatch('noise-filter',order_list))&& isempty(strmatch('stdev-filter',order_list)) && isempty(strmatch('derivative-filter',order_list))
            %% ask for confirmation
            button = questdlg('You did not choose any noise masking function! DO YOU REALLY WANT TO CONTINUE?','','CONTINUE','CANCEL','CANCEL');
            if strcmp(button, 'CANCEL')
                return
            end
        end
    end
    %% check autoseg version
    if strcmpi(get(handleStruct.autosegvers,'Enable'),'on')
        verlist = get(handleStruct.autosegvers,'String');
        autosegver = verlist{get(handleStruct.autosegvers,'Value')};
        if isempty(str2num(autosegver)) && get(handleStruct.autoseg,'Value') == 1
            msgbox('autoseg version is invalid. Choose another version.');
            return;
        end   
    end
    dataStruct = velomap_internal(handleStruct, dataStruct, 'database', 'saveFuncSetting',[]);
    dataStruct = velomap_internal(handleStruct, dataStruct, 'database', 'outToExcel',[]);
    act = 1; 
	
    %% chosen action will be performed
    if act        
        imaOP = local_get_image_operations(handleStruct,dataStruct);    
        dataStruct = velomap_internal(handleStruct, dataStruct, 'database', 'applyTodata',[]);
        dataStruct = velomap_internal(handleStruct, dataStruct, 'database', 'createMRstruct',  imaOP);
        %% MM: calculate pressure difference maps
        if get(handleStruct.pdCheck, 'Value')==1
            dataStruct = velomap_internal(handleStruct, dataStruct,'database', 'createPDmaps',[]);   
        end
                 
        if get(handleStruct.enSightCheck, 'Value')==1
            dataStruct = velomap_internal(handleStruct, dataStruct, 'database', 'enSight',[]);
            if (~isempty(get(handleStruct.AnatomDataList,'String')))&&(strcmp(dataStruct.status,'Conversion of data to EnSight-files was successful'))
               dataStruct = velomap_internal(handleStruct, dataStruct, 'database', 'appendToEnsight',[]); 
            end
        end
        if get(handleStruct.selMRstructCheck, 'Value')==1
            dataStruct = velomap_internal(handleStruct, dataStruct,'database', 'selStruct',[]);
        end
        if get(handleStruct.aviMovieCheck, 'Value') ==1
            dataStruct = velomap_internal(handleStruct, dataStruct,'database', 'avimovie',[]);
        end
        if get(handleStruct.pcmraToDICOMCheck, 'Value') ==1
            dataStruct = velomap_internal(handleStruct, dataStruct,'database', 'pcmraToDICOM',[]);
        end
           
        velomap_internal(handleStruct, dataStruct, 'gui', 'dummy',  []);
    end
%%  

elseif strcmp(entryStr,'outToExcelBut')
    dataStruct = velomap_internal(handleStruct, dataStruct, 'database', 'outToExcel',[]);
    velomap_internal(handleStruct, dataStruct, 'gui', 'dummy',  []);
%%    
  
% 20190729 Takashi modified to include autoseg
elseif (strcmp(entryStr,'enSightCheck') || strcmp(entryStr,'aviMovieCheck')||strcmp(entryStr,'selMRstructCheck')...
        ||strcmp(entryStr,'pcmraToDICOMCheck')) || strcmp(entryStr,'autoseg')
    if (get(handleStruct.enSightCheck, 'Value')==1 ||get(handleStruct.aviMovieCheck, 'Value')==1 ...
            ||get(handleStruct.selMRstructCheck, 'Value')==1 ||get(handleStruct.pcmraToDICOMCheck, 'Value')==1 ...
			||get(handleStruct.autoseg, 'Value')==1)        
        id = 1;
    else
        id = 0;
    end    
    dataStruct = velomap_internal(handleStruct, dataStruct, 'gui', 'applyEnable',id);
    
    if strcmp(entryStr,'selMRstructCheck')
       id = get(handleStruct.selMRstructCheck, 'Value');
       velomap_internal(handleStruct, dataStruct, 'gui', 'selMRstruct',id); 
    end
    if strcmp(entryStr,'enSightCheck')
       id = get(handleStruct.enSightCheck, 'Value');
       dataStruct = velomap_internal(handleStruct, dataStruct, 'gui', 'enSightCheck',id);
    end
    
    if strcmp(entryStr,'pcmraToDICOMCheck')
       id = get(handleStruct.pcmraToDICOMCheck, 'Value');
       if (id==1||get(handleStruct.pcmraCheck,'Value'))
           set(handleStruct.pcmraCheck,'Value',1);
           dataStruct = velomap_internal(handleStruct, dataStruct, 'gui', 'pcmraCheck',1);
       else
           dataStruct = velomap_internal(handleStruct, dataStruct, 'gui', 'pcmraCheck',0);
       end
    end
   
    
%%
elseif strcmp(entryStr,'nonuniformityCheck')
    id = get(handleStruct.nonuniformityCheck, 'Value');
    dataStruct = velomap_internal(handleStruct, dataStruct, 'gui', 'nonuniformityCheck',id);
    if id ==1
        [dataStruct,handleStruct] = velomap_internal(handleStruct, dataStruct, 'gui', 'updateTabGroupPreprocessing',  'nonuniformityTab');
    end
elseif strcmp(entryStr,'antiAliasCheck')   
    id = get(handleStruct.antiAliasCheck, 'Value');
    dataStruct = velomap_internal(handleStruct, dataStruct, 'gui', 'antiAliasCheck',id);
    if id
        [dataStruct,handleStruct] = velomap_internal(handleStruct, dataStruct, 'gui', 'updateTabGroupPreprocessing',  'antiAliasingTab');
    end
%%
elseif strcmp(entryStr,'filterCheck')   
    id=get(handleStruct.filterCheck, 'Value');
    dataStruct = velomap_internal(handleStruct, dataStruct, 'gui', 'filterCheck',id);
    if id
        [dataStruct,handleStruct] = velomap_internal(handleStruct, dataStruct, 'gui', 'updateTabGroupPreprocessing',  'noiseMaskingTab');
    end
elseif strcmp(entryStr,'radioORbut')
    set(handleStruct.radioANDBut, 'Value',0);
    set(handleStruct.radioORBut, 'Value',1);
elseif strcmp(entryStr,'radioANDbut')
    set(handleStruct.radioORBut, 'Value',0);
    set(handleStruct.radioANDBut, 'Value',1);
%%
elseif strcmp(entryStr,'appAnatomDataBut')
    dataStruct = velomap_internal(handleStruct, dataStruct, 'database', 'appAnatomData',[]);
%%
elseif strcmp(entryStr,'deleteAnatomDataBut')
    dataStruct = velomap_internal(handleStruct, dataStruct, 'database', 'deleteAnatomData',[]);
%%
elseif strcmp(entryStr,'AnatomDataList')
    id = double(get(handleStruct.figure,'CurrentCharacter'));
    if id==127
        dataStruct = velomap_internal(handleStruct, dataStruct, 'database', 'deleteAnatomData',[]);
    end
%%
elseif strcmp(entryStr,'queueList')
    id = double(get(handleStruct.figure,'CurrentCharacter'));
    if id==127
        dataStruct = velomap_internal(handleStruct, dataStruct, 'database', 'deleteQueueEntry',[]);
    end
%%% delete entries from functions list
elseif strcmp(entryStr,'orderList')
    id = double(get(handleStruct.figure,'CurrentCharacter'));
    if id==127 % delete entry
        n = get(handleStruct.orderList,'Value');
        func_list = get(handleStruct.orderList,'String');
        f = func_list{n};
        dataStruct = velomap_internal(handleStruct, dataStruct, 'gui', 'deleteOrderListEntry',[]);
        if strcmp(f,'anti-aliasing')            
            dataStruct = velomap_internal(handleStruct, dataStruct, 'gui', 'antiAliasCheck',0);        
            [dataStruct,handleStruct] = velomap_internal(handleStruct, dataStruct, 'gui', 'updateTabGroupPreprocessing',  'antiAliasingTab');        
        elseif strcmp(f,'eddy-current')
            dataStruct = velomap_internal(handleStruct, dataStruct, 'gui', 'eddyCurrentCheck',0);
            if ~(get(handleStruct.staticTissueCheck,'Value'))
                dataStruct = velomap_internal(handleStruct, dataStruct, 'gui', 'showStaticTissue',0);
            end
        elseif strcmp(f,'2nd-order-corr')     
            dataStruct = velomap_internal(handleStruct, dataStruct, 'gui', 'secOrderCorrCheck',0);
            if ~(get(handleStruct.staticTissueCheck,'Value'))
                dataStruct = velomap_internal(handleStruct, dataStruct, 'gui', 'showStaticTissue',0);
            end
        elseif strcmp(f,'-filter')||strcmp(f,'stdev-filter')||strcmp(f,'derivative-filter')||strcmp(f,'median-filter')
            mask_list = get(handleStruct.chooseFilterList,'String');
            func_list = get(handleStruct.orderList,'String');
            for m =1:size(mask_list,1)
                x(m,:) = strcmpi(func_list,mask_list{m});
            end
            dataStruct = velomap_internal(handleStruct, dataStruct, 'gui', 'filterCheck',any(x));
        end
    end
elseif strcmp(entryStr,'orderUpBut')||strcmp(entryStr,'orderDownBut')
    id = entryStr(6);    
    dataStruct = velomap_internal(handleStruct, dataStruct, 'gui', 'moveEntry',id);
%
elseif strcmp(entryStr,'backgroundSlide')
    dataValue = num2str(get(handleStruct.backgroundSlide, 'Value'));
    set(handleStruct.backgroundEdt,'String',dataValue);    

elseif strcmp(entryStr,'backgroundEdt')
    dataValue = str2double(get(handleStruct.backgroundEdt, 'String'));
    if (isnan(dataValue)|| dataValue<=0 || dataValue>=1)        
        dataValue = 0;
        set(handleStruct.backgroundEdt,'String','0');
    end
    set(handleStruct.backgroundSlide,'Value',dataValue);        
%
elseif strcmp(entryStr,'chooseFilter')
    index = get(handleStruct.chooseFilterList,'Value');
    id    = get(handleStruct.chooseFilterList, 'String');
    dataStruct = velomap_internal(handleStruct,  dataStruct, 'gui', 'chooseFilter',id{index});
%%
elseif strcmp(entryStr,'noiseValueSlide')
    set(handleStruct.noiseValueEdt,'String',num2str(get(handleStruct.noiseValueSlide, 'Value')));    
%
elseif strcmp(entryStr,'noiseValueEdt')
    dataValue = str2double(get(handleStruct.noiseValueEdt, 'String'));
    if (isnan(dataValue)|| dataValue<=0 || dataValue>=1)        
        dataValue = 0;
        set(handleStruct.noiseValueEdt,'String','0');
    end
    set(handleStruct.noiseValueSlide,'Value',dataValue);    
%
elseif strcmp(entryStr,'stdevNoiseValueSlide')
    set(handleStruct.stdevNoiseValueEdt,'String',num2str(get(handleStruct.stdevNoiseValueSlide, 'Value')));    
%%
elseif strcmp(entryStr,'stdevNoiseValueEdt')
    dataValue = str2double(get(handleStruct.stdevNoiseValueEdt, 'String'));
    if (isnan(dataValue)|| dataValue<=0 || dataValue>=1)        
        dataValue = 0;
        set(handleStruct.stdevNoiseValueEdt,'String','0');    
    end
    set(handleStruct.stdevNoiseValueSlide,'Value',dataValue);    
%% control the slider for derivative noise filter threshold
elseif strcmp(entryStr,'derivatNoiseValueSlide')
    set(handleStruct.derivatNoiseValueEdt,'String',num2str(get(handleStruct.derivatNoiseValueSlide, 'Value')));
%%  control the edit box for derivative noise filter threshold
elseif strcmp(entryStr,'derivatNoiseValueEdt')
    dataValue = str2double(get(handleStruct.derivatNoiseValueEdt, 'String'));
    if (isnan(dataValue)|| dataValue<=0 || dataValue>=1)        
        dataValue = 0;
        set(handleStruct.derivatNoiseValueEdt,'String','0');    
    end
    set(handleStruct.derivatNoiseValueSlide,'Value',dataValue); 
%%
elseif strcmp(entryStr,'eddyCurrentCheck')||strcmp(entryStr,'staticTissueCheck')||strcmp(entryStr,'secOrderCorrCheck')
    id(1,1) = get(handleStruct.staticTissueCheck, 'Value');
    id(1,2) = get(handleStruct.eddyCurrentCheck, 'Value');
    id(1,3) = get(handleStruct.secOrderCorrCheck, 'Value');
    
    dataStruct = velomap_internal(handleStruct, dataStruct, 'gui', 'showStaticTissue',any(id));
    if strcmp(entryStr,'eddyCurrentCheck')
        dataStruct = velomap_internal(handleStruct, dataStruct, 'gui', 'eddyCurrentCheck',id(1,2));
    end
    if strcmp(entryStr,'secOrderCorrCheck')
        dataStruct = velomap_internal(handleStruct, dataStruct, 'gui', 'secOrderCorrCheck',id(1,3));
    end
    if any(id)
        [dataStruct,handleStruct] = velomap_internal(handleStruct, dataStruct, 'gui', 'updateTabGroupPreprocessing',  'eddyCurrentCorrectionTab');
    end
%%
elseif strcmp(entryStr,'staticValueSlide')
    set(handleStruct.staticValueEdt,'String',num2str(get(handleStruct.staticValueSlide, 'Value')));    
%%
elseif strcmp(entryStr,'staticValueEdt')
    dataValue = str2double(get(handleStruct.staticValueEdt, 'String'));
    if (isnan(dataValue)|| dataValue<=0 || dataValue>=1)        
        dataValue = 0;
        set(handleStruct.staticValueEdt,'String','0');
    end
    set(handleStruct.staticValueSlide,'Value',dataValue);    
%%
elseif strcmp(entryStr,'pcmraCheck')
    id = get(handleStruct.pcmraCheck, 'Value');
    dataStruct = velomap_internal(handleStruct, dataStruct, 'gui', 'pcmraCheck',id);
    if id
        [dataStruct,handleStruct] = velomap_internal(handleStruct, dataStruct, 'gui', 'updateTabGroupPreprocessing',  'pcmraTab');
    end

%%
elseif strcmp(entryStr,'pdCheck')   
    id = get(handleStruct.pdCheck, 'Value');
    if id
        stat='on';
    else
        stat='off';
    end
    dataStruct = velomap_internal(handleStruct, dataStruct, 'gui', 'pdCheck',stat);
    if id
        [dataStruct,handleStruct] = velomap_internal(handleStruct, dataStruct, 'gui', 'updateTabGroupPreprocessing',  'pressureDiffTab');
    end

%% creating PD mask
elseif strcmp(entryStr,'previewPDMaskBut')
        dataStruct = velomap_internal(handleStruct, dataStruct, 'database', 'previewPDMask',[]);
    
%% end of my version
elseif strcmp(entryStr,'site')
        dataStruct = velomap_internal(handleStruct, dataStruct, 'gui','setVers', []);

%% close button
elseif strcmp(entryStr,'closeBut')
       delete(handleStruct.figure);
%%% End of: cancel button

%% destroy gui
elseif strcmp(entryStr,'gui_destroy')
   velomap_semaphore = 0;
   delete(handleStruct.figure);
   %dataStruct   =[];
   handleStruct =[];
%%% End of: destroy gui


%% detect invalid entry string
else
   warndlg(lasterr,'Sorry, entry switch ''%s'' was not recognized in ''local_handle_gui_event'', no action performed', entryStr);
end
%%% End of: detect invalid entry string

%%%%% End of: catch gui events
%%%%% 

%% local functions
function imaOP =local_get_image_operations(handleStruct,dataStruct)
        if (~isempty(get(handleStruct.orderList,'String'))|| strcmp(dataStruct.previewStatus,'on')||get(handleStruct.previewStatRegBut,'Value'))
            imaOP = 1;
        else
            imaOP = [];
        end
% Copyright (c) February 1th, 2006 by University of Freiburg, Dept. of Radiology, Section of Medical Physics