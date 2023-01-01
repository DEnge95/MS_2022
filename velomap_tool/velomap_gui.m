function varargout = velomap_gui(varargin)
%velomap_GUI M-file for velomap_gui.fig
%      velomap_GUI, by itself, creates a new velomap_GUI or raises the existing
%      singleton*.
%
%      H = velomap_GUI returns the handle to a new velomap_GUI or the handle to
%      the existing singleton*.
%
%      velomap_GUI('Property','Value',...) creates a new velomap_GUI using the
%      given property value pairs. Unrecognized properties are passed via
%      varargin to velomap_gui_OpeningFcn.  This calling syntax produces a
%      warning when there is an existing singleton*.
%
%      velomap_GUI('CALLBACK') and velomap_GUI('CALLBACK',hObject,...) call the
%      local function named CALLBACK in velomap_GUI.M with the given input
%      arguments.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help velomap_gui

% Last Modified by GUIDE v2.5 11-Jan-2011 14:16:06

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @velomap_gui_OpeningFcn, ...
                   'gui_OutputFcn',  @velomap_gui_OutputFcn, ...
                   'gui_LayoutFcn',  [], ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
   gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before velomap_gui is made visible.
function velomap_gui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   unrecognized PropertyName/PropertyValue pairs from the
%            command line (see VARARGIN)

% Choose default command line output for velomap_gui
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes velomap_gui wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = velomap_gui_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in velomap_outToExcel.
function velomap_outToExcel_Callback(hObject, eventdata, handles)
% hObject    handle to velomap_outToExcel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on selection change in velomap_numIterPop.
function velomap_numIterPop_Callback(hObject, eventdata, handles)
% hObject    handle to velomap_numIterPop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns velomap_numIterPop contents as cell array
%        contents{get(hObject,'Value')} returns selected item from velomap_numIterPop


% --- Executes during object creation, after setting all properties.
function velomap_numIterPop_CreateFcn(hObject, eventdata, handles)
% hObject    handle to velomap_numIterPop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in velomap_mraNormRadioBut.
function velomap_mraNormRadioBut_Callback(hObject, eventdata, handles)
% hObject    handle to velomap_mraNormRadioBut (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of velomap_mraNormRadioBut


% --- Executes on button press in velomap_mraSqrRadioBut.
function velomap_mraSqrRadioBut_Callback(hObject, eventdata, handles)
% hObject    handle to velomap_mraSqrRadioBut (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of velomap_mraSqrRadioBut



function velomap_lowlimSlicesEdt_Callback(hObject, eventdata, handles)
% hObject    handle to velomap_lowlimSlicesEdt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of velomap_lowlimSlicesEdt as text
%        str2double(get(hObject,'String')) returns contents of velomap_lowlimSlicesEdt as a double


% --- Executes during object creation, after setting all properties.
function velomap_lowlimSlicesEdt_CreateFcn(hObject, eventdata, handles)
% hObject    handle to velomap_lowlimSlicesEdt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function velomap_uplimSlicesEdt_Callback(hObject, eventdata, handles)
% hObject    handle to velomap_uplimSlicesEdt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of velomap_uplimSlicesEdt as text
%        str2double(get(hObject,'String')) returns contents of velomap_uplimSlicesEdt as a double


% --- Executes during object creation, after setting all properties.
function velomap_uplimSlicesEdt_CreateFcn(hObject, eventdata, handles)
% hObject    handle to velomap_uplimSlicesEdt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function velomap_lowlimPhasesEdt_Callback(hObject, eventdata, handles)
% hObject    handle to velomap_lowlimPhasesEdt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of velomap_lowlimPhasesEdt as text
%        str2double(get(hObject,'String')) returns contents of velomap_lowlimPhasesEdt as a double


% --- Executes during object creation, after setting all properties.
function velomap_lowlimPhasesEdt_CreateFcn(hObject, eventdata, handles)
% hObject    handle to velomap_lowlimPhasesEdt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function velomap_uplimPhasesEdt_Callback(hObject, eventdata, handles)
% hObject    handle to velomap_uplimPhasesEdt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of velomap_uplimPhasesEdt as text
%        str2double(get(hObject,'String')) returns contents of velomap_uplimPhasesEdt as a double


% --- Executes during object creation, after setting all properties.
function velomap_uplimPhasesEdt_CreateFcn(hObject, eventdata, handles)
% hObject    handle to velomap_uplimPhasesEdt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes when user attempts to close velomap_figure.
function velomap_figure_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to velomap_figure (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure
delete(hObject);


% --- Executes during object deletion, before destroying properties.
function velomap_previewPCmraBut_DeleteFcn(hObject, eventdata, handles)
% hObject    handle to velomap_previewPCmraBut (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --- Executes during object deletion, before destroying properties.
function velomap_maskPCmraBut_DeleteFcn(hObject, eventdata, handles)
% hObject    handle to velomap_previewPCmraBut (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --- Executes during object deletion, before destroying properties.
function velomap_colorbarAx_DeleteFcn(hObject, eventdata, handles)
% hObject    handle to velomap_colorbarAx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes during object creation, after setting all properties.
function velomap_figure_CreateFcn(hObject, eventdata, handles)
% hObject    handle to velomap_figure (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --------------------------------------------------------------------
function Untitled_1_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on selection change in velomap_changeViewPop.
function velomap_changeViewPop_Callback(hObject, eventdata, handles)
% hObject    handle to velomap_changeViewPop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns velomap_changeViewPop contents as cell array
%        contents{get(hObject,'Value')} returns selected item from velomap_changeViewPop


% --- Executes during object creation, after setting all properties.
function velomap_changeViewPop_CreateFcn(hObject, eventdata, handles)
% hObject    handle to velomap_changeViewPop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in velomap_appAnatomDataCheck.
function velomap_appAnatomDataCheck_Callback(hObject, eventdata, handles)
% hObject    handle to velomap_appAnatomDataCheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of velomap_appAnatomDataCheck


% --- Executes on button press in velomap_appAnatomDataBut.
function velomap_appAnatomDataBut_Callback(hObject, eventdata, handles)
% hObject    handle to velomap_appAnatomDataBut (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in velomap_deleteAnatomDataBut.
function velomap_deleteAnatomDataBut_Callback(hObject, eventdata, handles)
% hObject    handle to velomap_deleteAnatomDataBut (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on selection change in velomap_AnatomDataList.
function velomap_AnatomDataList_Callback(hObject, eventdata, handles)
% hObject    handle to velomap_AnatomDataList (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns velomap_AnatomDataList contents as cell array
%        contents{get(hObject,'Value')} returns selected item from velomap_AnatomDataList


% --- Executes during object creation, after setting all properties.
function velomap_AnatomDataList_CreateFcn(hObject, eventdata, handles)
% hObject    handle to velomap_AnatomDataList (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in velomap_queueList.
function velomap_queueList_Callback(hObject, eventdata, handles)
% hObject    handle to velomap_queueList (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns velomap_queueList contents as cell array
%        contents{get(hObject,'Value')} returns selected item from velomap_queueList


% --- Executes during object creation, after setting all properties.
function velomap_queueList_CreateFcn(hObject, eventdata, handles)
% hObject    handle to velomap_queueList (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in velomap_unwrapManualBut.
function velomap_unwrapManualBut_Callback(hObject, eventdata, handles)
% hObject    handle to velomap_unwrapManualBut (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of velomap_unwrapManualBut


% --- Executes on button press in radiobutton11.
function radiobutton11_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton11


% --- Executes on button press in velomap_pcmraToDICOMCheck.
function velomap_pcmraToDICOMCheck_Callback(hObject, eventdata, handles)
% hObject    handle to velomap_pcmraToDICOMCheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of velomap_pcmraToDICOMCheck


% --- Executes during object creation, after setting all properties.
function velomap_derivatNoiseValueSlide_CreateFcn(hObject, eventdata, handles)
% hObject    handle to velomap_derivatNoiseValueSlide (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes during object creation, after setting all properties.
function velomap_derivatNoiseValueEdt_CreateFcn(hObject, eventdata, handles)
% hObject    handle to velomap_derivatNoiseValueEdt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in velomap_secOrderCorrCheck.
function velomap_secOrderCorrCheck_Callback(hObject, eventdata, handles)
% hObject    handle to velomap_secOrderCorrCheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of velomap_secOrderCorrCheck


% --- Executes on selection change in velomap_interleavesPopup.
function velomap_interleavesPopup_Callback(hObject, eventdata, handles)
% hObject    handle to velomap_interleavesPopup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns velomap_interleavesPopup contents as cell array
%        contents{get(hObject,'Value')} returns selected item from velomap_interleavesPopup


% --- Executes during object creation, after setting all properties.
function velomap_interleavesPopup_CreateFcn(hObject, eventdata, handles)
% hObject    handle to velomap_interleavesPopup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in velomap_pdCheck.
function velomap_pdCheck_Callback(hObject, eventdata, handles)
% hObject    handle to velomap_pdCheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of velomap_pdCheck


% --- Executes on button press in velomap_pcmraSumSquaresCheck.
function velomap_pcmraSumSquaresCheck_Callback(hObject, eventdata, handles)
% hObject    handle to velomap_pcmraSumSquaresCheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of velomap_pcmraSumSquaresCheck


% --- Executes on button press in velomap_pcmraMeanAbsVelCheck.
function velomap_pcmraMeanAbsVelCheck_Callback(hObject, eventdata, handles)
% hObject    handle to velomap_pcmraMeanAbsVelCheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of velomap_pcmraMeanAbsVelCheck


% --- Executes on button press in velomap_pcmraPseudoComplDiffCheck.
function velomap_pcmraPseudoComplDiffCheck_Callback(hObject, eventdata, handles)
% hObject    handle to velomap_pcmraPseudoComplDiffCheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of velomap_pcmraPseudoComplDiffCheck

% --- Executes on button press in velomap_pcmraSumSquaresCheck.velomap_gui
function velomap_pcmraSqrtSumSquaresCheck_Callback(hObject, eventdata, handles)
% hObject    handle to velomap_pcmraSumSquaresCheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of velomap_pcmraSumSquaresCheck

% --- Executes on button press in velomap_pcmraSumSquaresCheck.velomap_gui
function velomap_pcmraVelCorrCheck_Callback(hObject, eventdata, handles)
% hObject    handle to velomap_pcmraSumSquaresCheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of velomap_pcmraSumSquaresCheck

% --- Executes on selection change in velomap_vencThresholdSlider.
function velomap_vencThresholdSlider_Callback(hObject, eventdata, handles)
% hObject    handle to velomap_vencThresholdSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns velomap_vencThresholdSlider contents as cell array
%        contents{get(hObject,'Value')} returns selected item from velomap_vencThresholdSlider


% --- Executes during object creation, after setting all properties.
function velomap_vencThresholdSlider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to velomap_vencThresholdSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit13_Callback(hObject, eventdata, handles)
% hObject    handle to edit13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit13 as text
%        str2double(get(hObject,'String')) returns contents of edit13 as a double


% --- Executes during object creation, after setting all properties.
function edit13_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function velomap_vencThresholdEdit_Callback(hObject, eventdata, handles)
% hObject    handle to velomap_vencThresholdEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of velomap_vencThresholdEdit as text
%        str2double(get(hObject,'String')) returns contents of velomap_vencThresholdEdit as a double


% --- Executes during object creation, after setting all properties.
function velomap_vencThresholdEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to velomap_vencThresholdEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function velomap_thresholdPDmaskEdit_Callback(hObject, eventdata, handles)
% hObject    handle to velomap_thresholdPDmaskEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of velomap_thresholdPDmaskEdit as text
%        str2double(get(hObject,'String')) returns contents of velomap_thresholdPDmaskEdit as a double


% --- Executes during object creation, after setting all properties.
function velomap_thresholdPDmaskEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to velomap_thresholdPDmaskEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in velomap_TabGroupPreprocessing.
function velomap_TabGroupPreprocessing_Callback(hObject, eventdata, handles)
% hObject    handle to velomap_TabGroupPreprocessing (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in velomap_EddyCurrentCorrTab.
function velomap_EddyCurrentCorrTab_Callback(hObject, eventdata, handles)
% hObject    handle to velomap_EddyCurrentCorrTab (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in velomap_PcmraTab.
function velomap_PcmraTab_Callback(hObject, eventdata, handles)
% hObject    handle to velomap_PcmraTab (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in velomap_PreprocessingTabFrame.
function velomap_PreprocessingTabFrame_Callback(hObject, eventdata, handles)
% hObject    handle to velomap_PreprocessingTabFrame (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in velomap_ConversionTabFrame.
function velomap_ConversionTabFrame_Callback(hObject, eventdata, handles)
% hObject    handle to velomap_ConversionTabFrame (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in velomap_TabGroupConversion.
function velomap_TabGroupConversion_Callback(hObject, eventdata, handles)
% hObject    handle to velomap_TabGroupConversion (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function uipanel24_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to uipanel24 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in checkbox19.
function checkbox19_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox19 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox19


% --- Executes on button press in velomap_shutdownCheck.
function velomap_shutdownCheck_Callback(hObject, eventdata, handles)
% hObject    handle to velomap_shutdownCheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of velomap_shutdownCheck


% --- Executes on button press in velomap_closeBut.
function velomap_closeBut_Callback(hObject, eventdata, handles)
% hObject    handle to velomap_closeBut (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes during object creation, after setting all properties.
function velomap_encodingPop_CreateFcn(hObject, eventdata, handles)
% hObject    handle to velomap_encodingPop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in velomap_encodingPop.
function velomap_encodingPop_Callback(hObject, eventdata, handles)
% hObject    handle to velomap_encodingPop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns velomap_encodingPop contents as cell array
%        contents{get(hObject,'Value')} returns selected item from velomap_encodingPop


% --- Executes on selection change in popupmenu15.
function popupmenu15_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu15 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popupmenu15 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu15


% --- Executes on selection change in velomap_encodingTypePop.
function velomap_encodingTypePop_Callback(hObject, eventdata, handles)
% hObject    handle to velomap_encodingTypePop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns velomap_encodingTypePop contents as cell array
%        contents{get(hObject,'Value')} returns selected item from velomap_encodingTypePop


% --- Executes on button press in velomap_loadUnwrapDatalBut.
function velomap_loadUnwrapDatalBut_Callback(hObject, eventdata, handles)
% hObject    handle to velomap_loadUnwrapDatalBut (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in velomap_loadUnwrapDataBut.
function velomap_loadUnwrapDataBut_Callback(hObject, eventdata, handles)
% hObject    handle to velomap_loadUnwrapDataBut (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in velomap_loadMaskBut.
function velomap_loadMaskBut_Callback(hObject, eventdata, handles)
% hObject    handle to velomap_loadMaskBut (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% --- Executes on key press with focus on velomap_orderList and none of its controls.
function velomap_orderList_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to velomap_orderList (see GCBO)
% eventdata  structure with the following fields (see UICONTROL)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in velomap_orderUpBut.
function velomap_orderUpBut_Callback(hObject, eventdata, handles)
% hObject    handle to velomap_orderUpBut (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


