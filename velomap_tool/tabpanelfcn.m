function handleStruct = tabpanelfcn(fcn, varargin)
%TABPANELFCN Helper function for tab panels constructed using GUIDE.
%
% handleStruct = tabpanelfcn('make_groups', group_name, panel_names, handleStruct, 1);
% - Puts a field into handleStruct for each entry in the cell array panel_names.
%   These fields contain a list of all component handleStruct with UserData matching
%   the name in panel_names.
% - Puts a field into the handleStruct structure called GROUP_NAME_all for easy 
%   visibility handling.
% 
% handleStruct = tabpanelfcn('tab_group_handler', hObject, handleStruct, group_name);
% - Hides components linked to group_name with user data not matching that of
%   hObject.
% - Shows components linked to group_name with the same user data as hObject
% - Raises the selected tab 2 pixels and returns the last_tab to its previous
%   position.
% - Updates the last_tab field.
%
% Note: UserData for panels must be unique for all tab groups in a figure.

switch fcn
 case 'tab_group_handler'
  handleStruct = tab_group_handler(varargin{:});
 case 'make_groups'
  handleStruct = make_groups(varargin{:});
 otherwise
  error('Unrecognized command option.  See HELP TABPANELFCN');
end


% --- helper function for tab buttons
% --- Makes controls for selected tab visible and hides others
% --- Note: this modifies the handleStruct structure - caller must store
function handleStruct = tab_group_handler(hObject, handleStruct, group_name)

if isempty(handleStruct.last_tab) | handleStruct.last_tab ~= hObject
    % define some constants for to improve cross platform look and feel
    if ispc
        edge_offset   = 1;
        bottom_offset = 1;
        width_offset  = 3;
        height        = 1;
    elseif isunix
        edge_offset   = 3;
        bottom_offset = 0;
        width_offset  = 6;
        height        = 2;
    end
    
    % create the blank text object if it has not been created yet
    if ~isfield(handleStruct, 'blank_patch')
        handleStruct.blank_patch = uicontrol('style','text', ...
            'units','pixels','Position',[0,0,1,height]);%,'backg','r');
    end
    
    % make only my children visible
    actTab = get(hObject,'UserData');
    my_kids = getfield(handleStruct, get(hObject,'UserData'));
    all_tabObjects = [group_name, '_all'];
    set(getfield(handleStruct, all_tabObjects),'Visible','off');
   
    set(my_kids,'Visible','on');
    
    if strcmp(actTab,'pcmra')
        set(handleStruct.staticTissueCheck,'Visible','on')
    end
    
    % make sure my kids are on top
    uistack(my_kids,'top');
    
    % get locations for patch and this button
    old_patch_position = get_position_using_pixels(handleStruct.blank_patch);
    old_this_position = get_position_using_pixels(hObject);

    % move the patch under this button
    new_patch_position = old_patch_position;
    % make sure bottom is correct
    new_patch_position(2) = old_this_position(2) + bottom_offset;
    % adjust left edge position
    new_patch_position(1) = old_this_position(1) + edge_offset;
    % adjust width
    new_patch_position(3) = old_this_position(3) - width_offset;
    % update position
    set_position_using_pixels(handleStruct.blank_patch, new_patch_position);
    
    % make this button 3 pixels taller and one pixel lower to give it a raised look
    new_this_position = old_this_position;
    new_this_position(2) = new_this_position(2) - 1;
    new_this_position(4) = new_this_position(4) + 2;
    set_position_using_pixels(hObject, new_this_position);
    
    % shrink previous button to original size and location
    if ~isempty(handleStruct.last_tab)%isfield(handleStruct,'last_tab')
        old_last_tab_position = get_position_using_pixels(handleStruct.last_tab);
        old_last_tab_position(4) = old_last_tab_position(4) - 2;
        old_last_tab_position(2) = old_last_tab_position(2) + 1;
        set_position_using_pixels(handleStruct.last_tab, old_last_tab_position);
    end
    
    % update last_tab
    handleStruct.last_tab= hObject;
    
end


% --- helper function for tab buttons
% --- Updates handleStruct structure with fields for tab_group_handler
function handleStruct = make_groups(group_name, panel_names, handleStruct, top_panel_index)

all_panel_children = [];
for i = 1:length(panel_names)
    % find all the components with the given UserData - this will include the
    % tab push button
     this_set = findobj('-regexp', 'UserData', panel_names{i});
%    this_set = findobj('-strcmp', 'UserData', panel_names{i});
    %% JB this_set = findobj(handleStruct.figure, 'UserData', panel_names{i});
    % remove the tab push button from the list
    this_set = setdiff(this_set, getfield(handleStruct, group_name));
    % update the handleStruct structure with a new field for these handleStruct - all
    % the components belonging to the panel with this UserData
    handleStruct = setfield(handleStruct, panel_names{i}, this_set);
    % keep a list of all panel children for access later
    all_panel_children = [all_panel_children; this_set];
end

% store all panel children for later showing/hiding
handleStruct = setfield(handleStruct, [group_name, '_all'], all_panel_children);

% bring top panel to front by faking a click on top panel's toggle button
top_panel_name = panel_names{top_panel_index};
toggle_button_handle = findobj('-regexp','Tag',group_name,'UserData',top_panel_name);
%toggle_button_handle = findobj(handleStruct,'Tag', group_name); % 'UserData', top_panel_name,                             
                               
handleStruct = tab_group_handler(toggle_button_handle, handleStruct, group_name);

% --- helper function for setting position
function set_position_using_pixels(hObject, position)
old_units = get(hObject, 'Units');
set(hObject, 'Units', 'pixels');
was_error = 0;
try
    set(hObject, 'Position', position)
catch
    was_error = 1;
end
set(hObject, 'Units', old_units);
if was_error
    %error(lasterr);
end

% --- helper function for getting position
function position = get_position_using_pixels(hObject)
old_units = get(hObject, 'Units');
set(hObject, 'Units', 'pixels');
position = get(hObject, 'Position');
set(hObject, 'Units', old_units);
