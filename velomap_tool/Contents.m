% VELOMAP_TOOL
% 20080331 updated (Jelena Bock)
%
% Files
%   calculate_pcmra     - function for calculating 3 different types of pc-mra from 3D velocity
%   createPDmask        - calculate pressure difference map from velocity data
%   fitParaboloid       - least square fit of a paraboloid to phase data
%   flow_anti_aliasing  - function for phase unwrapping
%   inifile             - Creates, reads, or writes data from/to a standard ini (ascii)
%   lsqf2               - least square fit of plane to mag weighted 2d pc data
%   magicwand2          - Simulation of Photoshop's magic wand tool.
%   my_ginput           - GINPUT Graphical input from mouse.
%   phase_unwrap_manual - function to unwrap phase manually
%   setWindowOnTop      - sets a figures Always On Top state on or off
%   tabpanelfcn         - Helper function for tab panels constructed using GUIDE.
%   velomap_gui         - velomap_GUI M-file for velomap_gui.fig
%   velomap_internal    - see velomap tool for more details
%   velomap_model       - see velomap_tool for a full description of the tool.
%   velomap_tool        - Import and visualization of dicom data(PC-MRI-data)
%   write_scaninfo      - get selected header information from dicom header of selected file
