function     bin_mask = magicwand2(im, tolerance, ylist, xlist)
%MAGICWAND2  Simulation of Photoshop's magic wand tool.
%            It allows selection of connected groups of pixels whose colors
%            are within a predefined tolerance of some reference pixels.
%
% SYNTAX
%    (1) bin_mask = magicwand2(im, tolerance, ylist, xlist);
%    (2) bin_mask = magicwand2(im, tolerance);
%    (3) bin_mask = magicwand2(im);
%    (4) bin_mask = magicwand2;
%
%    If xlist and ylist are omitted, the user is prompted to choose them
%    interactively from the current figure.
%    If tolerance is omitted or empty, it gets the default value (TOL=50)
%    If the image is omitted or empty, the user is prompted to get it
%    interactively.
%    The output bin_mask is automatically displayed in a different figure.
%
% INPUT
%   im:        input image RGB
%   tolerance: distance to reference pixels
%   ylist:     vector of row cordinates    (reference pixels)
%   xlist:     vector of column cordinates (reference pixels)
%
% OUTPUT
%   bin_mask: binary mask of selected regions
%

% VERSION HISTORY
% 1) Original version (magicwand.m) by Daniel Leo Lau, April 7, 1997
%    email: lau@ece.udel.edu  (mex file)
% 2) Updated for MATLAB 6.5 by Yoram Tal, June 30, 2003                                                            
%    email: yoram_tal@yahoo.com (mex file)                                                                         
% 3) Rewriten in MATLAB (no mex file) by Son Lam Phung, March 30, 2004 
%    email: s.phung@ecu.edu.au 
% 4) Current version by Yoram Tal, October 12 2004
%    email: yoram_tal@yahoo.com 
%     
%    The current version is based on Phung's version, but includes
%    some major changes:
%    a) It is much faster. (about 30 times for Phung's test image).
%    b) It is (optionally) interactive both for image and point selection. 
%    c) It works for both RGB and graylevel images.

% (C) Yoram Tal 2004

TOL = 50;         % default value
getPoints = 0;

if nargin == 3,
    error('Not enough input data');
elseif nargin < 3,
    getPoints = 1;
end

if nargin < 2 || isempty(tolerance),
    tolerance = TOL;
end
H = size(im, 1); % image height
W = size(im, 2); % image width
bin_mask = false(H,W);


% Get points interactively
if getPoints,
    
    %figpos = get (gcf,'Position');
    scnsize    = get(0,'ScreenSize');
    figpos(4)  = scnsize(4)-96;
    zoomfactor = figpos(4)/H;
    figpos(3)  = W*zoomfactor;
    
    nfig = figure;                  
    himage=imshow(im,[min(min(im)),max(max(im))]);

    %pos =[figpos(1) figpos(2) figpos(3) figpos(4)];
    pos =[8 40 figpos(3) figpos(4)];
    set(nfig,'Position',pos,'MenuBar','none','ToolBar','figure');
    iptsetpref('ImshowInitialMagnification','fit');
    set(gca,'position',[0 0 1 1])
    
    but   = 0;
    ii    = 0;
    xlist = [];
    ylist = [];
    hplot = [];
    hold on
        
    while but ~= 3,
        ii = ii + 1;        
        [x, y, but] = my_ginput(1,himage);
        xlist(ii)   = round(x);
        ylist(ii)   = round(y);
        hplot(ii)   = plot(x,y,'.');
    end
    xlist(end) = [];
    ylist(end) = [];
    delete(hplot);    
    hold off
end


% Check points validity
if isempty(xlist) || isempty(ylist),
    %('Point list is empty');
    return
end




k = ylist > 0 & ylist <= H;
k = k & xlist > 0 & xlist <= W;

if ~any(k),
    %error('Coordinates out of range');
    return
elseif ~all(k),
    disp('Warning: some coordinates out of range');
end

ylist = ylist(k);
xlist = xlist(k);

N = length(ylist); % Number of reference pixels


%Create the binary mask
color_mask = false(H, W);

if ndims(im) < 3,
    g = double(im);
    for i = 1:N,
        ref = double(im(ylist(i),xlist(i))); 
        color_mask = color_mask | (g - ref).^2 <= tolerance^2;
    end
elseif ndims(im) == 3,
    c_r = double(im(:, :, 1)); % Red channel
    c_g = double(im(:, :, 2)); % Green channel
    c_b = double(im(:, :, 3)); % Blue channel
    for i = 1:N,
        ref_r = double(im(ylist(i), xlist(i), 1));
        ref_g = double(im(ylist(i), xlist(i), 2));
        ref_b = double(im(ylist(i), xlist(i), 3));
        color_mask = color_mask | ...
            ((c_r - ref_r).^2 + (c_g - ref_g).^2 + (c_b - ref_b).^2)...
             <= tolerance^2;
    end
end

% Connected component labelling
[objects, count] = bwlabel(color_mask, 8); 


[y x v] = find(objects);
segList = [];

for i = 1:N,
    k = find(x == xlist(i) & y == ylist(i));
    segList = [segList; v(k)];
end

segList = unique(segList);


LUT = zeros(1,count+1);
LUT(segList+1) = 1;
bin_mask = LUT(objects+1);


