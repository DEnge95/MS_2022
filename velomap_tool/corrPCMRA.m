function [corrdata, mag] = corrPCMRA(vel,conn)
%CORRPCMRA Calculates a PCMRA by correlating the velocity data through time
%   This function first calculates the velocity magnitude at each voxel in
%   the image. Then is calculates a Pearson Correlation Coefficient for
%   each voxel. Briefly: think of each voxel in the [X,Y,Z] having a
%   velocity curve with a point at each time step. These curves are
%   correlated with the nearest neighbors. This results in a time-averaged
%   PCMRA after taking the median. The calculation evaluates the
%   correlation coefficient by shifting the mag matrix and calculating
%   values once for each nearest neighbor). This uses a helper function
%   PEARSONCC, which is much faster than the matlab built-in function
%   (which cannot evaluate the whole matrix along an arbitrary dimension).
%
%   Note that the normalization is done automatically. This normalization
%   accounts for and weighs edge voxels appropriately (so they are NOT
%   penalized for being on the edge). This behavior is intentional, since
%   edge voxels often contain ROI.
% Input:
%   - vel: Expects a 5D matrx: [X,Y,Z,T,V], where the velocity generally
%   has 3 components, this is typically the dataAy field from a vel_struct.
%   Should also work if passed the magnitude of the velocity. 
%   - conn: either 6 or 27, specifying the connectivity to check. 6 is
%   faster and tends to work better, and is the default if another value is
%   given to the function (in error)
% Outputs:
%   - corrdata: a matrix of size [X,Y,Z], which takes values [-1  1],
%   representing the average correlation coefficient at each voxel with its
%   six nearest neighbors. 
%   - mag: a matrix of size [X,Y,Z,T] that gives the velocity magnitude at
%   each voxel.
%
%
% Written by Mike Scott, NU

%% Calculate the velocity magnitude at each voxel / time
% Find the velocity dimension (ie dimension with size 3)
if size(vel,4) == 3
    mag = squeeze((sum(vel.^2,4)).^0.5);
elseif size(vel,5) == 3
    mag = squeeze((sum(vel.^2,5)).^0.5);
else
    error('No velocity dimension found')
end

% Make sure the conn value is numeric, if not set to 6 connectivity
if ~isnumeric(conn)
    conn = 6;
end

%% Correlate across time
% Pad the mag array
padmag = padarray(mag,[1 1 1 0 0],NaN,'both');

if conn == 27
    % Preallocate to store the 27 comparisons
    singleCorr = zeros(size(mag,1),size(mag,2),size(mag,3),27);
    % Set a counter to store the correlations
    counter = 1;
    % Move through neighbors
    for xx = -1:1
        for yy = -1:1
            for zz = -1:1
                if xx == 0 && yy == 0 && zz == 0
                    % Do nothing
                else
                    % Shift the matrix and correlate
                    singleCorr(:,:,:,counter) = pearsonCC(mag,padmag(2+xx:end-1+xx,2+yy:end-1+yy,2+zz:end-1+zz,:),4);
                    % Increment the counter
                    counter = counter + 1;
                end
            end
        end
    end    
else
    % Preallocate to store the 6 comparisons
    singleCorr = zeros(size(mag,1),size(mag,2),size(mag,3),6);
    % Set a counter to move through the correlations
    counter = 1;
    % Move through neighbors
    for xx = -1:1
        for yy = -1:1
            for zz = -1:1
                % Only use neighbors where one x,y,z is 1, others zero
                if (xx^2+yy^2+zz^2) ~= 1
                    % Do nothing
                else
                    % Shift the matrix and correlate
                    singleCorr(:,:,:,counter) = pearsonCC(mag,padmag(2+xx:end-1+xx,2+yy:end-1+yy,2+zz:end-1+zz,:),4);
                    counter = counter +1;
                end
            end
        end
    end    
end
% Take the median (better at edges than mean)
corrdata = nanmedian(singleCorr,4);
end

% Helper function
function r = pearsonCC(a,b,dim)
% Computes the Pearson correlation coefficient between elements in a and b
% along dim. For example, in a dataset of X,Y,Z,T: dim 4 will compute the
% correlation coefficient for each x,y,z voxel along the time dimension.
% Remove means
az = bsxfun(@minus, a, mean(a,dim));
bz = bsxfun(@minus, b, mean(b,dim));
% Standard Pearson correlation coefficient formula
a2 = az .^ 2;
b2 = bz .^ 2;
ab = az .* bz;
r = sum(ab, dim) ./ sqrt(sum(a2, dim) .* sum(b2, dim));
end
