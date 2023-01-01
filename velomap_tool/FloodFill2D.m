function pdstruct = FloodFill2D(pdstruct,row,col)%,,)
% [B] = FloodFill3D(A, slice);
% This program flood fills a 6-connected 3D region. The input matrix MUST
% be a binary image. The user will select a seed (point) in the matrix to
% initiate the flood fill. You must specify the matrix slice in which you
% wish to place the seed.
% 
% pdstruct: structure contains needed values for pressure maps calculation
% slice = a chosen slice in the matrix where you wish to place the seed.
%
% Enjoy,
% F. Dinath
if isstruct(pdstruct)
    A = single(pdstruct.MASK);
else % In case a logical matrix comes in.
    A = single(pdstruct);
end
A(1,:,:) = NaN;     % Pad the border of the matrix
A(end,:,:) = NaN;   % so the program doesn't attempt 
A(:,1,:) = NaN;     % to seek voxels outside the matrix
A(:,end,:) = NaN;   % boundry during the for loop below.


% figure;imagesc(A(:,:,slice));
% title('select seed on figure');
% 
% k = waitforbuttonpress;
% point = get(gca,'CurrentPoint'); % button down detected
% point = [fliplr(round(point(2,1:2))) slice]
% close;

point = [row, col];
if A(point(1), point(2))==1
    A(point(1), point(2)) = NaN;
    a{1} = sub2ind(size(A), point(1), point(2));

    i = 1;

    while 1

        i = i+1;
        a{i} = [];

        [x, y, z] = ind2sub(size(A), a{i-1});
        

        ob = nonzeros((A(sub2ind(size(A), x-1, y-1)) == 1).*sub2ind(size(A), x, y));
        A(ob) = NaN;
        a{i} = [a{i} ob'];

        ob = nonzeros((A(sub2ind(size(A), x-1, y+1)) == 1).*sub2ind(size(A), x, y));
        A(ob) = NaN;
        a{i} = [a{i} ob'];
        
        ob = nonzeros((A(sub2ind(size(A), x+1, y-1)) == 1).*sub2ind(size(A), x, y));
        A(ob) = NaN;
        a{i} = [a{i} ob'];

        ob = nonzeros((A(sub2ind(size(A), x+1, y+1)) == 1).*sub2ind(size(A), x, y));
        A(ob) = NaN;
        a{i} = [a{i} ob'];


        ob = nonzeros((A(sub2ind(size(A), x-1, y)) == 1).*sub2ind(size(A), x-1, y));
        A(ob) = NaN;
        a{i} = [a{i} ob'];

        ob = nonzeros((A(sub2ind(size(A), x+1, y)) == 1).*sub2ind(size(A), x+1, y));
        A(ob) = NaN;
        a{i} = [a{i} ob'];

        ob = nonzeros((A(sub2ind(size(A), x, y-1)) == 1).*sub2ind(size(A), x, y-1));
        A(ob) = NaN;
        a{i} = [a{i} ob'];

        ob = nonzeros((A(sub2ind(size(A), x, y+1)) == 1).*sub2ind(size(A), x, y+1));
        A(ob) = NaN;
        a{i} = [a{i} ob'];

        if isempty(a{i});            
            break;
        end
    end
end

b = cell2mat(a);
%b = sort(b,2);

B = false(size(A));
B(b) = true;
if isstruct(pdstruct)
    pdstruct.MASK = double(B);
    pdstruct.flood_points= b;
else
    pdstruct = double(B);
end

%imagesc(B(:,:,slice));