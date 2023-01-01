function maskdata = thresholding_percent(data,tvalue)
% get maximum
mx = max(data(:));
tvalue = mx*tvalue;

maskdata = (data<=tvalue);
dim =ndims(data);
if dim>2
    maskdata = all(maskdata,dim);
end
