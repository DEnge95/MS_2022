function dataStruct = timeseries_fft(dataStruct,staticFactor)

velocitymag = squeeze(sqrt(sum((dataStruct.imaFLOW/100).^2,4)));
velfft = abs(fft(velocitymag,[],3));
velfft = squeeze(velfft(:,:,1));
mx = max(velfft(:));            
mask1 = velfft<=(staticFactor*mx);
mask2 = velfft>0;
staticMask = cat(3,mask1, mask2);
staticMask = all(staticMask,3);
%%% remove single pixels and close holes
nloops = 1;
staticMask = bwmorph(staticMask,'majority',nloops);

dataStruct.statMask = staticMask;
