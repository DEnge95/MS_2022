function params = fitParaboloid(phiIn,maskIn,X,Y)
% least square fit of a paraboloid to phase data

% save this for offline test
% save myVars

% desired output vector
index = find(maskIn > 0);
outhat = phiIn(index);

% input vectors (all column vectors
u = X(index);
u2 = u.^2;
v = Y(index);
v2 = v.^2;
uv = u.*v;

% solve the system outhat = [u v u2 v2 uv 1] * beta

% calculate the pseudo-inverse matrix
pInvMatrix = pinv([u v u2 v2 uv ones(length(u),1)]);

% solve the system
params = pInvMatrix * outhat;

% plot the surface to show how good we are
% surf(params(1)*X + params(2)*Y + params(3) * (X.^2) + params(4) * (Y.^2) + params(5) * (X.*Y) + params(6))