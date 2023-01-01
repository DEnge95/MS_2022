%
% least square fit of plane to mag weighted 2d pc data
%

function [phi0,alpha,beta] = lsqf2(phiIn,magIn,X,Y)

magSq = magIn.^2;

M1(1,1) = sum(sum(magSq.*phiIn));
M1(2,1) = sum(sum(magSq.*phiIn.*X));
M1(3,1) = sum(sum(magSq.*phiIn.*Y));
M1(1,2) = sum(sum(magSq.*X));
M1(2,2) = sum(sum(magSq.*X.^2));
M1(3,2) = sum(sum(magSq.*X.*Y));
M1(1,3) = sum(sum(magSq.*Y));
M1(2,3) = sum(sum(magSq.*Y.*X));
M1(3,3) = sum(sum(magSq.*Y.^2));


M2(1,1) = sum(sum(magSq));
M2(2,1) = sum(sum(magSq.*X));
M2(3,1) = sum(sum(magSq.*Y));
M2(1,2) = sum(sum(magSq.*phiIn));
M2(2,2) = sum(sum(magSq.*phiIn.*X));
M2(3,2) = sum(sum(magSq.*phiIn.*Y));
M2(1,3) = sum(sum(magSq.*Y));
M2(2,3) = sum(sum(magSq.*Y.*X));
M2(3,3) = sum(sum(magSq.*Y.^2));


M3(1,1) = sum(sum(magSq));
M3(2,1) = sum(sum(magSq.*X));
M3(3,1) = sum(sum(magSq.*Y));
M3(1,2) = sum(sum(magSq.*X));
M3(2,2) = sum(sum(magSq.*X.^2));
M3(3,2) = sum(sum(magSq.*X.*Y));
M3(1,3) = sum(sum(magSq.*phiIn));
M3(2,3) = sum(sum(magSq.*phiIn.*X));
M3(3,3) = sum(sum(magSq.*phiIn.*Y));


Md(1,1) = sum(sum(magSq));
Md(2,1) = sum(sum(magSq.*X));
Md(3,1) = sum(sum(magSq.*Y));
Md(1,2) = sum(sum(magSq.*X));
Md(2,2) = sum(sum(magSq.*X.^2));
Md(3,2) = sum(sum(magSq.*X.*Y));
Md(1,3) = sum(sum(magSq.*Y));
Md(2,3) = sum(sum(magSq.*Y.*X));
Md(3,3) = sum(sum(magSq.*Y.^2));

phi0  = det(M1)/det(Md);
alpha = det(M2)/det(Md);
beta  = det(M3)/det(Md);
