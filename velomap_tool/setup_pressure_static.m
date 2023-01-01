function pdStruct = setup_pressure_static(pdStruct)

%%%%Fluid Properties
pdStruct.visc = 0/1000;   %% Viscosity in Cp divided by 1000 to get SI units
pdStruct.dens = 1060.0;      %% Fluid density in kg/m^3
%%%Calulation Parametrs
pdStruct.max_iter = 50;    %% Max pressure iterations
pdStruct.alpha = 0.2;        %% Under relaxation factor for iteration
pdStruct.poly_num = 2.0;     %% Order of polynomial fitting
pdStruct.max_error = 0.0001;   %% Stopping criteria for Iterations

%%%%Data Properties
% pdStruct.delX=1.0417;           %% Grid spacing in x
% pdStruct.delY=1.0417;           %% Grid spacing in y
% pdStruct.delZ=1.0;           %% Grid spacing in z
% pdStruct.tres=22;            %% Time spacing in s
pdStruct.tvals=3;            %% Number of points to use

%%%%Determine What Points Are Required For Calc
pdStruct.DIM = size(pdStruct.MASK);
pdStruct.plist = find( pdStruct.MASK );
pdStruct.npts  = length(pdStruct.plist);
[pdStruct.xlist pdStruct.ylist pdStruct.zlist] = ind2sub(size(pdStruct.MASK),pdStruct.plist);

%%%%%% Get Neighbors of Each Pixel %%%%%%%
if ndims(pdStruct.MASK)==3 
    pdStruct.nbs=zeros(pdStruct.npts,6);
else
    pdStruct.nbs=zeros(pdStruct.npts,4);
end
pdStruct.MASK(pdStruct.plist)=1:pdStruct.npts;
    
%%%%%%%%Mask has vales of indices %%%%%
for pos = 1:pdStruct.npts
%     if(pdStruct.verbose)
%        if 500*floor(pos/500)==pos
%             disp(['NBS:Point Number ',int2str(pos),' of ',int2str(pdStruct.npts)]);
%             drawnow;
%        end
%     end
    
    x0 = pdStruct.xlist(pos);
    y0 = pdStruct.ylist(pos);
    z0 = pdStruct.zlist(pos);

    if(pdStruct.MASK(x0+1,y0,z0) > 0)
        pdStruct.nbs(pos,1)=pdStruct.MASK(x0+1,y0,z0);
    end
    
    if(pdStruct.MASK(x0-1,y0,z0) > 0)
        pdStruct.nbs(pos,2)=pdStruct.MASK(x0-1,y0,z0);
    end
    
    if(pdStruct.MASK(x0,y0+1,z0) > 0)
        pdStruct.nbs(pos,3)=pdStruct.MASK(x0,y0+1,z0);
    end
    
    if(pdStruct.MASK(x0,y0-1,z0) > 0)
        pdStruct.nbs(pos,4)=pdStruct.MASK(x0,y0-1,z0);
    end
    if ndims(pdStruct.MASK)==3
        if(pdStruct.MASK(x0,y0,z0+1) > 0)
            pdStruct.nbs(pos,5)=pdStruct.MASK(x0,y0,z0+1);
        end

        if(pdStruct.MASK(x0,y0,z0-1) > 0)
            pdStruct.nbs(pos,6)=pdStruct.MASK(x0,y0,z0-1);
        end    
    end
end

pdStruct.MASK = (pdStruct.MASK >0);