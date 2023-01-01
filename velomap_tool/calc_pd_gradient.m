function pdStruct = calc_pd_gradient( pdStruct,time,encodingType)
    derivative_edge_type = 1;
    dataDim = ndims(pdStruct.MASK);
    sz =3;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%    Calculates the Pressure gradient and stores in
        %%%            pdstruct.GRADx(y/z)
        %%%%        teval is point in series to evaluate
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %%%DERIVATIVE CONVERSIONS
        conv_vel = 1; %mm/s to m/s
%         conv_vel2= 1/1000/1000;
%         conv_dx = 1/pdStruct.delX*1000; % mm to m
%         conv_dx2= 1/pdStruct.delX/pdStruct.delX*1000*1000;
%         conv_dy = 1/pdStruct.delY*1000;
%         conv_dy2= 1/pdStruct.delY/pdStruct.delY*1000*1000;
%         conv_dz = 1/pdStruct.delZ*1000;
%         conv_dz2= 1/pdStruct.delZ/pdStruct.delZ*1000*1000;
%         conv_dt = 1/pdStruct.tres*1000; %%ms to s


        %%%Loop Over Points
        for pos = 1:pdStruct.npts

            %     if( pdStruct.)
            %         if 50*floor(pos/50)==pos
            %             disp(['Point Number ',int2str(pos),' of ',int2str(pdStruct.npts)])
            %         end
            %     end

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%   Step 1 Derivative Space                        %%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


            %%%%Setup Derivative Space For Derivatives
            [x0 y0 z0] = ind2sub(pdStruct.DIM,pdStruct.plist(pos));


            if ( (x0 < 2 ) || (y0 < 2 ) || (x0 > pdStruct.DIM(1)-1) || (y0 > pdStruct.DIM(2)-1))
                %%KMJ These are uncomputable points
                GRADx(pos,:) = 0;
                GRADy(pos,:) = 0;
                if dataDim ==sz &&((z0 < 2)||(z0 > pdStruct.DIM(3)-1))
                    GRADz(pos,:) = 0;
                end
            else

                %%KMJ Time should =2
                time_plus1 = time+1;
                time_minus1= time-1;               

                %%% Either assume the vessel stops if the mask is zero or
                %%% assume that the velocity data outside the mask is good
                if derivative_edge_type == 1
                    mxp = pdStruct.MASK(x0+1,y0,z0);
                    myp = pdStruct.MASK(x0,y0+1,z0);                    
                    mxm = pdStruct.MASK(x0-1,y0,z0);
                    mym = pdStruct.MASK(x0,y0-1,z0);
                    if dataDim ==sz
                        mzm = pdStruct.MASK(x0,y0,z0-1);
                        mzp = pdStruct.MASK(x0,y0,z0+1);
                    end
                else
                    mxp = 1.0;
                    myp = 1.0;                    
                    mxm = 1.0;
                    mym = 1.0;
                    if dataDim ==sz
                        mzm = 1.0;
                        mzp = 1.0;
                    end
                end
                if strcmp(encodingType,'velocity')
                    %%%velocity Terms
                    vx = conv_vel*pdStruct.VELXt(x0,y0,z0,time);
                    vy = conv_vel*pdStruct.VELYt(x0,y0,z0,time);
                              
                    %%First Derivatives
                    dvxdt = conv_vel*( pdStruct.VELXt(x0,y0,z0,time_plus1) - pdStruct.VELXt(x0,y0,z0,time_minus1) )/(2*pdStruct.tres/1000);
                    dvydt = conv_vel*( pdStruct.VELYt(x0,y0,z0,time_plus1) - pdStruct.VELYt(x0,y0,z0,time_minus1) )/(2*pdStruct.tres/1000);
                    
                    dvxdx = conv_vel*( pdStruct.VELXt(x0+1,y0,z0,time)*mxp - pdStruct.VELXt(x0-1,y0,z0,time)*mxm )/(2*pdStruct.delX/1000);
                    dvxdy = conv_vel*( pdStruct.VELXt(x0,y0+1,z0,time)*myp - pdStruct.VELXt(x0,y0-1,z0,time)*mym )/(2*pdStruct.delY/1000);
                    
                    dvydx = conv_vel*( pdStruct.VELYt(x0+1,y0,z0,time)*mxp - pdStruct.VELYt(x0-1,y0,z0,time)*mxm )/(2*pdStruct.delX/1000);
                    dvydy = conv_vel*( pdStruct.VELYt(x0,y0+1,z0,time)*myp - pdStruct.VELYt(x0,y0-1,z0,time)*mym )/(2*pdStruct.delY/1000);
                    
                    dvzdx = conv_vel*( pdStruct.VELZt(x0+1,y0,z0,time)*mxp - pdStruct.VELZt(x0-1,y0,z0,time)*mxm )/(2*pdStruct.delX/1000);
                    dvzdy = conv_vel*( pdStruct.VELZt(x0,y0+1,z0,time)*myp - pdStruct.VELZt(x0,y0-1,z0,time)*mym )/(2*pdStruct.delY/1000);
                    
                    
                    if dataDim==sz
                        vz = conv_vel*pdStruct.VELZt(x0,y0,z0,time);  
                        dvzdt = conv_vel*( pdStruct.VELZt(x0,y0,z0,time_plus1) - pdStruct.VELZt(x0,y0,z0,time_minus1) )/(2*pdStruct.tres/1000);
                        
                        dvxdz = conv_vel*( pdStruct.VELXt(x0,y0,z0+1,time)*mzp - pdStruct.VELXt(x0,y0,z0-1,time)*mzm )/(2*pdStruct.delZ/1000);                    
                        dvydz = conv_vel*( pdStruct.VELYt(x0,y0,z0+1,time)*mzp - pdStruct.VELYt(x0,y0,z0-1,time)*mzm )/(2*pdStruct.delZ/1000);                    
                        dvzdz = conv_vel*( pdStruct.VELZt(x0,y0,z0+1,time)*mzp - pdStruct.VELZt(x0,y0,z0-1,time)*mzm )/(2*pdStruct.delZ/1000);
                    end

                    %%2nd Derivatives
                    dvxdx2 = conv_vel*( pdStruct.VELXt(x0+1,y0,z0,time)*mxp -2*pdStruct.VELXt(x0,y0,z0,time) + pdStruct.VELXt(x0-1,y0,z0,time)*mxm )/(pdStruct.delX/1000)^2;
                    dvxdy2 = conv_vel*( pdStruct.VELXt(x0,y0+1,z0,time)*myp -2*pdStruct.VELXt(x0,y0,z0,time) + pdStruct.VELXt(x0,y0-1,z0,time)*mym )/(pdStruct.delY/1000)^2;
                    
                    dvydx2 = conv_vel*( pdStruct.VELYt(x0+1,y0,z0,time)*mxp -2*pdStruct.VELYt(x0,y0,z0,time) + pdStruct.VELYt(x0-1,y0,z0,time)*mxm )/(pdStruct.delX/1000)^2;
                    dvydy2 = conv_vel*( pdStruct.VELYt(x0,y0+1,z0,time)*myp -2*pdStruct.VELYt(x0,y0,z0,time) + pdStruct.VELYt(x0,y0-1,z0,time)*mym )/(pdStruct.delY/1000)^2;
                    
                    dvzdx2 = conv_vel*( pdStruct.VELZt(x0+1,y0,z0,time)*mxp -2*pdStruct.VELZt(x0,y0,z0,time) + pdStruct.VELZt(x0-1,y0,z0,time)*mxm )/(pdStruct.delX/1000)^2;
                    dvzdy2 = conv_vel*( pdStruct.VELZt(x0,y0+1,z0,time)*myp -2*pdStruct.VELZt(x0,y0,z0,time) + pdStruct.VELZt(x0,y0-1,z0,time)*mym )/(pdStruct.delY/1000)^2;

                    if dataDim==sz
                        dvxdz2 = conv_vel*( pdStruct.VELXt(x0,y0,z0+1,time)*mzp -2*pdStruct.VELXt(x0,y0,z0,time) + pdStruct.VELXt(x0,y0,z0-1,time)*mzm )/(pdStruct.delZ/1000)^2;
                        dvydz2 = conv_vel*( pdStruct.VELYt(x0,y0,z0+1,time)*mzp -2*pdStruct.VELYt(x0,y0,z0,time) + pdStruct.VELYt(x0,y0,z0-1,time)*mzm )/(pdStruct.delZ/1000)^2;
                        dvzdz2 = conv_vel*( pdStruct.VELZt(x0,y0,z0+1,time)*mzp -2*pdStruct.VELZt(x0,y0,z0,time) + pdStruct.VELZt(x0,y0,z0-1,time)*mzm )/(pdStruct.delZ/1000)^2;
                    end
                    

                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %%%   Step 2 Navier-Stokes                           %%
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                    
                    if dataDim==sz
                        pdStruct.GRADx(pos,1) = single(-pdStruct.dens*dvxdt ...
                                                -pdStruct.dens*vx*dvxdx  ...
                                                -pdStruct.dens*vy*dvxdy  ...
                                                -pdStruct.dens*vz*dvxdz ...
                                                +pdStruct.visc*dvxdx2 ...
                                                +pdStruct.visc*dvxdy2 ...
                                                +pdStruct.visc*dvxdz2);
                        pdStruct.GRADy(pos,1) = single(-pdStruct.dens*dvydt  ...
                                                -pdStruct.dens*vx*dvydx ...
                                                -pdStruct.dens*vy*dvydy ...
                                                -pdStruct.dens*vz*dvydz ...
                                                +pdStruct.visc*dvydx2 ...
                                                +pdStruct.visc*dvydy2 ...
                                                +pdStruct.visc*dvydz2); 
                        pdStruct.GRADz(pos,1) = single(-pdStruct.dens*dvzdt ...
                                                -pdStruct.dens*vx*dvzdx ...
                                                -pdStruct.dens*vy*dvzdy ...
                                                -pdStruct.dens*vz*dvzdz  ...
                                                +pdStruct.visc*dvzdx2 ...
                                                +pdStruct.visc*dvzdy2 ...
                                                +pdStruct.visc*dvzdz2);
                    else
                        pdStruct.GRADx(pos,1) = single(-pdStruct.dens*dvxdt ...
                                                -pdStruct.dens*vx*dvxdx  ...
                                                -pdStruct.dens*vy*dvxdy  ...                                                
                                                +pdStruct.visc*dvxdx2 ...
                                                +pdStruct.visc*dvxdy2);
                                                
                        pdStruct.GRADy(pos,1) = single(-pdStruct.dens*dvydt  ...
                                                -pdStruct.dens*vx*dvydx ...
                                                -pdStruct.dens*vy*dvydy ...                                                
                                                +pdStruct.visc*dvydx2 ...
                                                +pdStruct.visc*dvydy2) ;
                    end
                    
                elseif strcmp(encodingType, 'acceleration')
                    %%% acceleration terms
                    ax = conv_vel*pdStruct.VELXt(x0,y0,z0,time);
                    ay = conv_vel*pdStruct.VELYt(x0,y0,z0,time);                    
                    
                    
%                     %%2nd Derivatives
%                     dvxdx2 = conv_vel*( pdStruct.VELXt(x0+1,y0,z0,time)*mxp -2*pdStruct.VELXt(x0,y0,z0,time) + pdStruct.VELXt(x0-1,y0,z0,time)*mxm )/(pdStruct.delX/1000)^2;
%                     dvxdy2 = conv_vel*( pdStruct.VELXt(x0,y0+1,z0,time)*myp -2*pdStruct.VELXt(x0,y0,z0,time) + pdStruct.VELXt(x0,y0-1,z0,time)*mym )/(pdStruct.delY/1000)^2;
%                     dvxdz2 = conv_vel*( pdStruct.VELXt(x0,y0,z0+1,time)*mzp -2*pdStruct.VELXt(x0,y0,z0,time) + pdStruct.VELXt(x0,y0,z0-1,time)*mzm )/(pdStruct.delZ/1000)^2;
% 
%                     dvydx2 = conv_vel*( pdStruct.VELYt(x0+1,y0,z0,time)*mxp -2*pdStruct.VELYt(x0,y0,z0,time) + pdStruct.VELYt(x0-1,y0,z0,time)*mxm )/(pdStruct.delX/1000)^2;
%                     dvydy2 = conv_vel*( pdStruct.VELYt(x0,y0+1,z0,time)*myp -2*pdStruct.VELYt(x0,y0,z0,time) + pdStruct.VELYt(x0,y0-1,z0,time)*mym )/(pdStruct.delY/1000)^2;
%                     dvydz2 = conv_vel*( pdStruct.VELYt(x0,y0,z0+1,time)*mzp -2*pdStruct.VELYt(x0,y0,z0,time) + pdStruct.VELYt(x0,y0,z0-1,time)*mzm )/(pdStruct.delZ/1000)^2;
% 
%                     dvzdx2 = conv_vel*( pdStruct.VELZt(x0+1,y0,z0,time)*mxp -2*pdStruct.VELZt(x0,y0,z0,time) + pdStruct.VELZt(x0-1,y0,z0,time)*mxm )/(pdStruct.delX/1000)^2;
%                     dvzdy2 = conv_vel*( pdStruct.VELZt(x0,y0+1,z0,time)*myp -2*pdStruct.VELZt(x0,y0,z0,time) + pdStruct.VELZt(x0,y0-1,z0,time)*mym )/(pdStruct.delY/1000)^2;
%                     dvzdz2 = conv_vel*( pdStruct.VELZt(x0,y0,z0+1,time)*mzp -2*pdStruct.VELZt(x0,y0,z0,time) + pdStruct.VELZt(x0,y0,z0-1,time)*mzm )/(pdStruct.delZ/1000)^2;
% 
%                     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                     %%%   Step 2 Navier-Stokes                           %%
%                     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                    pdStruct.GRADx(pos,1) =   -pdStruct.dens*ax; 
%                         +pdStruct.visc*dvxdx2 ...
%                         +pdStruct.visc*dvxdy2 ...
%                         +pdStruct.visc*dvxdz2;
                    pdStruct.GRADy(pos,1) =  -pdStruct.dens*ay;  
%                         +pdStruct.visc*dvydx2 ...
%                         +pdStruct.visc*dvydy2 ...
%                         +pdStruct.visc*dvydz2;
                    if dataDim==sz
                        az = conv_vel*pdStruct.VELZt(x0,y0,z0,time);
                        pdStruct.GRADz(pos,1) =  -pdStruct.dens*az;
%                             +pdStruct.visc*dvzdx2 ...
%                             +pdStruct.visc*dvzdy2 ...
%                             +pdStruct.visc*dvzdz2;
                    end
                end
            end
        end
