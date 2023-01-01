function pdStruct = update_pressure(pdStruct)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%    Initial Setup (Simple Integration)      %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
num_counted = 0;
counted = zeros(pdStruct.npts,1);
start_cnt = 1;
new_points =1;
pdStruct.PRESSURE = zeros(pdStruct.npts,1);
regions = 0;

dataDim = ndims(pdStruct.MASK);
sz =3;

signVx =  1;
signVy =  1;
signVz =  1;

while num_counted < pdStruct.npts     
    
    actpoint = pdStruct.flood_points(1,1);
    start_pos =  find(pdStruct.plist == actpoint);%find(counted==0,1);
    num_counted = num_counted + 1;
    counted(start_pos) = 1; 
    PRESSURE(start_pos)= 0;
    cnt_idx(num_counted)=start_pos;
    new_points = 100;
    
    regions = regions + 1;
    region_idx(regions) = start_pos;
%     disp(['Start Region Number ',int2str(regions)]);
    
    while new_points > 0

        num_countedp = num_counted;
        for pos = 1:num_counted

            idx = cnt_idx(pos);

            xp = pdStruct.nbs(idx,1);
            xm = pdStruct.nbs(idx,2);
            yp = pdStruct.nbs(idx,3);
            ym = pdStruct.nbs(idx,4);
            if dataDim==sz
                zp = pdStruct.nbs(idx,5);
                zm = pdStruct.nbs(idx,6);
            end

            if xp ~= 0 && counted(xp)==0
                counted(xp) = 1;
                num_counted = num_counted+1;
                cnt_idx(num_counted)=xp;
                 pdStruct.PRESSURE(xp) = pdStruct.PRESSURE(idx) + signVx*pdStruct.GRADx(idx)*pdStruct.delX/1000;
            end

            if xm ~= 0 && counted(xm)==0
                    counted(xm) = 1;
                    num_counted = num_counted+1;
                    cnt_idx(num_counted)=xm;
                    pdStruct.PRESSURE(xm) = pdStruct.PRESSURE(idx) - signVx*pdStruct.GRADx(idx)*pdStruct.delX/1000;
            end

            if yp ~= 0 &&  counted(yp)==0
                    counted(yp) = 1;
                    num_counted = num_counted+1;
                    cnt_idx(num_counted)=yp;
                    pdStruct.PRESSURE(yp) = pdStruct.PRESSURE(idx) + signVy*pdStruct.GRADy(idx)*pdStruct.delY/1000;
            end

            if ym ~= 0 && counted(ym)==0
                 counted(ym) = 1;
                 num_counted = num_counted+1;
                 cnt_idx(num_counted)=ym;
                    pdStruct.PRESSURE(ym) = pdStruct.PRESSURE(idx) - signVy*pdStruct.GRADy(idx)*pdStruct.delY/1000;
            end
            
            if dataDim==sz
                if zp ~= 0 && counted(zp)==0
                        counted(zp) = 1;
                        num_counted = num_counted+1;
                        cnt_idx(num_counted)=zp;
                        pdStruct.PRESSURE(zp) = pdStruct.PRESSURE(idx) + signVz*pdStruct.GRADz(idx)*pdStruct.delZ/1000;
                end

                if zm ~= 0 && counted(zm)==0
                        counted(zm) = 1;
                        num_counted = num_counted+1;
                        cnt_idx(num_counted)=zm;
                        pdStruct.PRESSURE(zm) = pdStruct.PRESSURE(idx) - signVz*pdStruct.GRADz(idx)*pdStruct.delZ/1000;
                end
            end
        end

        
        new_points = num_counted - num_countedp;
        
%           if(pdStruct.verbose)
%             disp(['Points Counted ',num2str(num_counted),' of ',num2str(pdStruct.npts)]);;
%           end
%         
    end
end
% % % % % % % % % %%
% % % % % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % % % %%%    Initial Setup (Simple Integration)      %%
% % % % % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % % % pdStruct.PRESSURE = zeros(pdStruct.npts,1);
% % % % % % %absGRAD = sqrt((pdStruct.GRADx).^2+(pdStruct.GRADy).^2+(pdStruct.GRADz).^2);
% % % % % % for pos=2:pdStruct.npts
% % % % % %     
% % % % % %     actpoint = pdStruct.flood_points(1,pos);
% % % % % %     idx = find(pdStruct.plist == actpoint);
% % % % % %     
% % % % % %     [actx acty actz]= ind2sub(pdStruct.DIM,actpoint);
% % % % % %     
% % % % % %     pdStruct.PRESSURE(idx,1) = pdStruct.GRADx(idx)*pdStruct.delX/1000+ pdStruct.GRADy(idx)*pdStruct.delY/1000+ pdStruct.GRADz(idx)*pdStruct.delZ/1000; 
% % % % % %     
% % % % % %     xp = sub2ind (pdStruct.DIM, actx+1, acty, actz);
% % % % % %     xp = find (pdStruct.plist == xp);
% % % % % %     if ~isempty(xp)
% % % % % %         pdStruct.PRESSURE(idx,1) = pdStruct.PRESSURE(idx,1)+ pdStruct.PRESSURE(xp,1);
% % % % % %     end
% % % % % %     
% % % % % %     xm = sub2ind (pdStruct.DIM, actx-1, acty, actz);
% % % % % %     xm = find (pdStruct.plist == xm);
% % % % % %     if ~isempty(xm)
% % % % % %         pdStruct.PRESSURE(idx,1) = pdStruct.PRESSURE(idx,1)+ pdStruct.PRESSURE(xm,1);
% % % % % %     end
% % % % % %     
% % % % % %     yp = sub2ind (pdStruct.DIM, actx, acty+1, actz);
% % % % % %     yp = find (pdStruct.plist == yp);
% % % % % %     if ~isempty(yp)
% % % % % %         pdStruct.PRESSURE(idx,1) = pdStruct.PRESSURE(idx,1)+ pdStruct.PRESSURE(yp,1);
% % % % % %     end
% % % % % %     
% % % % % %     ym = sub2ind (pdStruct.DIM, actx, acty-1, actz);
% % % % % %     ym = find (pdStruct.plist == ym);
% % % % % %     if ~isempty(ym)
% % % % % %         pdStruct.PRESSURE(idx,1) = pdStruct.PRESSURE(idx,1)+ pdStruct.PRESSURE(ym,1);
% % % % % %     end
% % % % % %     
% % % % % %     zp = sub2ind (pdStruct.DIM, actx, acty, actz+1);
% % % % % %     zp = find (pdStruct.plist == zp);
% % % % % %     if ~isempty(zp)
% % % % % %         pdStruct.PRESSURE(idx,1) = pdStruct.PRESSURE(idx,1)+ pdStruct.PRESSURE(zp,1);
% % % % % %     end
% % % % % %     
% % % % % %     zm = sub2ind (pdStruct.DIM, actx, acty, actz-1);
% % % % % %     zm = find (pdStruct.plist == zm);
% % % % % %     if ~isempty(zm)
% % % % % %         pdStruct.PRESSURE(idx,1) = pdStruct.PRESSURE(idx,1)+pdStruct.PRESSURE(zm,1);
% % % % % %     end 
% % % % % %     
% % % % % % end
%%        
% %%%%%%%%%%%Step 3 ITERATE TO CLEAN UP
 PRESSURE_old = pdStruct.PRESSURE;
 errord = 1e99;
iteration = 0;
pdStruct.max_error =0.001;
pdStruct.max_iter = 100;
 meandiff = 1e99;
while (errord > pdStruct.max_error) && (iteration < pdStruct.max_iter)
    
      iteration = iteration +1;    
    PRESSURE_old = pdStruct.PRESSURE;
      meandiff_old = meandiff;
    for pos = 1:pdStruct.npts
        
%         if sum(pos==region_idx)>0
%             pdStruct.PRESSURE(idx)=0.0;
%         else
        idx = find((pdStruct.plist) ==(pdStruct.flood_points(1,pos)));
        %start_pos =   == actpoint);
        
        xp = pdStruct.nbs(idx,1);
        xm = pdStruct.nbs(idx,2);
        yp = pdStruct.nbs(idx,3);
        ym = pdStruct.nbs(idx,4);
        if dataDim==sz
            zp = pdStruct.nbs(idx,5);
            zm = pdStruct.nbs(idx,6);
        end
%         idx = pos;
%         xp = pdStruct.nbs(pos,1);
%         xm = pdStruct.nbs(pos,2);
%         yp = pdStruct.nbs(pos,3);
%         ym = pdStruct.nbs(pos,4);
%         zp = pdStruct.nbs(pos,5);
%         zm = pdStruct.nbs(pos,6);
        
        GRAD_TERM = 0;
        nb_count = 0;
        
        if xp ~= 0
             GRAD_TERM = GRAD_TERM + PRESSURE_old(xp) - pdStruct.GRADx(idx)*pdStruct.delX/1000;
             nb_count = nb_count + 1;
        end

        if xm ~= 0
             GRAD_TERM =  GRAD_TERM + PRESSURE_old(xm) + pdStruct.GRADx(idx)*pdStruct.delX/1000;
             nb_count = nb_count + 1;
        end

        if yp ~= 0
             GRAD_TERM =  GRAD_TERM + PRESSURE_old(yp) - pdStruct.GRADy(idx)*pdStruct.delY/1000;
             nb_count = nb_count + 1;
        end

        if ym ~= 0
             GRAD_TERM =  GRAD_TERM + PRESSURE_old(ym) + pdStruct.GRADy(idx)*pdStruct.delY/1000;
             nb_count = nb_count + 1;
        end
        if dataDim==sz
            if zp ~= 0
                 GRAD_TERM =  GRAD_TERM + PRESSURE_old(zp) - pdStruct.GRADz(idx)*pdStruct.delZ/1000;
                 nb_count = nb_count + 1;
            end

            if zm ~= 0
                 GRAD_TERM =  GRAD_TERM + PRESSURE_old(zm) + pdStruct.GRADz(idx)*pdStruct.delZ/1000;
                 nb_count = nb_count + 1;
            end
        end
        
        if nb_count >0
            pdStruct.PRESSURE(idx) = (1- pdStruct.alpha)*PRESSURE_old(idx) + pdStruct.alpha*GRAD_TERM/nb_count;
        end
    end
% end
        
    %errord = sum( sqrt( (pdStruct.PRESSURE(:)- PRESSURE_old(:)).^2) )/sqrt(numel(pdStruct.PRESSURE));
    meandiff = mean(pdStruct.PRESSURE(:))-mean( PRESSURE_old(:));
    errord = abs(meandiff - meandiff_old);%.^2) )/sqrt(numel(pdStruct.PRESSURE));
    
    if(pdStruct.verbose)
    disp(['Error ',num2str(errord)])
    end
    pdStruct.tmp(1,1)=iteration;
    pdStruct.tmp(1,2)=errord;
 end

%  end

% actpoint = pdStruct.flood_points(1,1);
% start_pos =  find(pdStruct.plist == actpoint);
%pdStruct.PRESSURE(start_pos)= -13320;
pdStruct.PRESSURE_MAT = single(133320*ones(pdStruct.DIM));
%pdStruct.PRESSURE_MAT = NaN(pdStruct.DIM);
pdStruct.PRESSURE_MAT(pdStruct.plist) = single(pdStruct.PRESSURE);

if( pdStruct.GRAD_OUT == 1)
    pdStruct.GRADx_MAT = single(zeros(pdStruct.DIM));
    pdStruct.GRADy_MAT = single(zeros(pdStruct.DIM));
    pdStruct.GRADz_MAT = single(zeros(pdStruct.DIM));
    
    pdStruct.GRADx_MAT(pdStruct.plist) = single(pdStruct.GRADx);
    pdStruct.GRADy_MAT(pdStruct.plist) = single(pdStruct.GRADy);
    if dataDim==sz
        pdStruct.GRADz_MAT(pdStruct.plist) = single(pdStruct.GRADz);
    end
end