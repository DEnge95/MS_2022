%%  function for phase unwrapping
function dataStruct = flow_anti_aliasing(dataStruct,numOfIter)
     venc       = dataStruct.venc;
     
 for w=1:numOfIter    
    for v = 1:dataStruct.nFE      
      
         %% phase unwrapping in XY-plane
         for j=1:5
            tmpMX              = squeeze(dataStruct.imaFLOW(:,:,:,v));
            %% permuted array is needed for searching of phases with jumps
            flowMX             = permute(tmpMX,[3 1 2]);
            diffMX             = diff(flowMX,1,3);
            %% find if there are any phases with jumps > venc
            [phases,indXmx]    = find(diffMX <= -venc(v) |diffMX >= venc(v));
            clear ('diffMX','flowMX')
            if ~isempty(phases)
                %% calculate index of pixel to unwrap
                indXmx          = mod(indXmx,size(tmpMX,1));
                
                %% create array of coordinates [index of phase, X-index] of pixels to unwrap             
                indMX           = zeros(size(phases,1),2);
                indMX(:,1)      = phases;
                indMX(:,2)      = indXmx;
                [b,m,n]         = unique(indMX,'rows'); %% delete repeated pixel-coordinates
                indMX           = indMX(m,:);
                
                if ~isempty(indMX)
                   d = size(indMX,1); 
                   for i=1:d
                       indX       = indMX(i,2);
                       indPhase   = indMX(i,1);
                       
                       %% if pixel to unwrap is not in the first row of array                                             
                       if indX ~= 0
                         %% load in tmpVec the row with 'wraped-pixel' and the neighbour-row
                          tmpVec     = squeeze(dataStruct.imaFLOW(indX:indX+1,:,indPhase,v));
                         %% calculate if standard deviation is above a defined value                          
                          stdMX      = std(tmpVec,0,1);
                          [indRow indCol]   = find(stdMX >= venc(v)*2/3);
                         %% if the column components have the equal sign
                          equalSign   = find(abs(diff(sign(tmpVec(:,indCol)),1,1))==2);
                          
                          if ~isempty(equalSign)
                              diffVec    = diff(tmpVec(2,:),1,2);
                              diffVec    = find(diffVec<= -venc(v) |diffVec >= venc(v))+1;
                            %% if index of columns are equal 
                              indCol     = intersect(diffVec,indCol);
                              if ~isempty(indCol)
                                  %% change value of the pixel
                                  value    = squeeze(dataStruct.imaFLOW(indX,indCol,indPhase,v));
                                  value    = value - (sign(value)*venc(v)*2);
                                  dataStruct.imaFLOW(indX,indCol,indPhase,v)= value;
                              end
                          end
                       end
                   end
                end
            end
         end
         %% end of:  phase unwrapping in XY-plane

  
  
         %% phase unwrapping in temporal direction
  
         %% find if there are any phases with |jumps| > venc
         tmpMX              = squeeze(dataStruct.imaFLOW(:,:,:,v));         
         diffMX             = diff(tmpMX,1,3);
         [phases,x]         = find(diffMX <= -venc(v) |diffMX >= venc(v));

         %% correct jumps of one data pixel above all phases
         if ~isempty(phases)
            %% calculate difference along temporal direction 
            %% and find coordinates of pixels with |jumps| > venc              
            [indXmx,indYmx] = find(diffMX <= -venc(v) |diffMX >= venc(v));
            %% calculate y index 
            [x y z]         = size(tmpMX);
            indYmx          = mod(indYmx,y);
            index           = find(indYmx==0);
            indYmx(index)   = y;
            clear ('diffMX')
            %% create array of coordinates [X-index, Y-index] of
            %% pixels to unwrap
            indMX           = zeros(size(indXmx,1),2);
            indMX(:,1)      = indXmx;
            indMX(:,2)      = indYmx;
            [b,m,n]         = unique(indMX,'rows');  
            indMX           = flipud(indMX(m,:));
            
            if ~isempty(indMX)
                d = size(indMX,1);
                for i=1:d
                  indX            = indMX(i,1);
                  indY            = indMX(i,2);
                  %% if pixel coordinates are not on the edge
                  if indY>1 && indY<y && indX>1 && indX<x
                     %% load an array of voxel with wrapped phase and of 8 neighbour-voxels  
                     tmpVec          = squeeze(dataStruct.imaFLOW(indX-1:indX+1,indY-1:indY+1,:,v));
                     tmpVec          = reshape(tmpVec,9,[]);
                     %% calculate difference along temporal direction 
                     %% and find columns with jumps >= venc
                     diffVec         = diff(tmpVec,1,2);
                     [indRow,indCol] = find(diffVec<=-venc(v)|diffVec >= venc(v));                     
                     indCol          = unique(indCol);
                     
                     if ~isempty(indCol)
                     %% calculate whether the gradient(= vector of differences)has any sign-changes
                     %% which means there is an extremum
                        extrema     = diff(sign(diffVec),1,2);
                        extrema     = [zeros(9,1), extrema];
                     %% consider only the part from data 
                     %% from the index before first jump > venc until the index after the last jump 
                        permitInd   = union((indCol-1),(indCol+1));
                        if permitInd(1)~=0 && permitInd(length(permitInd))<((dataStruct.numPhases)-1)
                         extremaCut  = extrema(:,permitInd(1):permitInd(length(permitInd)));
                         %% find index of the curve with only one extremum 
                         refCurveInd = sum(abs(extremaCut),2);
                         refCurveInd = find(refCurveInd == 2);
                         if ~isempty(refCurveInd)
                            %%choice of reference curve 
                            absKurve = abs(tmpVec(refCurveInd,:));
                            [refIndX,refIndY] = find(absKurve == max(max(absKurve)));
                            %%% search for indices of phase to unwrap
                            index           = union(refCurveInd(refIndX),5);
                            phaseToChange   = tmpVec(index,:);
                            phaseToChange   = std(phaseToChange,0,1);
                            phaseToChange   = find(phaseToChange > venc(v));
                            %% unwrapping values
                            values          = tmpVec(5,phaseToChange);
                            for k=1:length(values)
                                values(k)    = values(k) - (sign(values(k))*venc(v)*2);
                            end
                            tmpVec(5,phaseToChange) = values;
                            %% overwrite data with corrected one
                            tmpVec = reshape(tmpVec,3,3,[]);
                            dataStruct.imaFLOW(indX-1:indX+1,indY-1:indY+1,:,v)= tmpVec;
                         end    
                        end
                     end
                  else
                    %% set pixels on the edge to zero
                     dataStruct.imaFLOW(indX,indY,:,v)= 0; 
                  end

                end
            end
         end
         %% end of :phase unwrapping in temporal direction
   
        %% phase unwrapping in XY-plane          
        tmpMX              = squeeze(dataStruct.imaFLOW(:,:,:,v));
        %% permuted array is needed for searching of phases with jumps
        flowMX             = permute(tmpMX,[3 2 1]);
        diffMX             = diff(flowMX,1,3);
        %% find if there are any phases with jumps > venc
        [phases,indYmx]    = find(diffMX <= -venc(v) |diffMX >= venc(v));
        clear ('diffMX','flowMX')
        if ~isempty(phases)            
            [x y z]         = size(tmpMX);
            %% calculate index of pixel to unwrap
            indYmx          = mod(indYmx,y);

            %% create array of coordinates [index of phase, X-index] of pixels to unwrap             
            indMX           = zeros(size(phases,1),2);
            indMX(:,1)      = phases;
            indMX(:,2)      = indYmx;
            [b,m,n]         = unique(indMX,'rows'); %% delete repeated pixel-coordinates
            indMX           = indMX(m,:);

            if ~isempty(indMX)
               d            = size(indMX,1); 
               for i=1:d
                   indY       = indMX(i,2);
                   indPhase   = indMX(i,1);

                   %% if pixel to unwrap is not in the first row of array                                             
                   if indY ~= 0
                     %% load in tmpVec the row with 'wraped-pixel' and the neighbour-row
                      tmpVec     = squeeze(dataStruct.imaFLOW(:,indY,indPhase,v));                        
                      diffVec    = diff(tmpVec,1,1);
                      indRow     = find(diffVec<= -venc(v) |diffVec >= venc(v));
                      
                      if ~isempty(indRow)
                        %% wenn Sprung hin und zurück stattfindet 
                        if length(indRow)==2 && (indRow(2)-indRow(1))<=3
                          for q=1:(indRow(2)-indRow(1))
                              indX =indRow(2)-q+1;
                              value    = squeeze(dataStruct.imaFLOW(indX,indY,indPhase,v));
                              value    = value - (sign(value)*venc(v)*2);
                              dataStruct.imaFLOW(indX,indY,indPhase,v)= value;
                          end
                        else
                          tmpVec     = squeeze(dataStruct.imaFLOW((indRow(1)):(indRow(1)+1),indY,:,v));                                  
                          stdVec     = std(tmpVec,0,2);
                          index      = find((ceil(stdVec))>= venc(v)*2/3);%nur eine oder beide Kurven haben Sprünge?
                          if ~isempty(index) 
                              for l=1:length(index)
                                  diffVec    = diff(tmpVec(l,:));
                                  [ind]      = find(diffVec<= -venc(v) |diffVec >= venc(v));
                                  if length(ind)==2 %% wenn es einen Sprung hin und zurück gibt
                                     indX = indRow(1)+index(l)-1; 
                                      for q=1:(ind(2)-ind(1))
                                          tmpIndPhase = ind(2)-q+1;
                                          value    = squeeze(dataStruct.imaFLOW(indX,indY,tmpIndPhase,v));
                                          value    = value - (sign(value)*venc(v)*2);
                                          dataStruct.imaFLOW(indX,indY,tmpIndPhase,v)= value;
                                      end
                                  end
                              end 
                          end
                        end
                      end
                   end
               end
            end
        end         
        %% end of:  phase unwrapping in XY-plane 
    end   
 end
%%% end of: function for phase unwrapping     