function dataStruct = phase_unwrap_manual(dataStruct, matr,actPhaseNum, ind, dirFE)
%%% function to unwrap phase manually

     %% change 
     dataStruct.imaFLOW (:,:,actPhaseNum,dirFE)= matr;

     [I, J]= ind2sub(size(matr),ind);
     for j=1:size(ind)
         tmpVec  = squeeze(dataStruct.imaFLOW(I(j),J(j),:,dirFE));
         refVec  = abs(tmpVec).*sign(matr(ind(j)));
         tmFrames = find(abs(tmpVec-refVec)>= dataStruct.venc(dirFE)); 

         if ~isempty(tmFrames)
             tmpVec(tmFrames)= tmpVec(tmFrames)-(sign(tmpVec(tmFrames))*(dataStruct.venc(dirFE))*2);
             dataStruct.imaFLOW(I(j),J(j),:,dirFE)= tmpVec;
         end
     end
     
        