function [blocks] = blockify( corr_mat,corr_threshold, min_size )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
mat=corr_mat>=corr_threshold;
blocks = false(size(corr_mat));
nT=length(mat);
for iT=1:nT
   %forward
   pot=false(1,nT);
   pot(iT)=1;
   iF=iT+1;
   while iF<=nT
       if mat(iT,iF)
           pot(iF)=true;
           iF=iF+1;
       else
           break
       end
   end
    %backward
    iB=iT-1;
    while iB>1
        if mat(iT,iB)
            pot(iB)=true;
            iB=iB-1;
        else
            break
        end
    end
    if nnz(pot)>=min_size
        start =iB+1;
        stop = iF-1;
        blocks(start:stop,start:stop)=true;
    end

end

