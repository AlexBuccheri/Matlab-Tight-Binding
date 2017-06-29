%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Function to generate kpoint grid.
%Inputs:
% HS1       Starting vector for each displacement vector
% HS2       Ending vector for each displacement vector
% nkpt      Number of k-points per displacement vector
%Output:
% kpoint    k-point grid of dim(Nkpt*Nkv,3)
% hk_index  Indices of high-symmetry k-points in continuous
%           k index (w.r.t. vectors of HS1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[kpoint,hk_index]=gen_kgrid(nkpt,HS1,HS2)

    %Number of k-vectors (4 for this path through the BZ)
    Nkv=size(HS1,1);
    %Number of strips/segments 
    Nstrips=nkpt-1;

    %Loop over k-vector
    cnt=0;
    for ikv=1:Nkv
    
        %Start vector
        ki=HS1(ikv,:);
        %End vector
        kf=HS2(ikv,:);
    
        %Displacement vector
        delk=(kf-ki);
   
        %k-spacing= displacement vector/NStrips
        dk=delk/double(Nstrips);
       
        %Index map for high-symmetry k-points
        hk_index(ikv)=cnt+1;
        
        %Iterate over sampling points on the displacement vector
        %nkpt-1 to prevent get high symmetry points included twice
        for ik=1:nkpt-1
            cnt=cnt+1;
            kpoint(cnt,:)=HS1(ikv,:)+(ik-1)*dk;         
        end
    end

   return
end