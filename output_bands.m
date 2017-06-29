%Generate band structure and output in format that gnuplot can plot
function[]=output_bands(HS1,HS2,eigen,name)


global nkpt compound ic N

%Number of k-vectors
Nkv=size(HS1,1);
%Number of strips/segments 
Nstrips=nkpt-1;

%File name
fname=strcat(compound(ic,:),'_eigenfull_',name,'.dat');      
%Open band structure file
fid=fopen(fname,'w');
           

%Initialise continuous k-value
kr=0.;
cnt=0;
%Loop over k-vector
for ikv=1:Nkv
    
    %Start vector
    ki=HS1(ikv,:);
    %End vector
    kf=HS2(ikv,:);
    %Displacement vector
    delk=(kf-ki);  
    %k-spacing= displacement vector/NStrips
    dk=delk/double(Nstrips);
   
    
    %Iterate over sampling points on the displacement vector
    %nkpt-1 to prevent get high symmetry points included twice
    for ik=1:nkpt-1
        
        %Continuous k index
        cnt=cnt+1;
        %k-vector 
        kvec=HS1(ikv,:)+(ik-1)*dk;
        %Magnitude of displacement k-vector(Euclidean norm)
        kmag=norm(dk);       
        %Continous k magnitude for plotting purposes
        kr=kr+ kmag;
           
        %Output continuous k-value and eigenvalues
        fprintf(fid, '%f\t',kr);
        for iprt=1:N
            fprintf(fid, '%f\t',eigen(cnt,iprt));
        end
        %New line
        fprintf(fid,'\n');
        
    end
end

%Close band structure file
fclose(fid);


end