%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Function to generate ideal eigenvalues for bulk system
%with SO coupling incorporated 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[eigen] = generate_eigenvalues(kpoint,hk_index,parameter,p_type)

global Eg_ideal N spin


%Declare eigenvalue array size: eigen(Nkpt,Nbands)
eigen=zeros(size(kpoint,1),N);


%Loop over k-points
for ik=1:size(kpoint,1)
  %Construct H and solve for eigenvalues, 
  H=ham(kpoint(ik,:),parameter);
  eigen(ik,:)=eig(H);
end


%For ideal eigenvalues, ensure the VBT is at 0 eV and the band-gap
%is equal to Eg. 
if strcmp(p_type,'ideal')==1
    %Highest VB index
    Nv=4*spin;
    %Gamma point in continuously-indexed k-grid
    igam=hk_index(2);
    %Valence band top at Gamma
    E_h1=eigen(igam,Nv);
    %Conduction band bottom at Gamma
    E_e1=eigen(igam,Nv+1);
    %Rigidly shift all valence bands, so VBT=0 eV
    eigen(:,1:Nv)=eigen(:,1:Nv)-E_h1;
    %Consistently modify position of all conduction states
    %so the bang-gap is equal to Eg
    E_h1=eigen(igam,Nv);
    dEc=Eg_ideal-(E_e1-E_h1);
    eigen(:,Nv+1:N)=eigen(:,Nv+1:N)+dEc;
    return
end


%For trial or final best-fit eigenvalues, one isn't required to do anything
 return    
end 
