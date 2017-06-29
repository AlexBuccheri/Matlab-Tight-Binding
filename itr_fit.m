
%%%%%%%%%%%%%%%%%%%%%%%%%%% Routine Description %%%%%%%%%%%%%%%%%%%%%%%%%%%
% Brute-force script to run the minimisation of a set of TB parameters 
% multiple times. Orthogonal NN sp3s* TB Model, with SO coupling.
% Note, SO parameters are not used in the fitting process
%
% Alex Buccheri. December 2016. University of Oxford. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%----------------------------
% User-defined Variables
%----------------------------
clear
global ic nkpt Eg_ideal Nsym b_weight parent_dir

%Specify material (index) of interest
ic=19;

%Number of k-sampling points (per displacement vector)
nkpt=25;

%Ideal band-gap at T=0K (eV)
Eg_ideal=1.85;

%Weighting of high-symmetry points in the fit w.r.t. all other k-points
Nsym=4.;

%Weightings for each band. VB=1-8 and CB=9-20
b_weight=[1., 1., 1., 1., 1., 1., 1., 1.,  ...
          1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.];

%Total Number of runs to make (vary seed parameters each time)
Nr=1;

%----------------------------
%Set Up Output directories
%----------------------------
%Get required material parameters
[Natoms,basis,spin,N,Np,Nc,db,al,compound,dv]=get_mat_param(ic);
        
%Directory to send outputs to
parent_dir=strcat(compound(ic,:),'_set1');

%Directory name per output
output_dir='output';

%Query existence of directory parent_dir
ios=exist(parent_dir,'dir');

%if directory doesn't exist, make it. 
if ios == 0
    %Make parent directory to contain the data from this run
    %Note, sprintf concantinates strings, preserving a space
    %between each string
    unix( sprintf('%s','mkdir ',parent_dir) )
%elseif ios==7
    %disp(sprintf('%s','Directory ',parent_dir,' exists.'))
    %unix( sprintf('%s','mkdir ',parent_dir,'_OVERWRITE') )
    %i=1;
end

% %Create a new directory if current one exists
% i=1;
% while ios==7
%     parent_dir=strcat(compound(ic,:),'_set',num2str(i));
%     ios=exist(parent_dir,'dir');
%     if ios ~=7
%         disp(sprintf('%s','Directory ',parent_dir,' exists.'));
%         unix( sprintf('%s','mkdir ',parent_dir) );
%         break
%     end
%     i=i+1;
% end
 
%Super useful for general formatting
%http://stackoverflow.com/questions/8805685/avoid-typing-conversion-specifier-for-every-column-in-large-table-in-textscan

%Output user settings
fname=strcat(parent_dir,'/parameters.dat');
fid=fopen(fname,'w+');
formatSpec = ['%s\t %d\t %.3f\t %d\n ' repmat('%f\t ', [1 19]) '%f\n'];
fprintf(fid,formatSpec,'Fitting Settings: nkpt, Eg_ideal, Nsym, b_weight ',...
            nkpt,Eg_ideal,Nsym,b_weight);  
fclose(fid);


    
%----------------------------
%Do fitting runs
%----------------------------
%Initialise timer for each iteration/fit
t = zeros(1,Nr);

%Do Nr fitting runs. 
for i=1:Nr
   %Time the routine per iteration
   tic;
    
   %Call min.b 
   [fval,m_eff,dH,dL,Eg_array,deltaSO_array,RMS,mxberr]=minb(i);
   
   %Fit value, m*, diff in HOMO and LUMO, Eg and delta_SO best-fit, 
   %RMS error and max error per band
   %band, per run
   fq(i,1)=fval;
   fq(i,2:3)=Eg_array(:);
   fq(i,4:5)=deltaSO_array(:); 
   fq(i,6:9)=m_eff(:);
   fq(i,10)=dH;
   fq(i,11)=dL;
   fq(i,12)=RMS(8);
   fq(i,13)=RMS(9);
   fq(i,14)=mxberr(8);
   fq(i,15)=mxberr(9);
     
   %Directory destination. 
   dir=strcat(parent_dir,'/',output_dir,num2str(i));
   unix(sprintf('%s','mkdir ',dir));
   %Move output files to 'dir'
   unix(sprintf('%s','mv compare.p ',dir));
   if i==1
       unix(sprintf('%s','mv CdSe_eigenfull_ideal.dat ',parent_dir));
   end
   unix(sprintf('%s','mv CdSe_eigenfull_fit.dat ',dir));
   unix(sprintf('%s','mv ',parent_dir,'/parameters.dat ',dir));
   
   %Save time for fit 'i'
   t(i)=toc;
    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Sort fq array w.r.t 1st column (fit value) only
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Produce dummy array B and index array I, which contains the indices
%of the original element positions 
[B,I]=sort(fq,1);

%Open file to write fq to
fname=strcat(parent_dir,'/results.dat');
fid=fopen(fname,'w+');
fprintf(fid,'%s\n','Ouput Index, Fit value, Eg_f (eV), Eg err (%), SO_f (eV),',...
     'SO err (%), me*, mhh*, mlh*, Diff: HOMO, LUMO. RMS: VBT, CBB. Max err: VBT, CBB');


%Write fq values to file
for i = 1:Nr
    %Mapped index, j
    j=I(i,1);
    %Print line to file
    fprintf(fid,'%d\t %f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\n',...
                 j,fq(j,:));
end
  
%Close fq values file
fclose(fid);


%Return total time for script
t_total=sum(t)

%Force Matlab to exit 
%error('Force Matlab to quit')
