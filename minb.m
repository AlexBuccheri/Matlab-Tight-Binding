%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Matlab function to minimise b
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function[fval,meff_array,dH,dL,Eg_array,deltaSO_array,band_RMS,mxbnd_err] = minb(itr)
    global ic nkpt N spin dv Natoms basis compound db al Nsym k_weight parent_dir

    %Clear working directory
    %unix('rm data_fits/*');

    %Get required material parameters
    [Natoms,basis,spin,N,Np,Nc,db,al,compound,dv]=get_mat_param(ic);
      
    %-----------------------------------
    %Generate k-point grid  
    %-----------------------------------
    %Starting High symmetry points
    HS1= [0.5,  0.5,  0.5;
          0.,   0.,   0. ;
          0.,   1.,   0. ;
          0.75, 0.75, 0. ];
           
    %Ending High symmetry points
    HS2= [0.,   0.,  0.  ;
          0.,   1.,  0.  ;
          0.25, 1.,  0.25;
          0.,   0.,  0.  ];
      
    [kpoint,hk_index]=gen_kgrid(nkpt,HS1,HS2);
    
    %----------------------------------------------------------------
    %k-point weights for use in fitting. High Symmetry points have 
    %weightings MSym times larger than all other k-points
    %----------------------------------------------------------------
    k_weight=zeros(size(kpoint,1));
    k_weight(:)=0.25;
    k_weight(hk_index(1))=Nsym*k_weight(hk_index(1));
    k_weight(hk_index(2))=Nsym*k_weight(hk_index(2));
    k_weight(hk_index(3))=Nsym*k_weight(hk_index(3));
    k_weight(hk_index(4))=Nsym*k_weight(hk_index(4));
    %where hk_index is defined w.r.t. HS1
    
    
   
    %-----------------------------------------------------------
    %Generate ideal eigenvalues, Eg, effective masses and Delta_SO
    %using orthogonal parameters. Correct band-gap is obtained by
    %a rigid shift of the eigenvalues (done in function call)
    %
    %Note, could be replaced by GW band structure or LDA with 
    %scissor-shift, for example. 
    %-----------------------------------------------------------
    
    %Read in orthogonal TB parameters
    fid=fopen('/Users/alexanderbuccheri/Documents/MATLAB/SO/parameters.dat','r');
    parameter=fscanf(fid,'%g',[Np,Nc]);
    fclose(fid);
    parameter=parameter';
    
    %Generate ideal eigenvalues, Eg and delta_SO
    p_type='ideal';
    eigen_ideal = generate_eigenvalues(kpoint,hk_index,parameter,p_type);
    %Output ideal band structure (once)
    if itr==1
        output_bands(HS1,HS2,eigen_ideal,'ideal');
    end
    
    %Ideal band-gap and bulk VB splitting
    VBT=4*spin; 
    igam=hk_index(2);
    Eg_ideal=eigen_ideal(igam,VBT+1)-eigen_ideal(igam,VBT);
    if spin == 1 
        deltaSO_ideal=0.;
    elseif spin == 2
        deltaSO_ideal=eigen_ideal(igam,5)-eigen_ideal(igam,3);
    end
      
    %Generate effective masses (these won't necessarily agree with 
    %experiment but can't do much better)
    %[m_e_ideal,m_hh_ideal,m_lh_ideal,m_so_ideal]=effective_mass(parameter(ic,:));
    
    
 
    %-----------------------------------------------------------
    %Initial guess at TB parameters 
    %-----------------------------------------------------------
    %Fix initial parameters guess using orthogonal CdSe parameters 
    b0= parameter(ic,:);
     
    
    %Create seeds which can vary up to ±100% of the orthogonal No SO values
%     b0= [-9.63, 1.47, 0.03, 4.73, 7.53, 5.72, -4.64, 2.64, 5.36, 4.57,...
%           5.54, 3.05, 2.49];  
%     for j=1:Np
%         %random number in range [-1:1]
%         r = -2.*rand(1,1) + 1;
%         Ar= r*b0(j);
%         b0(j)=b0(j)+Ar;
%     end
    
    %Output initial seed parameters
    fname=strcat(parent_dir,'/parameters.dat');
    fid=fopen(fname,'a+');
    formatSpec = ['%s\t ' repmat('%f\t ', [1 Np-1]) '%f\n'];
    fprintf(fid,formatSpec,'Seed Parameters: ',b0);  
    fclose(fid);
    
    %-----------------------------------------------------------   
    %Run minimisation of fval   
    %CONSIDER OTHER MINIMISATION FUNCTIONS INCLUDING ANNEALING
    %-----------------------------------------------------------
    %options = optimset('MaxFunEvals',60800)
    %@(b0) => b0 is the parameter to minimised. All following
    %function arguments do not vary with the minimisation procedure
    [b,fval]=fminsearch(@(b0)genf(b0,kpoint,hk_index,eigen_ideal),b0)
    
    %Output final best-fit parameters
    fname=strcat(parent_dir,'/parameters.dat');
    fid=fopen(fname,'a+');
    formatSpec = ['%s\t ' repmat('%f\t ', [1 Np-1]) '%f\n'];
    fprintf(fid,formatSpec,'Best-fit Parameters: ',b);  
    fclose(fid);
    
   
    
    %-----------------------------------------------------------
    %Assess the quality of the fit
    %-----------------------------------------------------------
    %Use best-fit parameters 'b' to compute final, best-fit eigenvalues.
    parameter(ic,:)=b(:); 
    p_type='fit';
    eigen_fit = generate_eigenvalues(kpoint,hk_index,parameter,p_type);
    
    %Output best-fit band structure
    output_bands(HS1,HS2,eigen_fit,'fit');
    
    %Calculate differences in HOMO and LUMO
    igam=hk_index(2);
    %where index 1= k-point (Gamma) and index 2 = band (highest VB or lowest CB)
    dH=eigen_ideal(igam,VBT)-eigen_fit(igam,VBT);
    dL=eigen_ideal(igam,VBT+1)-eigen_fit(igam,VBT+1);
    
    %Best-fit band-gap and bulk VB splitting
    Eg_fit=eigen_fit(igam,VBT+1)-eigen_fit(igam,VBT);
    Eg_array(1)=Eg_fit;
    Eg_array(2)=(abs(Eg_fit-Eg_ideal)/Eg_ideal)*100.;
    
    if spin==1
        deltaSO_fit=0.;
        deltaSO_array(1:2)=0.
    elseif spin ==2
        deltaSO_fit=eigen_fit(igam,5)-eigen_fit(igam,3);
        deltaSO_array(1)=deltaSO_fit;
        deltaSO_array(2)=(abs(deltaSO_fit-deltaSO_ideal)/deltaSO_ideal)*100.;
    end
    
    %Calculate effective masses and differences in them |m*ideal-m*fit|
    [m_e_fit,m_hh_fit,m_lh_fit,m_so_fit]=effective_mass(parameter(ic,:));
    meff_array(1)=m_e_fit;
    meff_array(2)=m_hh_fit;
    meff_array(3)=m_lh_fit;
    meff_array(4)=m_so_fit;
    
    %RMS and max error of each band (in eV)
    Nk=size(kpoint,1);
    for ib=1:N
        %Initialise max error per band and sum of the differences
        mxbnd_err(ib)=abs(eigen_ideal(1,ib)-eigen_fit(1,ib));
        sum_sq=0.;
        
        for ik=1:Nk
            %Evaluate energy diff per k-point
            if abs(eigen_ideal(ik,ib)-eigen_fit(ik,ib))>mxbnd_err(ib)
                mxbnd_err(ib)=abs(eigen_ideal(ik,ib)-eigen_fit(ik,ib));
            end
            %Sum differences squared
            sum_sq=sum_sq+(eigen_ideal(ik,ib)-eigen_fit(ik,ib))^2.;
        end
        band_RMS(ib)=sqrt(sum_sq/double(Nk));
    end
    
    
    %Generate GNU script to compare ideal and best-fit band structures
    fname='compare.p';
    gnu_script(fname);
    
    %Close any file handles that are open (shouldn't be an issue with this
    %script)
    fclose('all');

end

%END OF FUNCTION






