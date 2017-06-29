%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%genf function is a wrapper for 'generate_eigenvalues' that also calculates
%and returns a fit value. 
%Inputs 
%  b             Trial parameters for compound ic 
%  kpoint        k-point grid
%  hk_index      Index indentifying high-symmetry points in kpoint
%  eigen_ideal   Target/ideal eigenvalues, used in fit function definition
%Output
%  f             Fit value, to be minimised by varying b
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
function f=genf(b,kpoint,hk_index,eigen_ideal)  

    global ic parameter k_weight b_weight N

    %------------------------------------------------
    %Generate trial eigenvalues with trial parameters
    %------------------------------------------------
    %Set correct row in parameter array passed to generate eigenvalues
    %equal to trial parameters b for material ic
    parameter(ic,:)=b(:);
    
    %Generate trial eigenvalues with trial parameters 
    p_type='trial';
    eigen_trial = generate_eigenvalues(kpoint,hk_index,parameter,p_type);
    
    
    %-----------------------
    %Define the fit function
    %-----------------------
    %Could compute differences in effective masses
     
    %Band weights given in itr_fit, whereas k weights given in minb
    %Define the fit function value 
    f=0.;
    for ik=1:size(kpoint,1)
        for ib=1:N
            f= f+ k_weight(ik)*b_weight(ib)*(eigen_ideal(ik,ib)-eigen_trial(ik,ib))^2.; 
        end
    end
            
    return          
end 
  
 
