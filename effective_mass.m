%-------------------------------------------------------------------------
%Compute effective masses at Gamma for heavy, light and SO holes plus
%the lowest electron state
%-------------------------------------------------------------------------

function[meff_el,meff_hh,meff_lh,meff_so]=effective_mass(b)
    global al parameter ic

    %Define parameters equal to input parameters 'b'
    parameter(ic,:)=b(:);

    %Constants
    pi=3.14159265359;
    m_0=9.1093821545e-31;    %kg
    hbar=6.5821192815e-16;   %eV.s
    %hbar=1.0545718e-34;      %J.s
    qe=1.60217656535e-19;    %C
    
    
    %Set up a very fine, square k-grid centred at Gamma 3x3x3
    nkx=3;
    dk=(2*pi/al)*0.0001;
    k0=-dk;
    mid=int8(0.5*(nkx^3));  %Apparently should be 14 according to below

    
    %Solve for the energies of these k-points 
    cnt=0;
    for i=1:nkx
        for j=1:nkx
            for k=1:nkx
                cnt=cnt+1;
                %k-points(no need to store)
                kpoint(1)=k0+(i-1)*dk;
                kpoint(2)=k0+(j-1)*dk;
                kpoint(3)=k0+(k-1)*dk;

                %Check index for Gamma (should be centre of grid)
                %Indices required when computing forward difference
                if kpoint(1)==0. && kpoint(2)==0. && kpoint(3)==0.
                   m(1)=i;
                   m(2)=j;
                   m(3)=k;
                   %cnt
                   %kpoint(:)
                end
            
                %Construct H and solve for eigenvalues, 
                H=ham(kpoint,parameter);
                eigen=eig(H);
                %Store relevant bands (all bands doubly degenerate)
                E_el(i,j,k)=eigen(10);  %9 & 10
                E_hh(i,j,k)=eigen(7);  %7 & 8
                E_lh(i,j,k)=eigen(5);  %5 & 6 
                E_so(i,j,k)=eigen(3);  %3 & 4
                
            
            end
        end
    end
    
    
    %Compute elements of effective mass tensor
    for i=1:3
        for j=1:3    
            %Hessian matrices
            H_el(i,j)=ScndDeriv(E_el,i,j,m,dk);
            H_hh(i,j)=ScndDeriv(E_hh,i,j,m,dk);
            H_lh(i,j)=ScndDeriv(E_lh,i,j,m,dk);
            H_so(i,j)=ScndDeriv(E_so,i,j,m,dk);
              
            %Effective masses
            m_el(i,j)=((hbar)^2./H_el(i,j))*(qe*(1.e20)/m_0);
            m_hh(i,j)=((hbar)^2./H_hh(i,j))*(qe*(1.e20)/m_0);
            m_lh(i,j)=((hbar)^2./H_lh(i,j))*(qe*(1.e20)/m_0);
            m_so(i,j)=((hbar)^2./H_so(i,j))*(qe*(1.e20)/m_0);
        end
    end
    
    %Note, Matlab truncates outputs to 5 d.p. =>
    %if one has a mixture of very large and very small values, the small
    %values get truncated to zero. Hence use the 'format short e' command
    %Hessian diagonals should be ~10^1. If much smaller => should be zero.
    %Non-zero, small values will cause off-diagonal m* to blow up. 
    %
    %disp('Constant:')
    %con=((hbar^2)*qe*(1.e20)/m_0)
    %format short e
    %m_el
    %m_hh
    %m_lh
    %m_so
    
    %Take trace and normalise to get isotropic effective masses
    norm=size(m_el,1);
    meff_el=(1./norm)*trace(m_el);
    meff_hh=(1./norm)*trace(m_hh);
    meff_lh=(1./norm)*trace(m_lh);
    meff_so=(1./norm)*trace(m_so);

    return
end