%-----------------------------------------------------------------------
%Function to compute the 2nd order partialderivative of any 3D function
%defined on a linear grid
%-----------------------------------------------------------------------
%Inputs
% f      function f(x,y,z)
% v      Variable 1 = x, y or z (represented as an integer)
% u      Variable 2 = x, y or z (represented as an integer)
% m(1:3) Point at which 2nd derivative is evaluated (integer)
% dv     Grid spacing (uniform, hence a variable)
%-----------------------------------------------------------------------

function SD = ScndDeriv(f,u,v,m,dv)

    %----------------------
    %fxx, fyy and fzz
    %----------------------
    if u==v        
        %Initialise indices
        i=zeros(1,2);
        j=zeros(1,2);
        k=zeros(1,2);
        %Establish whether derivative is fxx, fyy or fzz
        %and assign appropriate indices for central difference
        if u==1 
            i=[1,-1];
        elseif u==2
            j=[1,-1];
        elseif u==3
            k=[1,-1];
        end
        
        %General expression for Central difference numerator 
        num = f(m(1)+i(1),m(2)+j(1),m(3)+k(1)) -2*f(m(1),m(2),m(3)) + ...
              f(m(1)+i(2),m(2)+j(2),m(3)+k(2));

        %Check numerator and demoninators behave, then return value
        if num==0.
            SD = 0.;
        elseif dv==0
            error('Grid spacing must be finite for 2nd derivative to work')
        else
            SD = num/dv^2.;
        end
        
        return
    end  
          
    %---------------------------------
    %Mixed Derivatives, like fxy etc
    %---------------------------------
    if u ~=v
        
        i=zeros(1,4);
        j=zeros(1,4);
        k=zeros(1,4);
           
       if u==1
           i=[1,1,-1,-1];
       elseif u==2
           j=[1,1,-1,-1];
       elseif u==3
           k=[1,1,-1,-1];
       end
       
        if v==1
            i=[1,-1,1,-1];
        elseif v==2
            j=[1,-1,1,-1];
        elseif v==3
            k=[1,-1,1,-1];
        end

       %Central or forward difference expression (not sure which)
       num = f(m(1)+i(1),m(2)+j(1),m(3)+k(1)) - f(m(1)+i(2),m(2)+j(2),m(3)+k(2)) -  ...
             f(m(1)+i(3),m(2)+j(3),m(3)+k(3)) + f(m(1)+i(4),m(2)+j(4),m(3)+k(4)); 
       denom=4.*dv*dv;

       %Output mixed, 2nd order partial derivative 
       if num==0.
           SD = 0.;
       elseif dv==0
           error('Grid spacing must be finite for 2nd derivative to work')
       else
           SD = num/denom;
       end    
       
       return    
    end
          
          
end