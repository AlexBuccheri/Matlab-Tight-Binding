%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Subroutine to construct Hamiltonian
%Input  = k-vector 
%Output = Hamiltonian
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [H]=ham(kvec,parameter)

%Declare global variables
global ic N spin dv  


%%%%%%%%%%%%%%%%%%%%%
%Assign parameters
%%%%%%%%%%%%%%%%%%%%%

%On-site
E_sa=parameter(ic,1);
E_pa=parameter(ic,2);
E_sc=parameter(ic,3);
E_pc=parameter(ic,4);
E_ssta=parameter(ic,5);
E_sstc=parameter(ic,6);

%Inter-atomic
Vss=parameter(ic,7);
Vxx=parameter(ic,8);
Vxy=parameter(ic,9);
Vsa_pc=parameter(ic,10);   
Vsc_pa=parameter(ic,11);   %Relation to Vpa_sc?
Vssta_pc=parameter(ic,12);
Vpa_sstc=parameter(ic,13);

%Vpa_sc
Vpa_sc=1.*Vsc_pa;

%Spin-orbit parameters (in eV) for CdS and CdSe, respectively
if ic == 17
    %Anion
    lma=0.03;
    %Cation
    lmc=0.0;
elseif ic ==19 
    lma=0.16;
    lmc=0.0755;
else
    error('No SO parameters for this material. Stop')
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Set up Hamiltonian
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Define Tight-binding parameters and periodic
%g-functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Initialise H
H=zeros(N);



%On-site, atomic energies

%Spins (up,up)
%Anion
H(1,1)=E_sa;
H(2,2)=E_pa;
H(3,3)=E_pa;
H(4,4)=E_pa;
H(5,5)=E_ssta;

%Cation
H(6,6)=E_sc;
H(7,7)=E_pc;
H(8,8)=E_pc;
H(9,9)=E_pc;
H(10,10)=E_sstc;

if spin==2 
   %Spins (down,down)
   %Anion
   H(11,11)=E_sa;
   H(12,12)=E_pa;
   H(13,13)=E_pa;
   H(14,14)=E_pa;
   H(15,15)=E_ssta;

   %Cation
   H(16,16)=E_sc;
   H(17,17)=E_pc;
   H(18,18)=E_pc;
   H(19,19)=E_pc;
   H(20,20)=E_sstc;
end

%Inter-atomic ia=anion, ic=cation using exponentials
H_ac=[    Vss*e0(kvec),  Vsa_pc*e1(kvec),   Vsa_pc*e2(kvec),  Vsa_pc*e3(kvec),   0.,       ; 
      -Vpa_sc*e1(kvec),     Vxx*e0(kvec),      Vxy*e3(kvec),     Vxy*e2(kvec),  -Vpa_sstc*e1(kvec), ; 
      -Vpa_sc*e2(kvec),     Vxy*e3(kvec),      Vxx*e0(kvec),     Vxy*e1(kvec),  -Vpa_sstc*e2(kvec), ; 
     -Vpa_sc*e3(kvec),     Vxy*e2(kvec),      Vxy*e1(kvec),     Vxx*e0(kvec),   -Vpa_sstc*e3(kvec), ; 
       0.,             Vssta_pc*e1(kvec), Vssta_pc*e2(kvec),Vssta_pc*e3(kvec),   0.          ];
   

%---------------------------
%Sub-matrix (up,up)
%---------------------------
%Put H_ac in H
H(1:5,6:10)=H_ac(:,:);

%Put H_ca in H, where H_ca= transpose(conjg(H_ac))
H(6:10,1:5)=H_ac';


%---------------------------
%Sub-matrix (down,down)
%---------------------------
if spin==2
   %Put H_ac in H
   H(11:15,16:20)=H_ac(:,:);

   %Put H_ca in H, where H_ca= transpose(conjg(H_ac))
   H(16:20,11:15)=H_ac';
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Spin-Orbit Hamiltonians
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Note, (down,down) and (down,up) generated from 
%symmetry relations

%SO Hamiltonian (up,up) or (1,1) -anion
H_aa_11=[ 0.,  0.,     0.,    0.,  0.,  ;
          0.,  0.,  -1i*lma,  0.,  0.,  ;
          0., 1i*lma,  0.,    0.,  0.,  ;  
          0.,  0.,     0.,    0.,  0.,  ;
          0.,  0.,     0.,    0.,  0.   ];
   
      
%SO Hamiltonian (up,up) or (1,1) -cation      
H_cc_11=[ 0.,  0.,     0.,    0.,  0.,  ;
          0.,  0.,  -1i*lmc,  0.,  0.,  ;
          0., 1i*lmc,  0.,    0.,  0.,  ;  
          0.,  0.,     0.,    0.,  0.,  ;
          0.,  0.,     0.,    0.,  0.   ];
      


%SO Hamiltonian (up,down) or (1,2) -anion
H_aa_12=[ 0.,  0.,     0.,    0.,     0.,  ;
          0.,  0.,     0.,    lma,    0.,  ;
          0.,  0.,     0.,  -1i*lma,  0.,  ;
          0.,  -lma, 1i*lma,   0.,    0.,  ;
          0.,  0.,     0.,     0.,    0.   ];



%SO Hamiltonian (up,down) or (1,2) -cation
H_cc_12=[ 0.,  0.,     0.,    0.,     0.,  ;
          0.,  0.,     0.,    lmc,    0.,  ;
          0.,  0.,     0.,  -1i*lmc,  0.,  ;
          0.,  -lmc, 1i*lmc,   0.,    0.,  ; 
          0.,  0.,     0.,     0.,    0.  ];

  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Add H_SO to full Hamiltonian
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if spin==2

   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
   %On-site Contributions
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %Spins (up,up)
   %Anion
   H(1:5,1:5)=H(1:5,1:5)+H_aa_11(:,:);
   %Cation
   H(6:10,6:10)=H(6:10,6:10)+H_cc_11(:,:);

   %Spins (down,down)
   %Anion
   H(11:15,11:15)=H(11:15,11:15)+conj(H_aa_11(:,:));
   %Cation
   H(16:20,16:20)=H(16:20,16:20)+conj(H_cc_11(:,:));

   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %Spin-mixing terms
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %Spins (up,down)
   %Anion
   H(1:5,11:15)=H_aa_12(:,:);
   %Cation
   H(6:10,16:20)=H_cc_12(:,:);

   %Spins (down,up)
   %Anion
   H(11:15,1:5)=H_aa_12';
   %Cation
   H(16:20,6:10)=H_cc_12';

end


%End of Hamiltonian subroutine
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%
%NN TB g-functions
%%%%%%%%%%%%%%%%%%%%%%%%

function g=g0(kvec)
g=cos((kvec(1)*pi)/2)*cos((kvec(2)*pi)/2.)*cos((kvec(3)*pi)/2.)...
  -1i*sin((kvec(1)*pi)/2.)*sin((kvec(2)*pi)/2.)*sin((kvec(3)*pi)/2.);
end

function g=g1(kvec)
g=-cos((kvec(1)*pi)/2)*sin((kvec(2)*pi)/2.)*sin((kvec(3)*pi)/2.)...
  +1i*sin((kvec(1)*pi)/2.)*cos((kvec(2)*pi)/2.)*cos((kvec(3)*pi)/2.);
end

function g=g2(kvec)
g=-sin((kvec(1)*pi)/2)*cos((kvec(2)*pi)/2.)*sin((kvec(3)*pi)/2.)...
  +1i*cos((kvec(1)*pi)/2.)*sin((kvec(2)*pi)/2.)*cos((kvec(3)*pi)/2.);
end

function g=g3(kvec)
g=-sin((kvec(1)*pi)/2)*sin((kvec(2)*pi)/2.)*cos((kvec(3)*pi)/2.)...
  +1i*cos((kvec(1)*pi)/2.)*cos((kvec(2)*pi)/2.)*sin((kvec(3)*pi)/2.);
end


%Exponential definitions
function e=e0(kvec)
global dv
e = 0.25*(exp(1i*0.5*pi*dot(kvec,dv(1,:)))+exp(1i*0.5*pi*dot(kvec,dv(2,:)))+...
          exp(1i*0.5*pi*dot(kvec,dv(3,:)))+exp(1i*0.5*pi*dot(kvec,dv(4,:))) );
end 

function e=e1(kvec)
global dv
e = 0.25*(exp(1i*0.5*pi*dot(kvec,dv(1,:)))+exp(1i*0.5*pi*dot(kvec,dv(2,:)))-...
          exp(1i*0.5*pi*dot(kvec,dv(3,:)))-exp(1i*0.5*pi*dot(kvec,dv(4,:))) );
end 

function e=e2(kvec)
global dv
e = 0.25*(exp(1i*0.5*pi*dot(kvec,dv(1,:)))-exp(1i*0.5*pi*dot(kvec,dv(2,:)))+...
          exp(1i*0.5*pi*dot(kvec,dv(3,:)))-exp(1i*0.5*pi*dot(kvec,dv(4,:))) );
end 

function e=e3(kvec)
global dv
e = 0.25*(exp(1i*0.5*pi*dot(kvec,dv(1,:)))-exp(1i*0.5*pi*dot(kvec,dv(2,:)))-...
          exp(1i*0.5*pi*dot(kvec,dv(3,:)))+exp(1i*0.5*pi*dot(kvec,dv(4,:))) );
end 