function[Natoms,basis,spin,N,Np,Nc,db,al,compound,dv]=get_mat_param(ic)
    %------------------------------
    %Material Parameters
    %------------------------------
    %Atoms in primitive unit cell
    Natoms=2;
    %Orbital basis
    basis=5;
    %Spin states (up and down)
    spin=1;
    %Hamiltonian Dimensions
    N=Natoms*basis*spin;
    %Number of parameters (sp3s* basis)
    Np=13;
    %Bond lengths (Ang)
    db = [1.54,2.35,2.45,2.81,1.88,2.36,2.45,2.66,2.36, ...
          2.45,2.64,2.54,2.62,2.81,2.45,2.64,2.53,2.35,2.62];
    %Lattice constant (Ang)
    al= (4./sqrt(3))*db(ic);
    %Compound names
    compound = ['C   ';'Si  ';'Ge  ';'Sn  ';'SiC ';'AlP ';'AlAs';
                'AlSb';'GaP ';'GaAs';'GaSb';'InP ';'InAs';'InSb';
                'ZnSe';'ZnTe';'CdS ';'ZnS ';'CdSe'];
    %Number of compounds
    Nc=size(compound,1);
    %NN vectors w.r.t anion, in reduced units
    dv= [ 1.,  1.,  1. ;
          1.  -1., -1. ;
         -1.,  1,  -1. ;
         -1., -1,   1. ];
     
end 