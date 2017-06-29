%Various approaches to generating parameter seeds 


%1. Choose all parameters at random within a given energy range.

    %Initial guess at parameter. Use random numbers
    %within the range [a,b]
    %Set random, positive values
    a=0;
    b=10;
    b0=a + (b-a).*rand(Np,1);

    %Set random, negative values
    a=-10;
    b=0;
    b0(1)=a + (b-a).*rand(1,1);
    b0(7)=a + (b-a).*rand(1,1);
    
    
%2. Fit inter-atomic only
    Requires Np=7 in minb.m

    %Initial guess at inter-atomic parameter. Use random numbers
    %within the range [a,b]
    %Set random, positive values
    a=0;
    b=10;
    b0=a + (b-a).*rand(Np,1);

    %Set random, negative value (Vss)
    a=-10;
    b=0;
    b0(1)=a + (b-a).*rand(1,1);
    
    
    
%2.5 Having fit inter-atomic parameters, keeping on-site fixed,
     %take orthogonal onsite and best-fit interatomic and allow
     %them all to vary
    
    Nr needs to be fixed within minb to the number of inter-atomic
    sets that will be used
     

    %Nr defined by user here and in itr_fit
    Nr=4
    
    %Fix initial parameters guess using orthogonal CdSe parameters 
    b0= [-9.63   1.47    0.03   4.73   7.53   5.72...
         -4.64   2.64    5.36   4.57   5.54   3.05   2.49]
    
    %Replace inter-atomic parameters with those read from file
    fid = fopen('CdSe_set10_seed.dat','r');
    %Skip header
    fgetl(fid);
    %Inter-atomic parameters
    p0=fscanf(fid,'%g',[7,Nr]);
    fclose(fid);
    p0=p0';
    b0(7:13)=p0(i,1:7)
 



     
%2.6 Could also use one set of best-fit inter-atomic parameters
     %and use random ons-ite seeds, for multiple runs, as a 3rd layer on the above. 
    
    %Initial guess at parameter. Use random numbers
    %within the range [a,b]
    %Set random, positive values
    a=0;
    b=10;
    b0=a + (b-a).*rand(Np,1);

    %Set random, negative values
    a=-10;
    b=0;
    b0(1)=a + (b-a).*rand(1,1);
    b0(7)=a + (b-a).*rand(1,1);
    
    
    %Replace inter-atomic parameters with those read from file
    fid = fopen('nor8_bestp.dat','r');
    fgetl(fid);
    p0=fscanf(fid,'%g',[7,Nr]);
    fclose(fid);
    p0=p0';
    b0(7:13)=p0(i,1:7)

    
     
%3. Use orthogonal parameters

    %CdSe
    b0= [-9.63   1.47    0.03   4.73   7.53   5.72...
     -4.64   2.64    5.36   4.57   5.54   3.05   2.49]
    
    %CdS
    b0=[-11.5300  0.5300  1.8300  5.8700  7.1300...  
        6.8700   -3.0700  1.7600  4.2300  2.1700  5.4800...   
        1.9900    3.0600]
    


    %Best-fit CdS, trying to refine VB
 
     b0 = [-11.793437 	0.625543 	6.248352 	5.682181 	14.287245...
           2.510075 	1.981246 	1.863334 	4.320821 	0.967844...
           2.918581 	-2.970730 	5.693129]; 
       

    
%4.  Keeping as reference
    %Take the best 4 non-orthogonal initial guesses from file
    %fid = fopen('nor8_bestp.dat','r');
    %fgetl(fid);
    %p0=fscanf(fid,'%g',[Np,Nr]);
    %fclose(fid);
    %p0=p0';
    %b0=p0(i,:)
    
