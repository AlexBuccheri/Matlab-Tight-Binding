
%Generate GNUPlot script to visualise ideal and best-fit eigenvalues
function[] = gnu_script(fname)

    global N

    %Open file 
    fid=fopen(fname,'w+');

    %Plot formatting 
    fprintf(fid,'%s\n','set term x11');
    fprintf(fid,'%s\n','set key off');
    fprintf(fid,'%s\n','set ylab''Energy (eV)''');
    fprintf(fid,'%s\n','set linestyle 1 pt 2 ps 1 lw 2 lc rgb''red''');
    fprintf(fid,'%s\n','set linestyle 2 lt 1 lw 1 lc rgb''blue''');

    %Ideal eigenvalue bands as red crosses
    str_plot='plot ''../CdSe_eigenfull_ideal.dat'' u 1:2 w p';
    fprintf(fid,'%s\n',str_plot);

    for i=3:N+1
       str_plot=sprintf('%s','replot ''../CdSe_eigenfull_ideal.dat'' u 1:',num2str(i),' w p ls 1');
         fprintf(fid,'%s\n',str_plot);
    end

    %Best-fit eigenvalue bands as blue lines
    for i=2:N+1
        str_plot=sprintf('%s','replot ''CdSe_eigenfull_fit.dat'' u 1:',num2str(i),' w l ls 2');
        fprintf(fid,'%s\n',str_plot);
    end
 
    %Shut file
    fclose(fid);

end