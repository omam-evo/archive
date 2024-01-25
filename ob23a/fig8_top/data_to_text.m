
f = fopen('data.txt','w');
fprintf(f,'N\tMU_RANGE\tTRIALS\tP_S\n')

for i=1:length(data)
    
    m = strjoin(cellstr(num2str(data(i).MU)),';')
    m = append('[',m,']');

    t = strjoin(cellstr(num2str(data(i).TR)),';')
    t = append('[',t,']');

    p = strjoin(cellstr(num2str(data(i).PS)),';')
    p = append('[',p,']');

    line = append(num2str(data(i).N),'\t',m,'\t',t,'\t',p,'\n');
    fprintf(f,line)


end

fclose(f);