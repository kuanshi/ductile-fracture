tartag = [5,6,7,8,188,189,190,191,192,193,194,195,215,216,217,218,219];
for i = 1:1:length(tartag)
%     ex = importdata(['./SpiralRaw/' num2str(i) '.txt']);
    ex = importdata(['./RecRaw/' num2str(tartag(i)) '.txt']);
    d_c = ex(:,1)*0.0393701;
    v_c = ex(:,2)*0.2248089;
    f = fopen(['./' num2str(tartag(i)) '.txt'],'wt');
    for j = 1:1:length(d_c)
        fprintf(f,'%6.4f\t',d_c(j));
        fprintf(f,'%6.4f\n',v_c(j));
    end
    fclose(f);
    figure;
    plot(d_c,v_c);
    clear d_c v_c;
end