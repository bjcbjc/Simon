

figure(1);
set(1, 'position', [200, 40, 1000, 650]);
set(1, 'paperpositionmode', 'auto');
set(gcf, 'renderer', 'painters');
%%
for chrmidx = 1:24
    clf
    varidx = curpatient.loc(:,1) == chrmidx;
    subplot(4,1,1); 
    %plot(curpatient.loc(varidx,2), log10(totalReadDNA(varidx,1)+1), '.')
    scatter(curpatient.loc(varidx,2), log10(curpatient.totalReadDNA(varidx,1)+1), 20, curpatient.variantAttr_mtx(varidx,7),'filled');
    ylabel('log10 reads', 'fontsize', 12);
    title(sprintf('DNA, normal, chrm %d',chrmidx), 'fontsize', 12);
    
    subplot(4,1,2); 
    %plot(curpatient.loc(varidx,2), log10(totalReadDNA(varidx,2)+1), '.')
    scatter(curpatient.loc(varidx,2), log10(curpatient.totalReadDNA(varidx,2)+1), 20, curpatient.variantAttr_mtx(varidx,7), 'filled');
    title('DNA, tumor', 'fontsize', 12);
    
    subplot(4,1,3); 
    %plot(curpatient.loc(varidx,2), log10(totalReadRNA(varidx,1)+1), '.')
    scatter(curpatient.loc(varidx,2), log10(curpatient.totalReadRNA(varidx,1)+1), 20, curpatient.variantAttr_mtx(varidx,7), 'filled');
    title('RNA, normal', 'fontsize', 12);
    
    subplot(4,1,4); 
    %plot(curpatient.loc(varidx,2), log10(totalReadRNA(varidx,2)+1), '.')
    scatter(curpatient.loc(varidx,2), log10(curpatient.totalReadRNA(varidx,2)+1), 20, curpatient.variantAttr_mtx(varidx,7), 'filled');
    title('RNA, tumor', 'fontsize', 12);
    
    
    saveas(1, sprintf('figures/chrmcoverage/chrm%02d, %s.png', chrmidx, curpatient.patient), 'png');
end

%%
for chrmidx = 1:24
    clf
    varidx = curpatient.loc(:,1)== chrmidx;
    cen = curpatient.distToCen(varidx) < 1e6;
    subplot(1,2,1)
    gscatter(log10(curpatient.totalReadDNA(varidx,1)+1), log10(curpatient.totalReadDNA(varidx,2)+1), cen);
    subplot(1,2,2)
    gscatter(log10(curpatient.totalReadRNA(varidx,1)+1), log10(curpatient.totalReadRNA(varidx,2)+1), cen);
    saveas(gcf, sprintf('figures/chrmcoverage/scatter chrm%02d, %s.png', chrmidx, curpatient.patient), 'png');
end
%%
for chrmidx = 1:24
    clf
    varidx = curpatient.loc(:,1)== chrmidx;    
    subplot(1,2,1)
    scatter(log10(curpatient.totalReadDNA(varidx,1)+1), log10(curpatient.totalReadDNA(varidx,2)+1), 30, curpatient.variantAttr_mtx(varidx,7));
    colorbar
    subplot(1,2,2)
    scatter(log10(curpatient.totalReadRNA(varidx,1)+1), log10(curpatient.totalReadRNA(varidx,2)+1),30, curpatient.variantAttr_mtx(varidx,7));
    colorbar
    saveas(gcf, sprintf('figures/chrmcoverage/MQ chrm%02d, %s.png', chrmidx, curpatient.patient), 'png');
end