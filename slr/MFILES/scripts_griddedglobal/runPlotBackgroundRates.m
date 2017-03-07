% Last updated by Robert Kopp, robert-dot-kopp-at-rutgers-dot-edu, Fri Jan 13 09:32:33 EST 2017


% plot GIA rates

clf;
scatter(mod(targsitecoords(:,2),360),targsitecoords(:,1),10,rateprojs,'filled');
hold on;
plotcont([],[],1);
colorbar;
colormap(jet);
caxis([-3 3]);
title('Median background rate (mm/yr)');
pdfwrite('bkgdrate_mean');

clf;
scatter(mod(targsitecoords(:,2),360),targsitecoords(:,1),10,rateprojssd,'filled');
hold on;
plotcont([],[],1);
colorbar;
colormap(jet);
caxis([0 2]);
title('Median background rate std. (mm/yr)');
pdfwrite('bkgdrate_std');

fid=fopen('bkgdrate.tsv','w');
fprintf(fid,'site\tID\tlat\tlong\tbkgd rate\tstd\n');
for iii=1:length(rateprojs)
    fprintf(fid,targregionnames{iii});
    fprintf(fid,'\t%0.0f',targregions(iii));
    fprintf(fid,'\t%0.2f',targsitecoords(iii,:));
    fprintf(fid,'\t%0.2f',[rateprojs(iii) rateprojssd(iii)]);
    fprintf(fid,'\n');
    
end

fclose(fid);
