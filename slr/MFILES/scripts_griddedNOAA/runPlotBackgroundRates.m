% Last updated by Robert Kopp, robert-dot-kopp-at-rutgers-dot-edu, Wed Aug 17 20:47:52 EDT 2016


% plot GIA rates

clf;
scatter(mod(targsitecoords(:,2),360),targsitecoords(:,1),10,rateprojs,'filled');
hold on;
plotcont([],[],1);
colorbar;
colormap(jet);
caxis([-3 3]);
xlim([150 320]); ylim([-60 73]);
title('Median background rate (mm/yr)');
pdfwrite('bkgdrate_mean');
xlim([220 300]); ylim([0 60]);
pdfwrite('bkgdrate_mean2');

clf;
scatter(mod(targsitecoords(:,2),360),targsitecoords(:,1),10,rateprojssd,'filled');
hold on;
plotcont([],[],1);
colorbar;
colormap(jet);
caxis([0 2]);
xlim([150 320]); ylim([-60 73]);
title('Median background rate std. (mm/yr)');
pdfwrite('bkgdrate_std');
xlim([220 300]); ylim([0 60]);
pdfwrite('bkgdrate_std2');

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
