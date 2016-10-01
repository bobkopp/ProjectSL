% Last updated by Robert Kopp, robert-dot-kopp-at-rutgers-dot-edu, Thu Sep 29 23:18:02 EDT 2016

% plot Oceanographic Changes

subt=find(OceanDynYears==2100);
subscen=1;

clf;
scatter(mod(targsitecoords(:,2),360),targsitecoords(:,1),10,OceanDynMean(subt,:,subscen),'filled');
hold on;
plotcont([],[],1);
colorbar;
colormap(jet);
caxis([-250 750]);
xlim([150 320]); ylim([-60 73]);
title(['Oceanographic mean, ' scens{subscen} ', 2100 (mm)']);
pdfwrite('Oceanographic_mean');
xlim([220 300]); ylim([0 60]);
pdfwrite('Oceanographic_mean2');

clf;
scatter(mod(targsitecoords(:,2),360),targsitecoords(:,1),10,OceanDynStd(subt,:,subscen),'filled');
hold on;
plotcont([],[],1);
colorbar;
colormap(jet);
crange=quantile(OceanDynStd(subt,:,subscen),[.01 .99]);
caxis(crange);
%caxis([0 2]);
xlim([150 320]); ylim([-60 73]);
title(['Ocean dynamic std, ' scens{subscen} ', 2100 (mm)']);
pdfwrite('Oceanographic_std');
xlim([220 300]); ylim([0 60]);
pdfwrite('Oceanographic_std2');

fid=fopen('Oceanographic.tsv','w');
fprintf(fid,['site\tID\tlat\tlong\t2100 ' scens{subscen} ' Ocean Dyn Mean\tstd\n']);
for iii=1:length(targregions)
    fprintf(fid,targregionnames{iii});
    fprintf(fid,'\t%0.0f',targregions(iii));
    fprintf(fid,'\t%0.2f',targsitecoords(iii,:));
    fprintf(fid,'\t%0.2f',[OceanDynMean(subt,iii,subscen), OceanDynStd(subt,iii,subscen),]);
    fprintf(fid,'\n');
    
end

fclose(fid);
