% Last updated by Robert Kopp, robert-dot-kopp-at-rutgers-dot-edu, Fri Jan 13 09:33:37 EST 2017

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
title(['Oceanographic mean, ' scens{subscen} ', 2100 (mm)']);
pdfwrite('Oceanographic_mean');

clf;
scatter(mod(targsitecoords(:,2),360),targsitecoords(:,1),10,OceanDynStd(subt,:,subscen),'filled');
hold on;
plotcont([],[],1);
colorbar;
colormap(jet);
crange=quantile(OceanDynStd(subt,:,subscen),[.01 .99]);
caxis(crange);
%caxis([0 2]);
title(['Ocean dynamic std, ' scens{subscen} ', 2100 (mm)']);
pdfwrite('Oceanographic_std');

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
