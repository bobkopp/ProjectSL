% Last updated by Robert Kopp, robert-dot-kopp-at-rutgers-dot-edu, Fri Sep 30 08:19:12 EDT 2016


% plot fingerprints
for iii=1:length(fplab)
    clf;
    scatter(mod(targsitecoords(:,2),360),targsitecoords(:,1),10,fpsite(:,iii),'filled');
    hold on;
    plotcont([],[],1);
    colorbar;
    colormap(jet);
    caxis([0 1.3]);
    xlim([150 320]); ylim([-60 73]);
    title(['fingerprint: ' fplab(iii)]);
    pdfwrite(['fprint_' fplab{iii} '_map1']);    
    xlim([220 300]); ylim([0 60]);
    pdfwrite(['fprint_' fplab{iii} '_map2']);    
end


fid=fopen('Fingerprints.tsv','w');
fprintf(fid,['site\tID\tlat\tlong']);
for jjj=1:length(fplab)
    fprintf(fid,['\t' fplab{jjj}]);
end
fprintf(fid,'\n');

for iii=1:length(targregions)
    fprintf(fid,targregionnames{iii});
    fprintf(fid,'\t%0.0f',targregions(iii));
    fprintf(fid,'\t%0.2f',targsitecoords(iii,:));
    fprintf(fid,'\t%0.2f',fpsite(iii,1:length(fplab)));
    fprintf(fid,'\n');
    
end

fclose(fid);
