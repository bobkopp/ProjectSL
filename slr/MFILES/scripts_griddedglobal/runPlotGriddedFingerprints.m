% Last updated by Robert Kopp, robert-dot-kopp-at-rutgers-dot-edu, Fri Jan 13 09:33:06 EST 2017


% plot fingerprints
for iii=1:length(fplab)
    clf;
    scatter(mod(targsitecoords(:,2),360),targsitecoords(:,1),10,fpsite(:,iii),'filled');
    hold on;
    plotcont([],[],1);
    colorbar;
    colormap(jet);
    caxis([0 1.3]);
    title(['fingerprint: ' fplab(iii)]);
    pdfwrite(['fprint_' fplab{iii} '_map1']);    
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
