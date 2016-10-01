% Last updated by Robert Kopp, robert-dot-kopp-at-rutgers-dot-edu, Fri Sep 30 09:29:22 EDT 2016

% plot grid

clf;
scatter(mod(targsitecoords(:,2),360),targsitecoords(:,1),10);
hold on;
plotcont;
pdfwrite('grid');


fid=fopen('grid.txt','w');
for iii=1:length(targsitecoords)
    fprintf(fid,'%0.0f',targregions(iii));
    fprintf(fid,['\t' targregionnames{iii}]);
    fprintf(fid,'\t%0.3f',targsitecoords(iii,:));
    fprintf(fid,'\n');
end
fclose(fid);