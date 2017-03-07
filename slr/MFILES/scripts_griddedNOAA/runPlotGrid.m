% Last updated by Robert Kopp, robert-dot-kopp-at-rutgers-dot-edu, Thu Oct 27 10:40:17 EDT 2016

% plot grid

clf;
%worldmap('world');
worldmap([-20 73],[140 -45]);
geoshow('landareas.shp','FaceColor','none');
scatterm(mod(targsitecoords(:,1),360),targsitecoords(:,2),5);
hold on;
pdfwrite('grid');


fid=fopen('grid.txt','w');
for iii=1:length(targsitecoords)
    fprintf(fid,'%0.0f',targregions(iii));
    fprintf(fid,['\t' targregionnames{iii}]);
    fprintf(fid,'\t%0.3f',targsitecoords(iii,:));
    fprintf(fid,'\n');
end
fclose(fid);