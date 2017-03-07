% Last updated by Robert Kopp, robert-dot-kopp-at-rutgers-dot-edu, Fri Jan 13 08:41:23 EST 2017
% plot grid

clf;
worldmap('world');
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