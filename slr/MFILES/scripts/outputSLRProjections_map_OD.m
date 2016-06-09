% outputSLRProjections_map_OD
% Last updated by Robert Kopp, robert-dot-kopp-at-rutgers-dot-edu, Fri May 23 15:26:23 EDT 2014

%%
% ocean dynamics map

%[jk,ia,ib]=intersect(OceanDynRegions,targregions);
%
%sub=1:length(ib);

%subyrs=find(OceanDynYears==2100);

%clf;
%plotcont; hold on;
%scatter(mod(sitecoords(ib(sub),2),360),sitecoords(ib(sub),1),15,OceanDynMean(subyrs,ia(sub))/1000,'filled');
%axis tight; box on;
%title('Oceanographic mean: RCP 8.5, 2100 (m)'); colorbar;
%pdfwrite('OD_mean');
%
%clf;
%plotcont; hold on;
%scatter(mod(sitecoords(ib(sub),2),360),sitecoords(ib(sub),1),15,OceanDynStd(subyrs,ia(sub))/1000,'filled');
%axis tight;
%title('Oceanographic std: RCP 8.5, 2100 (m)'); colorbar; box on;
%pdfwrite('OD_std');
%

subyear=find(targyears==2100);
kk=1;

nn=9;

clf;
worldmap('world');
geoshow('landareas.shp','facecolor','none');
setm(gca,'parallellabel','off','meridianlabel','off','flinewidth',1);

sb=find(quantlevs==0.5);
u=squeeze(quantlocscalefactors(sb,nn,subyear,:,kk));
hp=scatterm(sitecoords(:,1),sitecoords(:,2),15,u/10,'filled');
axis off; 
title({['Median ocean dynamics: RCP 8.5, 2100 (cm)']}); colorbar; 
pdfwrite(['OD_median']);

% $$$ plotcont; hold on;
% $$$ scatter(mod(sitecoords(:,2),360),sitecoords(:,1),15,u/10,'filled');
% $$$ axis tight; box on;
% $$$ title({['Median ocean dynamics: RCP 8.5, 2100 (cm)']}); colorbar; 
% $$$ pdfwrite(['OD_median']);

sb=find(round(quantlevs*100)==17);
sb(2)=find(round(quantlevs*100)==83);
delete(hp);
u=diff(squeeze(quantlocscalefactors(sb,nn,subyear,:,kk)),[],1);
hp=scatterm(sitecoords(:,1),sitecoords(:,2),15,u/10,'filled');
caxis([10 50]);
title({['Ocean dynamics range (17th-83rd %ile range): RCP 8.5, 2100 (cm)']}); colorbar; 
pdfwrite(['OD_range']);

sb=find(round(quantlevs*100)==17);
sb(2)=find(round(quantlevs*100)==83);
u=diff(squeeze(quantlocscalefactors(sb,nn,subyear,:,kk)),[],1);
delete(hp);
hp=scatterm(sitecoords(:,1),sitecoords(:,2),15,u,'filled');
caxis(quantile(u,[.05 .95]));
title({['Scale factor (17th-83rd range width): RCP 8.5, 2100']}); colorbar; 
pdfwrite(['locscalefactorrange' num2str(nn)]);

delete(hp);
sb=find(quantlevs==0.5);
u=squeeze(OceanDynN(subyear,:,kk));
hp=scatterm(sitecoords(:,1),sitecoords(:,2),15,u,'filled');
title({['Number of ocean dynamic models: RCP 8.5, 2100 (cm)']}); colorbar; 
caxis([min(u) max(u)]);
pdfwrite(['OD_N']);