% outputSLRProjections_map_GIAdomains
% Last updated by  Bob Kopp, robert-dot-kopp-at-rutgers-dot-edu, Wed 12 Feb 2014 08:41:45 EST

%%%

% domains of background rates
% must step through the identiifcation process in CalculateBackgroundRates to get targlong, targlat, nearest

[~,~,~,targcoord,~,~,~,nearest]=CalculateBackgroundRates([],[],[],0,PARAMDIR,IFILES);

if ~exist('sitecoords','var')
	[~,~,~,~,~,~,sitecoords]=ReadPSMSLData([],[],10,psmsldir,gslfile);
end

s=sitecoords(find(sitecoords(:,1)<1000),:);
[~,~,~,~,~,~,~,nearest2]=CalculateBackgroundRates([],s,[],0,PARAMDIR,IFILES);

figure;
worldmap('world');

geoshow('landareas.shp','facecolor','none');
setm(gca,'parallellabel','off','meridianlabel','off','flinewidth',1);

%plotcont; hold on; axis tight;
%s(:,2)=mod(s(:,2),360);
%plot(s(:,2),s(:,1),'g.');
%scatter(s(:,2),s(:,1),10,nearest2,'filled');

[targlong,ui]=unique(targcoord(:,2),'first');
TARGLAT=reshape(targcoord(:,1),length(targlong),[]);
targlat=TARGLAT(1,:);

NEAREST=reshape(nearest,length(targlong),length(targlat));
[sa,si]=sort(mod(targlong,360));
%[c,h]=contour(s,targlat,NEAREST(si,:)',-.5:.5:30); hold on
[c,h]=contourm(targlat,sa,NEAREST(si,:)',-.5:.5:30,'LineColor',[.4 .4 .4]); hold on
hp=scatterm(s(:,1),s(:,2),10,nearest2,'filled');
%title('Domains');
%box on; set(gca,'xtick',[],'ytick',[]);
colormap(prism);
pdfwrite('map_domains');
