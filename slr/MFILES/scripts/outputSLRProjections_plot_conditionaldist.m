% outputSLRProjections_plot_conditionaldist
% Last updated by  Bob Kopp, robert-dot-kopp-at-rutgers-dot-edu, 9 Jan 2014

% GSL
colrs='rgcb';

subcomp={colAIS,colGIS,colGIC,colTE};
labls={'AIS','GIS','GIC','TE'};
kk=1;
titl=['conditional distribution (' scens{kk} ')'];
rounder = 100;
subyear=find(targyears==2100);

wsamps=samps(:,:,subyear,kk);
wsampsgsl=round(sum(wsamps,2)/rounder)*rounder;
wqlevs=[.5 .167 .833];

u=unique(wsampsgsl);
clear wq wn;
for i=1:length(u)
	subgsl = find(wsampsgsl==u(i));
	wn(i)=length(subgsl);
	for j=1:length(subcomp)
		wq(i,:,j) = quantile(sum(wsamps(subgsl,subcomp{j}),2),wqlevs);
	end
end

clf;
clear hl hp;
for i=1:length(labls)
	sub=find(wn>10);
	hl(i)=plot(u(sub)/10,wq(sub,1,i)/10,colrs(i),'linew',2); hold on;
	plot(u(sub)/10,wq(sub,[2 3],i)/10,[colrs(i) '--'],'linew',1); hold on;
end
xlabel('GSL (cm)');
ylabel('Component (cm)');
legend(hl,labls,'Location','Northwest');
title(titl);
pdfwrite(['conditionaldist_' scens{kk} '_2100']);

% and locally


if colOD==colTE
	colTE2=[];
else
	colTE2=colTE;
end

subcomp={colAIS,colGIS,colGIC,[colTE2 colOD]};
labls={'AIS','GIS','GIC','Oc'};
kk=1;

colrs='rgcb';
sitesub=[12 299 396  188 161 10 405 155];


for jjj=1:length(sitesub)
	jj2=find(focussites==sitesub(jjj));
	if length(jj2)==1
		jj=find(targregions==sitesub(jjj));
		titl={targregionnames{jj},['conditional distribution (' scens{kk} ')']};

		wsamps=sampsregionfocus{jj2,kk}(:,:,subyear);
		wsampslsl=round(sum(wsamps,2)/rounder)*rounder;
		wqlevs=[.5 .167 .833];

		u=unique(wsampsgsl);
		clear wq wn;
		for i=1:length(u)
			sublsl = find(wsampslsl==u(i));
			wn(i)=length(sublsl);
			for j=1:length(subcomp)
				wq(i,:,j) = quantile(sum(wsamps(sublsl,subcomp{j}),2),wqlevs);
			end
		end
		
		clf;
		clear hl hp;
		for i=1:length(labls)
			sub=find(wn>10);
			hl(i)=plot(u(sub)/10,wq(sub,1,i)/10,colrs(i),'linew',2); hold on;
			plot(u(sub)/10,wq(sub,[2 3],i)/10,[colrs(i) '--'],'linew',1); hold on;
		end
		xlabel('LSL (cm)');
		ylabel('Component (cm)');
		legend(hl,labls,'Location','Northwest');
		title(titl);
		pdfwrite(['conditionaldist_' scens{kk} '_2100_' num2str(sitesub(jjj))]);
	end
end

