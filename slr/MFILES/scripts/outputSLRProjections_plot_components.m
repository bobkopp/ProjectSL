% outputSLRProjections_plot_components
% Last updated by  Bob Kopp, robert-dot-kopp-at-rutgers-dot-edu, Mon Dec 30 10:31:26 EST 2013

%% Component projections
quantsub=[.5 .167 .833 .05 .95 .005 .995];
[jk,ia,ib]=intersect(quantlevs,quantsub); [jk,ic]=sort(ib);
ia=ia(ic);

subcomp=[colGIS colAIStot colTE colGICtot ];
labls={'Greenland','Antarctica','Thermal expansion','Glaciers & ice caps'};

kk=1;

clf;
colrs='kkkkk';
linew=[2 1 1 1 1 1 1];
symb={'-','-','-','--','--',':',':'};
clear hl hp;
for i=1:length(labls)
	hp(i)=subplot(2,2,i);
	jj=1;
	hl(i)=plot([2000 targyears],[0 ; squeeze(quantcomponents(ia(jj),subcomp(i),:,kk))]/1000,[colrs(i) symb{jj}],'linew',linew(jj)); hold on;
	for jj=2:length(quantsub)
		plot([2000 targyears],[0 ; squeeze(quantcomponents(ia(jj),subcomp(i),:,kk))]/1000,[colrs(i) symb{jj}],'linew',linew(jj));
		if jj==5
			yl=get(gca,'ylim');
		end
	end
	title(labls{i});
	ylabel('m');
	xlim([2000 2200]);
	if i==1; ylim([0 0.6]);
	elseif i==2; ylim([0 1.5]);
	elseif i==3; ylim([-.5 1.5]);
	elseif i==4; ylim([0 1.5]);
	end
	longticks(gca);
end
pdfwrite(['components_' scens{kk}]);

set(hp,'xlim',[2000 2100]);
%set(hp(1),'ylim',[0 0.3]);
%set(hp(2),'ylim',[0 0.8]);
%set(hp(3),'ylim',[-.2 1]);
%set(hp(4),'ylim',[0 0.7]);
set(hp,'ylim',[-.15 .6]);
clear ht;
ltrs='abcd';
for ii=1:4
    axes(hp(ii)); ht(ii)=text(2005,0.6,ltrs(ii));
end
set(ht,'fontw','bold','fontsize',12,'verticalalignment','top','horizontalal','left')
pdfwrite(['components_' scens{kk} '_2100']);


subscens=[1 3 4];
linew=[2 1 1 1 1 1 1];
symb={'-','--','--'};
clf;
colrs='rcbg';
clear hl hp;
for i=1:length(labls)
	hp(i)=subplot(2,2,i);
	for kk=subscens
		for jj=2:3
			plot([2000 targyears],[0 ; squeeze(quantcomponents(ia(jj),subcomp(i),:,kk))]/1000,[colrs(kk) symb{jj}],'linew',linew(jj)); hold on;
		end
	end
	jj=1;
	for kk=subscens
		hl(kk)=plot([2000 targyears],[0 ; squeeze(quantcomponents(ia(jj),subcomp(i),:,kk))]/1000,[colrs(kk) symb{jj}],'linew',linew(jj)); hold on;
	end
	title(labls{i});
	ylabel('m');
	xlim([2000 2200]);
	if i==1; ylim([0 0.4]);
	elseif i==2; ylim([0 1]);
	elseif i==3; ylim([-.5 1]);
	elseif i==4; ylim([0 1.1]);
	end
	longticks(gca);
end
set(hp,'ylim',[-.5 1]);
clear ht;
ltrs='abcd';
for ii=1:4
    axes(hp(ii)); ht(ii)=text(2010,.95,ltrs(ii));
end
set(ht,'fontw','bold','fontsize',12,'verticalalignment','top','horizontalal','left')
pdfwrite(['components_allscens']);
