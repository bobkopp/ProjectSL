% outputSLRProjections_plot_LSL
% Last updated by  Bob Kopp, robert-dot-kopp-at-rutgers-dot-edu, Jan 8 2014

%%%

subscens=[1 3 4];
sitesub=[12 299 396  188 161 10 405 155];

% LSL projection plots
colrs='rcbg';

quantsub=[0.01 0.05 0.167 0.5 0.833 0.95 0.99];
[jk,ia,ib]=intersect(quantlevs,quantsub);

scenlab={};
for kk=subscens
	scenlab={scenlab{:},[upper(scens{kk}(1:3)) ' ' scens{kk}(4) '.' scens{kk}(5)]};
end

for jjj=1:length(sitesub)
	jj=find(targregions==sitesub(jjj));
	clf;
	subplot(2,2,1);

	for kk=subscens
		hp(kk)=plot([2000 targyears],[0  quanttotlocrise((ia(4)),:,jj,kk)/1e3],[colrs(kk)],'linew',2); hold on;
		plot([2000 targyears],[0  quanttotlocrise((ia([2])),:,jj,kk)/1e3],[colrs(kk),'--']);
		plot([2000 targyears],[0  quanttotlocrise((ia([6])),:,jj,kk)/1e3],[colrs(kk),'--']);
		plot([2000 targyears],[0  quanttotlocrise((ia([1 ])),:,jj,kk)/1e3],[colrs(kk),':']);
		plot([2000 targyears],[0 quanttotlocrise((ia([7])),:,jj,kk)/1e3],[colrs(kk),':']);
	end
	ylabel('m'); ylim([0 3]);
	title(targregionnames{jj})
	legend(hp(subscens),scenlab{:},'location','NorthWest'); xlim([2000 2200]);
	pdfwrite(['LSLrise_allscens_' num2str(sitesub(jjj))]);
end

