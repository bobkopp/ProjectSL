% outputSLRProjections_plot_GSL
% Last updated by  Bob Kopp, robert-dot-kopp-at-rutgers-dot-edu, Mon Dec 30 10:28:47 EST 2013

%%%

% GSL projection plots
colrs='rcbg';

quantsub=[0.01 0.05 0.167 0.5 0.833 0.95 0.99];
[jk,ia,ib]=intersect(quantlevs,quantsub);

scenlab={};
for kk=subscens
	scenlab={scenlab{:},[upper(scens{kk}(1:3)) ' ' scens{kk}(4) '.' scens{kk}(5)]};
end

clf; subplot(2,2,1); clear hp;
for kk=subscens
	hp(kk)=plot([2000 targyears],[0  quanttotrise((ia(4)),:,kk)/1e3],[colrs(kk)],'linew',2); hold on;
	plot([2000 targyears],[0  quanttotrise((ia([2])),:,kk)/1e3],[colrs(kk),'--']);
	plot([2000 targyears],[0  quanttotrise((ia([6])),:,kk)/1e3],[colrs(kk),'--']);
	plot([2000 targyears],[0  quanttotrise((ia([1 ])),:,kk)/1e3],[colrs(kk),':']);
	plot([2000 targyears],[0 quanttotrise((ia([7])),:,kk)/1e3],[colrs(kk),':']);

end
ylabel('m GSL rise'); ylim([0 3]);
legend(hp(subscens),scenlab{:},'location','NorthWest'); xlim([2000 2200]);
longticks(gca)
pdfwrite('GSLrise_allscens');