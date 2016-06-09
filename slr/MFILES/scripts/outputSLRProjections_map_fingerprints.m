% outputSLRProjections_plot_fingerprints
% Last updated by  Bob Kopp, robert-dot-kopp-at-rutgers-dot-edu, Mon Dec 30 10:31:56 EST 2013

icesheets={'gis','wais','eais'};

%% fingerprint maps

ii=10;
[LO,LA]=meshgrid(lo,la);
clear fptouse;
for i=1:length(fpmapperids)
	fptouse(i)=strmatch(fpmaps{i},fpname,'exact');
end

subql=find(quantlevs==0.5); subyear=find(targyears==2100); kk=1;
medianGIC = squeeze(quantcomponents(subql,[1:length(fptouse)],subyear,kk));
medianGICfpmap = sum(bsxfun(@times,medianGIC,reshape(fp(:,:,fptouse),length(lo)*length(la),length(fptouse))),2)/sum(medianGIC)*1000;
medianGICfpmap = reshape(medianGICfpmap,size(LO));

clf;
imagesc(lo,la,medianGICfpmap); plotcont;
colorbar; caxis([-.4 1.4]);
axis tight; set(gca,'ydir','norm');
title('Fingerprint: Glaciers');
pdfwrite('fp_medianGIC');

for i=1:length(icesheets);
	clf;
	fptouse=strmatch(icesheets{i},fpname,'exact');
	imagesc(lo,la,1000*fp(:,:,fptouse)); plotcont; colorbar; caxis([-.4 1.4]);
	axis tight; set(gca,'ydir','norm');
	title(['Fingerprint: ' upper(icesheets{i})]);
	pdfwrite(['fp_' icesheets{i}]);
end
