function [quanttotrise,quantcomponents]=SummarizeGSLProjections(samps,quantlevs,colGIC,colGIS,colAIS,colLS,colTE,coladdls,coladdlsorigin)

% Last updated by  Bob Kopp, robert-dot-kopp-at-rutgers-dot-edu, Tue Dec 31 13:53:07 EST 2013

colIS=union(colGIS,colAIS);

defval('coladdls',[]);
defval('coladdlsorigin',{[colGIC],[colAIS],[colIS],[colGIC colIS],[colGIC colIS colLS],[colGIC colIS colLS colTE]});

GSLsamps = squeeze(sum(samps,2));
quanttotrise = squeeze(quantile(GSLsamps,quantlevs,1));
quantcomponents = quantile([samps],quantlevs,1);
for i=1:length(coladdls)
	quantcomponents(:,coladdls(i),:,:) = quantile(sum(samps(:,coladdlsorigin{i},:,:),2),quantlevs,1);
end
