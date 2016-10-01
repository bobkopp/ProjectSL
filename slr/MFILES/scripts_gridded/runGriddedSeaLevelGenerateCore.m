% Last updated by Robert Kopp, robert-dot-kopp-at-rutgers-dot-edu, Thu Sep 29 13:25:08 EDT 2016

addpath(pwd)
workdir='workdir-160929';
if ~exist(workdir,'dir')
    mkdir(workdir)
end
cd(workdir);

%%%%%

runConfigureGriddedSLRProjections
runSetupSLRGrid
runPlotGrid

% GIA (slow)
[rateprojs,rateprojssd,rateprojs0,~,rateGIAproj,~,~,~,finescale]=CalculateBackgroundRates([],targsitecoords,[],[],PARAMDIR,IFILES);
runPlotBackgroundRates

% GIC

clear projGIC projGICse projGICyrs projGICmodel ;
for kk=1:length(scens)
    disp(scens{kk});
    [projGIC{kk},projGICse{kk},projGICyrs{kk},projGICmodel{kk},fpmapperids,fpmaps,GICnames]=readMarzeion(scens{kk},IFILES,PARAMDIR);

end
[fpsite,fp,fpname,lo,la]=AssignFingerprints(fpmapperids,fpmaps,targsitecoords,[],fullfile(IFILES,'FPRINT'));
fplab={GICnames{:},'gis','wais','eais'};
runPlotGriddedFingerprints

% now need to do oceanographic contribution (slow)

wSLRDIR='~/NASA/zosfinal';
filemode=1;
[ThermExpMean,ThermExpStd,ThermExpYears,ThermExpN,OceanDynMean,OceanDynStd,OceanDynYears,OceanDynRegions,OceanDynN,OceanDynTECorr,ZOSTOGA,sZOSTOGA,ZOSTOGAmodels,ZOSTOGAyrs] = CalculateOceanographicDist(scens,datayears,targregions,targsitecoords,mergeZOSZOSTOGA,IFILES,PARAMDIR,wSLRDIR,filemode);

runPlotGriddedOceanographicDist

%%%%%%

% seeds

seeds0=linspace(0,1,Nsamps+2); seeds0=seeds0(2:end-1);
seeds0=norminv(seeds0);
clear seeds
for i=1:Nseeds
	seeds(i,:) =seeds0(randperm(Nsamps));
end

lastseed=0;


% ICE SHEETS - BA and AR distributions

[thetEAIS,thetWAIS,thetGIS,ARthetAIS,ARthetGIS,ISLastDecade]=CalculateISDists(ratesmatrix2100,LastDecadeGt,ARIS2090);
seeds=seeds(lastseed+1:end,:); lastseed=0;

% run through alternative correlations
if iscell(BAcorrIS)
	for i=1:length(BAcorrIS)
		[ISaccelsamps{i},ARISaccelsamps{i},lastseed] = SampleISDists(thetEAIS,thetWAIS,thetGIS,ARthetAIS,ARthetGIS,ISLastDecade,BAcorrIS{i},ARcorrIS{i},seeds);
	end
else
	[ISaccelsamps{1},ARISaccelsamps{1},lastseed] = SampleISDists(thetEAIS,thetWAIS,thetGIS,ARthetAIS,ARthetGIS,ISLastDecade,BAcorrIS,ARcorrIS,seeds);
end


%% PROJECTION

quantlevs = [.001 .005 .01 .025 .05 .1 .167 .333 .5 .667 .833 .9 .95 .975 .99 .995 .999];
quantlevs=sort(union(quantlevs,0:.05:1));

seeds0=seeds;

[samps, colGIC, colIS, colGIS, colAIS, colLS, colTE, lastseed] = ProjectGSL(scens, targyears, projGICyrs, projGIC,projGICse, seeds, ISLastDecade, ISaccelsamps{1}, ARISaccelsamps{1}, ISmode,ThermExpYears, ThermExpMean,  ThermExpStd, ThermExpN, LWSpath, quantlevs);

%%%%

savefilecore='~/tmp/SLRProjections160930GRIDDEDcore';
save(savefilecore,'OceanDynMean','OceanDynN','OceanDynRegions','OceanDynStd','OceanDynTECorr','OceanDynYears','ThermExpMean','ThermExpStd','ThermExpYears','colAIS','colGIC','colGIS','colLS','colTE','fpsite','mergeZOSZOSTOGA','quantlevs','rateprojs','rateprojssd','samps','scens','seeds','targregionnames','targregions','targyears','nearestTG','targsitecoords');
