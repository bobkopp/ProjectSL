% Last updated by Robert Kopp, robert-dot-kopp-at-rutgers-dot-edu, Mon May 12 12:12:44 EDT 2014

configureSLRProjections;

%% PSMSL data
psmsldir=fullfile(IFILES,'rlr_annual');
gslfile=fullfile(IFILES,'CSIRO_Recons_gmsl_yr_2011.csv');

[~,~,~,~,regionsu,sitenames,sitecoords,sitelen,sitecoastline]=ReadPSMSLData([],[],10,psmsldir,gslfile);

sub=find(sitecoords(:,1)<1e3);

sitecoords=sitecoords(sub,:);
[targregions] = regionsu(sub); % tide gauges to work with
targregionnames = sitenames(sub);
targcoastlines=sitecoastline(sub);

% GIA (slow)

[rateprojs,rateprojssd,rateprojs0,~,rateGIAproj,~,~,~,finescale]=CalculateBackgroundRates([],sitecoords,[],[],PARAMDIR,IFILES);

% ocean dynamics (slow)
[ThermExpMean,ThermExpStd,ThermExpYears,ThermExpN,OceanDynMean,OceanDynStd,OceanDynYears,OceanDynRegions,OceanDynN,OceanDynTECorr,ZOSTOGA,sZOSTOGA,ZOSTOGAmodels,ZOSTOGAyrs] = CalculateOceanographicDist(scens,datayears,targregions,sitecoords,mergeZOSZOSTOGA,IFILES,PARAMDIR);

% Glaciers and ice caps

clear projGIC projGICse projGICyrs projGICmodel ;
for kk=1:length(scens)
	disp(scens{kk});
	[projGIC{kk},projGICse{kk},projGICyrs{kk},projGICmodel{kk},fpmapperids,fpmaps,GICnames]=readMarzeion(scens{kk},IFILES,PARAMDIR);

end

[fpsite,fp,fpname,lo,la]=AssignFingerprints(fpmapperids,fpmaps,sitecoords,[],fullfile(IFILES,'FPRINT'));

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

% run through alternative correlations
if iscell(BAcorrIS)
	for i=1:length(BAcorrIS)
		[ISaccelsamps{i},ARISaccelsamps{i},lastseed] = SampleISDists(thetEAIS,thetWAIS,thetGIS,ARthetAIS,ARthetGIS,ISLastDecade,BAcorrIS{i},ARcorrIS{i},seeds);
	end
else
	[ISaccelsamps{1},ARISaccelsamps{1},lastseed] = SampleISDists(thetEAIS,thetWAIS,thetGIS,ARthetAIS,ARthetGIS,ISLastDecade,BAcorrIS,ARcorrIS,seeds);
end

seeds=seeds(lastseed+1:end,:); lastseed=0;

if exist('savefile','var')
	save(savefile);
end

%%%%%%%%%%%%%%

%% PROJECTION

quantlevs = [.001 .005 .01 .025 .05 .1 .167 .333 .5 .667 .833 .9 .95 .975 .99 .995 .999];
quantlevs=sort(union(quantlevs,0:.05:1));

seeds0=seeds;

[samps, colGIC, colIS, colGIS, colAIS, colLS, colTE, lastseed] = ProjectGSL(scens, targyears, projGICyrs, projGIC,projGICse, seeds, ISLastDecade, ISaccelsamps{1}, ARISaccelsamps{1}, ISmode,ThermExpYears, ThermExpMean,  ThermExpStd, ThermExpN, LWSpath, quantlevs);

% REGIONALIZATION - PROJECTIONS

seeds=seeds(lastseed+1:end,:); lastseed=0;

focussites = [12 299 396  188 161 10 405 155 43 269 860 526 235 88 1];
sitesub=[
183 % Portland, ME
235 % Boston, MA
12 % New York, NY
351 % Newport, CT
180 % Atlantic City, NJ
135 % Philadelphia, PA
224 % Lewes, DE
148 % Baltimore, MD
360 % Washington, DC
299 % Sewell's Pt, Norfolk, VA
396 % Wilmington, NC
234 % Charleston, NC
395 % Fort Pulaski, Savannah, GA
112 % Ferdinanda Beach, FL
188 % Key West, FL
363 % Miami, FL
246 % Pensacola, FL
526 % Grand Isle, LA
161 % Galveston, TX
158 % San Diego, CA
10 % San Francisco, CA
265 % Astoria, OR
127 % Seattle, WA
405 % Juneau, AK
1067 % Anchorage, AK
155 % Honolulu, HI
78  % Stockholm
7  % Cuxhaven, Germany
438  % Kochi, India
134 % Kushimoto, Japan
499 % Valparaiso, Chile
42 % Sevastopol
]; sitesub=sitesub';
focussites=union(focussites,sitesub);

[quanttotlocrise, quantloccomponents,quantlocscalefactors,totlocalrisefocus,sampsregionfocus,colGIA,colOD,colGICtot,colAIStot,colIStot,colLItot,colLILStot,colLILSTEODtot,colLILSTEODGIAtot,lastseed] = ProjectLSL(scens,targregions,targregionnames,targyears,samps,seeds,targregions,OceanDynRegions,OceanDynYears,OceanDynMean,OceanDynStd,OceanDynN, ThermExpYears, ThermExpMean, ThermExpStd, OceanDynTECorr,  rateprojs,rateprojssd,mergeZOSZOSTOGA,fpsite,quantlevs,colGIC,colGIS,colAIS,colLS,colTE,focussites);


% summarize GSL

colLILSTEtot = colLILSTEODtot;
coladdls=[colGICtot colAIStot colIStot colLItot colLILStot colLILSTEtot];
coladdlsorigin = {[colGIC],[colAIS],[colIS],[colGIC colIS],[colGIC colIS colLS],[colGIC colIS colLS colTE]};

[quanttotrise,quantcomponents]=SummarizeGSLProjections(samps,quantlevs,colGIC,colGIS,colAIS,colLS,colTE,coladdls,coladdlsorigin);

% SENSITIVITY TESTS

[targregions2,ia,ib]=intersect(targregions,focussites);
[ib,ic]=sort(ib); ia=ia(ic); targregions2=targregions2(ic);
targregions2names=targregionnames(ia);

% default

sampsAlt=samps;
[quanttotlocriseDef,quantloccomponentsDef] =  ProjectLSL(scens,targregions(ia),targregionnames(ia),targyears,sampsAlt,seeds,targregions,OceanDynRegions,OceanDynYears,OceanDynMean,OceanDynStd,OceanDynN, ThermExpYears, ThermExpMean, ThermExpStd, OceanDynTECorr,  rateprojs,rateprojssd,mergeZOSZOSTOGA,fpsite,quantlevs,colGIC,colGIS,colAIS,colLS,colTE,focussites);

% - varying AR5 vs. BA vs. hybrid for ice sheets

[sampsAlt] = ProjectGSL(scens, targyears, projGICyrs, projGIC,projGICse, seeds0, ISLastDecade, ISaccelsamps{1}, ARISaccelsamps{1}, 'AR', ThermExpYears, ThermExpMean,  ThermExpStd, ThermExpN, LWSpath, quantlevs);
[quanttotriseAR,quantcomponentsAR]=SummarizeGSLProjections(sampsAlt,quantlevs,colGIC,colGIS,colAIS,colLS,colTE,coladdls,coladdlsorigin);
[quanttotlocriseAR] =  ProjectLSL(scens,targregions(ia),targregionnames(ia),targyears,sampsAlt,seeds,targregions,OceanDynRegions,OceanDynYears,OceanDynMean,OceanDynStd,OceanDynN, ThermExpYears, ThermExpMean, ThermExpStd, OceanDynTECorr,  rateprojs,rateprojssd,mergeZOSZOSTOGA,fpsite,quantlevs,colGIC,colGIS,colAIS,colLS,colTE,focussites);


[sampsAlt] = ProjectGSL(scens, targyears, projGICyrs, projGIC,projGICse, seeds0, ISLastDecade, ISaccelsamps{1}, ARISaccelsamps{1}, 'BA', ThermExpYears, ThermExpMean,  ThermExpStd, ThermExpN, LWSpath, quantlevs);
[quanttotriseBA,quantcomponentsBA]=SummarizeGSLProjections(sampsAlt,quantlevs,colGIC,colGIS,colAIS,colLS,colTE,coladdls,coladdlsorigin);
[quanttotlocriseBA] =  ProjectLSL(scens,targregions(ia),targregionnames(ia),targyears,sampsAlt,seeds,targregions,OceanDynRegions,OceanDynYears,OceanDynMean,OceanDynStd,OceanDynN, ThermExpYears, ThermExpMean, ThermExpStd, OceanDynTECorr,  rateprojs,rateprojssd,mergeZOSZOSTOGA,fpsite,quantlevs,colGIC,colGIS,colAIS,colLS,colTE,focussites);

% - varying ice sheet correlations

if iscell(BAcorrIS)
	for i=2:length(BAcorrIS)
		[sampsAlt] = ProjectGSL(scens, targyears, projGICyrs, projGIC,projGICse, seeds0, ISLastDecade, ISaccelsamps{i}, ARISaccelsamps{i}, ISmode, ThermExpYears, ThermExpMean,  ThermExpStd, ThermExpN, LWSpath, quantlevs);
		[quanttotriseAltCorr{i-1},quantcomponentsAltCorr{i-1}]=SummarizeGSLProjections(sampsAlt,quantlevs,colGIC,colGIS,colAIS,colLS,colTE,coladdls,coladdlsorigin);
		[quanttotlocriseAltCorr{i-1}] =  ProjectLSL(scens,targregions(ia),targregionnames(ia),targyears,sampsAlt,seeds,targregions,OceanDynRegions,OceanDynYears,OceanDynMean,OceanDynStd,OceanDynN, ThermExpYears, ThermExpMean, ThermExpStd, OceanDynTECorr,  rateprojs,rateprojssd,mergeZOSZOSTOGA,fpsite,quantlevs,colGIC,colGIS,colAIS,colLS,colTE,focussites);
	end
end

% - varying confidence in models

[sampsAlt] = ProjectGSL(scens, targyears, projGICyrs, projGIC,projGICse, seeds0, ISLastDecade, ISaccelsamps{1}, ARISaccelsamps{1}, ISmode, ThermExpYears, ThermExpMean,  ThermExpStd, ThermExpN, LWSpath, quantlevs,.95);
[quanttotriseHiGCMConf,quantcomponentsHiGCMConf]=SummarizeGSLProjections(sampsAlt,quantlevs,colGIC,colGIS,colAIS,colLS,colTE,coladdls,coladdlsorigin);
[quanttotlocriseHiGCMConf] =  ProjectLSL(scens,targregions(ia),targregionnames(ia),targyears,sampsAlt,seeds,targregions,OceanDynRegions,OceanDynYears,OceanDynMean,OceanDynStd,OceanDynN, ThermExpYears, ThermExpMean, ThermExpStd, OceanDynTECorr,  rateprojs,rateprojssd,mergeZOSZOSTOGA,fpsite,quantlevs,colGIC,colGIS,colAIS,colLS,colTE,focussites,.94);


[sampsAlt] = ProjectGSL(scens, targyears, projGICyrs, projGIC,projGICse, seeds0, ISLastDecade, ISaccelsamps{1}, ARISaccelsamps{1}, ISmode, ThermExpYears, ThermExpMean,  ThermExpStd, ThermExpN, LWSpath, quantlevs,[],5);
[quanttotriseRedDOF,quantcomponentsRedDOF]=SummarizeGSLProjections(sampsAlt,quantlevs,colGIC,colGIS,colAIS,colLS,colTE,coladdls,coladdlsorigin);
[quanttotlocriseRedDOF,quantloccomponentsRedDOF] =  ProjectLSL(scens,targregions(ia),targregionnames(ia),targyears,sampsAlt,seeds,targregions,OceanDynRegions,OceanDynYears,OceanDynMean,OceanDynStd,OceanDynN, ThermExpYears, ThermExpMean, ThermExpStd, OceanDynTECorr,  rateprojs,rateprojssd,mergeZOSZOSTOGA,fpsite,quantlevs,colGIC,colGIS,colAIS,colLS,colTE,focussites,[],5);


% - varying land water storage

[sampsAlt] = ProjectGSL(scens, targyears, projGICyrs, projGIC,projGICse, seeds0, ISLastDecade, ISaccelsamps{1}, ARISaccelsamps{1}, ISmode, ThermExpYears, ThermExpMean,  ThermExpStd, ThermExpN, LWSpathwPokhrel, quantlevs);
[quanttotrisePokhrel,quantcomponentsPokhrel]=SummarizeGSLProjections(sampsAlt,quantlevs,colGIC,colGIS,colAIS,colLS,colTE,coladdls,coladdlsorigin);
[quanttotlocrisePokhrel,quantloccomponentsPokhrel] =  ProjectLSL(scens,targregions(ia),targregionnames(ia),targyears,sampsAlt,seeds,targregions,OceanDynRegions,OceanDynYears,OceanDynMean,OceanDynStd,OceanDynN, ThermExpYears, ThermExpMean, ThermExpStd, OceanDynTECorr,  rateprojs,rateprojssd,mergeZOSZOSTOGA,fpsite,quantlevs,colGIC,colGIS,colAIS,colLS,colTE,focussites);
%sampsdiff=sampsregion(:,1:end-2,:)-samps;
%sampsdiff(:,end+1:end+2,:)=sampsregion(:,end-1:end,:);

if exist('savefile2','var')
save(savefile2,'-v7.3');	%save(savefile2,'quantlevs','quanttotrise','quantcomponents','quantlocscalefactors','quantloccomponents','quanttotlocrise','colGIC','colGIS','colAIS','colTE','colLS','colOD','colGIA','colGICtot','colAIStot','colIStot','colLItot','colLILStot','colLILSTEODtot','colLILSTEtot','totlocalrisefocus','focussites','sampsregionfocus');
end

%%%

if ~exist(outputdir,'dir')
	mkdir(outputdir);
end
cd(outputdir);

outputSLRProjections
