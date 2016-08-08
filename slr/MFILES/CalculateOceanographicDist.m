function [ThermExpMean,ThermExpStd,ThermExpYears,ThermExpN,OceanDynMean,OceanDynStd,OceanDynYears,OceanDynRegions,OceanDynN,OceanDynTECorr,ZOSTOGA,sZOSTOGA,ZOSTOGAmodels,ZOSTOGAyrs] = CalculateOceanographicDist(scens,datayears,targregions,sitecoords,mergeZOSZOSTOGA,IFILES,PARAMDIR)

% Last updated by Robert Kopp, robert-dot-kopp-at-rutgers-dot-edu, Mon Aug 08 18:24:20 EDT 2016

defval('mergeZOSZOSTOGA',1);
defval('ZOSpullregion',12); % New York
defval('IFILES','');
defval('PARAMDIR','');
defval('SLRDIR','SLR_ALL');

%% Thermal expansion

disp('Thermal Expansion...');
for kk=1:length(scens)
    disp(['   ' scens{kk}]);
    [ThermExpMean(:,kk),ThermExpStd(:,kk),ThermExpYears,ThermExpN(:,kk),ZOSTOGA{kk},sZOSTOGA{kk},CWdrift,histGICrate,ZOSTOGAmodels{kk},ZOSTOGAyrs{kk}]=readZOSTOGA(scens{kk},[],[],datayears,1,{'zostoga','zosga'},IFILES,PARAMDIR);
end

%% Dynamics

if nargout>4

    disp('Dynamics...');

    for kk=1:length(scens)
	disp(['   ' scens{kk}]);

	if length(ZOSTOGAmodels)>0
            doZSGAmodels=ZOSTOGAmodels{kk}; doZSGA=sZOSTOGA{kk}; doZSGAyrs=ZOSTOGAyrs{kk};
	else
            doZSGAmodels=[]; doZSGA=[]; doZSGAyrs=[];
 end

 sub=1:length(targregions);

 [OceanDynMean(:,:,kk),OceanDynStd(:,:,kk),OceanDynYears,OceanDynRegions,OceanDynN(:,:,kk),OceanDynTECorr(:,:,kk)]=readZOS(scens{kk},targregions(sub),sitecoords(sub,:),19,datayears,0,doZSGAmodels,doZSGA,doZSGAyrs,mergeZOSZOSTOGA,fullfile(IFILES,SLRDIR));

    end
end
