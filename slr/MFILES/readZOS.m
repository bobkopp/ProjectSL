function [OceanDynMean,OceanDynStd,OceanDynYears,OceanDynRegions,OceanDynN,OceanDynTECorr,ZOS,sZOS,modellist]=readZOS(scen,targregions,sitecoords,smoothwin,years,extrap,ZOSTOGAmodels,ZOSTOGA,ZOSTOGAyears,mergeZOSZOSTOGA,subdir,filemode)

% [OceanDynMean,OceanDynStd,OceanDynYears,OceanDynRegions,OceanDynN,OceanDynTECorr,ZOS,sZOS,modellist]=readZOS(scen,targregions,sitecoords,smoothwin,years,extrap,ZOSTOGAmodels,ZOSTOGA,ZOSTOGAyears,mergeZOSZOSTOGA,subdir)
%
% Read and smooth ZOS.
%
% INPUTS
% ------
% scen: string specifying scenario to be read (should be name of subdirectory)
% targregions: PSMSL IDs of sites to be retrieved
% sitecoords: [lat long] of sites
% smoothwin: window used for smoothing (default = 19 years)
% years: years to be retrieved
% extrap: number of years used for extrapolating beyond end of data (default = 0, meaning no extrapolation)
% ZOSTOGAmodels: names of ZOSTOGA models
% ZOSTOGA: ZOSTOGA values
% ZOSTOGAyears: years for ZOSTOGA
% mergeZOSZOSTOGA: add ZOSTOGA to ZOS? (default = 0, meaning no)
% subdir: path containing ZOS files
% filemode: use Kopp et al. (2014)-style text files (0) or GISS-processed regridded mat files (1)
%
% OUTPUTS
% -------
% OceanDynMean: Mean of smoothed ZOS (rows = years, columns = sites)
% OceanDynStd: Standard deviation of smoothed ZOS
% OceanDynYears: Years for rows of OceanDynMean and OceanDynStd
% OceanDynRegions: Region IDs for columns of OceanDynMean and OceanDynStd
% OceanDynN: Number of models contributing to each element of OceanDynMean
% OceanDynTECorr: correlation between smoothed ZOS and ZOSTOGA at each site
% ZOS: raw ZOS values (years x models x sites)
% sZOS: smoothed ZOS
% modellist: models corresponding to columns of ZOS
%
% Last updated by Robert Kopp, robert-dot-kopp-at-rutgers-dot-edu, Thu Sep 29 13:29:08 EDT 2016

defval('scen','rcp85');
defval('smoothwin',19);
defval('years',1860:2099);
defval('targregions',12);
defval('subdir','IFILES/slr/SLR_ALL');
defval('sitecoords',zeros(length(targregions),2))
defval('extrap',50);
defval('baseyear',2000);
defval('ZOSTOGA',[]);
defval('ZOSTOGAmodels',{});
defval('ZOSTOGAyears',[]);
defval('mergeZOSZOSTOGA',0);
defval('filemode',0);

if filemode==0
    [ZOSraw,ZOSmodels,scen,targregions,years]=readZOStgtxt(scen,targregions,years,subdir);
elseif filemode==1
    disp(['filemode 1 in ' subdir]);
    [ZOSraw,ZOSmodels,scen,~,years]=readZOSdabgrid(scen,sitecoords,years,subdir);
end

[OceanDynMean,OceanDynStd,OceanDynYears,OceanDynN,OceanDynTECorr,ZOS,sZOS,modellist]=processZOS(sitecoords,years,ZOSraw,ZOSmodels,smoothwin,extrap,ZOSTOGAmodels,ZOSTOGA,ZOSTOGAyears,mergeZOSZOSTOGA);
OceanDynRegions=targregions;
