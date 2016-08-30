
function [ZOS_final, modellist, scen, targsitecoords, years] = readZOSdabgrid(scen,targsitecoords,years,subdir)

% [ZOS,modellist]=dab_slrgridextract(scen,targsitecoords,years,subdir)
%
% Read ZOS files, saved in the format output by Dan Bader's post-processing
% of the CMIP5 archive.
%
% This routine expects files to be in the directory 'subdir' with subdirectories
% corresponding to different scenarios.
%
% INPUTS
% ------
% scen: string specifying scenario to be read (should be name of subdirectory)
% targsitecoords: [lat lon] matrix of sites to be loaded
% years: years to be retrieved
% subdir: path containing ZOS files
%
% OUTUTS
% ------
% ZOS: 3-D ZOS output array (years x models x sites)
% modellist: cell array of model names
%
% writted by Daniel Bader, dab2145-at-columbia-dot-edu, Tue Aug 30 13:24:00 EDT 2016 
% Last updated by Robert Kopp, robert-dot-kopp-at-rutgers-dot-edu, Tue Aug 30 18:32:05 EDT 2016

defval('scen','rcp85');
defval('years',1860:2099);
defval('targsitecoords',[40 286]);
defval('subdir','~/NASA/output');

angd = @(Lat0,Long0,lat,long) (180/pi)*(atan2( sqrt( (cosd(lat) .* sind(long-Long0)).^2 + (cosd(Lat0) .* sind(lat) - sind(Lat0) .* cosd(lat) .* cosd(long-Long0)).^2),(sind(Lat0) .* sind(lat) + cosd(Lat0) .* cosd(lat) .* cosd(long-Long0))));

pd=pwd;
cd(subdir)

%load in lat lon files (regridded domain) 
q=load('slrlat.mat');
lats=q.nlat(:);

q=load('slrlon.mat');
lons=q.nlon(:);

clear q;

%set variables 
Nyears = length(years); %
Nregions = size(targsitecoords,1); %number of points, assumes lat/lon comes in a N X 2 
Nmodels = 0; 

% move to directory with files 
cd(scen)
cd('zos') 

% count number of models 
files = dir('*.mat'); 
for ii = 1:length(files) 
    if files(ii).name(1)~='.' 
        Nmodels = Nmodels+1; 
    end 
end 

%output as years x models x points 
ZOS_final = NaN*ones(Nyears,Nmodels,Nregions); 

jj =1; 
modellist = {};

% identify nearest neighbors
disp('Finding nearest neighbors...');
goodsize=2*mean(diff(unique(lons)));
submatch=ones(size(targsitecoords,1),1)*NaN;
parfor kk = 1:size(targsitecoords,1)
    ad=angd(targsitecoords(kk,1),targsitecoords(kk,2),lats,lons);
    [m,mi]=min(ad);
    if m<goodsize
        submatch(kk)=mi;
    end
end

% select years 
years1 = years - 1859; %all files start 1860 and go to 2099

% load in files 
for ii = 1:length(files) 
    if files(ii).name(1)~='.' 
        
        % store model name 
        modelname = (files(ii).name);
        spot = min(find(modelname == '_'))-1; 
        curmodel = modelname(1:spot);
        modellist{jj} = lower(curmodel); 
        disp(curmodel) 
        
        % load in regridded model file 
        filname = strcat(curmodel,'_',scen,'_zos'); 
        dozos=load(filname); 
        
        % select years
        years_file=1859+[1:size(dozos.zosv2,3)];
        [doyears,doyearsi,doyearsj]=intersect(years,years_file);

        % setup file 
        zos_file = reshape(dozos.zosv2,[],1,length(years_file)); % changes to year x model x lat-lon
        zos_file=permute(zos_file,[3 2 1]);
        clear dozos;
        
        subgood=find(~isnan(submatch));
        ZOS_final(doyearsi,jj,subgood) = zos_file(doyearsj,1,submatch(subgood));
        jj = jj+1; 
    end
end

cd(pd);