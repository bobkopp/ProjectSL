% Last updated by Robert Kopp, robert-dot-kopp-at-rutgers-dot-edu, 2021-06-28 17:12:55 -0400

mf=mfilename('fullpath'); mfsl=strfind(mf,'/'); mypath=mf(1:mfsl(end)-1);
addpath(mypath,fullfile(mypath,'../lib/MFILES'),fullfile(mypath,'MFILES'),fullfile(mypath,'MFILES/scripts'));

scens={'rcp85','rcp60','rcp45','rcp26'};

Nsamps=10000; Nseeds=500;

datayears=1861:2299;
targyears = [2010:10:2300];

savefile='~/tmp/SLRProjections140523';
savefile2='~/tmp/SLRProjectionsquant140523';
outputdir='~/SLRprojections140523';
PARAMDIR=fullfile(mypath,'PARAMS/');
IFILES=fullfile(mypath,'../IFILES/slr');

% ICE SHEETS
% use Shepherd estimates for last decade (2000-2011), pooling APIS and WAIS
LastDecadeGt = [-211 -85-29 26 ];

% ice sheet distributions from Bamber & Aspinall (2013)

ratesmatrix2100 = [.8 1.0 1.2 2.4 5.8 ; .2 3.0 .3 1.5 11.8 ; -1.9 2.8 -1.5 0.2 10.2];
BAcorrIS{1} = eye(3);
BAcorrIS{2} =  [ 1 .7 -.2 ; .7 1 -.2 ; -.2 -.2 1 ];

% % ice sheet projections from AR5 -- these are 2081-2100 (i.e., 2090); third dimension corresponds to scens
% ARIS2090(:,:,4) = [0.06 0.04 0.10 ; 0.05 -0.03 0.14]*1000;
% ARIS2090(:,:,3) = [0.08 0.04 0.13 ; 0.05 -0.04 0.13]*1000;
% ARIS2090(:,:,2) = [0.08 0.04 0.13 ; 0.05 -0.04 0.13]*1000;
% ARIS2090(:,:,1) = [0.12 0.07 0.21 ; 0.04 -0.06 0.12]*1000;

% ice sheet projections from SROCC -- these are 2081-2100 (i.e., 2090); third dimension corresponds to scens
ARIS2090(:,:,4) = [0.06 0.04 0.10 ; 0.04 0.01 0.10]*1000;
ARIS2090(:,:,3) = [0.08 0.04 0.13 ; 0.05 0.01 0.13]*1000;
ARIS2090(:,:,2) = [0.08 0.04 0.13 ; 0.05 0.01 0.13]*1000;
ARIS2090(:,:,1) = [0.12 0.07 0.21 ; 0.10 0.02 0.23]*1000;


ARcorrIS{1} = eye(2);
ARcorrIS{2} = [1 .5 ; .5 1];

ISmode = 'hy'; % or 'AR' or 'BA'

% ZOS/ZOSTOGA
mergeZOSZOSTOGA=1;

% Land water storage
LWSpath.pophistfile=fullfile(IFILES,'LWS/UNWPP2012 population historical.csv');
LWSpath.damfile=fullfile(IFILES,'LWS/Chao2008 groundwater impoundment.csv');
LWSpath.GWDfiles={fullfile(IFILES,'LWS/Konikow2011 GWD.csv'),fullfile(IFILES,'LWS/Wada2012 GWD.csv')};
LWSpath.popscenfile=fullfile(IFILES,'LWS/UNWPP2012 population projections.csv');
LWSpath.dotriangular=0;

LWSpathwPokhrel=LWSpath;
LWSpathwPokhrel.GWDfiles={fullfile(IFILES,'LWS/Konikow2011 GWD.csv'),fullfile(IFILES,'LWS/Wada2012 GWD.csv'),fullfile(IFILES,'LWS/Pokhrel2012 GWD.csv')};
LWSpathwPokhrel.dotriangular=1;
