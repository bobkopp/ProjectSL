function [ZOS,modellist,scen,targregions,years]=readZOStgtxt(scen,targregions,years,subdir)

% [ZOS,modellist]=readZOStgtxt(scen,targregions,years,subdir)
%
% Read ZOS files, saved in the format used for Kopp et al. (2014), in which
% they have been processed into text files with rows corresponding to years and 
% columns to sites, labeled by PSMSL IDs.
%
% This routine expects files to be in the directory 'subdir' with subdirectories
% corresponding to different scenarios.
%
% INPUTS
% ------
% scen: string specifying scenario to be read (should be name of subdirectory)
% targregions: PSMSL IDs of sites to be retrieved
% years: years to be retrieved
% subdir: path containing ZOS files
%
% OUTUTS
% ------
% ZOS: 3-D ZOS output array (years x models x sites)
% modellist: cell array of model names
%
% Last updated by Robert Kopp, robert-dot-kopp-at-rutgers-dot-edu, Mon Aug 08 18:42:49 EDT 2016

defval('scen','rcp85');
defval('years',1860:2099);
defval('targregions',12);
defval('subdir','IFILES/slr/SLR_ALL');

Nyears=length(years);
Nregions=length(targregions);
Nmodels = 0;

pd=pwd;
cd(subdir)

% count number of models
cd(scen)
files=dir('*');
for ii=1:length(files)
    if (files(ii).isdir) && (files(ii).name(1)~='.')
        if exist([files(ii).name,'/zos'],'dir')
            curmodel=lower(files(ii).name);
            Nmodels=Nmodels+1;
        end
    end
end
cd ..;

% define matrix
ZOS = zeros(Nyears,Nmodels,Nregions);

% read files
cd(scen);
files=dir('*');
modellist={};
jj=1;
for ii=1:length(files)
    if (files(ii).isdir) && (files(ii).name(1)~='.')
        if exist([files(ii).name,'/zos'],'dir')
            curmodel=lower(files(ii).name);
            modellist{jj}=curmodel;
            disp([files(ii).name ' - zos']);
            cd(files(ii).name);
            cd zos
            files2=dir('*.txt');
            dat=importdata(files2(1).name,' ',2);			
            datids = str2num(dat.textdata{2}); dat=dat.data;
            ZOSyears=years;
            for ll=1:length(targregions)
                colselect = find(datids == targregions(ll));
                if length(colselect)>0
                    llA = colselect(1)+1;
                    ZOS(:,jj,ll) = interp1(dat(:,1),dat(:,llA),ZOSyears);
                end
            end
            jj=jj+1; cd ../..
        end
    end
end
cd(pd);
ZOS(find(ZOS<-999))=NaN;
