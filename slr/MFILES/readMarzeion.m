function [GIC,GICse,GICyrs,GICmodel,fpmapperids,fpmaps,GICnames]=readMarzeion(scen,IFILES,PARAMDIR,discardAntarctic)

% [GIC,GICse,GICyrs,GICmodel,fpmapperids,fpmaps,GICnames]=readMarzeion(scen,IFILES,PARAMDIR,discardAntarctic)
%
% Last updated by  Bob Kopp, robert-dot-kopp-at-rutgers-dot-edu, Tue Feb 11 23:16:05 EST 2014

defval('IFILES','IFILES/slr/');
defval('PARAMDIR','PARAMS/');
defval('subdir',fullfile(IFILES,'Marzeion2012supplement'));
defval('baseyear',2000);
defval('fpmapper','fingerprint_region_map.csv');
defval('discardAntarctic',1);

pd=pwd;
subKeep=[];
if (nargout>4)||(discardAntarctic)
	tmp = importdata(fullfile(PARAMDIR,fpmapper));
	fpmapperids=tmp.data;
	fpmaps=tmp.textdata(2:end,3);
	GICnames=tmp.textdata(2:end,1);
	if discardAntarctic
		subKeep = setdiff(1:length(GICnames),find(strcmpi('AntarcticGIC',GICnames)));
		fpmapperids=fpmapperids(subKeep);
		fp=fpmaps(subKeep);
		GICnames=GICnames(subKeep);
	end
end

cd(subdir);

files=dir(['dV_regional_accum*' scen '*.txt']);

for i=1:length(files)
	tmp=importdata(files(i).name);
	sub=find(tmp(:,1)==baseyear);
	tmp(:,2:end) = bsxfun(@minus,tmp(:,2:end),tmp(sub,2:end));
	ntmp = size(tmp,2)-1;
	if length(subKeep)==0
		subKeep=1:(ntmp/2);
	end
	GIC(:,:,i)=tmp(:,1+subKeep);
	GICse(:,:,i)=tmp(:,1+(ntmp/2)+subKeep);
	GICyrs(:,i) = tmp(:,1);
	s=regexp(files(i).name,'accum_(.+)_rcp','tokens'); GICmodel{i}=lower(s{1}{1});
end


cd(pd);