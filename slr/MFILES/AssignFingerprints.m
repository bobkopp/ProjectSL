function [fpsite,fp,fpname,lo,la]=AssignFingerprints(fpmapperids,fpmaps,sitecoords,Nuniform,subdir);

% Last updated by  Bob Kopp, robert-dot-kopp-at-rutgers-dot-edu, Sun Dec 29 19:47:49 EST 2013

defval('Nuniform',2);
defval('subdir','');

% fingerprints
[fp,fpname,lo,la] = readFingerprint(subdir);

clear fpsite;
for i=1:length(fpmapperids)
	fptouse=strmatch(fpmaps{i},fpname,'exact');
	fpsite(:,i) = 1000*interp2(lo,la,fp(:,:,fptouse),mod(sitecoords(:,2),360),sitecoords(:,1));
end

icesheets={'gis','wais','eais'};
for i=1:length(icesheets)
	fptouse=strmatch(icesheets{i},fpname,'exact');
	fpsite(:,length(fpmapperids)+i) = 1000*interp2(lo,la,fp(:,:,fptouse),mod(sitecoords(:,2),360),sitecoords(:,1));
end

for i=1:Nuniform
	fpsite(:,end+1) = 1;
end

sub=find(sitecoords(:,1)>1e3); fpsite(sub,:) = 1;
