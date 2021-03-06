function [X1,Y,dY,regions,regionsu,sitenames,sitecoords,sitelen,sitecoastline]=ReadPSMSLData(cl1,cl2,minlen,psmsldir,gslfile,addlsites,exclusions)

% [X1,Y,dY,regions,regionsu,sitenames,sitecoords,sitelen,sitecoastline]=ReadPSMSLData([cl1],[cl2],[minlen],[psmsldir],[gslfile],[addlsites],[exclusions])
%
% Last updated by  Bob Kopp, robert-dot-kopp-at-rutgers-dot-edu, 2021-03-06 15:28:45 -0500

defval('psmsldir','IFILES/slr/rlr_annual');
defval('gslfile','IFILES/slr/CSIRO_Recons_gmsl_yr_2011.csv')
defval('minlen',30);
defval('baserslerror',3);
defval('cl1',0);
defval('cl2',1000);
defval('addlsites',[]);
defval('exclusions',[]);

% first let's read in data

addpath(psmsldir);
data=readAnnual(psmsldir);

X1 =[]; regions=[]; Y=[]; dY=[];
regionsu=[];
sitenames={};
sitecoords=[];
sitecoastline=[];

% read in global data

if exist(gslfile)
	sitenames{1} = 'CSIRO GSL';
	fid=fopen(gslfile);
	C=textscan(fid,'%f%f%f','Delimiter',',','Headerlines',1);
	fclose(fid);
	
	years=[C{1}(:)];
	rsl = [C{2}(:)];
	rslunc = [C{3}(:)]/2;

	[m,mi] = min(abs(years-2005));
	rsl = rsl - rsl(mi);
	

	X1 = [1e6 * ones(length(years),2) years];
	regions = zeros(size(years));
	Y = [rsl];
	dY = [rslunc];
	sitelen = length(rsl);

	regionsu=0;
	sitenames={'CSIRO GSL'};
	sitecoords=[1e6 1e6];
	sitecoastline=[0];
else
	disp([gslfile ' does not exist -- not loading GSL estimate']);
end


sub=find(([data.coastline]>=cl1).*([data.coastline]<=cl2));
sub=union(sub,find(ismember([data.id],addlsites)));
sub=setdiff(sub,find(ismember([data.id],exclusions)));

tgsitenames={};
tgsitecoords=[];
for i=1:length(sub)
	tgsitenames{end+1} = data(sub(i)).name;
	tgsitecoords(end+1,:) = [data(sub(i)).latitude data(sub(i)).longitude];

    dataflag=data(sub(i)).dataflag;
    subgood=find(~dataflag);
	years=data(sub(i)).year(subgood);
	rsl=data(sub(i)).height(subgood);
	sub2 = find(~isnan(rsl));
 	years=years(sub2); rsl=rsl(sub2);
	
	[m,mi] = min(abs(years-2005));
	baselinesub=[max(1,(mi-5)):min((mi+5),length(years))];
	rsl = rsl - mean(rsl(baselinesub));
	
	rslunc = baserslerror*ones(size(rsl));
	
	if length(rsl)>minlen
		
		X1 = [X1 ; repmat(tgsitecoords(i,:),length(years),1) years];
		regions = [regions; data(sub(i)).id*ones(size(years))];
		regionsu = [regionsu ; data(sub(i)).id];
		Y = [Y ; rsl];
		dY = [dY ; rslunc];
		sitelen = [sitelen ; length(rsl)];
		sitenames={sitenames{:},tgsitenames{end}};
		sitecoords=[sitecoords ; tgsitecoords(end,:)];
		sitecoastline=[sitecoastline ; data(sub(i)).coastline];
	end
end

for i=1:length(sitenames)
	j=strfind(sitenames{i},' (');
	if length(j)>0
		sitenames{i}=sitenames{i}(1:(j(1)-1));
	end
	j=strfind(sitenames{i},',');
	if length(j)>0
		sitenames{i}=sitenames{i}(1:(j(1)-1));
	end
	j=strfind(sitenames{i},' / ');
	if length(j)>0
		sitenames{i}=sitenames{i}(1:(j(1)-1));
	end
end