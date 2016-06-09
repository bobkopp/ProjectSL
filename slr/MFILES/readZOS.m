function [OceanDynMean,OceanDynStd,OceanDynYears,OceanDynRegions,OceanDynN,OceanDynTECorr,ZOS,sZOS,modellist]=readZOS(scen,targregions,sitecoords,smoothwin,years,extrap,ZOSTOGAmodels,ZOSTOGA,ZOSTOGAyears,mergeZOSZOSTOGA,subdir)

% [OceanDynMean,OceanDynStd,OceanDynYears,OceanDynRegions,OceanDynN,ZOS,sZOS,modellist]=readZOS(scen,targregions,sitecoords,smoothwin,years)
%
% Last updated by  Bob Kopp, robert-dot-kopp-at-rutgers-dot-edu, Tue Feb 11 23:25:41 EST 2014

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

% read ZOS

clear ZOS;
clear sZOS;

Nyears=length(years);
Nregions=length(targregions);
Nmodels = 0;

pd=pwd;
cd(subdir)
cd(scen)

files=dir('*');
modellist={};
jj=1;
for ii=1:length(files)
	if (files(ii).isdir) && (files(ii).name(1)~='.')
		if exist([files(ii).name,'/zos'],'dir')
			curmodel=lower(files(ii).name);
			if length(ZOSTOGAmodels)==0
				Nmodels=Nmodels+1;
			elseif sum(strcmpi(curmodel,ZOSTOGAmodels))>0
				Nmodels = Nmodels+1;
			end
		end
	end
end

cd ..;

ZOS = zeros(Nyears,Nmodels,Nregions);

cd(scen);
files=dir('*');
modellist={};
jj=1;
for ii=1:length(files)
	if (files(ii).isdir) && (files(ii).name(1)~='.')
		if exist([files(ii).name,'/zos'],'dir')
			curmodel=lower(files(ii).name);
			s=find(strcmpi(curmodel,ZOSTOGAmodels));
			if length(s)>0
				s=s(1);
			else
				s=0;
			end
			if (length(ZOSTOGAmodels)==0) || (s>0)
				modellist{jj}=curmodel;
				disp([files(ii).name ' - zos']);
				cd(files(ii).name);
				cd zos
				files2=dir('*.txt');
				dat=importdata(files2(1).name,' ',2);			
				datids = str2num(dat.textdata{2}); dat=dat.data;
				ZOSyears=years;
				ZOSTOGAadj(:,jj)=zeros(length(ZOSyears),1);
				if (s>0)
					subnonan=find(~isnan(ZOSTOGA(:,s)));
					ZOSTOGAadj(:,jj)=interp1(ZOSTOGAyears(subnonan),ZOSTOGA(subnonan,s),ZOSyears);
				end
				for ll=1:length(targregions)
					colselect = find(datids == targregions(ll));
					if length(colselect)>0
						llA = colselect(1)+1;
						ZOS(:,jj,ll) = interp1(dat(:,1),dat(:,llA),ZOSyears);
						if mergeZOSZOSTOGA;	ZOS(:,jj,ll) = ZOS(:,jj,ll) + ZOSTOGAadj(:,jj); end;
					end
				end
				jj=jj+1; cd ../..
			end
		end
	end
end
cd(pd);
ZOS(find(ZOS<-999))=NaN;

sZOS = NaN*ZOS;
sZOSTOGAadj = NaN*ZOSTOGAadj;
for jj=1:size(ZOS,2);
%	disp([modellist{jj} ' - smoothing zos']);
	sub = find(~isnan(ZOS(:,jj,1)));
	for ll=1:size(ZOS,3)
		sZOS(sub,jj,ll) = smooth(ZOS(sub,jj,ll),smoothwin);
	end
	if length(ZOSTOGAmodels)>0
		sZOSTOGAadj(sub,jj) = smooth(ZOSTOGAadj(sub,jj),smoothwin);
	end
end


sZOS(:,:,:) = bsxfun(@minus,sZOS(:,:,:),sZOS(find(years==baseyear),:,:));
sZOSTOGAadj(:,:) = bsxfun(@minus,sZOSTOGAadj(:,:),sZOSTOGAadj(find(years==baseyear),:));


OceanDynMean=zeros(size(sZOS,1),size(ZOS,3));
OceanDynStd = OceanDynMean; OceanDynN = OceanDynMean;
OceanDynTECorr = OceanDynMean;

for ll=1:size(ZOS,3)
	if mod(ll,10)==0
		disp(['  Processing site ' num2str(ll) '/' num2str(size(ZOS,3))]);
	end
	sub1=find((years>2000).*(years<2100));
	sub=find(~isnan(sZOS(sub1(1),:,ll)));
        sub1=intersect(sub1,find(~isnan(sum(sZOS(:,sub,1),2))));
	sub=intersect(sub,find(std(sZOS(sub1,:,ll),[],1)>0));

        extremeness2=abs(sZOS(sub1(end),:,ll))/median(abs(sZOS(sub1(end),sub,ll)));
        sub=intersect(sub,find(extremeness2<10));

        extremeness=abs((sZOS(sub1(end),:,ll)-mean(sZOS(sub1(end),sub,ll)))/std(sZOS(sub1(end),sub,ll)));  % check for extreme outliers
        if std(sZOS(sub1(end),sub,ll))>.2
            sub=intersect(sub,find(extremeness<3));
        end
        
	% exclude MIROC and GISS at latitude > 55
	if abs(sitecoords(ll,1))>50 % at high latitudes, exclude miroc and giss
		sub=setdiff(sub,find(strncmp('miroc',modellist,length('miroc')))); 
		sub=setdiff(sub,find(strncmp('giss',modellist,length('giss')))); 
	end
	parfor jj=1:size(ZOS,1)
		sub2=intersect(sub,find(~isnan(sZOS(jj,:,ll))));
		OceanDynMean(jj,ll) = squeeze(mean(sZOS(jj,sub2,ll),2)*1000);
		OceanDynStd(jj,ll) = squeeze(std(sZOS(jj,sub2,ll),[],2)*1000);
		OceanDynN(jj,ll) = sum(~isnan(sZOS(jj,sub2,ll)),2);
		
		if (length(ZOSTOGAmodels)>0) && (length(sub2)>0)
			OceanDynTECorr(jj,ll) = squeeze(corr(sZOS(jj,sub2,ll)'*1000,sZOSTOGAadj(jj,sub2)'*1000));
		end

	end
	
	if extrap
		sub=find(OceanDynN(:,ll)>0);
		sub2=(sub(end)+1):size(OceanDynN,1);
		if length(sub2)>0
			slope1=(OceanDynMean(sub(end),ll)-OceanDynMean(sub(end)-extrap,ll))/extrap;
			slope2=(OceanDynStd(sub(end),ll)-OceanDynStd(sub(end)-extrap,ll))/extrap;
			OceanDynMean(sub2,ll)=OceanDynMean(sub(end),ll)+(sub2-sub(end))*slope1;
			OceanDynStd(sub2,ll)=OceanDynStd(sub(end),ll)+(sub2-sub(end))*slope2;
			if length(ZOSTOGAmodels)>0
				slope2=(OceanDynTECorr(sub(end),ll)-OceanDynTECorr(sub(end)-extrap,ll))/extrap;
				OceanDynTECorr(sub2,ll)=OceanDynTECorr(sub(end),ll)+(sub2-sub(end))*slope2;
			end

		end
	end
end

OceanDynYears = years;

OceanDynMean(end+1,:) = bsxfun(@plus,OceanDynMean(end,:),(OceanDynMean(end,:)-OceanDynMean(end-1,:)));
OceanDynStd(end+1,:) = bsxfun(@plus,OceanDynStd(end,:),(OceanDynStd(end,:)-OceanDynStd(end-1,:)));
if length(ZOSTOGAmodels)>0
	OceanDynTECorr(end+1,:) = bsxfun(@plus,OceanDynTECorr(end,:),(OceanDynTECorr(end,:)-OceanDynTECorr(end-1,:)));
end

OceanDynYears(end+1) = OceanDynYears(end) + (OceanDynYears(end)-OceanDynYears(end-1));
OceanDynN(end+1,:) = OceanDynN(end,:);

OceanDynRegions=targregions;
OceanDynTECorr=max(-1,min(1,OceanDynTECorr));
