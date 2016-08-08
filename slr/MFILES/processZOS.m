function [OceanDynMean,OceanDynStd,OceanDynYears,OceanDynN,OceanDynTECorr,ZOS,sZOS,modellist]=processZOS(sitecoords,years,ZOSraw,ZOSmodels,smoothwin,extrap,ZOSTOGAmodels,ZOSTOGA,ZOSTOGAyears,mergeZOSZOSTOGA)

% [OceanDynMean,OceanDynStd,OceanDynYears,OceanDynRegions,OceanDynN,OceanDynTECorr,ZOS,sZOS,modellist]=processZOS(sitecoords,years,ZOSraw,ZOSmodels,smoothwin,extrap,ZOSTOGAmodels,ZOSTOGA,ZOSTOGAyears,mergeZOSZOSTOGA)
%
% Smooth ZOS and calculate mean, standard deviation, and correlation with ZOSTOGA.
%
% INPUTS
% ------
% sitecoords: [lat long] of sites
% years: years of projections
% ZOSraw: raw ZOS values
% ZOSmodels: cell array of model names
% smoothwin: window used for smoothing (default = 19 years)
% extrap: number of years used for extrapolating beyond end of data (default = 0, meaning no extrapolation)
% ZOSTOGAmodels: names of ZOSTOGA models
% ZOSTOGA: ZOSTOGA values
% ZOSTOGAyears: years for ZOSTOGA
% mergeZOSZOSTOGA: add ZOSTOGA to ZOS? (default = 0, meaning no)
%
% OUTPUTS
% -------
% OceanDynMean: Mean of smoothed ZOS (rows = years, columns = sites)
% OceanDynStd: Standard deviation of smoothed ZOS
% OceanDynYears: Years for rows of OceanDynMean and OceanDynStd
% OceanDynN: Number of models contributing to each element of OceanDynMean
% OceanDynTECorr: correlation between smoothed ZOS and ZOSTOGA at each site
% ZOS: raw ZOS values (years x models x sites)
% sZOS: smoothed ZOS
% modellist: models corresponding to columns of ZOS
%
% Last updated by Robert Kopp, robert-dot-kopp-at-rutgers-dot-edu, Mon Aug 08 18:30:56 EDT 2016

defval('smoothwin',19); % smoothing window (years)
defval('sitecoords',zeros(size(ZOSraw,3),2)) % site coordinates (used for error checking)
defval('extrap',0); % number of years used for extrapolation 
defval('baseyear',2000); % base year for reference time period
defval('ZOSTOGA',[]); % values of ZOSTOGA (year x model)
defval('ZOSTOGAmodels',{}); % names of ZOSTOGA models
defval('ZOSTOGAyears',[]); % ZOSTOGA years
defval('mergeZOSZOSTOGA',0); % add ZOSTOGA to ZOS?

% read ZOS

modellist={};
jj=1;

if length(ZOSTOGAmodels)==0
    ZOS=ZOSraw;
    modellist=ZOSmodellist;  
else
    ll=1;
    for jj=1:length(ZOSTOGAmodels)
        curmodel=ZOSTOGAmodels{jj};
        s=find(strcmpi(curmodel,ZOSmodels));
        if (s>0)
            modellist{ll}=curmodel;
            %disp([modellist{ll} ' - zos']);
            ZOSTOGAadj(:,ll)=zeros(length(years),1);
            subnonan=find(~isnan(ZOSTOGA(:,jj)));
            ZOSTOGAadj(:,ll)=interp1(ZOSTOGAyears(subnonan),ZOSTOGA(subnonan,jj),years);
            if mergeZOSZOSTOGA
                ZOS(:,ll,:) = bsxfun(@plus,ZOSraw(:,s(1),:),ZOSTOGAadj(:,ll));    
            else
                ZOS(:,ll,:)=ZOSraw(:,s(1),:);
            end               
            ll=ll+1;
        end
    end
end

sZOS = NaN*ZOS;
sZOSTOGAadj = NaN*ZOSTOGAadj;
for jj=1:size(ZOS,2);
	for ll=1:size(ZOS,3)
            sub = find(~isnan(ZOS(:,jj,ll)));
            sZOS(sub,jj,ll) = smooth(ZOS(sub,jj,ll),smoothwin);
	end
	if length(ZOSTOGAmodels)>0
            sub = find(~isnan(ZOSTOGAadj(:,jj)));
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

OceanDynTECorr=max(-1,min(1,OceanDynTECorr));
