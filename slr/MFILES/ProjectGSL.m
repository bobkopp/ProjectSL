function [samps, colGIC, colIS, colGIS, colAIS, colLS, colTE, lastseed] = ProjectGSL(scens, targyears, projGICyrs, projGIC,projGICse, seeds, ISLastDecade, ISaccelsamps, ARISaccelsamps, ISmode, ThermExpYears, ThermExpMean,  ThermExpStd, ThermExpN, LWSpath,quantlevs,GCMprobscale,maxDOF);

% Last updated by Robert Kopp, robert-dot-kopp-at-rutgers-dot-edu, Fri May 23 08:11:51 EDT 2014

lastseed=0;
defval('GCMprobscale',.833);
defval('maxDOF',Inf);
defval('LWSpath',[]);

if isstruct(LWSpath)
    disp('Land water storage');
    [LWSsamps]=SampleLandWaterStorage(normcdf(seeds(lastseed+[1:4],:)),targyears,LWSpath.pophistfile,LWSpath.popscenfile,LWSpath.damfile,LWSpath.GWDfiles,LWSpath.dotriangular);
    lastseed=lastseed+4;
end

lastseed0=lastseed;
for kk=1:length(scens)
	disp(scens{kk});

        
            
	for ii=1:length(targyears)

		targyear = targyears(ii);
		disp(targyear);
		counter=0;
	
		% GLACIERS AND ICE CAPS
		% glacier and ice cap distribution from Marzeion et al. 2012

		GICyears=projGICyrs{kk}(:,1);

		if targyear<=2100
			sub = find((GICyears>=max(2000,targyear-4)).*(GICyears<=min(2100,targyear+4)));
		else
			sub = find((GICyears>=max(2000,targyear-4)).*(GICyears<=min(2300,targyear+4)));
		end		
		targ=mean(projGIC{kk}(sub,:,:),1);	tmp = squeeze(targ)';
		targSE=mean(projGICse{kk}(sub,:,:),1); tmpSE = squeeze(targSE)'; 
		sub=find(~isnan(sum(tmp,2)));
		meanGIC = mean(tmp(sub,:)); covGIC = cov(tmp(sub,:));
		tmpSE=mean(tmpSE(sub,:),1);
		covGIC = covGIC + diag(tmpSE).^2;
		NGIC=length(sub);

		GlacSD = sqrt(diag(covGIC));
		T = cholcov(covGIC);
		lastseed=lastseed0;
		
		if NGIC>0
			samps(:,counter+1:counter+length(meanGIC),ii,kk)=bsxfun(@plus,[tinv(normcdf(seeds(lastseed+[1:size(T,1)],:)),min(maxDOF,NGIC-1))'*T],meanGIC);
		else
			samps(:,counter+1:counter+length(meanGIC),ii,kk)=NaN;
		end
		lastseed=lastseed+length(meanGIC);
		counter=counter+length(meanGIC);
		colGIC = 1:length(meanGIC);

		% ICE SHEETS

		ISsampsBA = bsxfun(@plus,ISLastDecade * (targyear-2000),.5*ISaccelsamps*(targyear-2011).^2*(targyear>2011));
		ISsampsBAAIS = sum(ISsampsBA(:,2:3),2);
		BAq = quantile([ISsampsBA ISsampsBAAIS],[.167 .5 .833]);
		
		ARISLastDecade=[ISLastDecade(1) ISLastDecade(2)+ISLastDecade(3)];
		ISsampsAR = bsxfun(@plus,ARISLastDecade * (targyear-2000),.5*ARISaccelsamps(:,:,kk)*(targyear-2011).^2*(targyear>2011));
		ARq = quantile(ISsampsAR,[.167 .5 .833]);

		qslope=diff(ARq)./(diff(BAq(:,[1 4]))+eps);

		ISsampsHybrid = bsxfun(@minus,[ISsampsBA(:,1) ISsampsBAAIS],BAq(2,[1 4]));
		ISsampsHybrid = bsxfun(@plus, (ISsampsHybrid<0) .* bsxfun(@times,ISsampsHybrid,qslope(1,:)) + (ISsampsHybrid>0) .*  bsxfun(@times,ISsampsHybrid,qslope(2,:)), ARq(2,:));

		WAISsampsHybrid = ISsampsBA(:,2) - BAq(2,2);
		WAISsampsHybrid = (WAISsampsHybrid<0) .* (WAISsampsHybrid*qslope(1,2)) + (WAISsampsHybrid>0) .*  (WAISsampsHybrid*qslope(2,2)) + BAq(2,2) * ARq(2,2)/BAq(2,4);
		WAISsampsHybrid=min(WAISsampsHybrid,5000);
		
		ISsampsHybrid=[ISsampsHybrid(:,1) WAISsampsHybrid ISsampsHybrid(:,2)-WAISsampsHybrid];

		ISsampsHybrid(:,1) = min(ISsampsHybrid(:,1),7000);
		ISsampsHybrid(:,2) = min(ISsampsHybrid(:,2),5000);

		%samps(:,counter+1:counter+3,ii)=bsxfun(@plus,ISLastDecade * (targyear-2000),.5*ISaccelsamps*(targyear-2011).^2*(targyear>2011));
	
		if strcmpi(ISmode,'AR')
			samps(:,counter+1:counter+3,ii,kk)=[ISsampsAR(:,1) WAISsampsHybrid ISsampsAR(:,2)-WAISsampsHybrid];
		elseif strcmpi(ISmode,'BA')
			samps(:,counter+1:counter+3,ii,kk)=ISsampsBA;
		else
			samps(:,counter+1:counter+3,ii,kk)=ISsampsHybrid;
		end
		counter = counter+3;
		
		colIS = colGIC(end)+[1:3];
		colGIS = colIS(1); colAIS = colIS(2:3);

                % land water storage
                if ~isstruct(LWSpath)
                    % land water storage from AR5
                    meanLS = 40*(targyear-2000)/95;
                    stdLS = 90/(2*norminv(.833))*(targyear-2000)/95;
                    samps(:,counter+1,ii,kk) = [seeds(lastseed+1,:)'*stdLS] + meanLS;
                    lastseed = lastseed+1;
                else
                    samps(:,counter+1,ii,kk) = LWSsamps(ii,:)';
                end
                
		colLS=counter+1;
		counter = counter+1;

		% global thermal expansion
		
		sub1=find(ThermExpYears==targyear);
		samps(:,counter+1,ii,kk) = [norminv(.95)/norminv(GCMprobscale)*tinv(normcdf(seeds(lastseed+1,:)),min(maxDOF,ThermExpN(ii,kk)-1))'*ThermExpStd(sub1,kk)] + ThermExpMean(sub1,kk); % scale per AR5
			lastseed = lastseed+1;
		colTE = counter+1;
		counter=counter+1;


	end
end
