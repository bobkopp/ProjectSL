function [quanttotlocrise, quantloccomponents,quantscalefactors,totlocalrisefocus,sampsregionfocus,colGIA,colOD,colGICtot,colAIStot,colIStot,colLItot,colLILStot,colLILSTEODtot,colLILSTEODGIAtot,lastseed] = ProjectLSL(scens,targregions,targregionnames,targyears,samps,seeds,regionsu,OceanDynRegions,OceanDynYears,OceanDynMean,OceanDynStd,OceanDynN,ThermExpYears, ThermExpMean, ThermExpStd,OceanDynTECorr,rateprojs,rateprojssd,mergeZOSZOSTOGA,fpsite,quantlevs,colGIC,colGIS,colAIS,colLS,colTE,focussites,GCMprobscale,maxDOF);

% Last updated by  Bob Kopp, robert-dot-kopp-at-rutgers-dot-edu, Jan 12 2014

colIS=union(colAIS,colGIS);
defval('GCMprobscale',.833);
defval('maxDOF',Inf);
defval('regionsu',targregions);
lastseed0=0;

for kk=1:length(scens)
	clear sampsregion;
	for jj=1:length(targregions)
		lastseed = lastseed0;
	
		% fingerprints
		sub=find(regionsu==targregions(jj));
		subOD = find(OceanDynRegions==targregions(jj));
		n0 = size(samps,2);	
		for ii=1:length(targyears)
			targyear=targyears(ii);
			lastseed=lastseed0;
			sampsregion(:,1:n0,ii) = bsxfun(@times,samps(:,:,ii,kk),fpsite(sub,:));

			%GIA
			colGIA = n0+1;
			GIAproj = rateprojs(sub)*(targyear-2000); GIAprojsd = rateprojssd(sub)*(targyear-2000);
			sampsregion(:,colGIA,ii) = GIAproj + [seeds(lastseed+1,:)'*GIAprojsd];
			lastseed = lastseed+1;

			% ocean dynamics
			if mergeZOSZOSTOGA
				colOD=colTE;
			else
				colOD = n0+2;
			end
			if length(subOD)>0
		
				sub1=find(OceanDynYears==targyear);
				sub2=find(ThermExpYears==targyear);
				
				condmean = OceanDynMean(sub1(1),subOD(1),kk) +  OceanDynStd(sub1(1),subOD(1),kk) * OceanDynTECorr(sub1(1),subOD(1),kk) * (samps(:,colTE,ii,kk)-ThermExpMean(sub2(1),kk))/ThermExpStd(sub2(1),kk);
				
				condstd = (norminv(.95)/norminv(GCMprobscale) * OceanDynStd(sub1,subOD(1),kk)) * sqrt(1 - OceanDynTECorr(sub1(1),subOD(1),kk)^2);
				
				
%				sampsregion(:,colOD,ii) = [norminv(.95)/norminv(GCMprobscale)*tinv(normcdf(seeds(lastseed+1,:)),min(OceanDynN(ii,subOD(1),kk)-1,maxDOF))'*OceanDynStd(sub1,subOD(1),kk)] + OceanDynMean(sub1,subOD(1),kk);

				sampsregion(:,colOD,ii) = [tinv(normcdf(seeds(lastseed+1,:)),min(OceanDynN(ii,subOD(1),kk)-1,maxDOF))'*condstd] + condmean;

			end

		end	
		
		totlocalrise = sum(sampsregion,2);
                subfocussite = find(focussites==targregions(jj));
		if length(subfocussite)>0
		     totlocalrisefocus{subfocussite(1),kk} = squeeze(totlocalrise);
             sampsregionfocus{subfocussite(1),kk} = sampsregion;
		end

		if colOD==colTE
			colTE2=[];
		else
			colTE2=colTE;
		end

		clear scalefactors;
		scalefactors(:,1,:) = sum(sampsregion,2)./sum(samps(:,:,:,kk),2);
		scalefactors(:,2,:) = sum(sampsregion,2)-sum(samps(:,:,:,kk),2);;
		scalefactors(:,3,:) = sum(sampsregion(:,[colGIC colIS colTE2 colOD colLS],:),2)./sum(samps(:,:,:,kk),2);
		scalefactors(:,4,:) = sum(sampsregion(:,[colGIC colIS colTE2 colOD],:),2)./sum(samps(:,[colGIC colIS colTE],:,kk),2);
		scalefactors(:,5,:) = sum(sampsregion(:,[colGIC colIS colTE2],:),2)./sum(samps(:,[colGIC colIS colTE],:,kk),2);
		scalefactors(:,6,:) = sum(sampsregion(:,[colGIC colIS],:),2)./sum(samps(:,[colGIC colIS ],:,kk),2);
		scalefactors(:,7,:) = sum(sampsregion(:,[colGIC],:),2)./sum(samps(:,[colGIC],:,kk),2);
		scalefactors(:,8,:) = sum(sampsregion(:,[colTE2 colOD],:),2)./sum(samps(:,[colTE],:,kk),2);
		scalefactors(:,9,:) = sum(sampsregion(:,[colTE2 colOD],:),2)-sum(samps(:,[colTE],:,kk),2);

		quantloccomponents(:,:,:,jj,kk) = quantile([sampsregion sum(sampsregion(:,colGIC,:),2) sum(sampsregion(:,colAIS,:),2) sum(sampsregion(:,colIS,:),2) sum(sampsregion(:,[colGIC colIS],:),2) sum(sampsregion(:,[colGIC colIS colLS],:),2) sum(sampsregion(:,[colGIC colIS colLS colTE2 colOD],:),2) sum(sampsregion(:,[colGIC colIS colLS colTE2 colOD colGIA],:),2)],quantlevs,1);
		quanttotlocrise(:,:,jj,kk) = squeeze(quantile(totlocalrise,quantlevs,1));
		
		quantscalefactors(:,:,:,jj,kk) = quantile(scalefactors,quantlevs,1);

		colGICtot = size(sampsregion,2)+1;
		colAIStot = size(sampsregion,2)+2;
		colIStot = size(sampsregion,2)+3;
		colLItot = size(sampsregion,2)+4;
		colLILStot = size(sampsregion,2)+5;
		colLILSTEODtot = size(sampsregion,2)+6;
		colLILSTEODGIAtot = size(sampsregion,2)+7;
		
		sub2100=find(targyears==2100);
		disp(sprintf(['Region %0.0f (' targregionnames{jj} ')- median in ' scens{kk} ' 2100 of %0.2f mm'],[targregions(jj) median(totlocalrise(:,:,sub2100))]));

	end
end
