% outputSLRProjections_table_LSLfactors
% Last updated by  Bob Kopp, robert-dot-kopp-at-rutgers-dot-edu, Sun Jan 26 12:22:52 EST 2014


% local factors

factorfiles={'full','goodbkgdrates','subset'};
factorsubsets={targregions,targregions(find(abs(rateprojs)>2*rateprojssd)),union(focussites,[235 12 363 520 526 440 333 1034 1366 269 43 860 1099])};

for nnn=1:length(factorfiles)

	icesheets={'GIS','WAIS','EAIS'};

	fid=fopen(['LSLfactors_' factorfiles{nnn} '.tsv'],'w');
	fprintf(fid,'Site\tID\tLat\tLong');
	fprintf(fid,['\tBkgd (mm/y)\t+/-2s\tMedian ' scens{1} ' dynamics in 2100 (mm)\t17th\t83rd\t5th\t95th' '\tMedian climatic ' scens{1} ' scaler in 2100\t17th\t83rd\t5th\t95th']);
	for ii=1:length(icesheets)
		fprintf(fid,['\tFP: ',upper(icesheets{ii})]);
	end
	fprintf(fid,'\tFP: Median GIC');
	for ii=1:length(GICnames)
		fprintf(fid,['\tFP: ',GICnames{ii}]);
	end

	u=sort(unique(targcoastlines));
	for mmm=1:length(u)
		subcl=find(targcoastlines==u(mmm));
		subcl=intersect(subcl,find(ismember(targregions,factorsubsets{nnn})));
		if length(subcl)>0
			fprintf(fid,'\nCOASTLINE %0.0f',u(mmm));
			[s,si]=sort(targregions(subcl));
			for jjj=s(:)';
				jj=find(targregions==jjj);
				subOD = find(OceanDynRegions==targregions(jj));
				subODyear = find(OceanDynYears==2100);
				subyear = find(targyears==2100);
				fprintf(fid,['\n' targregionnames{jj}]);
				fprintf(fid,'\t%0.0f',targregions(jj));
				fprintf(fid,'\t%0.2f',sitecoords(jj,:));
				fprintf(fid,'\t%0.2f\t%0.2f',[rateprojs(jj) rateprojssd(jj)*2]);

				subql=find(quantlevs==0.5);
				subql(2)=find(quantlevs==0.167);
				subql(3)=find(quantlevs==0.833);
				subql(4)=find(quantlevs==0.05);
				subql(5)=find(quantlevs==0.95);

				fprintf(fid,'\t%0.2f',quantlocscalefactors(subql,9,subyear,jj,1));
				fprintf(fid,'\t%0.2f',quantlocscalefactors(subql,4,subyear,jj,1));
			
				fprintf(fid,'\t%0.2f',fpsite(jj,length(GICnames)+[1:3]));
				subql=find(quantlevs==0.5);
				medianfp = quantloccomponents(subql,27,subyear,jj,1)/quantcomponents(subql,27,subyear,1);
				fprintf(fid,'\t%0.2f',medianfp);
				fprintf(fid,'\t%0.2f',fpsite(jj,1:length(GICnames)));
			end
		end
	end
	fclose(fid);
end
