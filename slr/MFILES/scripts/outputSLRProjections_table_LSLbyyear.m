% outputSLRProjections_table_LSLbyyear
% Last updated by Robert Kopp, robert-dot-kopp-at-rutgers-dot-edu, Fri May 23 14:50:49 EDT 2014

%%%

subscens=[1 3 4];
defval('targyears2',[2030 2050 2100 2150 2200]);
quantsub=[.5 .167 .833 .05 .95 .005 .995 .999];
roundafter=2100;

% Table: Local projections

[jk,ia,ib]=intersect(quantlevs,quantsub); [jk,ic]=sort(ib);
ia=ia(ic);

fid=fopen('LSLproj.tsv','w');

fprintf(fid,'cm ');
for kk=subscens
	fprintf(fid,[' \t ' scens{kk}]);
	if kk~=subscens(end)
		fprintf(fid,'\t \t \t \t');
	end
end
fprintf(fid,'\n ');

for kk=subscens
	fprintf(fid,['\t 50 \t 17--83 \t 5--95 \t 0.5--99.5 \t 99.9 ' ]);
end

sitesub=[12 299 188 161 10 405 155 7 78 438 134 499];
for jjj=1:length(sitesub)
	jj=find(targregions==sitesub(jjj));
	fprintf(fid,['\n' targregionnames{jj} ' (Bkgd: %0.2f +/- %0.2f mm/y)'], [rateprojs(jj) rateprojssd(jj)*2]);	
	for mm=1:length(targyears2)
		subyear = find(targyears==targyears2(mm));
		fprintf(fid,'\n%4.0f',targyears2(mm)); 
		for kk=subscens
			fprintf(fid,' \t ');
                        if targyears2(mm)>roundafter
                            fprintf(fid,'%0.0f \t',10*round(quanttotlocrise(ia(1),subyear,jj,kk)/10/10));
                            fprintf(fid,'%0.0f--%0.0f \t ',10*round(quanttotlocrise(ia(2:7),subyear,jj,kk)/10/10));
                            fprintf(fid,'<%0.0f',10*round(quanttotlocrise(ia(8),subyear,jj,kk)/10/10));
                        else
                            fprintf(fid,'%0.0f \t',quanttotlocrise(ia(1),subyear,jj,kk)/10);
                            fprintf(fid,'%0.0f--%0.0f \t ',quanttotlocrise(ia(2:7),subyear,jj,kk)/10);
                            fprintf(fid,'<%0.0f',5*round(quanttotlocrise(ia(8),subyear,jj,kk)/10/5));
                        end
                        
		end
	end
end


fclose(fid);


if mergeZOSZOSTOGA
	subcomp=[colGICtot colGIS colAIStot colOD colLS colGIA];
	complabl={'GIC','GIS','AIS','Oceanographic','Global Water Storage','GIA/Tectonics'};
else
	subcomp=[colGICtot colGIS colAIStot colTE colLS colOD colGIA];
	complabl={'GIC','GIS','AIS','Thermal Exp.','Water Storage','Ocean Dynamics','GIA/Tectonics'};
end

subyear = find(targyears==2100);
fid=fopen(['LSLprojcomp_' num2str(targyears(subyear)) '.tsv'],'w');

fprintf(fid,'cm ');
for kk=subscens
	fprintf(fid,[' \t ' scens{kk}]);
	if kk~=subscens(end)
		fprintf(fid,'\t \t \t \t');
	end
end
fprintf(fid,'\n ');

for kk=subscens
	fprintf(fid,['\t 50 \t 17--83 \t 5--95 \t 0.5--99.5 \t 99.9 ' ]);
end

sitesub=[12 299 188 161 10 405 155 7 78 438 134 499];
for jjj=1:length(sitesub)
	jj=find(targregions==sitesub(jjj));
	fprintf(fid,['\n' targregionnames{jj}]);
	for ii=1:length(subcomp)	
		fprintf(fid,['\n' complabl{ii}]); 
		for kk=subscens
			fprintf(fid,' \t ');
			fprintf(fid,'%0.0f \t',quantloccomponents(ia(1),subcomp(ii),subyear,jj,kk)/10);
			fprintf(fid,'%0.0f--%0.0f \t ',quantloccomponents(ia(2:7),subcomp(ii),subyear,jj,kk)/10);
			fprintf(fid,'<%0.0f',5*round(quantloccomponents(ia(8),subcomp(ii),subyear,jj,kk)/10/5));
		end
	end
	fprintf(fid,['\nTOTAL']); 
	for kk=subscens
		fprintf(fid,' \t ');
		fprintf(fid,'%0.0f \t',quanttotlocrise(ia(1),subyear,jj,kk)/10);
		fprintf(fid,'%0.0f--%0.0f \t ',quanttotlocrise(ia(2:7),subyear,jj,kk)/10);
		fprintf(fid,'<%0.0f',5*round(quanttotlocrise(ia(8),subyear,jj,kk)/10/5));
	end
end

fclose(fid);

% local projections for all sites

[jk,ia,ib]=intersect(quantlevs,quantsub); [jk,ic]=sort(ib);
ia=ia(ic);

fid=fopen('LSLproj_full.tsv','w');

fprintf(fid,'cm ');
for kk=subscens
	fprintf(fid,[' \t ' scens{kk}]);
	if kk~=subscens(end)
		fprintf(fid,'\t \t \t \t');
	end
end
fprintf(fid,'\n ');

for kk=subscens
	fprintf(fid,['\t 50 \t 17--83 \t 5--95 \t 0.5--99.5 \t 99.9 ' ]);
end

u=sort(unique(targcoastlines));
for mmm=1:length(u)
	subcl=find(targcoastlines==u(mmm)); subcl=subcl(:)';
	fprintf(fid,'\nCOASTLINE %0.0f',u(mmm));
	[s,si]=sort(targregions(subcl));
	for jjj=s(:)';
		jj=find(targregions==jjj);
		fprintf(fid,['\n' targregionnames{jj} ' [%0.0f] (Bkgd: %0.2f +/- %0.2f mm/y)'], [targregions(jj) rateprojs(jj) rateprojssd(jj)*2]);	
		for mm=1:length(targyears2)
			subyear = find(targyears==targyears2(mm));
			fprintf(fid,'\n%4.0f',targyears2(mm)); 
			for kk=subscens
				fprintf(fid,' \t ');
                                if targyears2(mm)>roundafter
                                    fprintf(fid,'%0.0f \t',10*round(quanttotlocrise(ia(1),subyear,jj,kk)/10/10));
                                    fprintf(fid,'%0.0f--%0.0f \t ',10*round(quanttotlocrise(ia(2:7),subyear,jj,kk)/10/10));
                                    fprintf(fid,'<%0.0f',10*round(quanttotlocrise(ia(8),subyear,jj,kk)/10/10));
                                else
                                    fprintf(fid,'%0.0f \t',quanttotlocrise(ia(1),subyear,jj,kk)/10);
                                    fprintf(fid,'%0.0f--%0.0f \t ',quanttotlocrise(ia(2:7),subyear,jj,kk)/10);
                                    fprintf(fid,'<%0.0f',5*round(quanttotlocrise(ia(8),subyear,jj,kk)/10/5));
                                end
   end
		end
	end
end


fclose(fid);


% local projections for all sites and quantiles

fid=fopen('LSLproj_full_allquants.tsv','w');

fprintf(fid,'cm \t');
for kk=subscens
	fprintf(fid,[ scens{kk}]);
	if kk~=subscens(end)
		fprintf(fid,repmat('\t ',1,length(quantlevs)));
	end
end
fprintf(fid,'\n ');

for kk=subscens
	fprintf(fid,['\t %0.1f'],quantlevs*100);
end

u=sort(unique(targcoastlines));
for mmm=1:length(u)
	subcl=find(targcoastlines==u(mmm)); subcl=subcl(:)';
	fprintf(fid,'\nCOASTLINE %0.0f',u(mmm));
	[s,si]=sort(targregions(subcl));
	for jjj=s(:)';
		jj=find(targregions==jjj);
		fprintf(fid,['\n' targregionnames{jj} ' [%0.0f] (Bkgd: %0.2f +/- %0.2f mm/y)'], [targregions(jj) rateprojs(jj) rateprojssd(jj)*2]);	
		for mm=1:length(targyears2)
			subyear = find(targyears==targyears2(mm));
			fprintf(fid,'\n%4.0f',targyears2(mm)); 
			for kk=subscens
				fprintf(fid,'\t %0.0f',quanttotlocrise(:,subyear,jj,kk)/10);
			end
		end
	end
end


fclose(fid);
