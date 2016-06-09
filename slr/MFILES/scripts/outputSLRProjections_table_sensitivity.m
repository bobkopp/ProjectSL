% outputSLRProjections_table_sensitivity
% Last updated by Robert Kopp, robert-dot-kopp-at-rutgers-dot-edu, Fri May 23 15:00:09 EDT 2014

targyears2=2100;
subscens=1;
quantsub=[.5 .167 .833 .05 .95 .005 .995 .999];
[jk,ia,ib]=intersect(quantlevs,quantsub); [jk,ic]=sort(ib);
ia=ia(ic);

sensquantsGSL={quanttotrise,quanttotriseAR,quanttotriseBA,quanttotriseAltCorr{1},quanttotriseHiGCMConf,quanttotriseRedDOF,quanttotrisePokhrel};

sensquants={quanttotlocriseDef,quanttotlocriseAR,quanttotlocriseBA,quanttotlocriseAltCorr{1},quanttotlocriseHiGCMConf,quanttotlocriseRedDOF};
senslabels={'Default','AR5','BA','Alt. Corr.','High GCM Conf.','Red. DOF','LWS'};

fid=fopen('SLRproj_sensitivity.tsv','w');

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
subyear=find(targyears==targyears2(1));


fprintf(fid,'\nGSL');
for pp=1:length(sensquantsGSL)
	fprintf(fid,'\n');
	fprintf(fid,senslabels{pp}); 
	for kk=subscens
		fprintf(fid,' \t ');
		fprintf(fid,'%0.0f \t',sensquantsGSL{pp}(ia(1),subyear,kk)/10);
		fprintf(fid,'%0.0f--%0.0f \t ',sensquantsGSL{pp}(ia(2:7),subyear,kk)/10);
		fprintf(fid,'<%0.0f',5*round(sensquantsGSL{pp}(ia(8),subyear,kk)/10/5));
	end
end

sitesub=focussites;
for jjj=1:length(sitesub)
	jj=find(targregions==sitesub(jjj));
	fprintf(fid,['\n' targregionnames{jj} ' (Bkgd: %0.2f +/- %0.2f mm/y)'], [rateprojs(jj) rateprojssd(jj)*2]);	
	jj=find(targregions2==sitesub(jjj));
	for pp=1:length(sensquants)
		fprintf(fid,'\n');
		fprintf(fid,senslabels{pp}); 
		for kk=subscens
			fprintf(fid,' \t ');
			fprintf(fid,'%0.0f \t',sensquants{pp}(ia(1),subyear,jj,kk)/10);
			fprintf(fid,'%0.0f--%0.0f \t ',sensquants{pp}(ia(2:7),subyear,jj,kk)/10);
			fprintf(fid,'<%0.0f',5*round(sensquants{pp}(ia(8),subyear,jj,kk)/10/5));
		end
	end
end

fprintf(fid,'\n');
fclose(fid);

%%


sensquants={quantcomponents,quantcomponentsAR,quantcomponentsBA,quantcomponentsAltCorr{1}};
senslabels={'Default','AR5','BA','Alt. Corr.'};

colselect=[colGIS colAIStot];
colselectlabels={'GIS','AIS'};

kk=1;

fid=fopen(['ISproj_sensitivity_' scens{kk} '.tsv'],'w');

fprintf(fid,'cm ');
for pp=1:length(colselect)
	fprintf(fid,[' \t ' colselectlabels{pp}]);
	if pp~=length(colselect)
		fprintf(fid,'\t \t \t \t');
	end
end
fprintf(fid,'\n ');
for pp=colselect
	fprintf(fid,['\t 50 \t 17--83 \t 5--95 \t 0.5--99.5 \t 99.9 ' ]);
end
subyear=find(targyears==targyears2(1));


for pp=1:length(sensquants)
	fprintf(fid,'\n');
	fprintf(fid,senslabels{pp}); 
	for qq=colselect
		fprintf(fid,' \t ');
		fprintf(fid,'%0.0f \t',sensquants{pp}(ia(1),qq,subyear,kk)/10);
		fprintf(fid,'%0.0f--%0.0f \t ',sensquants{pp}(ia(2:7),qq,subyear,kk)/10);
		fprintf(fid,'<%0.0f',5*round(sensquants{pp}(ia(8),qq,subyear,kk)/10/5));
	end
end



