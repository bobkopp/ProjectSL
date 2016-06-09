% outputSLRProjections_table_GSLcomponents
% Last updated by Robert Kopp, robert-dot-kopp-at-rutgers-dot-edu, Fri May 23 14:43:11 EDT 2014

% Table: GSL components

roundafter=2100;

subcomp=[colGICtot colGIS colAIStot colTE colLS];
complabl={'GIC','GIS','AIS','Thermal Exp.','Water Storage'};

subyear = find(targyears==2100);
fid=fopen(['GSLprojcomp_' num2str(targyears(subyear)) '.tsv'],'w');

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

for ii=1:length(subcomp)	
	fprintf(fid,['\n' complabl{ii}]); 
	for kk=subscens
		fprintf(fid,' \t ');
		fprintf(fid,'%0.0f \t',quantcomponents(ia(1),subcomp(ii),subyear,kk)/10);
		fprintf(fid,'%0.0f--%0.0f \t ',quantcomponents(ia(2:7),subcomp(ii),subyear,kk)/10);
		fprintf(fid,'<%0.0f',5*round(quantcomponents(ia(8),subcomp(ii),subyear,kk)/10/5));
	end
end
fprintf(fid,['\nTOTAL']); 
for kk=subscens
	fprintf(fid,' \t ');
	fprintf(fid,'%0.0f \t',quanttotrise(ia(1),subyear,kk)/10);
	fprintf(fid,'%0.0f--%0.0f \t ',quanttotrise(ia(2:7),subyear,kk)/10);
	fprintf(fid,'<%0.0f',5*round(quanttotrise(ia(8),subyear,kk)/10/5));
end

fclose(fid);

% and repeat for all key years

targyears2=[2030 2050 2100 2150 2200];

subcomp=[colGICtot colGIS colAIStot colTE colLS];
complabl={'GIC','GIS','AIS','Thermal Exp.','Water Storage'};

fid=fopen(['GSLprojcomp_allyrs.tsv'],'w');
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

for nn=1:length(targyears2)
	subyear = find(targyears==targyears2(nn));

	fprintf(fid,'\n%4.0f',targyears2(nn));

	for ii=1:length(subcomp)	
		fprintf(fid,['\n' complabl{ii}]); 
		for kk=subscens
			fprintf(fid,' \t ');
                        if targyears2(nn)>roundafter
                            fprintf(fid,'%0.0f \t',10*round(quantcomponents(ia(1),subcomp(ii),subyear,kk)/10/10));
                            fprintf(fid,'%0.0f--%0.0f \t ',10*round(quantcomponents(ia(2:7),subcomp(ii),subyear,kk)/10/10));
                            fprintf(fid,'<%0.0f',10*round(quantcomponents(ia(8),subcomp(ii),subyear,kk)/10/10));
                        else
                            
                            fprintf(fid,'%0.0f \t',quantcomponents(ia(1),subcomp(ii),subyear,kk)/10);
                            fprintf(fid,'%0.0f--%0.0f \t ',quantcomponents(ia(2:7),subcomp(ii),subyear,kk)/10);
                            fprintf(fid,'<%0.0f',5*round(quantcomponents(ia(8),subcomp(ii),subyear,kk)/10/5));
                        end
                        
		end
	end
	fprintf(fid,['\nTOTAL']); 
	for kk=subscens
		fprintf(fid,' \t ');
                if targyears2(nn)>roundafter
                    fprintf(fid,'%0.0f \t',10*round(quanttotrise(ia(1),subyear,kk)/10/10));
                    fprintf(fid,'%0.0f--%0.0f \t ',10*round(quanttotrise(ia(2:7),subyear,kk)/10/10));
                    fprintf(fid,'<%0.0f',10*round(quanttotrise(ia(8),subyear,kk)/10/10));
                else
                    fprintf(fid,'%0.0f \t',quanttotrise(ia(1),subyear,kk)/10);
                    fprintf(fid,'%0.0f--%0.0f \t ',quanttotrise(ia(2:7),subyear,kk)/10);
                    fprintf(fid,'<%0.0f',5*round(quanttotrise(ia(8),subyear,kk)/10/5));
                end
                
	end
end

fclose(fid);


% and repeat for all key years and quantiles

targyears2=[2030 2050 2100 2150 2200];

subcomp=[colGICtot colGIS colAIStot colTE colLS];
complabl={'GIC','GIS','AIS','Thermal Exp.','Water Storage'};

fid=fopen(['GSLprojcomp_allyrs_allquants.tsv'],'w');
fprintf(fid,'cm \t');
for kk=subscens
	fprintf(fid,[ scens{kk}]);
	if kk~=subscens(end)
		fprintf(fid,repmat('\t ',1,length(quantlevs)));
	end
end
fprintf(fid,'\n ');
for kk=subscens
	fprintf(fid,['\t %0.1f' ],quantlevs*100);
end

for nn=1:length(targyears2)
	subyear = find(targyears==targyears2(nn));

	fprintf(fid,'\n%4.0f',targyears2(nn));

	for ii=1:length(subcomp)	
		fprintf(fid,['\n' complabl{ii}]); 
		for kk=subscens
			fprintf(fid,'\t %0.0f',quantcomponents(:,subcomp(ii),subyear,kk)/10);
		end
	end
	fprintf(fid,['\nTOTAL']); 
	for kk=subscens
		fprintf(fid,'\t %0.0f',quanttotrise(:,subyear,kk)/10);
	end
end

fclose(fid);
