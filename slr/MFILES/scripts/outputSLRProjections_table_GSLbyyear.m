% outputSLRProjections_table_GSLbyyear
% Last updated by Robert Kopp, robert-dot-kopp-at-rutgers-dot-edu, Fri May 23 14:39:07 EDT 2014

% Table: GSL projections

subscens=[1 3 4];
defval('targyears2',[2030 2050 2100 2150 2200]);
roundafter=2100;

quantsub=[.5 .167 .833 .05 .95 .005 .995 .999];
[jk,ia,ib]=intersect(quantlevs,quantsub); [jk,ic]=sort(ib);
ia=ia(ic);

fid=fopen('GSLproj.tsv','w');

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
for mm=1:length(targyears2)
	subyear = find(targyears==targyears2(mm));
	fprintf(fid,'\n%4.0f',targyears2(mm)); 
	for kk=subscens
		fprintf(fid,' \t ');
                if targyears2(mm)>roundafter
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

fprintf(fid,'\n');
fprintf(fid,'AR5 - 2100');
fprintf(fid,'\t 73 \t 53--97 \t \t \t'); 
fprintf(fid,'\t 52 \t 35--70 \t \t \t'); 
fprintf(fid,'\t 43 \t 28--60 \t \t \t'); 

fprintf(fid,'\n');
fprintf(fid,'H14 - 2100');
fprintf(fid,'\t  \t 70--120 \t 50--150 \t \t'); 
fprintf(fid,'\t \t  \t \t \t'); 
fprintf(fid,'\t \t 40--60 \t 25--70 \t \t'); 

fprintf(fid,'\n');
fprintf(fid,'J12 - 2100');
fprintf(fid,'\t 110  \t  \t 81--165 \t \t'); 
fprintf(fid,'\t 75 \t  \t 52--110 \t \t'); 
fprintf(fid,'\t 57 \t \t 36--83 \t \t'); 

fprintf(fid,'\n');
fprintf(fid,'S12 - 2100');
fprintf(fid,'\t  \t  \t \t \t'); 
fprintf(fid,'\t 90 \t  \t 64--121 \t \t'); 
fprintf(fid,'\t 75 \t \t 52--96 \t \t'); 


fprintf(fid,'\n');

fclose(fid);