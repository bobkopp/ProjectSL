% outputSLRProjections_table_LSLNorthAmerica
% Last updated by Robert Kopp, robert-dot-kopp-at-rutgers-dot-edu, Fri May 23 15:51:16 EDT 2014

% North America all years

subscens=[1 3 4];
targyears2 = targyears(find(targyears<=2200));
quantsub=[.001 .005 .01 .025 .05:.05:.95 .975 .99 .995 .999];
[jk,ia,ib]=intersect(round(quantlevs*10000),round(quantsub*10000)); [jk,ic]=sort(ib);
ia=ia(ic);

fid=fopen('NorthAmerica_LSLproj_full.tsv','w');

fprintf(fid,'cm \t');
for kk=subscens
	fprintf(fid,[ scens{kk}]);
	if kk~=subscens(end)
		fprintf(fid,repmat('\t ',1,length(quantlevs(ia))));
	end
end
fprintf(fid,'\n ');

for kk=subscens
	fprintf(fid,['\t %0.1f'],quantlevs(ia)*100);
end

u=[760 820 821 822 823 940 960 970];
cllabel={'Hawaii','Aleutian Islands','Alaska','Pacific Canada','Pacific US','Gulf of Mexico US','Atlantic US','Atlantic Canada'};
for mmm=1:length(u)
	subcl=find(targcoastlines==u(mmm)); subcl=subcl(:)';
	if length(subcl)>0
%	fprintf(fid,'\n\nCOASTLINE %0.0f',u(mmm));
        fprintf(fid,['\n\n' cllabel{mmm}]);
        [s,si]=sort(targregions(subcl));
        for jjj=s(:)';
            jj=find(targregions==jjj);
            fprintf(fid,['\n\n' targregionnames{jj} ]);
            fprintf(fid,'\t[%0.0f]',targregions(jj));
            fprintf(fid,'\t%0.2f',sitecoords(jj,:));
            for mm=1:length(targyears2)
                subyear = find(targyears==targyears2(mm));
                fprintf(fid,'\n%4.0f',targyears2(mm)); 
                for kk=subscens
                    fprintf(fid,'\t %0.0f',quanttotlocrise(ia,subyear,jj,kk)/10);
                end
            end
        end
    end
end


fclose(fid);

% look at NA in 2100, downselecting projections

subscens=[1 3 4];
quantsub=[.01 .025 .05 .1:.1:.9 .95 .975 .99 .995 .999];
overallquantlevs = [.005 .01 .05 .167 .333 .5 .667 .833 .95 .99 .995 .999];
quantgrid=0:.001:1;
[jk,ia,ib]=intersect(round(quantlevs*10000),round(quantsub*10000)); [jk,ic]=sort(ib);
ia=ia(ic);
subyear = find(targyears==2100);

clear ulevs uquants;

u0=[760 820 821 822 823 940 960 970];
for mmm=1:length(u0)
	subcl=find(targcoastlines==u0(mmm)); subcl=subcl(:)';
	[s,si]=sort(targregions(subcl));
	for jjj=s(:)';
		jj=find(targregions==jjj);

		for tt=1:length(targyears)

			clear u;
			for kkk=1:length(subscens)
				kk=subscens(kkk);
				u(:,kkk) = interp1([quantlevs],quanttotlocrise(:,tt,jj,kk)/10,quantgrid);
			end
			u=u(find(~isnan(u)));
			uquants(:,tt,jj) = quantile(u,overallquantlevs);
		end

		for kk=subscens
			sub=find(~isnan(quanttotlocrise(:,subyear,jj,kk)));
			if length(sub)>1
				[s1,si]=sort(quanttotlocrise(sub,subyear,jj,kk)/10);
				s2=(quantlevs(sub(si)));
				ulevs(:,kk,jj)=interp1(s1,s2,uquants(:,subyear,jj));
			end
		end

	end
end

fid=fopen('NorthAmerica_LSLproj_weighted.tsv','w');

fprintf(fid,'cm');
fprintf(fid,'\n ');

u0=[760 820 821 822 823 940 960 970];
targyears2=targyears(find(targyears<=2100));
for mmm=1:length(u0)
	subcl=find(targcoastlines==u0(mmm)); subcl=subcl(:)';
%	fprintf(fid,'\n\nCOASTLINE %0.0f',u(mmm));
	fprintf(fid,['\n\n' cllabel{mmm}]);
	[s,si]=sort(targregions(subcl));
	for jjj=s(:)';
		jj=find(targregions==jjj);
		fprintf(fid,['\n\n' targregionnames{jj} ]);
		fprintf(fid,'\t[%0.0f]',targregions(jj));
		fprintf(fid,'\t%0.2f',sitecoords(jj,:));
		fprintf(fid,'\n');
		fprintf(fid,'Overall');
		fprintf(fid,'\t %0.2f',overallquantlevs*100);
		fprintf(fid,'\n');		
		for kk=subscens
			fprintf(fid,[scens{kk} ' %%ile']);
			fprintf(fid,'\t %0.2f',ulevs(:,kk,jj)*100);
			fprintf(fid,'\n');
		end
		for mm=1:length(targyears2)
			subyear = find(targyears==targyears2(mm));
			fprintf(fid,'\n%4.0f',targyears2(mm)); 
			fprintf(fid,'\t%0.0f',uquants(:,mm,jj));
		end
	end
end

fclose(fid);

%%%

% Table: full local MC for select sites


subscens=[1 3 4];

smoothwin=19;
%sitesub=[235 12 299 396 363  188 520 526 161 10 127 158 405 155];
sitesub=[
183 % Portland, ME
235 % Boston, MA
12 % New York, NY
351 % Newport, CT
180 % Atlantic City, NJ
135 % Philadelphia, PA
224 % Lewes, DE
148 % Baltimore, MD
360 % Washington, DC
299 % Sewell's Pt, Norfolk, VA
396 % Wilmington, NC
234 % Charleston, NC
395 % Fort Pulaski, Savannah, GA
112 % Ferdinanda Beach, FL
188 % Key West, FL
363 % Miami, FL
246 % Pensacola, FL
526 % Grand Isle, LA
161 % Galveston, TX
158 % San Diego, CA
10 % San Francisco, CA
265 % Astoria, OR
127 % Seattle, WA
405 % Juneau, AK
1067 % Anchorage, AK
155 % Honolulu, HI
78  % Stockholm
7  % Cuxhaven, Germany
438  % Kochi, India
134 % Kushimoto, Japan
499 % Valparaiso, Chile
42 % Sevastopol
]; sitesub=sitesub';

sitesub=sitesub';

[histX1,histY,histdY,histregions,histregionsu]=ReadPSMSLData(-1,-1,[],psmsldir,gslfile,sitesub);

testt=1890:2012;
clear histY1 histY1s;
for jj=1:length(sitesub)
	sub=find(histregions==sitesub(jj));
	histY1(:,jj)=interp1(histX1(sub,3),histY(sub),testt,'linear');
	histY1s(:,jj) = smooth(histY1(:,jj),smoothwin);
	subanchor=find(testt==2000);
	histY1s(:,jj) = histY1s(:,jj)-histY1s(subanchor,jj);
	histY1(:,jj) = histY1(:,jj)-histY1s(subanchor,jj);
	sub=find(isnan(histY1(:,jj))); histY1s(sub,jj)=NaN;
end
subtimes=find((testt>(testt(1)+smoothwin/2)).*(testt<(testt(end)-smoothwin/2)));
testt=testt(subtimes); histY1=histY1(subtimes,:); histY1s=histY1s(subtimes,:);

for jjj=1:length(sitesub)
	jj2=find(focussites==sitesub(jjj));
	if length(jj2)==1
		jj=find(targregions==sitesub(jjj));
		for kk=subscens

			fid=fopen(['LSLprojMC_' num2str(sitesub(jjj)) '_' scens{kk} '.tsv'],'w');
			fprintf(fid,[targregionnames{jj} '\t' scens{kk} '\t(mm)\n']);
			
			sub = find(~isnan(histY1s(:,jjj)));
			sub=intersect(sub,find(mod(testt,10)==0));
			for mm=1:length(sub)
				fprintf(fid,'%0.0f',testt(sub(mm)));
				fprintf(fid,'\t%4.0f',repmat(histY1s(sub(mm),jjj),1,size(totlocalrisefocus{jj2,kk},1)));
				fprintf(fid,'\n');
			end
	
			for mm=1:length(targyears)
				fprintf(fid,'%0.0f',targyears(mm));
				fprintf(fid,'\t%4.0f',totlocalrisefocus{jj2,kk}(:,mm));
				fprintf(fid,'\n');
			end
			fclose(fid);

	
		end
	end
end

fid=fopen('TideGaugeIndex.txt','w');
for jjj=1:length(sitesub)
	jj=find(targregions==sitesub(jjj));
	fprintf(fid,['%0.0f\t' targregionnames{jj} '\n'],sitesub(jjj));
end
fclose(fid);



%%%%

% local projections for all sites and quantiles

targyears2 = [2010:10:2200];
quantsub=[.001 .005 .01 .025 .05:.05:.95 .975 .99 .995 .999];

[jk,ia,ib]=intersect(round(quantlevs*10000),round(quantsub*10000)); [jk,ic]=sort(ib);
ia=ia(ic);

for doall=0:1

    if doall
        fid=fopen('NorthAmerica_LSLproj_full_allquants.tsv','w');
    else
        fid=fopen('NorthAmerica_LSLproj_subset_allquants.tsv','w');
        sitesub=[235 12 299 396 363  188 520 526 161 10 127 158 405 155];
    end
    
    fprintf(fid,'cm \t');
    for kk=subscens
        fprintf(fid,[ scens{kk}]);
        if kk~=subscens(end)
            fprintf(fid,repmat('\t ',1,length(quantlevs(ia))));
        end
    end
    fprintf(fid,'\n ');

    for kk=subscens
        fprintf(fid,['\t %0.1f'],quantlevs(ia)*100);
    end

    fprintf(fid,'\nGlobal');
    for mm=1:length(targyears2)
        subyear = find(targyears==targyears2(mm));
        fprintf(fid,'\n%4.0f',targyears2(mm)); 
        for kk=subscens
            fprintf(fid,'\t %0.0f',quanttotrise(ia,subyear,kk)/10);
        end
    end

    u=[760 820 821 822 823 940 960 970];
    cllabel={'Hawaii','Aleutian Islands','Alaska','Pacific Canada','Pacific US','Gulf of Mexico US','Atlantic US','Atlantic Canada'};

    for mmm=1:length(u)
        subcl=find(targcoastlines==u(mmm)); subcl=subcl(:)';
        if ~doall
            subcl=intersect(subcl,find(ismember(targregions,sitesub)));
        end
    %	fprintf(fid,'\n\nCOASTLINE %0.0f',u(mmm));
        if length(subcl)>0
            fprintf(fid,['\n\n' cllabel{mmm}]);
            [s,si]=sort(targregions(subcl));
            for jjj=s(:)';
                jj=find(targregions==jjj);
                fprintf(fid,['\n' targregionnames{jj} ' [%0.0f / %0.2f, %0.2f]'], [targregions(jj) sitecoords(jj,:)]);	
                for mm=1:length(targyears2)
                    subyear = find(targyears==targyears2(mm));
                    fprintf(fid,'\n%4.0f',targyears2(mm)); 
                    for kk=subscens
                        fprintf(fid,'\t %0.0f',quanttotlocrise(ia,subyear,jj,kk)/10);
                    end
                end
            end
        end
    end


    fclose(fid);
end

% pick special levels
slev=[0.05 2050 ; 0.5 2050 ; 0.05 2100 ; 0.5 2100 ;  0.95 2100 ; 0.05 2200 ; 0.5 2200];
subscens=[1 3 4];
targyears2=[2030 2050 2090 2100 2200];

fid=fopen('NorthAmerica_LSLproj_heightlevels.tsv','w')

fprintf(fid,'cm above 2000\t');
fprintf(fid,'Probability inundation exceeds');
fprintf(fid,'\n\t');
for ii=1:length(subscens)
    fprintf(fid,scenlab{ii});
    for jj=1:length(targyears2);
        fprintf(fid,'\t');
    end
end
fprintf(fid,'\n');
for ii=1:length(subscens)
    for jj=1:length(targyears2);
        fprintf(fid,'\t%0.0f',targyears2(jj));
    end
end

for ii=1:length(sitesub)
    jj = find((targregions==sitesub(ii)));
    fprintf(fid,['\n\n' targregionnames{jj} ' [%0.0f / %0.2f, %0.2f]'], [targregions(jj) sitecoords(jj,:)]);  
    kk=1;
    clear sheight sheightq;
    for jjj=1:size(slev,1)
        ia = find(round(quantlevs*10000)==slev(jjj,1)*10000);
        subyear = find((targyears==slev(jjj,2)));
        sheight(jjj) = quanttotlocrise(ia,subyear,jj,kk)/10;
    end
    sheight=sort(sheight);
    sub=[1 find(diff(sheight)>10)+1];
    sheight=sheight(sub);
    for jjj=1:length(sheight)
        fprintf(fid,'\n%0.0f',sheight(jjj));
        for kkk=1:length(subscens)
            for mmm=1:length(targyears2)
                subyear = find((targyears==targyears2(mmm)));
                kk=subscens(kkk);
                sheightq(jjj,kkk,mmm) = interp1(quanttotlocrise(:,subyear,jj,kk)/10,quantlevs,sheight(jjj));
                if isnan(sheightq(jjj,kkk,mmm))
                    if sheight(jjj)>(quanttotlocrise(end,subyear,jj,kk)/10)
                        sheightq(jjj,kkk,mmm) = 1;
                    elseif sheight(jjj)<(quanttotlocrise(1,subyear,jj,kk)/10)
                        sheightq(jjj,kkk,mmm) = 0;
                    end
                end
                fprintf(fid,'\t%0.1f%%',100-sheightq(jjj,kkk,mmm)*100);
            end 
        end
    end
 
 end

fclose(fid);  

% median projections table

targyears2 = [2100];
quantsub=[.5];

[jk,ia,ib]=intersect(round(quantlevs*10000),round(quantsub*10000)); [jk,ic]=sort(ib);
ia=ia(ic);

fid=fopen('NorthAmerica_LSLproj_full_median.tsv','w');   
fprintf(fid,'cm \t\t\t\t\t\t');
for kk=subscens
    fprintf(fid,[ scens{kk}]);
    if kk~=subscens(end)
        fprintf(fid,repmat('\t ',1,length(quantlevs(ia))));
    end
end
fprintf(fid,'\n ');

fprintf(fid,'Coast\tSite\tSiteID\tlat\tlong\tyear');

for kk=subscens
    fprintf(fid,['\t' scens{kk} '_%0.1f'],quantlevs(ia)*100);
end

for mm=1:length(targyears2)
    subyear = find(targyears==targyears2(mm));
    fprintf(fid,['\nGlobal\tGlobal\t0\tNaN\tNaN']);
    fprintf(fid,'\t%4.0f',targyears2(mm)); 
    for kk=subscens
        fprintf(fid,'\t %0.0f',quanttotrise(ia,subyear,kk)/10);
    end
end

u=[760 820 821 822 823 940 960 970];
cllabel={'Hawaii','Aleutian Islands','Alaska','Pacific Canada','Pacific US','Gulf of Mexico US','Atlantic US','Atlantic Canada'};

for mmm=1:length(u)
    subcl=find(targcoastlines==u(mmm)); subcl=subcl(:)';
    if length(subcl)>0
        [s,si]=sort(targregions(subcl));
        for jjj=s(:)';
            jj=find(targregions==jjj);
            for mm=1:length(targyears2)
                subyear = find(targyears==targyears2(mm));
                fprintf(fid,['\n' cllabel{mmm} '\t' targregionnames{jj} '\t%0.0f\t%0.2f\t%0.2f'], [targregions(jj) sitecoords(jj,:)]);	
                fprintf(fid,'\t%4.0f',targyears2(mm)); 
                for kk=subscens
                    fprintf(fid,'\t %0.0f',quanttotlocrise(ia,subyear,jj,kk)/10);
                end
            end
        end
    end
end


fclose(fid);

