% outputSLRProjections_table_CMIP5models
% Last updated by Robert Kopp, robert-dot-kopp-at-rutgers-dot-edu, Fri May 23 16:05:10 EDT 2014

%%
% table of CMIP5 models

subscens=[1 3 4];
fullmodels=union(projGICmodel{1},ZOSTOGAmodels{1});

if exist('ZOSmodels','var')
	clear ZOSmodellast;
	for kk=subscens
		workmodels=ZOSmodels{kk};
		for ii=1:length(workmodels)
			u=find(strcmp(workmodels{ii},fullmodels));
			if length(u)==0
				fullmodels={fullmodels{:},workmodels{ii}};
				u=length(fullmodels);
			end
			 q=find(~isnan(ZOS{kk}(:,ii)));
			 if length(q)>0
				ZOSmodellast(u,kk) = OceanDynYears(q(end));
			end
		end
	end
end

clear ZOSTOGAmodellast;
for kk=subscens
	workmodels=ZOSTOGAmodels{kk};
	for ii=1:length(workmodels)
		u=find(strcmp(workmodels{ii},fullmodels));
		if length(u)==0
			fullmodels={fullmodels{:},workmodels{ii}};
			u=length(fullmodels);
		end
		 q=find(~isnan(ZOSTOGA{kk}(:,ii)));
		 if length(q)>0
			ZOSTOGAmodellast(u,kk) = ThermExpYears(q(end));
		end
	end
end

clear GICmodellast
for kk=subscens
	workmodels=projGICmodel{kk};
	for ii=1:length(workmodels)
		u=find(strcmp(workmodels{ii},fullmodels));
		if length(u)==0
			fullmodels={fullmodels{:},workmodels{ii}};
			u=length(fullmodels);
		end
		 q=find(~isnan(projGIC{kk}(:,1,ii)));
		 if length(q)>0
			GICmodellast(u,kk) = projGICyrs{kk}(q(end),ii);
		end
	end
end

[s,si]=sort(fullmodels);
if exist('ZOSmodels','var')
	modellast={ZOSTOGAmodellast,ZOSmodellast,GICmodellast};
else
	modellast={ZOSTOGAmodellast,GICmodellast};
end

fid=fopen('CMIP5models.tsv','w');
spacer=repmat('\t',1,length(subscens));
fprintf(fid,['Model\tOceanographic' spacer]);
if exist('ZOSmodels','var')
	fprintf(fid,[ 'ZOS' spacer]);
end
fprintf(fid,['GIC' spacer '\n']);
fprintf(fid,['\t']);
Ntodo=2+exist('ZOSmodels','var');
for i=1:Ntodo
	for j=subscens
		fprintf(fid,[scens{j} '\t']);
	end
end
fprintf(fid,'\n');
for i=si
	fprintf(fid,[fullmodels{i}]);
	for jj=1:length(modellast)
		for kk=subscens
			fprintf(fid,'\t');
			mrkr='';
			if size(modellast{jj},1)>=i
				if modellast{jj}(i,kk)>2290
					mrkr='3';
				elseif modellast{jj}(i,kk)>=2099
					mrkr='1';
				end
			end
			fprintf(fid,[mrkr]);
		end
	end
	fprintf(fid,'\n');
end
fclose(fid);