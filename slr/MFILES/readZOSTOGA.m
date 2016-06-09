function [ThermExpMean,ThermExpStd,ThermExpYears,ThermExpN,ZOSTOGA,sZOSTOGAd2,CWdrift,histGICrate,modellist,years]=readZOSTOGA(scen,smoothwin,projGIC85,years,driftcorr,varname,IFILES,PARAMDIR)

% Last updated by  Bob Kopp, robert-dot-kopp-at-rutgers-dot-edu, Tue Feb 11 22:53:08 EST 2014


defval('scen','rcp85');
defval('smoothwin',19);
defval('doplot',0);
defval('IFILES','IFILES/slr');
defval('PARAMDIR','PARAMS/');
defval('subdir',fullfile(IFILES,'/SLR_ALL'));
defval('years',1860:2099);
defval('projGIC85',[]);
defval('baseyear',2000);
defval('varname','zostoga');

defval('driftcorr',1);

ZOSTOGAbaseyear=[baseyear-10 baseyear+10];

% read ZOSTOGA

if iscell(varname)
	Nvars=length(varname);
else
	Nvars=1;
	varname={varname};
end

clear ZOSTOGA;
clear sZOSTOGA;

pd=pwd;
cd(subdir);
cd(scen);
files=dir('*');
modellist={};
jj=1;

for ii=1:length(files)
	if (files(ii).isdir) && (files(ii).name(1)~='.')
		curmodel=lower(files(ii).name);
		foundvar=0;
		for mm=1:Nvars
			if (exist([files(ii).name,'/' varname{mm}],'dir'))&&~foundvar
				foundvar=1;
				cd(files(ii).name);
				cd(varname{mm})
				files2=dir('*.txt');
				dat=importdata(files2(1).name); dat=dat.data;
				incorporate=1;
				% ignore if total change is too small
				subnonan=find(~isnan(dat(:,2)));
				totalchange=sum((diff(dat(subnonan(1):subnonan(end),2))))/(dat(subnonan(end),1)-dat(subnonan(1),1));
				if totalchange < 2e-4
					incorporate=0;
					disp([curmodel ' - DISCARDED']);
				end
				
				% ignore if doesn't have data near 2000
				if min(abs(dat(subnonan,1)-2000))>2
					incorporate=0;
				end
				
				if incorporate
					modellist{jj}=curmodel;
					disp([curmodel ' - ' varname{mm}]);
					subnonan=find(~isnan(dat(:,2)));
					ZOSTOGA(:,jj)=interp1(dat(subnonan,1),dat(subnonan,2),years);
					jj=jj+1;
				end
				cd ../..
			end
		end
	end
end
cd(pd);
sub=find((years<=ZOSTOGAbaseyear(2)).*(years>=ZOSTOGAbaseyear(1)));
ZOSTOGA0=ZOSTOGA;
sZOSTOGA = NaN*ZOSTOGA;

for jj=1:size(ZOSTOGA,2);
	sub1=sub(find(~isnan(ZOSTOGA(sub,jj))));
	if length(sub1)>0
		ZOSTOGA(:,jj)=ZOSTOGA(:,jj)-mean(ZOSTOGA(sub1,jj));
	end

	% fix ZOSTOGA suturing problem
	sub2=find(years==2007);
	diffa = (ZOSTOGA(sub2,jj)-ZOSTOGA(sub2-1,jj));
	diffb = (ZOSTOGA(sub2-2,jj)-ZOSTOGA(sub2-1,jj));
	if abs(diffa)>20*abs(diffb)
		offset=diffa-diffb;
		ZOSTOGA(sub2:end,jj)=ZOSTOGA(sub2:end,jj)-offset;
	end

	sub = find(~isnan(ZOSTOGA(:,jj)));
	sZOSTOGA(sub,jj) = smooth(ZOSTOGA(sub,jj),smoothwin);
end


sZOSTOGA(:,:) = bsxfun(@minus,sZOSTOGA(:,:),sZOSTOGA(find(years==baseyear),:));

if driftcorr

	% read GSL curve
	psmsldir=fullfile(IFILES,'rlr_annual');
    gslfile=fullfile(IFILES,'CSIRO_Recons_gmsl_yr_2011.csv');

	[GSLx,GSLy]=ReadDenoisedPSMSLData(0,0,[],[],1,[],psmsldir,gslfile);
	GSLy = GSLy/1000;
	sub=find(GSLx(:,3)<=1900);
	CWdrift = (GSLy(sub(end))-GSLy(sub(1)))/(GSLx(sub(end),3)-GSLx(sub(1),3));
	GSLy = GSLy-mean(GSLy(find(floor(GSLx(:,3))==baseyear)));

	selectyears=[1861 floor(GSLx(sub(end),3))];

	% read GIC projections for drift

	if length(projGIC85)==0
		[projGIC85,projGIC85se,projGIC85yrs,projGIC85model] = readMarzeion('rcp85',IFILES,PARAMDIR);
	end

	GICyears=projGIC85yrs(:,1);
	clear GICmodelmap;
	sub1=find(GICyears==selectyears(1));
	sub2=find(GICyears==selectyears(2));
	histGIC=squeeze(sum(projGIC85(:,:,:),2));
	dhistGIC=histGIC(sub2,:)-histGIC(sub1,:);
	histGICrate = dhistGIC/diff(selectyears)/1000;

	[jk,subi,subj]=intersect(floor(GSLx(:,3)),GICyears);
	GSLlessGIC = GSLy(subi)-mean(sum(projGIC85(subj,:,:),2),3)/1000;
	GSLlessGICa = GSLy(subi)-min(sum(projGIC85(subj,:,:),2),[],3)/1000;
	GSLlessGICb = GSLy(subi)-max(sum(projGIC85(subj,:,:),2),[],3)/1000;
	sub=find(jk==baseyear); GSLlessGIC=GSLlessGIC-GSLlessGIC(sub);
	GSLlessGICa=GSLlessGICa-GSLlessGICa(sub);GSLlessGICb=GSLlessGICb-GSLlessGICb(sub);

	for i=1:length(modellist)
		 s=find(strcmpi(modellist{i},projGIC85model));
		 if length(s)>0
			GICmodelmap(i) = s(1);
		 else
			GICmodelmap(i) = 0;
		 end
	end
	sub = find(GICmodelmap>0);
	histresrate(sub) = CWdrift - histGICrate(GICmodelmap(sub));
	sub = find(GICmodelmap==0);
	histresrate(sub) = CWdrift - mean(histGICrate);
	histresratea = CWdrift-mean(histGICrate);

	for jj=1:size(ZOSTOGA,2);
		sub1=find(years==selectyears(1)); sub2=find(years==selectyears(2));
		drift(jj) = (sZOSTOGA(sub2,jj)-sZOSTOGA(sub1,jj))/diff(selectyears);
		sZOSTOGAd(:,jj) = sZOSTOGA(:,jj)-(drift(jj))*(years'-selectyears(1));
	%		sZOSTOGAd2(:,jj) = sZOSTOGAd(:,jj)+histresrate(jj)*(years'-selectyears(1));
		sZOSTOGAd2(:,jj) = sZOSTOGAd(:,jj)+histresratea*(years'-selectyears(1));
		sub=find(years==baseyear);
		sZOSTOGAd(:,jj)=sZOSTOGAd(:,jj)-mean(sZOSTOGAd(sub,jj));
	%		sZOSTOGAd2(:,jj)=sZOSTOGAd2(:,jj)-mean(sZOSTOGAd2(sub,jj));
		sZOSTOGAd2(:,jj)=sZOSTOGAd2(:,jj)-mean(sZOSTOGAd2(sub,jj));
	end
else
	sZOSTOGAd2 = sZOSTOGA;
	CWdrift=NaN; histGICrate=NaN;
end

ThermExpYears = years;
for jj=1:size(ZOSTOGA,1)
	sub=find(~isnan(sZOSTOGAd2(jj,:)));
	ThermExpMean(jj) = squeeze(mean(sZOSTOGAd2(jj,sub),2))*1000;
	ThermExpStd(jj)= squeeze(std(sZOSTOGAd2(jj,sub),[],2))*1000;
	ThermExpN(jj) = length(sub);
	if driftcorr
		ThermExpStd(jj) = sqrt(ThermExpStd(jj).^2 + (std(histGICrate)*(ThermExpYears(jj)-selectyears(1))).^2); % add error to account for uncertainty in histGICrate
	end
end

ThermExpMean(end+1) = ThermExpMean(end) + (ThermExpMean(end)-ThermExpMean(end-1));
ThermExpStd(end+1) = ThermExpStd(end) + (ThermExpStd(end)-ThermExpStd(end-1));
ThermExpYears(end+1) = ThermExpYears(end) + (ThermExpYears(end)-ThermExpYears(end-1));
ThermExpN(end+1)=ThermExpN(end);

ThermExpMean=ThermExpMean(:); ThermExpStd=ThermExpStd(:); ThermExpYears=ThermExpYears(:); ThermExpN=ThermExpN(:);

if doplot
	clf;
	plot(GSLx(:,3),GSLy,'k','linew',2);
	hold on
	plot(GSLx(subi,3),GSLlessGIC(subi),'k--','linew',2);
	plot(years,sZOSTOGAd);
	legend('GSL','GSL - [mean/max/min] GIC proj','Location','Northwest');
	plot(GSLx(:,3),GSLy,'k','linew',2);
	plot(GSLx(subi,3),GSLlessGIC(subi),'k--','linew',2);
	plot(GSLx(subi,3),GSLlessGICa(subi),'k--','linew',2);
	plot(GSLx(subi,3),GSLlessGICb(subi),'k--','linew',2);
	xlim([1860 2100]);  ylim([-.2 .5]);
	title('smoothed zostoga after drift removal');
	pdfwrite('zostoga_driftremoved');

	clf;
	plot(GSLx(:,3),GSLy,'k','linew',2);
	hold on
	plot(GSLx(subi,3),GSLlessGIC(subi),'k--','linew',2);
	plot(years,sZOSTOGA);
	legend('GSL','GSL - [mean/max/min] GIC proj','Location','Northwest');
	plot(GSLx(:,3),GSLy,'k','linew',2);
	plot(GSLx(subi,3),GSLlessGIC(subi),'k--','linew',2);
	plot(GSLx(subi,3),GSLlessGICa(subi),'k--','linew',2);
	plot(GSLx(subi,3),GSLlessGICb(subi),'k--','linew',2);
	xlim([1860 2100]); ylim([-.2 .50]);
	title('smoothed zostoga before drift removal');
	pdfwrite('zostoga');

	clf;
	plot(GSLx(:,3),GSLy,'k','linew',2);
	hold on
	plot(GSLx(subi,3),GSLlessGIC(subi),'k--','linew',2);
	plot(years,sZOSTOGAd2);
	legend('GSL','GSL - [mean/max/min] GIC proj','Location','Northwest');
	plot(GSLx(:,3),GSLy,'k','linew',2);
	plot(GSLx(subi,3),GSLlessGIC(subi),'k--','linew',2);
	plot(GSLx(subi,3),GSLlessGICa(subi),'k--','linew',2);
	plot(GSLx(subi,3),GSLlessGICb(subi),'k--','linew',2);
	xlim([1860 2100]); ylim([-.2 .50]);
	title('smoothed zostoga with drift substitution');
	pdfwrite('zostoga_driftsubst');

	clf;
	plot(GSLx(:,3),GSLy,'k','linew',2);
	hold on
	plot(GSLx(subi,3),GSLlessGIC(subi),'k--','linew',2);
	plot(years,sZOSTOGAd2);
	legend('GSL','GSL - [mean/max/min] GIC proj','Location','Northwest');
	plot(GSLx(:,3),GSLy,'k','linew',2);
	plot(GSLx(subi,3),GSLlessGIC(subi),'k--','linew',2);
	plot(GSLx(subi,3),GSLlessGICa(subi),'k--','linew',2);
	plot(GSLx(subi,3),GSLlessGICb(subi),'k--','linew',2);
	xlim([1860 2100]); ylim([-.2 .50]);
	title('smoothed zostoga with drift substitution (mean GIC)');
	pdfwrite('zostoga_driftsubstA');
end