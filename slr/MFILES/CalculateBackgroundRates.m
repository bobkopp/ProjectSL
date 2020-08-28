function [rateprojs,rateprojssd,rateprojs0,targcoord,rateGIAproj,priorsd,thetGLR,nearest,finescale]=CalculateBackgroundRates(coastlineset,targcoord,finescale,doregression,PARAMDIR,IFILES,regionalonly,yearrange,psmsldir,gslfile)

% [rateprojs,rateprojssd,rateprojs0,targcoord,rateGIAproj,priorsd,thetGLR,nearest,finescale]=CalculateBackgroundRates(coastlineset,targcoord,finescale,doregression,ROOTDIR,[regionalonly],[yearrange])
%
% set targcoord to -1 to use tide gauge locations; otherwise provide a list of sites (default will be 1 degree grid)
%
% set regionalonly to 1 to drop local component
% set regionalonly to 2 to not incorporate global sea-level curve
% set regionalonly to 3 to not mask out global sea level
%
% Last updated by Robert Kopp, robert-dot-kopp-at-rutgers-dot-edu, 2020-08-28 12:25:31 -0400

defval('IFILES','IFILES/');
defval('PARAMDIR','PARAMS/');
defval('COASTLINESFILE',fullfile(PARAMDIR,'coastlines.txt'));
coastlines=importdata(COASTLINESFILE);

defval('coastlineset',1:size(coastlines.data,1));
defval('GIAFILE',fullfile(IFILES,'dsea250.1grid.ICE5Gv1.3_VM2_L90_2012.nc'));
defval('WORKDIR',PARAMDIR);
defval('targcoord',[]);
defval('finescale',[]);
defval('doregression',1);
defval('regionalonly',0);
defval('yearrange',[1900 2005]);
if length(targcoord)==1
	doself=(targcoord==-1);
else
	doself=0;
end

rateGIAproj=[];
priorsd=[];
thetGLR=[];

% initialize

GIAanchoryear=yearrange(2);
GPSLDefineCovFuncs;

if length(targcoord)==0
	targlat=-80:90;
	targlong=-179:180;
	[TARGLAT,TARGLONG]=meshgrid(targlat,targlong);
	targcoord=[TARGLAT(:) TARGLONG(:)];
end

giamodel.gia=ncread(GIAFILE,'Dsea_250');
giamodel.lat=ncread(GIAFILE,'Lat');
giamodel.long=ncread(GIAFILE,'Lon');

ICE5Gin=giamodel.gia;
ICE5Glat=giamodel.lat;
ICE5Glon=giamodel.long;
sub=find(ICE5Glon>180); ICE5Glon(sub)=ICE5Glon(sub)-360;
[ICE5Glon,si]=sort(ICE5Glon); ICE5Gin=ICE5Gin(si,:);

% first load all the coordinates into memory
coastlinenames=coastlines.textdata(2:end);
coastlineset=unique(coastlineset);
for i=coastlineset
	disp(['Loading ' coastlinenames{i} '...']);

    defval('psmsldir',fullfile(IFILES,'rlr_annual'));
    defval('gslfile',fullfile(IFILES,'CSIRO_Recons_gmsl_yr_2011.csv'));
	
    [X1{i},Y{i},dY{i},regions{i},regionsu{i},sitenames{i},sitecoords{i},sitelen{i}]=ReadPSMSLData(coastlines.data(i,1),coastlines.data(i,2),15,psmsldir,gslfile,[],coastlines.data(i,3:end));

    if regionalonly==2
        sub=find(regionsu{i}~=0);
        regionsu{i}=regionsu{i}(sub);
        sitenames{i}=sitenames{i}(sub);
        sitecoords{i}=sitecoords{i}(sub,:);
        sitelen{i}=sitelen{i}(sub);
        
        sub=find(regions{i}~=0);
        regions{i}=regions{i}(sub);
        X1{i}=X1{i}(sub,:);
        Y{i}=Y{i}(sub);
        dY{i}=dY{i}(sub);
    end
    
	sub = find(regionsu{i}>0);
	if length(sub)>0
		if length(sub)>5
			minlength(i)=quantile(sitelen{i}(sub),.5);
		else
			minlength(i)=15;
		end
	else
		minlength(i)=0;
	end
	Nsites(i)=length(sub);
	
	
	
	% subtract GIA
	Y0=Y;
    [Y2{i},GIAproju{i},GIAproj{i}] = SubtractGIAfromTG(Y{i},ICE5Glat,ICE5Glon,ICE5Gin,sitecoords{i},regionsu{i},X1{i},regions{i},yearrange(2));

	
end

targregions=[];
if targcoord==-1
	targcoord=[];
	for i=coastlineset
		targcoord=[targcoord ; sitecoords{i}];
		targregions=[targregions; regionsu{i}];
	end
end

for jj=1:length(targcoord)
	if mod(jj,500)==0
		disp(['Computing distances: ' num2str(jj) '/' num2str(length(targcoord))]);
	end
	for kk=1:length(coastlineset)
		diststemp(kk) = min(dDist(targcoord(jj,:),sitecoords{coastlineset(kk)}));
	end
	m=min(diststemp);
	mi=find(diststemp==m);
	[mx,mxi]=max(Nsites(mi));
	mi=mi(mxi);
	nearest(jj)=coastlineset(mi);
end

%figure;
%scatter(mod(targcoord(:,2),360),targcoord(:,1),4,nearest,'filled');
%plotcont;
%title('Domains');
%pdfwrite('map_domains');

rateprojs0 = NaN*ones(size(targcoord,1),1);
rateprojs = NaN*ones(size(targcoord,1),1);
rateprojssd = NaN*ones(size(targcoord,1),1);

if doregression
	clear thetGLR priorsd;
	
	jjj=1;
	for i=coastlineset
		% Regress for slopes from tide gauge data
	    if sum(nearest==i)
            disp(['Regressing ' coastlinenames{i} '...']);
        
            fid=fopen(fullfile(WORKDIR,['SLModel_' coastlines.textdata{i+1} '.txt']),'r');
            thetGLRtemp=fscanf(fid,'%f');
            thetGLR(i,:)=thetGLRtemp';
            fclose(fid); 
    
            priorsd(i)=thetGLR(i,5);

            noiseMasks = ones(1,size(thetGLR,2));
            if regionalonly ~= 3
                noiseMasks(1,[1 2 8]) = 0; % no red noise or global
            elseif regionalonly == 3
                noiseMasks(1,[2 8])=0; % no red noise
            end
            

            % subtract GIA model

            testyears = [yearrange]';
            testlat = repmat(sitecoords{i}(:,1)',length(testyears),1);
            testlong = repmat(sitecoords{i}(:,2)',length(testyears),1);
            testX = [testlat(:) testlong(:) repmat(testyears,size(sitecoords{i},1),1)];
            testt=testX(:,3);
            testRegions = repmat(regionsu{i}',length(testyears),1); testRegions = testRegions(:);

            testGIAproj=zeros(size(testX,1),1);
            for ii=1:length(regionsu{i})
                sub=find(testRegions==regionsu{i}(ii));
                testGIAproj(sub)=GIAproju{i}(ii).*(testX(sub,3)-yearrange(2));
            end
    
            clear dx1x1 dx1x2 dx2x2 dt1t1 dt1t2 dt2t2;
            for jj=1:size(X1{i},1)
                if mod(jj,100)==0
                    disp(['Computing distance matrix 1: ' num2str(jj) '/' num2str(size(X1{i},1))]);
                end
                dx1x1(:,jj)=dDist(X1{i}(jj,:),X1{i});
                dx1x2(:,jj)=dDist(X1{i}(jj,:),testX);
                dt1t1(:,jj)=dYears(X1{i}(jj,3),X1{i}(:,3));
                dt1t2(:,jj)=dYears(X1{i}(jj,3),testt);
            end
            for jj=1:size(testX,1)
                if mod(jj,100)==0
                    disp(['Computing distance matrix 2: ' num2str(jj) '/' num2str(size(testX,1))]);
                end
                dx2x2(:,jj)=dDist(testX(jj,:),testX);
                dt2t2(:,jj)=dYears(testX(jj,3),testt);
            end

            traincv = @(thetas) cvfuncGLR(X1{i}(:,3),X1{i}(:,3),thetas,dx1x1,dt1t1) + diag(dY{i}.^2);
            testcv = @(thetas) cvfuncGLR(X1{i}(:,3),testX(:,3),thetas,dx1x2,dt1t2);
            testcv2 = @(thetas) cvfuncGLR(testX(:,3),testX(:,3),thetas,dx2x2,dt2t2);
            trcv = traincv(thetGLR(i,:));

            [f2,V2,logp,~,~,~,invcv] = GaussianProcessRegression([],Y2{i},[],trcv,testcv(thetGLR(i,:).*noiseMasks(1,:))',testcv2(thetGLR(i,:).*noiseMasks(1,:)));
            sd2=sqrt(diag(V2));
    
            clear rate ratesd;
            for jj=1:length(regionsu{i})
                sub=find(testRegions==regionsu{i}(jj));
                M = zeros(1,length(sub)); M(1)=-1; M(end)=1; 
                M=M/(testX(sub(end),3)-testX(sub(1),3));
                rate(jj) = M*f2(sub) + GIAproju{i}(jj);
                ratesd(jj)=sqrt(M*V2(sub,sub)*M');
            end
            rate=rate(:); ratesd=ratesd(:);

            % Find length scale of finer component



            trainsub=find(sitecoords{i}(:,1)<1e3);
            dx0x0 = dDist(sitecoords{i}(trainsub,:),sitecoords{i}(trainsub,:));
            traincvS = @(thetas) cvfuncS(dx0x0,thetas) + diag(ratesd(trainsub).^2);
            testcvS = @(thetas) cvfuncS(dx0x0,thetas);
            thetS = [thetGLR(i,5:7) thetGLR(i,6)*.1];
            lbS = .01; ubS = thetGLR(i,6);
            if length(finescale)<jjj
                disp('Finding length scale of finer component...');
                if thetS(3)>1e-4
                    thetS(end) = SLGPOptimize(rate(trainsub),@(x) traincvS([thetS(1:3) x]),thetS(end),lbS,ubS,1);
                end
                finescale(jjj)=thetS(end);
            else
                thetS(end)=finescale(jjj);
            end
    
            % Regress for slopes everywhere
            trcv = traincvS(thetS);

            [~,~,~,~,~,~,invcv] = GaussianProcessRegression([],rate(trainsub)-GIAproju{i}(trainsub),[],trcv,testcvS(thetS)',testcvS(thetS));

            sub=find(nearest==i);
            clear dx0x2 dx2x2 rateproj rateprojsd rateGIAproj;
            docoords=targcoord;
            docoords(:,1)=min(max(ICE5Glat),max(min(ICE5Glat),docoords(:,1)));
            docoords(:,2)=min(max(ICE5Glon),max(min(ICE5Glon),docoords(:,2)));
            for jj=1:length(sub)
                if mod(jj,500)==0
                    disp(['Regressing slopes for ' coastlinenames{i} ': ' num2str(jj) '/' num2str(length(sub))]);
                end
                dx0x2 = dDist(sitecoords{i}(trainsub,:),targcoord(sub(jj),:));
                dx2x2 = 0;
                nm = ones(size(thetS));
                if regionalonly==1
                    nm(3)=0;
                    nm(1)=sqrt(1-thetS(3));
                end
                testcvS = cvfuncS(dx0x2,thetS.*nm);
                testcv2S = cvfuncS(dx2x2,thetS.*nm);
                [rateproj(jj),rateprojV] = GaussianProcessRegression([],rate(trainsub)-GIAproju{i}(trainsub),[],trcv,testcvS',testcv2S,invcv);
                rateprojsd(jj)=sqrt(rateprojV);
                rateGIAproj(jj)=interp2(ICE5Glat,ICE5Glon,ICE5Gin,docoords(sub(jj),1),docoords(sub(jj),2));

            end
            rateprojs0(sub) = rateproj;
            rateprojs(sub) = rateproj + rateGIAproj ;
            rateprojssd(sub) = rateprojsd;
		end
		jjj=jjj+1;
	end

	if doself
		nearest=targregions;
	end
end

% plot rates where standard deviation is less than 2/3 prior
%
%rateprojs2 = rateprojs;
%for i=coastlineset
%	sub=find((nearest(:)==i).*(rateprojssd>(.67*thetGLR(i,5))));
%	rateprojs2(sub)=NaN;
%end
%figure;
%sub=find(~isnan(rateprojs2));
%scatter(mod(targcoord(sub,2),360),targcoord(sub,1),2,rateprojs2(sub),'filled');
%[~,handl]=plotcont; hold on; set(handl,'linew',.5);
%axis tight;
%colorbar;
%caxis([-5 5]);
%title('Background linear rates, data-constrained (mm/y)');
%pdfwrite('map_linrates');
%
%clf;
%sub=find(~isnan(rateprojs));
%scatter(mod(targcoord(sub,2),360),targcoord(sub,1),2,rateprojs(sub),'filled');
%[~,handl]=plotcont; hold on; set(handl,'linew',.5);
%axis tight;
%colorbar;
%caxis([-5 5]);
%title('Background linear rates (mm/y)');
%pdfwrite('map_linrates2');
%saveas(gcf,'map_linrates2','png');
%
%clf;
%scatter(mod(targcoord(sub,2),360),targcoord(sub,1),2,rateprojs2(sub) - (rateprojs(sub)-rateprojs0(sub)),'filled');
%[~,handl]=plotcont; hold on; set(handl,'linew',.5);
%axis tight;
%colorbar;
%caxis([-3 3]);
%title('Deviation from ICE-5G VM2-90 (mm/y)');
%pdfwrite('map_ice5gdev');
%%%