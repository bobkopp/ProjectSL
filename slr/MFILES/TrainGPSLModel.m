function [thetGLR,cvfuncGLR,thetGLRA] = TrainGPSLModel(X1,Y,dY,regions,regionsu,sitenames,sitecoords,sitelen,minlength,giamodel,thetG,dospatial,Nlongest)

% [thetGLR,cvfuncGLR,thetGLRA] = TrainGPSLModel(X1,Y,dY,regions,regionsu,sitenames,sitecoords,sitelen,minlength,giamodel,thetG,dospatial,Nlongest)
%
% Last updated by  Bob Kopp, robert-dot-kopp-at-rutgers-dot-edu, Mon Feb 24 10:42:04 EST 2014

    defval('minlength',50);
    defval('timesc',[]);
    defval('thetG',[]);
    defval('dospatial',1);
    defval('Nlongest',[15 5]);
    
    if length(Nlongest)<2
        Nlongest=[Nlongest Nlongest];
    end

    GPSLDefineCovFuncs;
	angdists = dDist(X1,X1);
	dt1t2 = dYears(X1(:,3),X1(:,3));
	
    % subtract out GIA model

    ICE5Gin=giamodel.gia;
    ICE5Glat=giamodel.lat;
    ICE5Glon=giamodel.long;
    sub=find(ICE5Glon>180); ICE5Glon(sub)=ICE5Glon(sub)-360;
    [ICE5Glon,si]=sort(ICE5Glon); ICE5Gin=ICE5Gin(si,:);

    Y0=Y;
    [Y,GIAproju,GIAproj] = SubtractGIAfromTG(Y,ICE5Glat,ICE5Glon,ICE5Gin,sitecoords,regionsu,X1,regions,GIAanchoryear);

    for i=1:length(regionsu)
        sub=find(regions==regionsu(i));
        [m,mi] = min(abs(X1(sub,3)-2005));
        baselinesub=[max(1,(mi-5)):min((mi+5),length(sub))];
        Y(sub) = Y(sub) - mean(Y(sub(baselinesub)));
    end

    %set up model starting point and boundaries

	% for slopes
	
	traincvS = @(ad,thetas,errs) cvfuncS0(ad,thetas) + diag(errs.^2);
	
	tluS = [
	4   0    100 % amplitude
	2   .1   15  % spatial scale
	.5  0     1 % nugget fraction
	];

	thetS=tluS(:,1)';
	lbS=tluS(:,2)';
	ubS=tluS(:,3)';
	
	% local only

	traincvL = @(thetas,X1,errs) cvfuncL(X1,X1,dYears(X1,X1),thetas) + diag(errs.^2);

	tluL = [
	1   0    100 % DP

	5   0   1e3   % MatG
	1   .05  300
	1.5  .5  5.5

	];

	thetL=tluL(:,1)';
	lbL=tluL(:,2)';
	ubL=tluL(:,3)';

	% global plus regional plus local
	
	traincvGLR = @(thetas,X1,errs,ad,dt1t2) cvfuncGLR(X1,X1,thetas,ad,dt1t2) + diag(errs.^2);

	tluGLR = [
	1   0    100 % DP global

	5   0   1e3   % MatG global
	1   .05  300
	1.5  .5  5.5

	5   .5    100; % DP
	2   1    8; % DP spatial scale
	.5  0    1; % DP nugget fraction
	
	50    .5   1e4; % MatG
	1    .5  300
	1.5  .5  5.5
    
  	.5  0    1; % oceanographic nugget fraction
	2   1    8; % oceanographic spatial scale
	
  	];

	thetGLR=tluGLR(:,1)';
	lbGLR=tluGLR(:,2)';
	ubGLR=tluGLR(:,3)';

	%%% PART 3: OPTIMIZE HYPERPARAMETERS

	% first optimize for global
	
	if length(thetG)<4
        trainsub1=find(regions==0);
        if length(trainsub1>0)
            disp('Optimizing GSL terms...')
            thetG = SLGPOptimize(Y(trainsub1),@(x) traincvL(x,X1(trainsub1,3),dY(trainsub1)),thetL,lbL,ubL,1);
            disp(['Optimized GSL terms:' sprintf('%0.3f ',thetG)]);
        end
    else
        thetG = thetG(1:4);
    end

	% now optimizing temporally

	disp('Optimizing temporally...');
	if length(thetG)>=4
    	thetGLR(1:4) = thetG;
    	subfixed=1:4;
    else
        subfixed=[];
    end
	thetGLR(7) = 1;
	thetGLR(11)=1;
	subfixed=[subfixed  6:7 11:12];
	
	subnotfixed=setdiff(1:length(thetGLR),subfixed);
	Mfixed=sparse(length(subnotfixed),length(thetGLR));
	for i=1:length(subnotfixed)
		Mfixed(i,subnotfixed(i))=1;
	end
	fixedvect = 0*thetGLR; fixedvect(subfixed)=thetGLR(subfixed);

	thetGLR0=thetGLR;
	
	[s,si]=sort(sitelen(:)',2,'descend');

	% fifteen longest records
	trainsub=find(ismember(regions,[0 ; regionsu(si(1:min(Nlongest(1),length(s))))]));

	for globoptim=[0 1]
		[thetGLR(subnotfixed)] = SLGPOptimize(Y(trainsub),@(x) traincvGLR(x*Mfixed+fixedvect,X1(trainsub,3),dY(trainsub),angdists(trainsub,trainsub),dt1t2(trainsub,trainsub)),thetGLR(subnotfixed),lbGLR(subnotfixed),ubGLR(subnotfixed),globoptim);
	end

    if dospatial
	
    	disp(['Optimized local terms, first iteration:' sprintf('%0.3f ',thetGLR(subnotfixed))]);
	
	% now optimizing linear term spatially

        disp('Optimizing linear term spatially...');

        noiseMask = ones(size(thetGLR));
        noiseMask([1 2 8])=0;
        trcv=traincvGLR(thetGLR,X1(:,3),dY(:),angdists(:,:),dt1t2);	
        testcv=traincvGLR(thetGLR.*noiseMask,X1(:,3),0*dY(:),angdists(:,:),dt1t2);
        [floclin,Vloclin] = GaussianProcessRegression([],Y,[],trcv,testcv,testcv);
        clear fslopelin sdslopelin;
        for i=1:length(regionsu)
            Mslope=zeros(1,length(floclin));
            sub=find(regions==regionsu(i));
            dt=X1(sub(end),3)-X1(sub(end-1),3);
            Mslope(sub(end))=1/dt; Mslope(sub(end-1))=-1/dt;
            sub2=find(abs(Mslope)>0);
            fslopelin(i)=Mslope(sub2)*floclin(sub2)+GIAproju(i);
            sdslopelin(i)=sqrt(Mslope(sub2)*Vloclin(sub2,sub2)*Mslope(sub2)');
        end
    
        adu = dDist(sitecoords,sitecoords);
        sub=find(regionsu>0);
        fslopelin=fslopelin(:);
        [thetS] = SLGPOptimize(fslopelin(sub),@(x) traincvS(adu(sub,sub),x,sdslopelin(sub)),thetS,lbS,ubS,1);

        disp(['Spatially optimized linear term:' sprintf('%0.3f ',thetS)]);
    

        % now optimize Mat term locally
    
        disp('Optimizing Matern term locally...');
    
        thetGLR(5:7) = thetS;
        thetGLR(11)=1;
        subfixed=[1:4 5:7 11:12];
    
        subnotfixed=setdiff(1:length(thetGLR),subfixed);
        Mfixed=sparse(length(subnotfixed),length(thetGLR));
        for i=1:length(subnotfixed)
            Mfixed(i,subnotfixed(i))=1;
        end
        fixedvect = 0*thetGLR; fixedvect(subfixed)=thetGLR(subfixed);

        thetGLR0=thetGLR;
    
        [s,si]=sort(sitelen(:)',2,'descend');
        trainsub=find(ismember(regions,[0 ; regionsu(si(1:min(Nlongest(1),length(s))))]));

        for globoptim=[0 1]
            [thetGLR(subnotfixed)] = SLGPOptimize(Y(trainsub),@(x) traincvGLR(x*Mfixed+fixedvect,X1(trainsub,3),dY(trainsub),angdists(trainsub,trainsub),dt1t2(trainsub,trainsub)),thetGLR(subnotfixed),lbGLR(subnotfixed),ubGLR(subnotfixed),globoptim);
        end	


        disp(['Temporally optimized Matern term:' sprintf('%0.3f ',thetGLR(subnotfixed))]);

        % now optimizing Mat spatially

        disp('Optimizing Matern term spatially...');

        subfixed=[1:4 5:7 8:10];
    
        thetGLR(11)=.5;
        subnotfixed=setdiff(1:length(thetGLR),subfixed);
        Mfixed=sparse(length(subnotfixed),length(thetGLR));
        for i=1:length(subnotfixed)
            Mfixed(i,subnotfixed(i))=1;
        end
        fixedvect = 0*thetGLR; fixedvect(subfixed)=thetGLR(subfixed);

        thetGLR0=thetGLR;
    
        [s,si]=sort(sitelen(:)',2,'descend');
        trainsub=find(ismember(regions,[0 ; regionsu(si(1:min(Nlongest(2),length(s))))]));

        for globoptim=[0 1]
            [thetGLR(subnotfixed)] = SLGPOptimize(Y(trainsub),@(x) traincvGLR(x*Mfixed+fixedvect,X1(trainsub,3),dY(trainsub),angdists(trainsub,trainsub),dt1t2(trainsub,trainsub)),thetGLR(subnotfixed),lbGLR(subnotfixed),ubGLR(subnotfixed),globoptim);
        end

        if length(regions)<500
            Nruns=1;
        else
            Nruns=20;
        end

        for n=1:Nruns

            disp(n)
            thetGLR=thetGLR0;

            [s,si]=sort(sitelen(:)',2,'descend');
            trainsub=find(ismember(regions,[0 ; regionsu(si(1:min(Nlongest(2),length(s))))]));

            subexcl = setdiff(1:size(X1,1),find(X1(:,1)>1000));
            subexcl = setdiff(subexcl,trainsub);

            timedists = dYears(X1(:,3),X1(:,3));

            blocklength=30;
            Nnewsamps=500;
            trainsub=ConstructSubset(trainsub,blocklength,Nnewsamps,regions,subexcl,X1,angdists,timedists);

            if n==1
                globoptim=0;
            else
                globoptim=0;
            end

            [thetGLR(subnotfixed)] = SLGPOptimize(Y(trainsub),@(x) traincvGLR(x*Mfixed+fixedvect,X1(trainsub,3),dY(trainsub),angdists(trainsub,trainsub),dt1t2(trainsub,trainsub)),thetGLR(subnotfixed),lbGLR(subnotfixed),ubGLR(subnotfixed),globoptim);

            thetGLRA(n,:) = thetGLR;
            disp(sprintf('%0.3f ',thetGLR));

        end
        thetGLR=median(thetGLRA,1);

        disp(['Spatially optimized Matern term:' sprintf('%0.3f ',thetGLR(subnotfixed))]);
    else
        thetGLRA=thetGLR;
    end
    disp(['theta = ' sprintf('%0.3f ',thetGLR)]);
end

%%%%

function trainsubnew=ConstructSubset(trainsub,blocklength,Nnewsamps,regions,subexcl,X1,angdists,timedists)

	defval('blocklength',50);
	defval('Nnewsamps',750);

	trainsubnew=[trainsub];
	while length(trainsubnew)<min((length(trainsub)+Nnewsamps),length(regions))
		l0 = length(trainsubnew)-length(trainsub);
		%disp(l0)
		inactiveset=setdiff(subexcl,trainsubnew);

		% find geographically most distant
		[u,ui]=unique(X1(trainsubnew,2)*1000+X1(trainsubnew,1));

		[u2,ui2]=unique(X1(inactiveset,2)*1000+X1(inactiveset,1));

		[c,ci]=setdiff(u2,u);
		if length(ci)>0
			ui2=ui2(ci);
		end

		metric=mean(angdists(trainsubnew(ui),inactiveset(ui2)));
		metric=cumsum(metric)/sum(metric);
		mi = find(rand<=metric);mi=mi(1);

		%[m,mi] = max(mean(angdists(trainsubnew(ui),inactiveset(ui2)),1));
		sub=intersect(find(regions==regions(inactiveset(ui2(mi)))),inactiveset);

		if length(sub)<=blocklength
			trainsubnew = [trainsubnew(:) ; sub(:)];
		else
			% now find temporally most distance 
			%sub2 = sub(ceil(1+(blocklength-1)/2):floor(length(sub)-(blocklength-1)/2));

			sub2=sub;
			metric=mean(timedists(trainsubnew,sub2),1);
			metric=cumsum(metric)/sum(metric);
			mi = find(rand<=metric);mi=mi(1);
			%[m,mi]=max(mean(timedists(trainsubnew,sub2),1));
			q=min(length(sub2),ceil(mi+(blocklength-1)/2));
			sub3=sub2(max(1,q-(blocklength-1)):q);
			trainsubnew=[trainsubnew(:) ; sub3(:)];
		end
		trainsubnew=unique(trainsubnew);

		if length(trainsubnew)==l0
			inactiveset=setdiff(subexcl,trainsubnew);
			trainsubnew=[trainsubnew ; inactiveset(randsample(length(inactiveset),1))];
		end	
	end
	trainsub=unique(trainsubnew);
end
