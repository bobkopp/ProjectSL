function [TGtargcoords,fTG,VTG,sdTG,TGtargid,TGsiteid,sitenames,TGsitecoords,sitelen,cvfuncTG] = ReadDenoisedPSMSLData(cl1,cl2,thetTG,bedrocksiteids,thinyrs,minlen,psmsldir,gslfile,addlsites,exclusions, firstyear,targsites)

%
%
% Last updated by  Bob Kopp, robert-dot-kopp-at-rutgers-dot-edu, Sat Oct 12 18:48:20 EDT 2013

defval('cl1',960);
defval('cl2',960);
defval('minlen',[]);
defval('psmsldir',[]);
defval('gslfile',[]);
defval('bedrocksiteids',[]);
defval('thinyrs',5);
defval('addlsites',[]);
defval('exclusions',[]);
defval('firstyear',[]);
defval('targsites',[]);

% % 18 Apr: tide gauge
defval('thetTG',[1.6585 36.1312 1e3 0.05 0 1.0217 0.5970   28.6848   11.9728    1.2888    5.7005   26.6414    4.1376   28.2188    7.3826]);

disp('Reading data...');

[TGcoords,TGrsl,TGrslunc,TGid,TGsiteid,sitenames,TGsitecoords,sitelen]=ReadPSMSLData(cl1,cl2,minlen,psmsldir,gslfile,addlsites);
TGlat=TGcoords(:,1);
TGlong=TGcoords(:,2);
TGyears=TGcoords(:,3);
TGsitelat = TGsitecoords(:,1);

GPSLDefineCovFuncs;

bedrockMask = @(r1,r2) [repmat(sparse(~ismember(r1,bedrocksiteids)),1,length(r2)).*(repmat(sparse(~ismember(r2,bedrocksiteids))',length(r1),1))]';

cvfuncTG = @(t1,t2,dt1t2,thetas,ad,bedmask) kDP(t1,t2,1) .* (thetas(1).^2 + kDELTAG(ad,thetas(7)).*bedmask + kGEOG(ad,thetas(10:11))) + kRQ(dt1t2,[1 thetas(3:4)]) .* (thetas(2).^2 + kDELTAG(ad,thetas(8)) + kGEOG(ad,thetas(12:13))) + kMat1(dt1t2,[1 thetas(6)]).* (thetas(5).^2 + kDELTAG(ad,thetas(9)) + kGEOG(ad,thetas(14:15)));


%%%

sub=[];
TGtargcoords = [];
TGtargyears = [];
TGtargid = [];
if length(targsites)==0
	targsites=TGsiteid;
end
for i=1:length(targsites)
	sub1=find(TGid==targsites(i));
	if length(firstyear)==0
		trg = TGyears(sub1(1)):thinyrs:TGyears(sub1(end));
	else
		trg = firstyear:thinyrs:TGyears(sub1(end));
	end		
	
	sub2=find(TGsiteid==targsites(i)); sub2=sub2(1);
	TGtargcoords = [TGtargcoords ; repmat(TGsitecoords(sub2,:),length(trg),1)];
	TGtargyears = [TGtargyears ; trg(:)];
	TGtargid = [TGtargid ; repmat(targsites(i),length(trg),1)];
end
TGtargcoords = [TGtargcoords TGtargyears];

traincv=@(thet) cvfuncTG(TGyears,TGyears,dYears(TGyears,TGyears),thet,dDist(TGcoords,TGcoords),bedrockMask(TGid,TGid)) + diag(TGrslunc.^2);
testcv= @(thet) cvfuncTG(TGyears,TGtargyears,dYears(TGyears,TGtargyears),thet,dDist(TGcoords,TGtargcoords),bedrockMask(TGid,TGtargid));
testcv2= @(thet) cvfuncTG(TGtargyears,TGtargyears,dYears(TGtargyears,TGtargyears),thet,dDist(TGtargcoords,TGtargcoords),bedrockMask(TGtargid,TGtargid));

noiseMask=ones(size(thetTG));
noiseMask([5 9 14])=0; % no red noise

disp('Gaussian process smoothing...')

[fTG,VTG,logp] = GaussianProcessRegression([],TGrsl,[],traincv(thetTG),testcv(thetTG.*noiseMask)',testcv2(thetTG.*noiseMask));
sdTG=sqrt(diag(VTG));
