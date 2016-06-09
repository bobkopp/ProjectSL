function [ISaccelsamps,ARISaccelsamps,lastseed] = SampleISDists(thetEAIS,thetWAIS,thetGIS,ARthetAIS,ARthetGIS,ISLastDecade,BAcorrIS,ARcorrIS,seeds)

% Last updated by  Bob Kopp, robert-dot-kopp-at-rutgers-dot-edu, Sun Dec 29 19:38:51 EST 2013

lastseed=0;

% sample ice sheet distribution from Bamber & Aspinall (2013)

sigmas=[thetGIS(3)  thetWAIS(3) thetEAIS(3)];
mus = [thetGIS(2) thetWAIS(2) thetEAIS(2)];
offsets = [thetGIS(1) thetWAIS(1) thetEAIS(1)];
covIS = diag(sigmas)*BAcorrIS*diag(sigmas);

T=cholcov(covIS);
sampeps = [seeds(lastseed+[1:length(mus)],:)'*T];
sampISrates2100 = bsxfun(@plus,offsets,exp(bsxfun(@plus,sampeps,mus)));
lastseed=lastseed+length(mus);

ISaccelsamps = bsxfun(@minus,sampISrates2100,ISLastDecade)/(2100-2011);

% sample AR5 ice sheet distribution
ARISLastDecade=[ISLastDecade(1) ISLastDecade(2)+ISLastDecade(3)];

for kk=1:length(ARthetAIS)
	ARsigmas=[ARthetGIS{kk}(3) ARthetAIS{kk}(3)];
	ARmus=[ARthetGIS{kk}(2) ARthetAIS{kk}(2)];
	ARoffsets=[ARthetGIS{kk}(1) ARthetAIS{kk}(1)];
	ARcovIS = diag(ARsigmas)*ARcorrIS*diag(ARsigmas);

	T=cholcov(ARcovIS);
	sampeps = [seeds(lastseed+[1:length(ARmus)],:)'*T];
	ARsampISrates2100 = bsxfun(@plus,ARoffsets,exp(bsxfun(@plus,sampeps,ARmus)));
	ARISaccelsamps(:,:,kk) = bsxfun(@minus,ARsampISrates2100,ARISLastDecade)/(2100-2011);
end
lastseed=lastseed+length(ARmus);