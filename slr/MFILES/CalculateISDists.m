function [thetEAIS,thetWAIS,thetGIS,ARthetAIS,ARthetGIS,ISLastDecade]=CalculateISDists(ratesmatrix2100,LastDecadeGt,ARIS2090);

% Last updated by  Bob Kopp, robert-dot-kopp-at-rutgers-dot-edu, Sun Dec 29 19:37:02 EST 2013

% use Shepherd estimates for last decade (2000-2011), pooling APIS and WAIS
surfarea = 3.61e14; %m^2
densitysw = 1020; %kg/m^3
Gt2mm = -1e12/surfarea/densitysw * 1000;
ISLastDecade = LastDecadeGt * Gt2mm;

% fit ice sheet distributions from Bamber & Aspinall (2013)

i=3; [thetEAIS,fval]=FitLNDistributionQuantiles(ratesmatrix2100(i,2),ratesmatrix2100(i,3:5),ratesmatrix2100(i,1),[.05 .5 .95],1);
i=2; [thetWAIS,fval]=FitLNDistributionQuantiles(ratesmatrix2100(i,2),ratesmatrix2100(i,3:5),ratesmatrix2100(i,1),[.05 .5 .95],1);
i=1; [thetGIS,fval]=FitLNDistributionQuantiles(ratesmatrix2100(i,2),ratesmatrix2100(i,3:5),ratesmatrix2100(i,1),[.05 .5 .95],1);

% fit AR5 ice sheet distributions  -- these are 2081-2100 (i.e., 2090); third dimension corresponds to scens

ARISLastDecade=[ISLastDecade(1) ISLastDecade(2)+ISLastDecade(3)];
ARbase = (2090-1995)*ARISLastDecade;
ARISaddl0 = bsxfun(@minus,ARIS2090,ARbase');
ARISaccel0 = 2*ARISaddl0/(2090-2011)^2;

ARratesmatrix2100 = bsxfun(@plus,ARISaccel0*(2100-2011),ARISLastDecade');


for kk=1:size(ARratesmatrix2100,3)
	i=2; [ARthetAIS{kk},fval]=FitLNDistributionQuantiles(ARratesmatrix2100(i,1,kk),ARratesmatrix2100(i,:,kk),-3,[.5 .167 .833],1);
	i=1; [ARthetGIS{kk},fval]=FitLNDistributionQuantiles(ARratesmatrix2100(i,1,kk),ARratesmatrix2100(i,:,kk),0,[.5 .167 .833],1);
end

