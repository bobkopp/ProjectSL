% Last updated by  Bob Kopp, robert-dot-kopp-at-rutgers-dot-edu, Mon Feb 10 22:32:11 EST 2014

defval('GIAanchoryear',2005);

% define covariance functions we will use

dYears = @(years1,years2) abs(bsxfun(@minus,years1',years2));
dYears0 = @(years1,years2) (bsxfun(@minus,years1',years2));

angd = @(Lat0,Long0,lat,long) (180/pi)*(atan2( sqrt( (cosd(lat) .* sind(long-Long0)).^2 + (cosd(Lat0) .* sind(lat) - sind(Lat0) .* cosd(lat) .* cosd(long-Long0)).^2),(sind(Lat0) .* sind(lat) + cosd(Lat0) .* cosd(lat) .* cosd(long-Long0))));

dDist = @(x1,x2) angd(repmat(x1(:,1),1,size(x2,1)),repmat(x1(:,2),1,size(x2,1)),repmat(x2(:,1)',size(x1,1),1),repmat(x2(:,2)',size(x1,1),1))' + 1e6*(bsxfun(@plus,x1(:,1)',x2(:,1))>1000);

kMat1 = @(dx,thetas) thetas(1).^2 .* (1).*exp(-dx/thetas(2));

kMat3 = @(dx,thetas) thetas(1).^2 .* (1 + sqrt(3)*dx/thetas(2)).*exp(-sqrt(3)*dx/thetas(2));

kMat5 = @(dx,thetas) thetas(1).^2 .* (1 + (sqrt(5)*dx/thetas(2)).*(1 + sqrt(5)*dx/thetas(2)/3)).*exp(-sqrt(5)*dx/thetas(2));

kSE = @(dx,thetas) thetas(1).^2 * exp(-(dx.^2)/(2*thetas(2).^2));

kDELTA = @(dx,thetas) thetas(1).^2 .* (dx==0);

kDP = @(years1,years2,thetas) thetas(1).^2 * bsxfun(@times,(years1-GIAanchoryear)',(years2-GIAanchoryear));

kMatG = @(dx,thetas) thetas(1).^2 .* 2.^(1-thetas(3))./gamma(thetas(3)) .* (sqrt(2*thetas(3))*(dx+eps)/thetas(2)).^thetas(3) .* besselk(thetas(3),sqrt(2*thetas(3))*(dx+eps)/thetas(2));

kRQ = @(dx,thetas) thetas(1).^2 * (1 + dx.^2/(2*thetas(2)*thetas(3))).^-thetas(3);

kDELTAG = @(ad,thetas) kDELTA(ad,thetas) .* (ad<360);

kGEOG = @(ad,thetas) kMat5(ad,thetas) .* (ad<360);

%% Define composite covariance functions


cvfuncL =  @(t1,t2,dt1t2,thetas) kDP(t1,t2,thetas(1)) + kMatG(dt1t2,thetas(2:4));

cvfuncGLR = @(t1,t2,thetas,ad,dt1t2) cvfuncL(t1,t2,dt1t2,thetas(1:4)) + kDP(t1,t2,thetas(5)) .* (kDELTAG(ad,sqrt(thetas(7))) + kGEOG(ad,[sqrt(1-thetas(7)) thetas(6)])) + kMatG(dt1t2,thetas(8:10)) .* (kDELTAG(ad,sqrt(thetas(11))) + kGEOG(ad,[sqrt(1-thetas(11)) thetas(12)])) + kDELTAG(ad,sqrt(thetas(1)^2+thetas(5)^2)*50);

cvfuncS0 = @(ad,thetas) thetas(1).^2 * (kDELTAG(ad,sqrt(thetas(3))) + kGEOG(ad,[sqrt(1-thetas(3)) thetas(2)]));
	
cvfuncS = @(ad,thetas) thetas(1).^2 * (kGEOG(ad,[sqrt(thetas(3)) thetas(4)]) + kGEOG(ad,[sqrt(1-thetas(3)) thetas(2)]));
