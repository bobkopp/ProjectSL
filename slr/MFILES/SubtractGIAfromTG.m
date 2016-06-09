function [Y2,GIAproju,GIAproj] = SubtractGIAfromTG(Y,ICE5Glat,ICE5Glon,ICE5Gin,sitecoords,regionsu,X1,regions,GIAanchoryear)

% [Y2,GIAproju,GIAproj] = SubtractGIAfromTG(Y,ICE5Glat,ICE5Glon,ICE5Gin,sitecoords,regionsu,X1,regions,GIAanchoryear)
%
% Last updated by  Bob Kopp, robert-dot-kopp-at-rutgers-dot-edu, Mon Feb 10 21:46:34 EST 2014

GIAproju=zeros(size(regionsu));
GIAproj=zeros(size(Y));
docoords=sitecoords;
docoords(:,1)=min(max(ICE5Glat),max(min(ICE5Glat),docoords(:,1)));
docoords(:,2)=min(max(ICE5Glon),max(min(ICE5Glon),docoords(:,2)));

for i=1:length(GIAproju)
	if regionsu(i)>0
		GIAproju(i)=interp2(ICE5Glat,ICE5Glon,ICE5Gin,docoords(i,1),docoords(i,2));
		sub=find(regions==regionsu(i));
		GIAproj(sub)=GIAproju(i).*(X1(sub,3)-GIAanchoryear);
	end
end
Y2=Y-GIAproj;