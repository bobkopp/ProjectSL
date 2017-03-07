% Last updated by Robert Kopp, robert-dot-kopp-at-rutgers-dot-edu, Fri Jan 13 09:11:34 EST 2017

%%%%%%%%%%%%%%%%%%

angd = @(Lat0,Long0,lat,long) (180/pi)*(atan2( sqrt( (cosd(lat) .* sind(long-Long0)).^2 + (cosd(Lat0) .* sind(lat) - sind(Lat0) .* cosd(lat) .* cosd(long-Long0)).^2),(sind(Lat0) .* sind(lat) + cosd(Lat0) .* cosd(lat) .* cosd(long-Long0))));


% set up grid

GX2=[]; GY2=[]; siteids2=[]; sitenames2={};
TGregions=[]; TGsitenames=[]; TGsitecoords=[]; TGsitelen=[];

[X1,Y,dY,regions,regionsu,sitenames,sitecoords,sitelen]=ReadPSMSLData(0,10000,15,psmsldir,gslfile,[],[]);

c=load('coastlines');
[u,ui]=unique(1e4*round(c.coastlat/gridsz)+round(c.coastlon/gridsz));
GY=round(c.coastlat(ui)/gridsz)*gridsz;
GX=round(c.coastlon(ui)/gridsz)*gridsz;


     sub=find(sitecoords(:,1)<1000);

    [Dx1,Dx2]=meshgrid(GX(:),sitecoords(sub,2));
    [Dy1,Dy2]=meshgrid(GY(:),sitecoords(sub,1));
    dist=angd(Dy1,Dx1,Dy2,Dx2)';
    clear Dx1 Dx2 Dy1 Dy2;
    mdist=min(dist,[],2);
    subdo=find(mdist<distcutoff);

    GXa=[GX(subdo);sitecoords(sub,2)];
    GYa=[GY(subdo);sitecoords(sub,1)];
    siteids2=[siteids2 ; 1e9+1e4*[round((90-GY(subdo))*10)]+round(mod(GX(subdo),360)*10) ; regionsu(sub)];
    for lll=1:length(subdo)
        sitenames2{end+1}=sprintf('grid_%0.1f_%0.1f',[GY(subdo(lll)) mod(GX(subdo(lll)),360)]);
    end
    sitenames2={sitenames2{:},sitenames{sub}};

    GX2=[GX2; GXa]; GY2=[GY2; GYa];
    
    TGregions=[TGregions ; regionsu]; TGsitenames=[TGsitenames sitenames]; TGsitecoords=[TGsitecoords ; sitecoords]; TGsitelen=[TGsitelen ; sitelen];

%%

[targregions,ui]=unique(siteids2);
targsitecoords=[GY2(ui) GX2(ui)];
targregionnames=sitenames2(ui);

[Dx1,Dx2]=meshgrid(targsitecoords(:,2),TGsitecoords(:,2));
[Dy1,Dy2]=meshgrid(targsitecoords(:,1),TGsitecoords(:,1));
dist=angd(Dy1,Dx1,Dy2,Dx2)';
clear Dx1 Dx2 Dy1 Dy2;
[m,mi]=min(dist,[],2);
nearestTG=TGregions(mi);
[~,nearestTGmap]=ismember(nearestTG,TGregions);

