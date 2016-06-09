% outputSLRProjections_map_LSL
% Last updated by Robert Kopp, robert-dot-kopp-at-rutgers-dot-edu, Tue May 27 19:28:53 EDT 2014

%% 

% maps
kk=1;
subyear=find(targyears==2100);


for doNA=[0 1 2]
    clear subsites;
    if doNA==1
        labl2='NA';
        subsites={[1:length(targcoastlines)]};
    elseif doNA==2
        labl2='US';
        uscl={[823 940 960],[820 821],[760]};
        for ppp=1:length(uscl);
            subsites{ppp}=find(ismember(targcoastlines,uscl{ppp}));
        end
    else
        labl2='';
        subsites={[1:length(targcoastlines)]};
    end
    
    clf;
    if doNA==1
        worldmap('united states');
        setm(gca,'parallellabel','off','meridianlabel','off','flinewidth',1);
        ax=gca;
        geoshow('landareas.shp','facecolor','none');
    elseif doNA==2
        ax = usamap('all');
        set(ax, 'Visible', 'off')
        states = shaperead('usastatelo', 'UseGeoCoords', true);
        names = {states.Name};
        indexHawaii = strcmp('Hawaii',names);
        indexAlaska = strcmp('Alaska',names);
        indexConus = 1:numel(states);
        indexConus(indexHawaii|indexAlaska) = []; 
        stateColor = 'none';
        geoshow(ax(1), states(indexConus),  'FaceColor', stateColor)
        geoshow(ax(2), states(indexAlaska), 'FaceColor', stateColor)
        geoshow(ax(3), states(indexHawaii), 'FaceColor', stateColor)
        for k = 1:3
            setm(ax(k), 'Frame', 'off', 'Grid', 'off',...
                        'ParallelLabel', 'off', 'MeridianLabel', 'off')
        end
        posa=get(ax(2),'Position');
        posb1(1:2)=posa(1:2)+posa(3:4)/2.7;
        posb1(3:4)=posa(3:4)*.5;
        set(ax(2),'Position',posb1);
        posa=get(ax(3),'Position');
        %        posb2(1:2)=posa(1:2)+posa(3:4)/2.5;
        posb2(1:2)=posb1(1:2);
        posb2(1)=posb2(1)+.8*posb1(3);
        posb2(2)=posb2(2)+.3*posb1(4);
        posb2(3:4)=posa(3:4)*.7;
        set(ax(3),'Position',posb2);
        
    else
        worldmap('world');
        setm(gca,'parallellabel','off','meridianlabel','off','flinewidth',1);
        ax=gca;
        geoshow('landareas.shp','facecolor','none');
    end


    sb=find(quantlevs==0.5);
    for ppp=1:length(ax)
        hp=scatterm(ax(ppp),sitecoords(subsites{ppp},1),sitecoords(subsites{ppp},2),15,squeeze(quanttotlocrise(sb,subyear,subsites{ppp},kk))/1000,'filled');
        caxis(ax(ppp),[0 1.5]);
    end
    
    title({'Median projection: RCP 8.5',['GSL = ' sprintf('%0.2f',quanttotrise(sb,subyear)/1000) ' m']});
    hcl=colorbar('peer',ax(1)); set(hcl,'ytick',0:.25:1.5,'ytickl',{'<0','0.25','0.50','0.75','1.00','1.25','>1.50'});
    pdfwrite(['medianprojections' labl2]);
    
    % version of median projections in feet
    if doNA==2
        delete(hp); clear hp;
        sb=find(quantlevs==0.5);
        for ppp=1:length(ax)
            hp=scatterm(ax(ppp),sitecoords(subsites{ppp},1),sitecoords(subsites{ppp},2),15,squeeze(quanttotlocrise(sb,subyear,subsites{ppp},kk))/25.4/12,'filled');
            caxis(ax(ppp),[0 5]);
        end
        
        title({'Median projection: RCP 8.5, 2100',['GSL = ' sprintf('%0.1f',quanttotrise(sb,subyear)/25.4/12) ' ft']});
        hcl=colorbar('peer',ax(1)); set(hcl,'ytick',0:0.5:5,'ytickl',{'<0.0','0.5','1.0','1.5','2.0','2.5','3.0','3.5','4.0','4.5','>5.0'});
        ylabel(hcl,'feet');
        pdfwrite(['medianprojections_feet' labl2]);    
    end        

    for ql=[0.05 0.95]
	sb=find(quantlevs==ql);
	delete(hp); clear hp;
        for ppp=1:length(ax)
            hp(ppp)=scatterm(ax(ppp),sitecoords(subsites{ppp},1),sitecoords(subsites{ppp},2),15,squeeze(quanttotlocrise(sb,subyear,subsites{ppp},kk))/1000,'filled');
        end
        
	title({['Projection (' sprintf('%0.0f',ql*100) '%): RCP 8.5'],['GSL = ' sprintf('%0.2f',quanttotrise(sb,subyear,kk)/1000) ' m']}); colorbar('peer',ax(1)); 
	pdfwrite(['quantprojections' labl2 sprintf('%0.0f',ql*100)]);
    end


    delete(hp); clear hp;
    
    ql=[0.167 0.833];
    sb=find(quantlevs==ql(1));
    sb(2)=find(quantlevs==ql(2));
    rngs=[];
    for ppp=1:length(ax)
        rng=squeeze(quanttotlocrise(sb(2),subyear,subsites{ppp},kk)-quanttotlocrise(sb(1),subyear,subsites{ppp},kk))/1000;
        hp(ppp)=scatterm(ax(ppp),sitecoords(subsites{ppp},1),sitecoords(subsites{ppp},2),15,rng,'filled');
        rngs=union(rngs,rng);
    end
    
    title({['Projection (' sprintf('%0.0f%%-%0.0f%%',ql*100) '): RCP 8.5']}); colorbar('peer',ax(1)); 
    for ppp=1:length(ax)
        caxis(ax(ppp),quantile(rngs,[.02 .98]));
    end
    
    pdfwrite(['rangeprojections' labl2]);
    
    delete(hp);
    clear hp hp2;
    for ppp=1:length(ax)
        sub1=find(2*rateprojssd<abs(rateprojs));
        sub2=find(2*rateprojssd>=abs(rateprojs));
        sub1=intersect(sub1,subsites{ppp}); sub2=intersect(sub2,subsites{ppp});
        hp(ppp)=scatterm(ax(ppp),sitecoords(sub2,1),sitecoords(sub2,2),8,rateprojs(sub2)); hold on;
        hp2(ppp)=scatterm(ax(ppp),sitecoords(sub1,1),sitecoords(sub1,2),12,rateprojs(sub1),'filled'); hold on;
        caxis(ax(ppp),[-6 6]);
    end
    
    title('Background rate (mm/y)'); hc=colorbar('peer',ax(1));
    set(hc,'Ytick',[-6:2:6],'ytickl',{'<-6','-4','-2','0','2','4','>6'});longticks(hc,2);
    pdfwrite(['GIAmap' labl2]);
    
    delete(hp);
    delete(hp2);

    clear hp hp2;
    for ppp=1:length(ax)
        sub1=find(2*rateprojssd<abs(rateprojs));
        sub2=find(2*rateprojssd>=abs(rateprojs));
        sub1=intersect(sub1,subsites{ppp}); sub2=intersect(sub2,subsites{ppp});
        hp(ppp)=scatterm(ax(ppp),sitecoords(sub2,1),sitecoords(sub2,2),8,rateprojssd(sub2)); hold on;
        hp2(ppp)=scatterm(ax(ppp),sitecoords(sub1,1),sitecoords(sub1,2),12,rateprojssd(sub1),'filled'); hold on;
        caxis(ax(ppp),[0 2]);
    end
    
    title('Std. Dev. of Background rate (mm/y)'); hc=colorbar('peer',ax(1)); 
    set(hc,'Ytick',[0:2],'ytickl',{'0','1','>2'});longticks(hc,2);
    pdfwrite(['GIAsdmap' labl2]);

    delete(hp2); clear hp2;
    
    for nn=4
        delete(hp); clear hp;
        
	sb=find(quantlevs==0.5);
        us=[];
        for ppp=1:length(ax)
            u=squeeze(quantlocscalefactors(sb,nn,subyear,subsites{ppp},kk));
            hp(ppp)=scatterm(ax(ppp),sitecoords(subsites{ppp},1),sitecoords(subsites{ppp},2),15,u,'filled');
            us=union(us,u);
        end
        
	title({['Median scale factor: RCP 8.5, 2100']}); colorbar('peer',ax(1)); 
        for ppp=1:length(ax)
            caxis(ax(ppp),quantile(us,[.05 .95]));
        end
        
	pdfwrite(['locscalefactormedian' num2str(nn) labl2]);

        delete(hp); clear hp;
        
	sb=find(round(quantlevs*100)==17);
	sb(2)=find(round(quantlevs*100)==83);
        us=[];
        for ppp=1:length(ax)
            u=diff(squeeze(quantlocscalefactors(sb,nn,subyear,subsites{ppp},kk)),[],1);
            hp(ppp)=scatterm(ax(ppp),sitecoords(subsites{ppp},1),sitecoords(subsites{ppp},2),15,u,'filled');
            us=union(u,us);
        end
        
        for ppp=1:length(ax)
            caxis(ax(ppp),quantile(us,[.05 .95]));
        end
        
	title({['Scale factor (17th-83rd %ile range): RCP 8.5, 2100']}); colorbar('peer',ax(1)); 
	pdfwrite(['locscalefactorrange' num2str(nn) labl2]);

	for ql=[0.05 0.95]
            delete(hp); clear hp;
            sb=find(quantlevs==ql);
            us=[];
            for ppp=1:length(ax)
                u=squeeze(quantlocscalefactors(sb,nn,subyear,subsites{ppp},kk));
                hp(ppp)=scatterm(ax(ppp),sitecoords(subsites{ppp},1),sitecoords(subsites{ppp},2),15,u,'filled');
                us=union(u,us);
            end
            
            for ppp=1:length(ax)
                caxis(ax(ppp),quantile(u,[.05 .95]));
            end
            
            title({['Scale factor (' sprintf('%0.0f',ql*100) '%): RCP 8.5, 2100']}); colorbar('peer',ax(1)); 
            pdfwrite(['locscalefactorquant' num2str(nn) labl2 '_' sprintf('%0.0f',ql*100)]);
         end
 

    end
end



