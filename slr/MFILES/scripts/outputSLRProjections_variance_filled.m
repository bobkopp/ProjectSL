% outputSLRProjections_variance_filled
% Last updated by Robert Kopp, robert-dot-kopp-at-rutgers-dot-edu, Tue May 27 16:04:34 EDT 2014

% variance plots

dopooledscens=0;
if dopooledscens
    colrs='rcbgmyk';
else
    colrs='rcbgmyk';
end
scenlabpretty={'RCP 8.5','RCP 6.0','RCP 4.5','RCP 2.6'};
subyears=find(targyears<=2100);

for kk=[1 4]

    % GSL
    linew=[1 1 1 1 1 1];

    subcomp={colAIS,colIS,[colIS colTE],[colGIC colIS colTE],1:size(samps,2)};
    %    labls={'AIS','IS','IS+TE','LI+TE',scenlabpretty{kk}};
    labls={'AIS','GIS','TE','GIC','LWS'};

    titl=['GSL variance (' scenlabpretty{kk} ')'];

    clf;
    clear hl hp;
    hp(1)=subplot(2,2,1);
    yrs=[targyears(subyears)];
    vlast=zeros(1,length(yrs));
    for i=1:length(labls)
        u=squeeze(sum(samps(:,subcomp{i},subyears,kk),2));
        vcur=[var(u)]/1e6;
        %        hl(i)=plot(yrs,[0 var(u)]/1e6,[colrs(i) '-'],'linew',linew(i)); hold on;
        hl(i)=patch([yrs yrs(end:-1:1)],[vcur vlast(end:-1:1)],colrs(i)); hold on;
        vlast=vcur;
    end

    clear vcur;
    if dopooledscens
        u0A=squeeze(sum(samps(:,subcomp{end},subyears,:),2));
        u0A=reshape(permute(u0A,[1 3 2]),size(u0A,1)*size(u0A,3),[]);
        for jj=1:size(u0A,2)
             sub=find(~isnan(u0A(:,jj)));
                 varu0A(jj)=var(u0A(sub,jj));
        end
        vcur=[varu0A]/1e6;
        hl(i+1)=patch([yrs yrs(end:-1:1)],[vcur vlast(end:-1:1)],colrs(i+1)); hold on;
    end
    
    title(titl);
    ylabel('m^2');
    if dopooledscens
        legend(hl(end:-1:1),'All RCPs',labls{end:-1:1},'Location','Northwest');
    else
        legend(hl(end:-1:1),labls{end:-1:1},'Location','Northwest');
    end
    longticks(gca);
    
    % fraction of variance plots
    % GSL

    %clf;
    hp(2)=subplot(2,2,2);
    clear hl;
    
    if ~dopooledscens
        varu0A=var(squeeze(sum(samps(:,subcomp{end},subyears,kk),2)));
    end
    yrs=[targyears(subyears)];
    vlast=zeros(1,length(yrs));
    for i=1:length(labls)
        u=squeeze(sum(samps(:,subcomp{i},subyears,kk),2));
        vcur=[var(u)./varu0A];
        hl(i)=patch([yrs yrs(end:-1:1)],[vcur vlast(end:-1:1)],colrs(i)); hold on;
        vlast=vcur;
    end
    if dopooledscens
        vcur=[varu0A./varu0A];
        hl(i)=patch([yrs yrs(end:-1:1)],[vcur vlast(end:-1:1)],colrs(i)); hold on;
        vlast=vcur;
    end
    
    title(titl);
    ylabel('Fraction of variance');
    %legend(labls{:},'All RCPs','Location','Northwest');
    longticks(gca);
    set(hp,'xlim',[2010 2100]);
    pdfwrite(['componentvar_' scens{kk}]);
 

    % and locally

    if colOD==colTE
        colTE2=[];
    else
        colTE2=colTE;
    end

    subcomp={colAIS,colIS,[colIS colTE2 colOD],[colGIC colIS colTE2 colOD],[colGIC colIS colLS colTE2 colOD],[colGIC colIS colLS colTE2 colOD colGIA]};
    labls={'AIS','GIS','Ocean','GIC','LWS','Bkgd'};
    titl=['variance (' scenlabpretty{kk} ')'];

    linew=[1 1 1 1 1 1 1];
    sitesub=[12 299 396  188 161 10 405 155 78 7 134 499];

    for jjj=1:length(sitesub)
        jj2=find(focussites==sitesub(jjj));
        if length(jj2)==1
            jj=find(targregions==sitesub(jjj));
            titl={targregionnames{jj},['variance (' scenlabpretty{kk} ')']};

            clf;
            clear hl hp;
            hp(1)=subplot(2,2,1);
            yrs=[targyears(subyears)];
            vlast=zeros(1,length(yrs));

            for i=1:length(labls)		
                u=squeeze(sum(sampsregionfocus{jj2,kk}(:,subcomp{i},subyears),2));
                vcur=[var(u)]/1e6;
                hl(i)=patch([yrs yrs(end:-1:1)],[vcur vlast(end:-1:1)],colrs(i)); hold on;
                %                hl(i)=plot(yrs,[0 var(u)]/1e6,[colrs(i) '-'],'linew',linew(i)); hold on;
                vlast=vcur;
            end

            if dopooledscens
                u0A = [];
                for kk2=1:size(sampsregionfocus,2)
                    u0A=[u0A ; squeeze(sum(sampsregionfocus{jj2,kk2}(:,subcomp{end},subyears),2))];
                end

                for jj=1:size(u0A,2)
                         sub=find(~isnan(u0A(:,jj)));
                         varu0A(jj)=var(u0A(sub,jj));
                end
                vcur=[varu0A]/1e6;
                %                hl(i+1)=plot(yrs,[0 varu0A]/1e6,[colrs(i+1) '-'],'linew',linew(i+1)); hold on;
                hl(i+1)=patch([yrs yrs(end:-1:1)],[vcur vlast(end:-1:1)],colrs(i+1)); hold on;
            end

            title(titl);
            ylabel('m^2');
            if dopooledscens
                legend(hl(end:-1:1),'Scen',labls{end:-1:1},'Location','Northwest');
            else
                legend(hl(end:-1:1),labls{end:-1:1},'Location','Northwest');
            end
             longticks(gca);
                    
            %
            %clf;
            hp(2)=subplot(2,2,2);
            yrs=[targyears(subyears)];
            vlast=zeros(1,length(yrs));
            clear hl;
            
            if ~dopooledscens
                varu0A=var(squeeze(sum(sampsregionfocus{jj2,kk}(:,subcomp{end},subyears),2)));
            end
            
            for i=1:length(labls)		
                u=squeeze(sum(sampsregionfocus{jj2,kk}(:,subcomp{i},subyears),2));
                vcur=[var(u)./varu0A];
                hl(i)=patch([yrs yrs(end:-1:1)],[vcur vlast(end:-1:1)],colrs(i)); hold on;
                %                hl(i)=plot([ targyears],,[colrs(i) '-'],'linew',linew(i)); hold on;
                vlast=vcur;
            end
            if dopooledscens
                vcur=[varu0A./varu0A];
                hl(i+1)=patch([yrs yrs(end:-1:1)],[vcur vlast(end:-1:1)],colrs(i+1)); hold on;
                %                hl(i+1)=plot([ targyears],[varu0A./varu0A],[colrs(i+1) '-'],'linew',linew(i+1)); hold on;
            end
            title(titl);
            ylabel('Fraction of variance');
            ylim([0 1]);
    %		xlim([2010 2100]); %ylim([0 1]);
    %		legend(labls{:},'All RCPs','Location','Northwest');
            longticks(gca);
            set(hp,'xlim',[2010 2100]);
            pdfwrite(['componentvar_' scens{kk} '_' num2str(sitesub(jjj))]);
 
        
        end
    end
end