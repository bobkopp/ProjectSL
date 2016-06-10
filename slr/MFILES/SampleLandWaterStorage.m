function [LWSsamps,GWDsamps,damsamps,thetdam0,mGWD,stdGWD,mGWDs,hp]=SampleLandWaterStorage(Nsamps,yrs,pophistfile,popscenfile,damfile,GWDfiles,dotriangular,doplot)

% [LWSsamps,GWDsamps,damsamps]=SampleLandWaterStorage([Nsamps],[yrs],[pophistfile],[popscenfile],[damfile],[GWDfiles],[doplot])
%
% Last updated by Robert Kopp, robert-dot-kopp-at-rutgers-dot-edu, Thu Jun 09 16:44:29 EDT 2016

defval('dotriangular',0);
defval('doplot',0);
defval('Nsamps',100);
defval('pophistfile','UNWPP2012 population historical.csv');
defval('damfile','Chao2008 groundwater impoundment.csv');
defval('GWDfiles',{'Konikow2011 GWD.csv','Wada2012 GWD.csv','Pokhrel2012 GWD.csv'});
defval('popscenfile','UNWPP2012 population projections.csv');
defval('hp',[]);

defval('GWDslopepcterr',.25);
defval('dampcterror',.25);
defval('yrs',2000:10:2200);
if length(GWDfiles)<3
    dotriangular=0;
end



% import historical population change
if ~isnumeric(pophistfile)
    dat=importdata(pophistfile);
    t=dat.data(:,1);
    pop=dat.data(:,2);
else
    t=pophistfile(:,1);
    pop=pophistfile(:,2);
end

t0=t; pop0=pop;
t=t(1:5:end);
pop=pop(1:5:end);

t2=(t(1:end-1)+t(2:end))/2;
pop2=(pop(1:end-1)+pop(2:end))/2;

% import dams

if ~isnumeric(damfile)
    dat=importdata(damfile);
    dams=interp1(dat.data(:,1),dat.data(:,2),t0);

    fitfunc = @(theta,x) theta(1)*erf((x-theta(2))/theta(3))+theta(4);
    fitter = @(theta) sum((dams-fitfunc(theta,pop0/1e6)).^2);
    tlu = [
        max(dams)/2 0 1000*max(dams)
        1 0 1e2
        1 .1 1e2
        max(dams)/2 -max(dams) max(dams)
          ];
    thet0=tlu(:,1);
    lb=tlu(:,2);
    ub=tlu(:,3);
    [thetdam0,fval]=fmincon(fitter,thet0,[],[],[],[],lb,ub,[],optimset('maxfunevals',2000,'display','off'));
else
    thetdam0 = damfile;
end

% import alternative groundwater depletion

if ~isnumeric(GWDfiles)
    for jj=1:length(GWDfiles)
        dat=importdata(GWDfiles{jj});
        [u,ui]=unique(dat.data(:,1));
        GWD=interp1(dat.data(ui,1),dat.data(ui,2),t);
        sub=find(~isnan(GWD));
        dGWD{jj}=diff(GWD(sub))./diff(t(sub));
        pop2GWD{jj}=(pop(sub(1:end-1))+pop(sub(2:end)))/2;
        mGWDs(jj) =pop2GWD{jj}\dGWD{jj};
        [a,b,c]=fileparts(GWDfiles{jj});
        GWDlabel{jj}=b(1:min(length(GWDfiles{jj}),4));
    end
    % assume errors of +/- 25% on slopes
    if ~dotriangular
        mGWD = mean(mGWDs);
        stdGWD = sqrt((std(mGWDs).^2+(GWDslopepcterr*mGWD).^2));
    else
        mGWD = median(mGWDs);
        stdGWD = [min(mGWDs) max(mGWDs)];
    end
else
    mGWD=GWDfiles(1);
    stdGWD=GWDfiles(2:end);
end


% now import population projections
if ~isnumeric(popscenfile)
    dat = importdata(popscenfile);
    popdata=dat.data;
else
    popdata=popscenfile;
end

popscenyr=popdata(:,1);
popscen=popdata(:,2:end);
[s,si]=sort(popscen(end,:));
popscen=popscen(:,si);
popscenids=linspace(0,1,size(popscen,2));

popscenyr0=popscenyr;
if min(popscenyr)>2000
    popscenyr=[ [2000:(min(popscenyr)-1)]' ; popscenyr];
    popscen=[repmat(interp1(t0,pop0,[2000:(min(popscenyr0)-1)])',1,3) ; popscen];
end
if max(popscenyr)==2100
    popscengrowthUN=[ % from UN(2008)
                    -0.75 -0.15 0.35
                    -0.60 -0.11 0.4
                    -0.48 -0.03 0.45
                    -0.38 0.03 0.51 ]/100+1;
    popscengrowth=interp1([0 .5 1],popscengrowthUN',popscenids)';
    for ii=1:size(popscengrowth,1)
        popscenyr=[popscenyr ; popscenyr(end)+[1:25]'];
        popscen=[popscen ;  bsxfun(@times,popscen(end,:),bsxfun(@power,popscengrowth(ii,:),[1:25]'))];
    end
end

if max(yrs)>max(popscenyr)
    popscenyr(end+1)=max(yrs);
    popscen(end+1,:)=popscen(end,:);
end


popdraw = @(x) interp1(popscenids,popscen',x)';
if ~dotriangular
    GWDdraw = @(x1,x2,x3) interp1(popscenyr,cumsum(popdraw(x1) * (mGWD+norminv(x2)*stdGWD)),yrs');
else
    GWDdraw = @(x1,x2,x3) interp1(popscenyr,cumsum(popdraw(x1) * interp1([0 .5 1],[stdGWD(1) mGWD stdGWD(2)],x2))*(1+norminv(x3)*GWDslopepcterr),yrs');
end


pop2000=pop0(find(t0==2000)); 
damdraw = @(x1,x2) interp1(popscenyr,-(fitfunc(thetdam0,bsxfun(@max,popdraw(x1),pop2000)/1e6)-fitfunc(thetdam0,pop2000/1e6))*(1+norminv(x2)*dampcterror),yrs');

if length(Nsamps)>1
    seeds=Nsamps;
    seeds0=Nsamps;
else
    seeds0=linspace(0,1,Nsamps+2); seeds0=seeds0(2:end-1);
end
if prod(size(seeds0))==length(seeds0)
    clear seeds;
    for jj=1:4
        seeds(jj,:)=seeds0(randperm(length(seeds0)));
    end
end

parfor ii=1:size(seeds,2)
    GWDsamps(:,ii) = GWDdraw(seeds(1,ii),seeds(2,ii),seeds(3,ii));
    damsamps(:,ii) = damdraw(seeds(1,ii),seeds(4,ii));
end

LWSsamps=GWDsamps+damsamps;

if doplot
        clf;
        
        hp(1)=subplot(2,2,1);
        plot(pop0/1e6,dams,'b','linew',2); hold on
        plot(pop0/1e6,fitfunc(thetdam0,pop0/1e6),'r');
        xlabel('Population (billions)');
        ylabel({'Impoundment (mm esl)'});
        
        hp(2)=subplot(2,2,2);
        colrs={'b','g','r','c','m','y'};
        for ii=1:length(pop2GWD)
            hpg(ii)=plot(pop2GWD{ii}/1e6,dGWD{ii},colrs{ii},'linew',2); hold on  
        end
        for ii=1:length(pop2GWD)
            plot(pop0/1e6,pop0*mGWDs(ii),'k--'); hold on  
        end        
        legend(hpg,GWDlabel,'Location','Northwest');
        xlabel('Population (billions)');
        ylabel('GWD (mm/y esl)');
        
        hp(3)=subplot(2,2,3);
        plot(popscenyr,popdraw(.5)/1e6,'linew',2); hold on
        plot(popscenyr,popdraw(.167)/1e6); hold on  
        plot(popscenyr,popdraw(.833)/1e6); hold on
        plot(popscenyr,popdraw(.05)/1e6,'--'); hold on  
        plot(popscenyr,popdraw(.95)/1e6,'--'); hold on
        plot(popscenyr,popdraw(.005)/1e6,':'); hold on  
        plot(popscenyr,popdraw(.995)/1e6,':'); hold on
        xlabel('Year');
        ylabel('Population (billions)');
        
        hp(4)=subplot(2,2,4);
        plot(yrs,quantile(LWSsamps',.5),'linew',2); hold on
        plot(yrs,quantile(LWSsamps',.167)); hold on  
        plot(yrs,quantile(LWSsamps',.833)); hold on
        plot(yrs,quantile(LWSsamps',.05),'--'); hold on  
        plot(yrs,quantile(LWSsamps',.95),'--'); hold on
        plot(yrs,quantile(LWSsamps',.005),':'); hold on  
        plot(yrs,quantile(LWSsamps',.995),':'); hold on
        xlabel('Year');
        ylabel('LWS (mm esl)');        
    end
end
