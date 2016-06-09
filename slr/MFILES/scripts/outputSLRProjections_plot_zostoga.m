% outputSLRProjections_plot_zostoga
% Last updated by Robert Kopp, robert-dot-kopp-at-rutgers-dot-edu, Mon May 12 11:54:42 EDT 2014


%%
% figure: drift correction of thermal expansion


kk=1;
baseyear=2000;
[GSLx,GSLy,~,GSLysd]=ReadDenoisedPSMSLData(0,0,[],[],1,[],psmsldir,gslfile);
GSLy = GSLy/1000;
GSLysd=GSLysd/1000;
sub=find(GSLx(:,3)<=1900);
GSLy = GSLy-mean(GSLy(find(floor(GSLx(:,3))==baseyear)));

selectyears=[1861 floor(GSLx(sub(end),3))];

GICyears=projGICyrs{1}(:,1);
clear GICmodelmap;
sub1=find(GICyears==selectyears(1));
sub2=find(GICyears==selectyears(2));
histGIC=squeeze(sum(projGIC{1}(:,:,:),2));
dhistGIC=histGIC(sub2,:)-histGIC(sub1,:);
histGICrate = dhistGIC/diff(selectyears)/1000;

[jk,subi,subj]=intersect(floor(GSLx(:,3)),GICyears);
GSLlessGIC = GSLy(subi)-mean(sum(projGIC{1}(subj,:,:),2),3)/1000;
sub=find(jk==baseyear); GSLlessGIC=GSLlessGIC-GSLlessGIC(sub);
GSLlessGICa = GSLlessGIC + sqrt(GSLysd(subi).^2 + (std(sum(projGIC{1}(subj,:,:),2),[],3)/1000).^2)*2;
GSLlessGICb = GSLlessGIC - sqrt(GSLysd(subi).^2 + (std(sum(projGIC{1}(subj,:,:),2),[],3)/1000).^2)*2;

clf;

subplot(2,2,1)
plot(GSLx(:,3),GSLy,'k','linew',2);
hold on
plot(GSLx(subi,3),GSLlessGIC(subi),'k--','linew',2);
plot(ZOSTOGAyrs{kk},ZOSTOGA{kk});
%legend('GSL','GSL - [mean/max/min] GIC proj','Location','Northwest');
plot(GSLx(:,3),GSLy,'k','linew',2);
plot(GSLx(subi,3),GSLlessGIC(subi),'k--','linew',2);
plot(GSLx(subi,3),GSLlessGICa(subi),'k--','linew',2);
plot(GSLx(subi,3),GSLlessGICb(subi),'k--','linew',2);
xlim([1860 2200]); ylim([-.2 .50]);
longticks;
ylabel('m');
%title('smoothed zostoga before drift removal');
%pdfwrite('zostoga');

subplot(2,2,2);
plot(GSLx(:,3),GSLy,'k','linew',2);
hold on
plot(GSLx(subi,3),GSLlessGIC(subi),'k--','linew',2);
plot(ThermExpYears(1:end-1,kk),sZOSTOGA{kk});
%legend('GSL','GSL - GIC','Location','Northwest');
plot(GSLx(:,3),GSLy,'k','linew',2);
plot(GSLx(subi,3),GSLlessGIC(subi),'k--','linew',2);
plot(GSLx(subi,3),GSLlessGICa(subi),'k--','linew',2);
plot(GSLx(subi,3),GSLlessGICb(subi),'k--','linew',2);
xlim([1860 2200]); ylim([-.2 .50]);
longticks;
ylabel('m');
%title('smoothed zostoga with drift substitution (mean GIC)');
%pdfwrite('zostoga_driftsubst');
pdfwrite('zostoga');

