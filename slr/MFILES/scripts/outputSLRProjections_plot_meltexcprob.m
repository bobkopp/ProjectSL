% outputSLRProjections_plot_meltexcprob
% Last updated by  Bob Kopp, robert-dot-kopp-at-rutgers-dot-edu, Fri Jan 3 01:35:27 EST 2014

%%%

% AR5, BA and hybrid AIS and GIS projections
 targyear=2100; kk=1;


		ISsampsBA = bsxfun(@plus,ISLastDecade * (targyear-2000),.5*ISaccelsamps{1}*(targyear-2011).^2*(targyear>2011));
		ISsampsBAAIS = sum(ISsampsBA(:,2:3),2);
		BAq = quantile([ISsampsBA ISsampsBAAIS],[.167 .5 .833]);
		
		ARISLastDecade=[ISLastDecade(1) ISLastDecade(2)+ISLastDecade(3)];
		ISsampsAR = bsxfun(@plus,ARISLastDecade * (targyear-2000),.5*ARISaccelsamps{1}(:,:,kk)*(targyear-2011).^2*(targyear>2011));
		ARq = quantile(ISsampsAR,[.167 .5 .833]);

		qslope=diff(ARq)./(eps+diff(BAq(:,[1 4])));

		ISsampsHybrid = bsxfun(@minus,[ISsampsBA(:,1) ISsampsBAAIS],BAq(2,[1 4]));
		ISsampsHybrid = bsxfun(@plus, (ISsampsHybrid<0) .* bsxfun(@times,ISsampsHybrid,qslope(1,:)) + (ISsampsHybrid>0) .*  bsxfun(@times,ISsampsHybrid,qslope(2,:)), ARq(2,:));

		WAISsampsHybrid = ISsampsBA(:,2) - BAq(2,2);
		WAISsampsHybrid = (WAISsampsHybrid<0) .* (WAISsampsHybrid*qslope(1,2)) + (WAISsampsHybrid>0) .*  (WAISsampsHybrid*qslope(2,2)) + BAq(2,2) * ARq(2,2)/BAq(2,4);
		WAISsampsHybrid=min(WAISsampsHybrid,5000);
		
		ISsampsHybrid=[ISsampsHybrid(:,1) WAISsampsHybrid ISsampsHybrid(:,2)-WAISsampsHybrid];

		ISsampsHybrid(:,1) = min(ISsampsHybrid(:,1),7000);
		ISsampsHybrid(:,2) = min(ISsampsHybrid(:,2),5000);

		%samps(:,counter+1:counter+3,ii)=bsxfun(@plus,ISLastDecade * (targyear-2000),.5*ISaccelsamps{1}*(targyear-2011).^2*(targyear>2011));


testq=10.^[-3:.1:-.3];
testq=union(testq,1-testq);

clf; clear hp;

hp(1)=subplot(1,2,1);

plot(quantile(ISsampsBA(:,1),testq)/1000,1-testq,'g','linew',2); hold on
plot(quantile(ISsampsAR(:,1),testq)/1000,1-testq,'b','linew',2); hold on
plot(quantile(ISsampsHybrid(:,1),testq)/1000,1-testq,'r','linew',1)
ylabel('Exceedance probability');
legend('BA','AR5','Hybrid','Location','Northeast');
xlabel('m GIS mass loss');
xlim([0 1]);
longticks(gca);

hp(2)=subplot(1,2,2);

plot(quantile(ISsampsBA(:,2)+ISsampsBA(:,3),testq)/1000,1-testq,'g','linew',2); hold on
plot(quantile(ISsampsAR(:,2),testq)/1000,1-testq,'b','linew',2); hold on
plot(quantile(ISsampsHybrid(:,2)+ISsampsHybrid(:,3),testq)/1000,1-testq,'r','linew',1)
ylabel('Exceedance probability');
xlabel('m AIS mass loss');
xlim([-.3 2]);
longticks(gca);

pdfwrite('meltexcprob');
