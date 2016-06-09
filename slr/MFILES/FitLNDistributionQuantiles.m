function [thet,fval]=FitLNDistributionQuantiles(targetmeanval,targetquantvals,minval,quantlevs,doglob)

% [thet,fval]=FitLNDistributionQuantiles(targetmeanval,targetquantvals,minval,[quantlevs],[doglob])
%
% Fit log-normal distributions base on quantiles.
% targetquants are used for fit target, other inputs for initial point.
% 
% INPUT:
%   targetmeanval (just used for estimating starting value)
%   targetquantvals [5% 50% 95%]
%   minval
%
% OUTPUT:
%   thet: [minval mu sigma]
%   fval: residual sum of squares
%
% Last updated by  Bob Kopp, robert-dot-kopp-at-rutgers-dot-edu, Thu Jul 18 13:38:23 EDT 2013

% log normal distribution with offset

defval('targetmeanval',1);
defval('targetquantvals',[-1.5 .2 10.2]);
defval('minval',-1.9);
defval('quantlevs',[.05 .5 .95]);
defval('doglob',0);

targetmean = @(x0,mu,sigma) x0 + exp(mu) .* exp(sigma.^2/2);
targetquants = @(x0,mu,sigma) x0 + exp(mu + norminv(quantlevs).*sigma);
targetf=@(x0,mu,sigma) [targetmean(x0,mu,sigma) targetquants(x0,mu,sigma)];

sigmacond = @(thet,targ) sqrt(log((targ-thet(1))./exp(thet(2)))*2);

muguess=log(targetquantvals(2)-minval);
sigmaguess = sigmacond([minval muguess],targetmeanval);
mintarget = @(thet) sum((targetquants(thet(1),thet(2),thet(3))-targetquantvals).^2);

thet=[minval muguess sigmaguess];

if doglob
	problem = createOptimProblem('fmincon','x0',thet,'objective',mintarget,'lb',[-100e6 -100e6 1e-9],'ub',[100e6 100e6 100],'options',optimset('Display','iter','maxfunevals',5000));
	[thet fval eflag output] = fmincon(problem);
	problem = createOptimProblem('fmincon','x0',thet,'objective',mintarget,'lb',[-100e6 -100e6 1e-9],'ub',[100e6 100e6 100],'options',optimset('Display','iter','maxfunevals',5000));
	gs=GlobalSearch('MaxTime',3000,'Display','iter');
	[thet,optm1.fval,optm1.exitflag,optm1.output,optm1.manymin] = run(gs,problem);
else
	[thet,fval]=fmincon(mintarget,thet,[],[],[],[],[-100e6 -100e6 1e-9],[100e6 100e6 100],[],optimset('Display','iter','maxfunevals',5000));
end