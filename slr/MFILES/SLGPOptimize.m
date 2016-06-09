function [thet]=SLGPOptimize(y0,traincv,thet0,lb,ub,globl)

% Last updated by  Bob Kopp, robert-dot-kopp-at-rutgers-dot-edu, Tue Feb 11 23:37:39 EST 2014

    defval('globl',0)

    lb=max(1e-12,lb);
    ub=max(1e-12,ub);
    fitoptions=optimset('Display','iter','MaxFunEval',8000,'Algorithm','sqp','UseParallel','always','TolX',1e-3);
    hessin=[];

    if globl==1

        problem = createOptimProblem('fmincon','x0',log(thet0),'objective',(@(x) -logprob(y0,traincv,exp(x))),'lb',log(lb),'ub',log(ub),'options',fitoptions);
        [x fval eflag output] = fmincon(problem);
        %ms=MultiStart('TolX',1e-2,'MaxTime',3000,'Display','iter');
        gs=GlobalSearch('TolX',1e-2,'MaxTime',3000,'Display','iter');
        [optm1.coeffs,optm1.fval,optm1.exitflag,optm1.output,optm1.manymin] = run(gs,problem);
        if nargout>1
            hessin = optm1.manymin;
        end

        if nargout>2
            [optm1.coeffs,optm1.fval,optm1.exitflag,optm1.output,optm1.lambda,optm1.grad,hessin2] = fmincon(@(x) -logprob(y0,traincv,exp(x)),optm1.coeffs,[],[],[],[],log(lb),log(ub),[],fitoptions);
        end
    elseif globl==2
            rng(10,'twister') % for reproducibility

            fitoptions=gaoptimset('Display','iter','useparallel','always');
            [optm1.coeffs,optm1.fval] = ga(@(x) -logprob(y0,traincv,exp(x)),length(thet0),[],[],[],[],log(lb),log(ub),[],fitoptions);	
    elseif globl==3
            rng(10,'twister') % for reproducibility
            fitoptions=saoptimset('Display','iter','MaxFunEval',8000,'TolFun',2e-2,'temperaturefcn',@temperaturefast,'TimeLimit',3000);
            [optm1.coeffs,optm1.fval] = simulannealbnd(@(x) -logprob(y0,traincv,exp(x)),log(thet0),log(lb),log(ub),fitoptions);
    else
        if nargout>1
            [optm1.coeffs,optm1.fval,optm1.exitflag,optm1.output,optm1.lambda,optm1.grad,hessin] = fmincon(@(x) -logprob(y0,traincv,exp(x)),log(thet0),[],[],[],[],log(lb),log(ub),[],fitoptions);

        else
            [optm1.coeffs,optm1.fval] = fmincon(@(x) -logprob(y0,traincv,exp(x)),log(thet0),[],[],[],[],log(lb),log(ub),[],fitoptions);

        end	
    end
    thet = exp(optm1.coeffs);
    logp=optm1.fval;
end

function logp=logprob(y0,traincv,x)

	errorflags=0;
	
	tcv=traincv(x);
	try
			L=chol((tcv),'lower');
			alfa=L'\(L\y0);
			doChol=1;
	catch
		disp('Not positive definite!')
		errorflags=-1;
		doChol = 0;
        [m,n] = size(tcv);
        [U,S,V] = svd(tcv,0);
        s = diag(S);
        tol = max(m,n) * eps(max(s));
        r = sum(s > tol);
        invtraincv = V(:,1:r)*diag(s(1:r).^-1)*U(:,1:r)';      
        alfa = invtraincv * y0;
	end

    logpterms(1) = -.5*abs(y0'*alfa);
    logpterms(3) = -.5*length(y0)*log(2*pi);
    % logpterms(2) = -.5 * log(det(traincv));
    if doChol
        logpterms(2) = - sum(log(diag(L)));
    else
        logpterms(2) = -.5 * sum(log((s(1:r))));
    end
    logp=sum(logpterms);
    if errorflags == -1
        logp = -1e20;
    end
end
