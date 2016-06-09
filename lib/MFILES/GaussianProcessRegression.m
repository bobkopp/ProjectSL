function [f,V,logp,alfa,errorflags,invtraincv,invcv] = GaussianProcessRegression(x0,y0,x,traincv,cvfunc,varargin)

% [f,V,logp,alfa,errorflags,invtraincv,invcv] = GaussianProcessRegression(x0,y0,x,traincv,cvfunc,[testcv2],[invcv])
%
% Modeled on Rasmussen & Williams (2006)
%
% INPUT
%
% 	x0			x at training points
%	y0			y at training points (demeaned)
%	x			x at which to approximate
%	traincv		covariance matrix among covariance points
%	cvfunc		EITHER a handle for a function that takes x and x0 as inputs
%				OR the covariance matrix between x and x0
%	testcv2		IF cvfunc is passed as a matrix, then the covariance matrix among x
%   invcv       structure with invcv.alfa, invcv.invtraincv, invcv.s, invcv.r
% 
% Assumes a zero-mean Gaussian process; add/subtract the means separately otherwise.
%
% Last updated by  Bob Kopp, robert-dot-kopp-at-rutgers-dot-edu, Fri Aug 30 15:29:10 EDT 2013

%%%%%

	errorflags=0;
	defval('tol',1e-6);
	defval('doChol',0);
	defval('invcv',[]);

	if nargin == 5
		testcv = feval(cvfunc,x0,x);
		testcv2 = feval(cvfunc,x,x);
	else
		testcv = cvfunc;
		testcv2 = varargin{1};
	end

	if nargin>6
		invcv=varargin{2};
	end
	
	if length(invcv)==0
		if doChol
			try
				L=chol(traincv,'lower');
				if nargout>5
					invtraincv = L'\(L\eye(size(L)));
				end
				alfa=L'\(L\y0);
			catch
				disp('Not positive definite!')
				errorflags=-1;
				[alfa,invtraincv,s,r]=svdinv(traincv,y0)
				doChol = 0;
			end
		else
			[alfa,invtraincv,s,r]=svdinv(traincv,y0);
		end

		if nargout>6
			invcv.alfa=alfa;
			invcv.invtraincv = invtraincv;
			invcv.s=s; invcv.r = r;
		end
	else
		alfa=invcv.alfa;
		invtraincv=invcv.invtraincv;
		s = invcv.s; r=invcv.r;
		clear invcv;
	end
	f=testcv'*alfa;

	if nargout>1

		logpterms(1) = -.5*abs(y0'*alfa);
		logpterms(3) = -.5*length(y0)*log(2*pi);
		% logpterms(2) = -.5 * log(det(traincv));
		if doChol
			v = L\testcv;
			V=testcv2-v'*v;
			logpterms(2) = - sum(log(diag(L)));
		else
			V=testcv2-testcv'*invtraincv*testcv;
			logpterms(2) = -.5 * sum(log((s(1:r))));
		end
		V=.5*(V+V');
		logp=sum(logpterms);

		if min(diag(V))<0
			disp('Negative posterior variances!')
		end

	%	logp = -.5*y0'*alfa - .5*log(prod(s(1:r)))- .5*length(y0)*log(2*pi);
	
	end

end

function [alfa,invtraincv,s,r]=svdinv(traincv,y0)
	[m,n] = size(traincv);
	[U,S,V] = svd(traincv,0);
	s = diag(S);
	tol = max(m,n) * eps(max(s));
	r = sum(s > tol);
	invtraincv = V(:,1:r)*diag(s(1:r).^-1)*U(:,1:r)';      
	alfa = invtraincv * y0;
end