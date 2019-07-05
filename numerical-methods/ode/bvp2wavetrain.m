% Wed  1 Nov 11:40:08 CET 2017
% Karl Kastner, Berlin
%
%% solve second order boundary value problem by repeated integration
%
%function [x z1 q1 Q1 Z QQ F] = wavetrainq(cfun,X,z10,omega,hfun,wfun,opt)
function [x, y, cflag] = bvp2wavetrain(odefun,bcfun,X,opt)
	
	nout = nargout();
	
	if (nargin() < 7)
		opt = struct();
	end
	if (~isfield(opt,'nx'))
		nx = 1024;
	else
		nx = opt.nx;
	end
	if (~isfield(opt,'kmax'))
		kmax = 100;
	else
		kmax = opt.k;
	end
	
	% discretise domain
	% TODO, use numerically more stable version p1 X1 + p2 X2
	dx = (X(2)-X(1))/(nx-1);
	x  = X(1) + dx*(0:nx-1)';
	% w  = wfun(x);

	% initial estimate
	% incase of river tides, this is Q1
	% TODO, better initial condition with is z1 w sqrt(gh)
	y = ones(size(x));

%	Q1 = wavetrainq_(Q1);
	[y cflag] = picard(@wavetrainq_,y);
	%[Q cflag] = picard(@[x z q F] = wavetrainq_(cfun,X,omega,opt,hfun),Q1); %,sopt);	

	if (nout > 3)
		QQ = fliplr(cumsum(F.').');
		% TODO devide by w
		ZZ = -1./(1i*omega)*cdiff(QQ)/dx;	
	end

function y = wavetrainq_(y)

	% coefficients of the ode
	c = odefun(x,y);
	% normalize
	c = bsxfun(@times,c,1./c(:,1));

	% roots of the characteristic polynomial
	c_ = c;
	r = roots2(c);
	r = [r(:,2),r(:,1)];

	c=c(:,3);
	%c = conj(c);
	%c = -conj(c);
	
	% convergence term
	dcdx = cdiff(c)./dx;
	
	% reflection/transmission coefficient
	R    = -1/4*dcdx./c;
	
	% integral of the transmitted and damped travelling wave
	iL   = cumintL(-sqrt(-c) + R, dx);
	iL_  = cumintL(-sqrt(-c), dx);
	iR   = cumintR(-sqrt(-c) - R, dx);
	iR_  = cumintR(-sqrt(-c) , dx);
	
	f0   = -1;
	
	% incoming wave
	f  = f0*exp(iL);
	f_ = f0*((c(1)./c).^(1/4)).*exp(iL_);
	f__ = exp(cumintL(r(:,1),x));

	if (nout>3)	
		F  = zeros(nx,k);
		F_ = zeros(nx,k);
		F(:,1)  = f;
		F_(:,1) = f_;
	end
	fs  = f;
	fs_ = f_;
	fs__ = f__;
	idx = 1;
	fs__old = 0;
	fsold = 0;
	while (true)
		if (idx>=kmax)
			disp('stopping prematurely');
			break;
		end
		idx=idx+1;
		% waves reflected an odd number of times
		f  = exp(iR).*cumintR(R.*f.*exp(-iR),dx);
		f_ = ((1./c).^(1/4)).*exp(iR_) .* ...
                     cumintR(R.*f_.*(c).^(1/4).*exp(-iR_),dx);

		% compute the residual
		%res = c_(:,1).*cdiff(f__,2)./dx.^2 + c_(:,2).*cdiff(f__)./dx + c_(:,3).*f__;
		%res = c_(:,1).*cdiff(fs__,2)./dx.^2 + c_(:,2).*cdiff(fs__)./dx + c_(:,3).*fs__;

%		figure(1e4)		
%		plot(abs(res)./abs(fs__))
%		pause(1)

		res = cdiff(r(:,1))./cdiff(x).*f__;
		%res =  cdiff(1i./f__.*cdiff(f__)./cdiff(x))./cdiff(x).*f;
			%c_(:,1).*cdiff(f__,2)./dx.^2 + c_(:,2).*cdiff(f__)./dx + c_(:,3).*f__;
		f__0 = f__;
		%f__ = exp(-cumintR(r(:,1),dx)).*cumintR(1./(2*r(:,1)).*res.*exp(cumintR(r(:,1),dx)),dx);
		f__ = exp(cumintR(-r(:,1),dx)).*cumintR(1./(2*r(:,1)).*res.*exp(cumintR(r(:,1),dx)),dx);

		%f_ = ((c(end)./c).^(1/4)).*exp(iR_) .* ...
                %     cumintR(R.*f_.*(c./c(end)).^(1/4).*exp(-iR_),dx);
		fs = fs+f;
		fs_ = fs_+f_;
		fs__ = fs__+f__;
		if (nout>3)
			F(:,idx) = f;
			F_(:,idx) = f;
		end

		if (idx>=kmax)
			disp('stopping prematurely');
			break;
		end
		idx = idx+1;

		% waves reflected an even number of times
		f  = exp(iL).*cumintL(-R.*f.*exp(-iL),dx);
		f_ = (1./c).^(1/4).*exp(iL_) .* ...
		     cumintL(R.*f_.*(c).^(1/4).*exp(-iL_),dx); % -R ?
		%f_ = (c(1)./c).^(1/4).*exp(iL_) .* ...
		%     cumintL(-R.*f_.*(c./c(1)).^(1/4).*exp(-iL_),dx);

		% compute the residual
		%res = -(c_(:,1).*cdiff(f__,2)./dx.^2 + c_(:,2).*cdiff(f__)./dx + c_(:,3).*f__);
		% i k f = df/dx -> k = i/f df/dx -> dk/dx = d/dx(i/f df/dx) 
		res = cdiff(r(:,1))./cdiff(x).*f__;
		%res =  cdiff(1i./f__.*cdiff(f__)./cdiff(x))./cdiff(x).*f;
		f__ = exp(cumintL(-r(:,2),dx)).*cumintL(1./(2*r(:,2)).*res.*exp(cumintL(r(:,2),dx)),dx);

		fs = fs+f;
		fs_ = fs_+f_;
		fs__ = fs__+f__;
		if (nout>3)
			F(:,idx) = f;
			F_(:,idx) = f_;
		end % if nargout

	%if (max(abs(fs__ - fs__old)) < sqrt(eps))
	if (max(abs(fs - fsold)) < sqrt(eps))
		fprintf('Converged in %d iterations\n',idx)
		break;
	end
	fs__old = fs__;
	fsold = fs;

	end % while

	y = fs;

%	q1 = Q1./w;

	% compute surface elevation
%	z1 = -1./(1i*omega)*cdiff(q1)/dx;
	
	% apply initial value
%	scale = z10/z1(1);
%	Q1    = scale*Q1;
%	q1    = scale*q1;
%	z1    = scale*z1;

	% apply boundary condition at left end
	% (the current implementation works only with zero value at right end)
	% TODO extend for general case and make solution a linear combination
	% of the left and right-going waves
	[bv bp] = bcfun(X(1),y(1));
	scale   = bv/(bp(1)*y(1) + bp(2)*(y(2)-y(1))/dx);
	y       = scale*y;
	
end % wavetrainq_

end % wavetrainq

