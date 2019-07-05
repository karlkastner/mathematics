% 29 2013-08-18 15:33:17
% Sat Nov  9 20:35:43 UTC 2013
% Karl Kästner, Berlin
%% convolutions with padding
function A = conv1_man(f,A,mode)
	if (nargin() < 3)
		mode = '';
	end
	switch (mode)
	case {'paddconst'}
		n = length(f);
		A = [A(1)*ones(size(f)); A; A(end)*ones(size(f))];
		A = conv(A,f,'same');
		A = A(n+1:end-n);
	case {'paddzero'}
		A = conv(A,f,'same');
	case {'truncate'}
		% TODO, only works for odd nf
		% TODO, this is problematic, when l(A) < l(f)
		t = isrow(A);
		n = length(A);
		f = rvec(f);
		A = cvec(A);
		B = zeros(n,1);
		nf = length(f);
		nfh = 0.5*(nf+1);
		sf = sum(f);
		% left
		for idx=1:nfh-1
			% sub-filter
			f_ = f(nfh-idx+1:nf);
			% normalise
			f_ = sf/sum(f_)*f_;
			B(idx) = f_*A(1:nfh-1+idx);
		end
		% centre
		B(nfh:end-nfh+1) = conv(A,fliplr(f),'valid');
		%for idx=nfh:n-nfh+1
		%	B(idx) = f*A(idx-nfh+1:idx+nfh-1);
		%end
		% right
		for idx=1:nfh-1
			% sub filter
			f_ = f(1:nf-idx);
			f_ = sf/sum(f_)*f_;
			B(n-nfh+idx+1) = f_*A(n-nf+idx+1:n);
		end
		A = B;
		if (t) A = A'; end
	otherwise % case {'extrapolate'}
		t = isrow(A);
		A = cvec(A);
		f = rvec(f);
		n = length(A);
		B = zeros(n,1);
		nf = length(f);
		nfh = 0.5*(nf+1);
		% linear extrapolation of end points
		for idx=1:nfh-1
			a = A(1:nfh-1+idx);
			R = [ones(nfh+idx-1,1), (1:nfh-1+idx)'];
			b = R\a; % todo weighing
			% etrapolate
			R = [ones(nfh-idx,1), (-nfh+idx+1:0)'];
			a = [R*b; a];
			% filter
			B(idx) = f*a; 
		end
		% centre
		B(nfh:end-nfh+1) = conv(A,fliplr(f),'valid');
		%for idx=nfh:n-nfh+1
		%	B(idx) = f*A(idx-nfh+1:idx+nfh-1);
		%end
		% right
		for idx=1:nfh-1
			a = A(n-nf+idx+1:n);
			% TODO shift by n to improve conditioning
			R = [ones(nf-idx,1), (n-nf+idx+1:n)'];
			b = R\a;
			R = [ones(idx,1), (n+1:n+idx)'];
			a = [a; R*b];
			B(n-nfh+idx+1) = f*a;
		end
		A = B;
		if (t) A = A'; end
	end % switch
end

