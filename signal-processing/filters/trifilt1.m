% Mo 30. Mär 12:10:42 CEST 2015
% Karl Kastner, Berlin
%
%% filter with triangular window
%% trifilt1 is ident to twice applying rectfilt1 (meanfilt1) with half the domain size
%% note : inifnitely many convolution yield a gaussian
%
% function [Y f] = trifilt1(X,nf)
%
function [Y f] = trifilt1(X,nf)

%	if (isvector(X))
%		Y = conv(X,f,'same');
%	else
	f =[];
	if (1)
	nf_ = floor(nf/2);
	X = meanfilt1(X,nf_);
	Y = meanfilt1(X,nf-nf_);

	else
	n1 = floor(nf/2);
	nc = ceil(nf/2);
	nx = size(X,1);
	f = [(1:n1) (nf-n1):-1:1];
	f = f/sum(f);
		Y = zeros(size(X));
		for idx=1:size(X,2);
			% centre
			Y(:,idx) = conv(X(:,idx),f,'same');
			% phase in
			for jdx=1:nc-1
				f_         = rvec(f(n1+1-jdx+1:end));
				f_	   = f_/sum(f_);
				Y(jdx,idx) = f_*X(1:nc+jdx-1,idx);
			end % for jdx
			% phase out
			for jdx=1:nc-1
				f_ 	   = rvec(f(1:nc+jdx-1));
				f_	   = f_/sum(f_);
				Y(end-jdx+1,idx) = f_*X(end-nc-jdx+2:end,idx);
%			end
		end % for idx (each colum)
	end % if not a vectpr
	end
end % trifilt1

