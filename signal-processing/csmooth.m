% 2017-05-08 13:26:24.568005932 +0200
% note: smoothing radius ~ sqrt(n)
%% smooth recursively with [1,2,1]/4 kernel
%% function x = csmooth(x,n,p,circ)
function x = csmooth(x,n,p,bndrule)
	if (nargin()<2 || isempty(n))
		n = 1;
	end
	if (nargin() < 3 || isempty(p))
		p = 1;
	end
	if (nargin() < 4 || isempty(bndrule))
		bndrule = 'keep';
	end
	for idx=1:n
		switch (bndrule)
		case {'keep'}
		x = (1-p)*x ...
                     + p*[(x(1,:));
	                  0.25*(x(1:end-2,:) + 2*x(2:end-1,:) + x(3:end,:));
                          (x(end,:))];
		case {'circ'}
		x = (1-p)*x ...
                     + p*[0.25*x(end,:) + 0.5*x(1,:) + 0.25*x(2,:);
	                  0.25*(x(1:end-2,:) + 2*x(2:end-1,:) + x(3:end,:));
                          0.25*x(end-1,:)+0.5*x(end,:)+0.25*x(1,:)];
		case {'extrap'}
		x = (1-p)*x ...
                     + p*[0.5*(x(1,:) + x(2,:));
	                  0.25*(x(1:end-2,:) + 2*x(2:end-1,:) + x(3:end,:));
                          0.5*(x(end-1,:)+x(end,:))];
		end
	end % for idx
end % csmooth

