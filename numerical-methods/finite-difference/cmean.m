% 2016-11-22 16:11:13.487438288 +0100
% Karl Kastner, Berlin
%% single gaussian smoothing step with kernel 1/4*[1,2,1]
% function x = cmean(x,dim,p,nk,keepends)
function x = cmean(x,dim,p,n,nk,keepends)
	if (nargin()<2||isempty(dim))
		% dimension
		dim = 1;
	end
	if (nargin()<3 || isempty(p))
		% relaxation
		p = 0.5;
	end
	if (nargin()<4||isempty(n))
		n = 1;
	end
	if (nargin()<5||isempty(nk))
		% kernel width
		nk = 3;
	end
	if (nargin()<6||isempty(keepends))
		% boundary treatment
		keepends = false;
	end
	if (2 == dim)
		x  = x.';
	end
	for idx=1:n
	switch (nk)
	case {3} % 3 points
		switch (keepends)
		case {0}
		x = (1-p)*x + p*[    (2*x(2,:)      - x(3,:))
	                         0.5*(  x(1:end-2,:)+x(3:end,:))
		                     (2*x(end-1,:)  - x(end-2,:))
				];
		case {1}
	             x(2:end-1,:)  = (1-p)*x(2:end-1,:) + ...
				   p*0.5*(x(1:end-2,:)+x(3:end,:));
		otherwise
			error('keepends has to be true or false');
		end
	case {5} % 5 points
		switch (keepends)
		case {0}
			x = (1-p)*x + p*[     (2*x(2,:) - x(3,:));
					  0.5*(x(1,:)   + x(3,:));
				          p/6*(   -x(1:end-4,:)  +  4*x(2:end-3,:) ...
					       + 4*x(4:end-1,:)  -    x(5:end,:)   );
					  0.5*(x(end-2,:)    + x(end,:));
		                              (2*x(end-1,:)  - x(end-2,:))
					];
		case {1}
			x(2:end-1,:) = (1-p)*x(2:end-1,:) ...
				+ p*[	  0.5*(x(1,:)   + x(3,:));
				          p/6*(   -x(1:end-4,:)  +  4*x(2:end-3,:) ...
					       + 4*x(4:end-1,:)  -    x(5:end,:)   );
					  0.5*(x(end-2,:)    + x(end,:))
				    ];
		otherwise
			error('keepends has to be true or false');
		end
	otherwise
		error('kernel width has to be 3 or 5');
	end
	end % for idx
	if (2 == dim)
		x  = x.';
	end
end

