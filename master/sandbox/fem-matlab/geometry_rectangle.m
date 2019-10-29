function [x y] = geometry_rectangle(bs,s,varargin)
L0 = 10;
%function [x y] = geometry_rectangle(L0, bs,s)
	% boundary structure
	d=[
	  0 0 0 0 % start parameter value
	  1 1 1 1 % end parameter value
	  1 1 1 1 % left hand region
	  0 0 0 0 % right hand region
	];

	switch (nargin())
		case {0} % no input: number of boundary-segments
			x=4;
		case {1} % one input argument: structure of boundary
			x=d(:,bs(:)');
		case {2}
			% two input arguments: points on boundary

			% midpoint and radius of circle
			mcx=0; mcy=0; rc=4;

			% midpoint and radius of hole
			m = max(s);
			x=zeros(size(s));
			y=zeros(size(s));

			% expand bs
			if numel(bs)==1
				bs=bs*ones(size(s));
			end

			if ~isempty(s)
			  k=find(bs==1); % boundary segment 1
				x(k) = L0*s(k);
				y(k) = zeros(size(s(k)));
			  k=find(bs==2);% boundary segment 2
				x(k) = L0*ones(size(s(k)));
				y(k) = L0*s(k);
			  k=find(bs==3);% boundary segment 3
				x(k) = L0*s(flipud(k));
				y(k) = L0*ones(size(s(k)));
			  k=find(bs==4);% boundary segment 4
				x(k) = zeros(size(s(k)));
				y(k) = L0*s(flipud(k));
			end % if
		otherwise
			error('geometry','here');
	end % switch
end
