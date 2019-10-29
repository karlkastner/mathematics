function [x,y] = geometry_circle_with_hole(bs,s)

	% boundary structure
	d=[
	  0 0 0 0 % start parameter value
	  1 1 1 1 % end parameter value
	  1 1 0 0 % left hand region
	  0 0 1 1 % right hand region
	];

	switch (nargin)
		case {0} % no input: number of boundary-segments
			x=4;
		case {1} % one input argument: structure of boundary
			d
			bs(:)
			x=d(:,bs(:)');
		case {2}
			% two input arguments: points on boundary

			% midpoint and radius of circle
			mcx=0; mcy=0; rc=4;

			% midpoint and radius of hole
			mhx=0; mhy=0; rh=1;
			x=zeros(size(s));
			y=zeros(size(s));

			% expand bs
			if numel(bs)==1
				bs=bs*ones(size(s));
			end

			if ~isempty(s)
			  k=find(bs==1);% boundary segment 1: circle
				x(k)=mcx+rc*cos((pi)*s(k));
				y(k)=mcy+rc*sin((pi)*s(k));
			  k=find(bs==2);% boundary segment 2: circle
			  	x(k)=mcx+rc*cos((pi)*s(k)+pi);
				y(k)=mcy+rc*sin((pi)*s(k)+pi)
			  k=find(bs==3);% boundary segment 3: hole
				x(k)=mhx+rh*cos((pi)*s(k));
				y(k)=mhy+rh*sin((pi)*s(k));
			  k=find(bs==4);% boundary segment 4: hole
				x(k)=mhx+rh*cos((pi)*s(k)+pi);
				y(k)=mhy+rh*sin((pi)*s(k)+pi)
		end
		otherwise
			error('geometry','here');
	end % switch
end % circle_with_hole

