% Sun 18 Dec 21:13:26 CET 2016
% 
%% mean angle
% TODO this is inefficient, simply comute mean of sin and cos and back-transform
% TODO weighing
function am = meanangle2(a,dim)
	if (nargin()<2)
		dim = 1;
		if (isvector(alpha))
			alpha = cvec(alpha);
		end
	end
	n = size(alpha,1);

	a   = wrapToPi(a);
	a   = sort(a);
	am  = zeros(1,size(a,2));
	for idx=1:size(a,2)
		ldx = 1;
		rdx = size(a,1); 
		while (ldx ~= rdx)
			ac = 0.5*(a(ldx,idx)+a(rdx,idx));
			ac = binsearch(a(:,idx),ac,ldx,rdx);
			if ( cdx-ldx > rdx-cdx )
				% more points in the left half, choose left half
				rdx = cdx;
			else
				% more points in right half, choose right half
				ldx = cdx;
			end % else of if
		end % while
		am(idx) = a(cdx,idx);
	end % for
end % function
