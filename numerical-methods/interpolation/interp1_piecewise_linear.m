% Thu 16 Apr 21:13:23 +08 2020
function yi = interp1_piecewise_linear(cx,cy,xi)
	siz = size(xi);
	xi  = flat(xi);

	cx = rvec(cx);
	cy = rvec(cy);

	dy_dx = diff(cy)./diff(cx); 
	%for idx=1:length(cy)-1

	% constant exrapolation left end and right end
	yi      = ( cy(1)*(xi<cx(1)) ...
		  + sum(  (xi>=cx(1:end-1) & xi<cx(2:end) ) ...
		       .* (cy(1:end-1) + (xi-cx(1:end-1)).*dy_dx), 2 ) ...
		  + cy(end).*(xi>cx(end)) ...
		  );

	yi = reshape(yi,siz);
end

