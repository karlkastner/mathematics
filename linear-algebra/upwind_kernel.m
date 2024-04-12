% Sun  3 Mar 18:47:25 CET 2024
function k = upwind_kernel(a)
	p = (cvec(a)<0);
	k = cvec(a).*(p*[0,-1,1] + (1-p)*[-1,1,0]);
end

