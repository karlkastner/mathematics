% Mon  4 Oct 15:11:38 CEST 2021
% 
% note that this effectively low-pass filters the random phase (underestimates the variation)
function [d] = determine_phase_shift(y,p)
	n  = length(y);
	% number of segments spanning about 2 wavelengths
	m  = round(1/p)
	% number of taps per wavelength
	ni = round(n/m);
	nj = 2*ni;
	% window
	%w = tukeywin(n);
% TODO wrap around or not wrap
	y(end+(1:ni)) = 0;
	%y = [y;y];

	% number of segment pairs
	na  = n*(n-1)/2;
%	A   = zeros(na,m);
	A = [];
	rhs = []; %zeros(na,m);
	
	% for each segment
	k = 0;
	for idx=1:m-1
	 for jdx=idx+1:m
		% current segment
		id = (idx-1)*ni+(1:ni);
		yi = y(id);
		% next segment
		% TODO, needs overlap
		jd = (jdx-1)*ni+(1:nj)+(ni-nj)/2;
%limits(jd) 
		yj = y(jd);
		% determine cross correlation
		xc = xcorr(yj,yi);
		% get shift of indices
		[mx,mdx] = max(xc);
		centre = (nj+ni)/2;
%subplot(m,m,(idx-1)*m+jdx)
%cla
%plot(jd-jd(1),yj)
%hold on
%plot(id-id(1)+centre,yi)
%plot(id-id(1)-centre+mdx+(ni-nj)/2,yi,'--');
		k        = k+1;
		A(k,idx) = -1;
		A(k,jdx) =  1;
		% TODO
%[mdx-centre,ni,nj,(ni-nj)]
		rhs(k,1)   = (jdx-idx)*ni + mdx - centre + (ni-nj);
	 end
	end
	% get the relative indices with respect to the first segment
	% A*b is : 
	%     [0, b-a, c-a
        %     [ ,   0, c-b  
	% delete first column of A, to make determinant of A'A nonzero

	% distance from start
	d = [0;A(:,2:end) \ rhs]
%	dd = [rhs,A*d]
%	dd(:,2)-dd(:,1)
%pause

	% correct from segment start
	d = d - (0:m-1)'*ni;
end

