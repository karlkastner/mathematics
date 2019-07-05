% Di 22. Dez 13:49:04 CET 2015
% Karl Kastner, Berlin
%
%% interpolation along a parametric curve with variable step width
%
function [Xi Yi Hi] = resample1(X,Y,H)
	X = cvec(X);
	Y = cvec(Y);
	H = cvec(H);
mode = 'geometric';
switch (mode)
case {'geometric'}
	% geometric scaling
	dS = hypot(diff(X),diff(Y));
%	dN = abs(log(H(1:end-1)./H(2:end)));
	dN = dS./sqrt(H(1:end-1).*H(2:end));
	N  = [0; cumsum(dN)];
	% round to full integer
	% ceil is required, otherwise N may be all zeros
	N = (ceil(N(end))/N(end)) * N;
	Ni = (0:N(end))';
	% TODO, steps in N are still linear
	XYHi = interp1(N,[X Y log(H)],Ni,'linear');
	Hi = exp(XYHi(:,3));
case {'arithmetic'}
	dS   = hypot(diff(X),diff(Y));
	Hc   = 0.5*(H(1:end-1)+H(2:end));
	dN   = dS./Hc;
	N    = [0; dN];

	% round to full integer
	% ceil is required, otherwise N may be all zeros
	N = (ceil(N(end))/N(end)) * N;

	Ni = (0:N(end))';

	% TODO, averaging X and Y is better than sampling
	XYHi = interp1(N,[X Y H],Ni,'linear');
	Hi = XYHi(:,3);
otherwise
	n   = length(X);
	S   = [0; cumsum(hypot(diff(cvec(X)),diff(cvec(Y))))];
	Si  = 0;
	Hi  = H(1);
	idx = 1;
	jdx = 1;
	while (idx < n)
		Si(jdx+1)   = Si(jdx)+H(idx);
		jdx = jdx+1;
		while (idx < n && S(idx) < Si(jdx))
			idx=idx+1;
		end
	end
	% rescale Si
	Si = (S(end)/Si(end))*Si;

	% interpolate X and Y
	[XYHi] = interp1(S,[X Y H],Si);
	Hi = XYHi(:,3);
end % witch mode
	Xi = XYHi(:,1);
	Yi = XYHi(:,2);
end % resample1

