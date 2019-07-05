% 2015-01-30 21:28:59.350605293 +0100
% Karl Kastner, Berlin
%
%% expand a continous fourier series at times t
%
% 
% TODO make this a class, where c is an matraix horizontally concatenated vectors
%
% function [val] = fourier_predict(T, t0, c, t)
%
function [val,F] = fourier_predict(T, t0, c, t)
%	nt = length(t);
%	no = length(T);
%	t = cvec(t);
%	c = reshape(c,2*no+1,[]);

	t   = t-t0(1);
	F   = fourier_matrix(T, t);
	val = F*c;

%	A = zeros(nt,2*no+1);
%	A(:,1) = 1;
%	for idx=1:no
%		A(:,2*idx)   = sin(2*pi*t/T(idx));
%		A(:,2*idx+1) = cos(2*pi*t/T(idx));
%	end
%	val = A*c; %cvec(c);
	% TODO prediction error
end

