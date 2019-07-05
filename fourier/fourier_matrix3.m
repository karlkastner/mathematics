% Sat  6 May 10:38:42 CEST 2017
% Karl Kastner, Berlin
%
%% transformation matrix for the continous fourier transform
%% this is a matrix with (2*n+1) real columns
%
% function A = fouriermtx(t,T)
%
function A = fouriermtx(t,T)
	n = length(t);
	% set up regression matrix
	A      = zeros(n,1+2*length(T));
	A(:,1) = 1;
	t_ = 2*pi*t;
	for idx=1:length(T)
		A(:,2*idx)   = cos(t_/T(idx));
                A(:,2*idx+1) = sin(t_/T(idx));
	end
end

