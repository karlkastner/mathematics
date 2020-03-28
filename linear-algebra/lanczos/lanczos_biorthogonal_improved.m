% Oct 23 00:47
% Karl Kästner, Berlin


% Computational electromagnetics  By A. Bondeson, Thomas Rylander, Pär Ingelström
% Iterative methods for sparselinear systems  By Y. Saad
% V*T*W' - A
function [T V W] = lanczos_biorthogonal(A)
	n = size(A,1);
	v = zeros(n,1); v(1) = 1;
	w = v;

	V(:,1) = v;
	W(:,1) = w;

	v_old = 0.*v;
	w_old = 0.*w;
	b = 0; d = 0;
	for idx=1:n
%{
		v_next = A*v - b*v_old;
		w_next = A*w - d*w_old;
		a = v_next'*w; AA(idx) = a;
		v_next = v_next - a*w;
		w_next = w_next - a*v_next;
%}
		a = (A*v)'*w;
		v_next = A*v - a*v - b*v_old; 
		w_next = A'*w - a*w - d*w_old;
		%w_next = A'*v - a*w - d*w_old; % typo in Bondenson

		d = sqrt(abs(v_next.'*w_next));
		b = v_next'*w_next/d;

		v_old = v;
		w_old = w;
		
		v = v_next/d;
		w = w_next/b;
		% store results
 		AA(idx) = a;
		D(idx) = d;
		B(idx) = b;
		if (idx < n)
			V(:,idx+1) = v;
			W(:,idx+1) = w;
		end
	end % for idx
	T = diag(B(1:end-1),+1) + diag(AA) + diag(D(1:end-1),-1);
end % function l3

