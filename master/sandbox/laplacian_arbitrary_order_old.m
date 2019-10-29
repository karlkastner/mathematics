% Tue Nov  1 03:02:24 MSK 2011
% Karl Kästner, Berlin

% TODO: use integers and later apply 1/h^2, otherwise huge cancellation error
% even for small matrices!
% a = int64([2^30 2^31 2^32]), b=int64([2 2 2]), a.*b, b.*a

% variable grid laplacian
function L = laplacian(X, order)
	n = length(X);
	% X = int32(X);
	% X = int64(X); % not supported :(
	% X = vpa(X);
	Xl = X(1:end-2);
	Xc = X(2:end-1);
	Xr = X(3:end);
	if (order <= 2)
		D_num = 2*[ Xc-Xl Xl-Xr Xr-Xc ];
		D_den = (((Xr-Xl).*(Xr-Xc).*(Xc-Xl))*[1 1 1]);
		D = double(D_num)./double(D_den);
		L = spdiags(D,-1:1,n-2,n-2);
	else
%	D = 
%	[ 2*(Xl*(Xr-2Xc+Xrr) + Xc*( -2*Xr - 2*Xrr +3*Xc) + Xr*Xrr),
% 	 -2*( Xll*(Xr-2*Xc+Xrr) +Xc*(-2*Xr-2*Xrr+3*Xc)),
%	  2*( Xll*(Xl-3*Xc+Xr+Xrr) + Xl*(-3*Xc+Xr+Xrr) + Xc*(-3*Xr-3*Xrr+6*Xc) +Xr*Xrr)
%	 -2*(Xll*(Xl
%	 -2*( - 2*Xc*Xll - 2*Xc*Xrr - 2*Xc*Xl + Xl*Xrr + Xll*Xrr + 3*Xc^2), 
%	 2*(Xl*Xll - 2*Xc*Xll - 2*Xc*Xr - 2*Xc*Xl + Xl*Xr + Xll*Xr + 3*Xc^2)]
	% linear extrapolation of X
		Xll = [ 2*X(1)-X(2); X(1:end-3)];
		Xrr = [ X(4:end); 2*X(end)-X(end-1) ];
		D_num = [ 2.*Xl.*Xr - 4.*Xc.*Xr - 4.*Xc.*Xrr - 4.*Xc.*Xl + 2.*Xl.*Xrr + 2.*Xr.*Xrr + 6.*Xc.*Xc, ...
			4.*Xc.*Xll + 4.*Xc.*Xr + 4.*Xc.*Xrr - 2.*Xll.*Xr - 2.*Xll.*Xrr - 2.*Xr.*Xrr - 6.*Xc.*Xc, ...
			2.*Xl.*Xll - 6.*Xc.*Xll - 6.*Xc.*Xr - 6.*Xc.*Xrr - 6.*Xc.*Xl + 2.*Xl.*Xr + 2.*Xll.*Xr + 2.*Xl.*Xrr + 2.*Xll.*Xrr + 2.*Xr.*Xrr + 12.*Xc.*Xc, ...
			4.*Xc.*Xl + 4.*Xc.*Xll + 4.*Xc.*Xrr - 2.*Xl.*Xll - 2.*Xl.*Xrr - 2.*Xll.*Xrr - 6.*Xc.*Xc, ...
			2.*Xl.*Xll - 4.*Xc.*Xll - 4.*Xc.*Xr -	4.*Xc.*Xl + 2.*Xl.*Xr + 2.*Xll.*Xr + 6.*Xc.*Xc];
		D_den = [ (Xc - Xll).*(Xl - Xll).*(Xll - Xr).*(Xll - Xrr), ...
			  (Xc - Xl).*(Xl - Xll).*(Xl - Xr).*(Xl - Xrr), ...
		          (Xc - Xl).*(Xc - Xll).*(Xc - Xr).*(Xc - Xrr), ...
        	          (Xc - Xr).*(Xl - Xr).*(Xll - Xr).*(Xr - Xrr), ...
	        	  (Xc - Xrr).*(Xl - Xrr).*(Xll - Xrr).*(Xr - Xrr)];
		D = double(D_num) ./ double(D_den);
		L = spdiags(D,-2:2,n-2,n-2);
	% interpolation kernel
		h1 = (X(1) - X(2));
		h2 = (X(1) - X(3));
		h3 = (X(1) - X(4));
		h4 = (X(1) - X(5));
		ik_den = 12*h1^2;
		ik_num = [((h1 + h2)*(h1 + h3)*(h1 + h4))/((h1 - h2)*(h1 - h3)*(h1 - h4)), ...
			   -(2*h1^2*(h1 + h3)*(h1 + h4))/(h2*(h1 - h2)*(h2 - h3)*(h2 - h4)), ...
			    (2*h1^2*(h1 + h2)*(h1 + h4))/(h3*(h1 - h3)*(h2 - h3)*(h3 - h4)), ...
			   -(2*h1^2*(h1 + h2)*(h1 + h3))/(h4*(h1 - h4)*(h2 - h4)*(h3 - h4)) ];
		L(1,1:4) = L(1,1:4) - double(ik_num)/double(ik_den);

		h1 = (X(end) - X(end-1));
		h2 = (X(end) - X(end-2));
		h3 = (X(end) - X(end-3));
		h4 = (X(end) - X(end-4));
		ik_den = 12*h1^2;
		ik_num = [((h1 + h2)*(h1 + h3)*(h1 + h4))/((h1 - h2)*(h1 - h3)*(h1 - h4)), ...
			   -(2*h1^2*(h1 + h3)*(h1 + h4))/(h2*(h1 - h2)*(h2 - h3)*(h2 - h4)), ...
			    (2*h1^2*(h1 + h2)*(h1 + h4))/(h3*(h1 - h3)*(h2 - h3)*(h3 - h4)), ...
			   -(2*h1^2*(h1 + h2)*(h1 + h3))/(h4*(h1 - h4)*(h2 - h4)*(h3 - h4)) ];
		ik_num = fliplr(ik_num);
		L(end,end-3:end) = L(end,end-3:end) - double(ik_num)/(ik_den);
	end
end % laplacian

