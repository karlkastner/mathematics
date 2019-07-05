% Tue 12 Dec 14:19:11 CET 2017
% Karl Kastner, Berlin

%% second order odes 
%% transmittded and reflected wave
%
% y : value of left and right traveling wave
% l : left and right eigenvalue
% 
% y(x)    = y1 exp(-l1 x) + y2 exp(l2 x)
% d/dx y = (d/dx y1 - l1 y1) exp(-l1 x) + ...
%
% TODO, it would be better, to define take the sign into both roots
function dydx_ = ode2characteristics(x,l,y)
	dl_dx  = cdiff(l)./(cdiff(x)*[1, 1]);
	d      = l(:,1)-l(:,2);
	dydx_(:,1) = (l(:,1) ...
                      - 1./d.*dldx(:,1)).*y(:,1) ...
                      - 1./d.*dldx(:,2).*y(:,2);

	% TODO there is something wrong, paratheses are different
	dydx_(:,2) = l(:,2) ...
                     + 1./d.*dldx(:,2).*y(:,2) ...
                     + 1./d.*dldx(:,1).*y(:,1);
end

