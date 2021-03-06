% Do 9. Jul 14:44:43 CEST 2015
% Karl Kastner, Berlin
%
%% variance correction factor for an autocorrelated finite process
%% n   : [1 .. inf] population size
%% m   : [1 .. n]   samples size
%% rho : [ -1 < rho < 1 (for convergence) ] correlation of samples
function f = ar1_var_factor(rho,n,m)
	f = (n-m)./(n*m)*(1+rho)/(1-rho) ...
            - 2*rho/(1-rho)^2 ...
	       * (  1./m.^2.*( 1 - rho.^m ) ...
                  - 1./(n*m).*(1 - rho.^m + rho.^(n-m) - rho^n) ...
	          + 1/n^2*(1 - rho^n) ...
                 );
end % ar1_var_factor

