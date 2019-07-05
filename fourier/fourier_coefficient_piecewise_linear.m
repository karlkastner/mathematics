% Tue 14 Aug 10:48:53 CEST 2018
% Karl Kastner, Berlin
%
%% fourier series coefficients of a piecewise linear function
%% (not coefficient of discrete fourier transform)
%% function can be discontinuous between intervals
%% scales domain length to 2pi
%%
%% input :
%% l,r : end points of piecewise linear function
%% lval, rval : values at end points
%% L : length of domain
%% n : number of samples/highest frequency   
%%
%% output :
%% a, b : coefficients for frequency components
%
function [a,b] = fourier_coefficient_piecewise_linear(l,r,lval,rval,L,n)
	a  = 2/L*sum(-(1.0./pi.^2.*1.0./n.^2.*(L.^2.*(lval-rval).*(cos((pi.*l.*n.*2.0)./L) ...
		 -cos((pi.*n.*r.*2.0)./L)).*(1.0./4.0)+L.*pi.*n.*(lval.*sin((pi.*l.*n.*2.0)./L)...
		 -rval.*sin((pi.*n.*r.*2.0)./L)).*(l-r).*(1.0./2.0)))./(l-r));
	b = 2/L*sum(-(1.0./pi.^2.*1.0./n.^2.*(L.^2.*(sin((pi.*l.*n.*2.0)./L)...
		 -sin((pi.*n.*r.*2.0)./L)).*(lval-rval).*(1.0./4.0)...
		 -L.*pi.*n.*(lval.*cos((pi.*l.*n.*2.0)./L) ...
		 -rval.*cos((pi.*n.*r.*2.0)./L)).*(l-r).*(1.0./2.0)))./(l-r));
	a(n==0) = 2/L*sum(-((lval + rval).*(l - r))/4);
	b(n==0) = 0;
end

