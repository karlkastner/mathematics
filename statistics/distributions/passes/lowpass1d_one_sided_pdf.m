% Tue 30 Nov 16:57:02 CET 2021
% Karl Kastner, Berlin
%
% note : this is the discrete one-sided filter,
%        the equation is only valid for p = 1
function S = lowpass1d_one_sided_pdf(j,r,p,n)
	dk = 2*pi/n;
	S = (  1 - r*cos(dk*j) - r^n*(1 - r*cos(dk*j))) ...                                
            ./((1 - 2*r*cos(dk*j) + r^2));
	S = S.^p;
end
