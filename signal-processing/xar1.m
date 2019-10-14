% mi 7. Okt 16:52:02 CEST 2015
% Karl Kastner, Berlin

function s = xar1(s1,s2,r1,r2,r12,n1,n2,m)
	s = r12*s1*s2*( ...
		(n1-m)/(n1*m)*(1-r1*r2)/((1-r1)*(1-r2)) ...
		- r1/(1-r1)^2*( (1-r1^m)/m^2    +       (1-r1^n1)/(n1*n2) ...
	                      - (1-r1^m)/(m*n2) - r1^n1*(r1^-m-1)/(m*n1) ) ...
		- r2/(1-r2)^2*( (1-r2^m)/m^2    + r2^n2*(r2^-n1-1)/(n1*n2) ...
                              - (1-r2^m)/(m*n1) - r2^n2*(r2^-m-1)/(m*n2) ) );
%	s1_ = ar1_var_range(s1,r1,n1,m);
%	s2_ = ar1_var_range(s2,r2,n2,m);
%	f1 = ar1_var_factor(r1,n1,m);
%	f2 = ar1_var_factor(r2,n2,m);
	%s = r12*s1*s2*(((0.5*(f1 + f2)) ));
	%s = r12*s1*s2*(((0.5*(n1*f1 + n2*f2)/sqrt(n1*n2)) ));
%f1
%f2
%m^2*f1
%m^2*f2
%	s = r12*s1*s2*(((0.5*(m^2*f1 + m^2*f2 - 2*m)/m^2) ));
%	 s = xar1_first_term(s1,s2,r1,r2,r12,m) ...
%	     - xar1_mid_term(s1,s2,r1,r2,r12,n1,n2,m) ...
%	     - xar1_mid_term(s2,s1,r2,r1,r12,n2,n1,m) ...
%	     + xar1_last_term(s1,s2,r1,r2,r12,n1,n2);
end


