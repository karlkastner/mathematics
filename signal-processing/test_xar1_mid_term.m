% mi 7. Okt 16:52:02 CEST 2015
% Karl Kastner, Berlin

function test_xar1_mid_term()
n1 = 20;
n2 = 30;
m  = 10;


s1 = exp(1);
s2 = pi;

rho1  = 0.8; %1/2;
rho2  = 0.9; %1/3;
rho12 = 1/2;

%rho1 = 0;
%s1 = 1;
%n2 = n1;
%rho2 = rho1;
%s2 = s1;

S = []


S(1,1) = xar1(s1,s2,rho1,rho2,rho12,n1,n2,m);
S(1,2) = xar1_(s1,s2,rho1,rho2,rho12,n1,n2,m);

S(2,1) = xar1(s2,s1,rho2,rho1,rho12,n2,n1,m);
S(2,2) = xar1_(s2,s1,rho2,rho1,rho12,n2,n1,m);

S(3,1) = xar1_first_term(s1,s2,rho1,rho2,rho12,m);
S(3,2) = xar1_first_term_(s1,s2,rho1,rho2,rho12,m);

S(4,1) = xar1_first_term(s2,s1,rho2,rho1,rho12,m);
S(4,2) = xar1_first_term_(s2,s1,rho2,rho1,rho12,m);

S(5,1) = xar1_last_term(s1,s2,rho1,rho2,rho12,n1,n2);
S(5,2) = xar1_last_term_(s1,s2,rho1,rho2,rho12,n1,n2);

%S(1,1) = xar1_first_term(s1,s2,rho1,rho2,rho12,m,n2);
S(6,1) = xar1_mid_term(s1,s2,rho1,rho2,rho12,n1,n2,m);
S(6,2) = xar1_mid_term_(s1,s2,rho1,rho2,rho12,n1,n2,m);

S(7,1) = xar1_mid_term(s2,s1,rho2,rho1,rho12,n2,n1,m);
S(7,2) = xar1_mid_term_(s2,s1,rho2,rho1,rho12,n2,n1,m);

S(8,1) = xar1_first_last(s1,s2,rho1,rho2,rho12,n1,n2,m);
S(8,2) = xar1_first_last_(s1,s2,rho1,rho2,rho12,n1,n2,m);

S(9,1) = xar1_mid_mid(s1,s2,rho1,rho2,rho12,n1,n2,m);
S(9,2) = xar1_mid_mid_(s1,s2,rho1,rho2,rho12,n1,n2,m);
[S S(:,2)-S(:,1)]

%for idx=1:1e6
[X Y] = randar1_dual(s1,s2,rho1,rho2,rho12,max(n1,n2),1e6);
% finite population correction
eX = bsxfun(@minus,X,mean(X(1:n1,:)));
eY = bsxfun(@minus,Y,mean(Y(1:n2,:)));
E = [mean(eX(1:m,:))'.^2 mean(eX(1:m,:))'.*mean(eY(1:m,:))' mean(eY(1:m,:))'.^2];
mean(E)
serr(E)
[ar1_var_range(s1,rho1,n1,m) xar1(s1,s2,rho1,rho2,rho12,n1,n2,m) ar1_var_range(s2,rho2,n2,m)]

end

function s = xar1_last_term(s1,s2,r1,r2,r12,n1,n2)
	if (n2<n1)
		s = xar1_last_term(s2,s1,r2,r1,r12,n2,n1);
	else	
	s = r12*s1*s2/(n1*n2)*( ...
	    n1*(1-r1*r2)/((1-r1)*(1-r2)) ...
	  - r1/(1-r1)^2*(1-r1^n1) ...
	  - r2^(n2+1)/(1-r2)^2*(r2^-n1-1) );
%	s = r12*s1*s2/(n1*n2)*( ...
%	    1/(1-r1)*(n1 - r1*(1-r1^n1)/(1-r1)) ...
%	  + r2/(1-r2)*(n1 - r2^n2*1/(1-r2)*(r2^-n1-1)) );
	end
end

function s = xar1_first_term(s1,s2,r1,r2,r12,m)
	s = r12*s1*s2/m^2*( ...
	    m*(1-r1*r2)/((1-r1)*(1-r2)) ...
	    -r1/(1-r1)^2*(1-r1^m) ...
	    -r2/(1-r2)^2*(1-r2^m));
%	s = r12*s1*s2/m^2*( ...
%	    1/(1-r1)*(m - r1*(1-r1^m)/(1-r1)) ...
%	  + 1/(1-r2)*(r2*m - r2*(1-r2^m)/(1-r2)) );
%[1/(1-r2)*(r2*m - r2*(1-r2^m)/(1-r2)) 1/(1-r2)*(m - r2*(1-r2^m)/(1-r2))]
end

function s = xar1_first_last(s1,s2,r1,r2,r12,n1,n2,m)
	s = r12*s1*s2*( ... 
	(1-r1*r2)/((1-r1)*(1-r2))*(1/m+1/n2) ...
	- r1/(1-r1)^2*(       (1-r1^n1)/(n1*n2)  + (1-r1^m)/m^2 ) ...
	- r2/(1-r2)^2*( r2^n2*(r2^-n1-1)/(n1*n2) + (1-r2^m)/m^2 ) );
end

function s = xar1_mid_mid(s1,s2,r1,r2,r12,n1,n2,m)
	s = -r12*s1*s2*( ...
		-(1/n1+1/n2)*(1-r1*r2)/((1-r1)*(1-r2)) ...
		+ r1/(1-r1)^2*((1-r1^m)/(m*n2) + r1^n1*(r1^-m - 1)/(m*n1) ) ...
		+ r2/(1-r2)^2*((1-r2^m)/(m*n1) + r2^n2*(r2^-m - 1)/(m*n2) ) );
end





function s = xar1_(s1,s2,r1,r2,r12,n1,n2,m)
	 s = xar1_first_term_(s1,s2,r1,r2,r12,m) ...
	     - xar1_mid_term_(s1,s2,r1,r2,r12,n1,n2,m) ...
	     - xar1_mid_term_(s2,s1,r2,r1,r12,n2,n1,m) ...
	     + xar1_last_term_(s1,s2,r1,r2,r12,n1,n2);
end

function s = xar1_last_term_(s1,s2,r1,r2,r12,n1,n2)
	I  = repmat((1:n1)',1,n2);
	J  = repmat(1:n2,n1,1);
	D  = I-J;
	R1 =  r1.^abs(D);
	R2 =  r2.^abs(D);
	s = r12*s1*s2/(n1*n2)*(sum(sum(R1.*(I>J)))+sum(sum(R2.*(I<=J))));
end

function s = xar1_first_term_(s1,s2,r1,r2,r12,m)
	s = xar1_last_term_(s1,s2,r1,r2,r12,m,m);
end

function s = xar1_first_last_(s1,s2,r1,r2,r12,n1,n2,m)
	s =   xar1_first_term_(s1,s2,r1,r2,r12,m) ...
	    + xar1_last_term_(s1,s2,r1,r2,r12,n1,n2);
end



function s = xar1_mid_term_(s1,s2,rho1,rho2,rho12,n1,n2,m)
	I  = repmat((1:m)',1,n2);
	J  = repmat(1:n2,m,1);
	D  = I-J;
	R1 =  rho1.^abs(D);
	R2 =  rho2.^abs(D);
	s = rho12*s1*s2/(n2*m)*(sum(sum(R1.*(I>J)))+sum(sum(R2.*(I<=J))));
end

function s = xar1_mid_mid_(s1,s2,r1,r2,r12,n1,n2,m)
	s =   xar1_mid_term_(s1,s2,r1,r2,r12,n1,n2,m) ...
	    + xar1_mid_term_(s2,s1,r2,r1,r12,n2,n1,m);
end

function s = xar1_mid_term(s1,s2,r1,r2,r12,n1,n2,m)
	s = r12*s1*s2/(m*n2)*( ...
                  m*(1-r1*r2)/((1-r1)*(1-r2)) ...
		- r1/(1-r1)^2*(1-r1^m) ...
		- r2^n2/(1-r2)^2*(r2^-m-1));
%	s = r12*s1*s2/(m*n2)*(m*(1/(1-r1) + r2/(1-r2)) ...
%		- r1/(1-r1)^2*(1-r1^m) ...
%		- r2/(1-r2)^2*(r2^(n2-m)-r2^n2));
end


