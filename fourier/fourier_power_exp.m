% Thu 23 Mar 10:51:00 CET 2017
% Karl Kastner, Berlin
%
%% powers of the continuous fourier series 
%% 	a^p = (ur + u1 sin(ot) + u2 sin(ot+dp))^p
%% phase of first component assumed 0
%% 
%% higher orders than 2 ignored input
%% higher order than 3 not computed in output
%%
%% y = a_0 + sum (a_j sin(jot) + b_j cos(jot))
%%   = Real(sum_{i=0}^inf c_i exp(1i*omega), c_i = a_i + b_i
%%
%% NOT the alternative sum_{i=-inf}^inf \tilde c_i, tile c_j = 1/2 a_j + 1/2i b_j
%
% TODO use cos omega t !!!
function cp = fourier_exp_power(c)
	if (length(c)<3)
		% dummy
		c(3) = 0;
	end
	if (length(c)>3)
		warning('only implemented for two frequency components');
	end
	c0 = c(1);
	c1 = c(2);
	c2 = c(3);

	% dimensions of output cp
	% first dimension  : frequency (mean, 1, 2, 3)
	% second dimension : power (1, 2, 3)
	
	%
	% first power a^1
	%
	% a^1 zero frequency term (mean)
	cp(1,1) = c0;

	% a^1 diurnal
	cp(2,1) = c1;

	% a^1 semidiurnal
	cp(3,1) = c2;

	% a^1 terdiurnal (zero, as input only has two species)
	cp(4,1) = 0;

	%
	% second power a^2
	%

	% a^2 mean (zero frequency)
	% (u0 + u1 sin(x) + u2 sin(2x+dp) )^2  
	%                     = u0^2 + 2 u0 u1 sin(x) + 2u0 u2 sin(2x+dp) + 2 u1 u2 sin(x)sin(2x dp) + u1^2 sin(x)^2 + u2^2 sin(2x dp)^2
	%		      = ur^2 + 2 ur ut sin(a) + ut^2 (1/2 - 1/2 cos(2a))
	% note: there is no interaction in u2 with respect to the mean between D1 and D2,
	% if the chebycheff expansion is only computed until third order
	cp(1,2) = (c0^2 + 1/2*(abs(c1)^2 + abs(c2)^2); % 1/2,1,1

	% a^2 diurnal component
	% a^2_1 = 2*a0*a1*sin(x) + a1*a2*cos(x + dp)
	%      = 2*a0*a1*sin(x) + a1*a2*(-sin(x)*sin(dp) + cos(x)*cos(dp))
	cp(2,2) = 2*c0*c1 + conj(c1)*c2;	% 1 1 1

	% a^2 semidiurnal component
	cp(3,2) = 2*c0*c2 + 1/2*c1^2;  % 1,1,1

	% a^2 terdiurnal component
	cp(4,2) = c1*c2;	% 1,1,1

	% a^2 quaterdiurunal component
	cp(5,2) = 1/2*c2^2;	% 1,1,1

	%
	% third power c^3
	%

	% mean
%        cp(1,3) = c0^3 + (3*c0*abs(c1)^2)/2 + (3*c0*abs(c2)^2)/2 + real((c2*conj(c1)^2)/4 + (c1^2*conj(c2))/2);
%	the real is missing in the derivation, how come?
	cp(1,3) = c0*(c0^2 + 3/2*(abs(c1)^2 + abs(c2)^2)) + 3/8*(c1^2*conj(c2) + conj(c1)^2*c2);
	% diurnal
	cp(2,3) = c1*(3*c0^2 + 3/4*abs(c1)^2 + 3/2*abs(c2)^2) + 3*c0*c2*conj(c1)
	% semi
        cp(3,3) = c2*(3*c0^2 + 3/2*abs(c1)^2 + 3/4*abs(c2)^2) + 3/2*c0*c1^2;
	% ter
        cp(4,3) = 1/4*c1^3 + 3/4*c2^2*conj(c1) + 3*c0*c1*c2
	% quarter
        cp(5,3) = 3/4*c2*(c1^2 + 2*c0*c2);
	% quin
        cp(6,3) = 3/4*c1*c2^2;
	% seis
	cp(7,3) =  1/4*c2^3;

%	c(:,1) = [cp(1,1) cp(2,1) cp(3,1) cp(4,1) ];
%	c(:,2) = [cp(1,2) cp(2,2) cp(3,2) cp(4,2) ];
%	c(:,3) = [cp(1,3) cp(2,3) cp(3,3) cp(4,3) ];
%	c(:,4) = [cp(1,4) cp(2,4) cp(2,4) cp(3,4) ];
end % fourier_exp_power

