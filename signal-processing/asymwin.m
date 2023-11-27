% 2016-03-03 19:10:36.799831064 +0100
% 2016-03-03 18:26:54.259194140 +0100
% 2016-03-03 18:04:12.326400388 +0100
% 2016-03-03 18:01:23.219263480 +0100
%
%  This program is free software: you can redistribute it and/or modify
%  it under the terms of the GNU General Public License as published by
%  the Free Software Foundation, either version 3 of the License, or
%  (at your option) any later version.
%
%  This program is distributed in the hope that it will be useful,
%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU General Public License for more details.
%
%  You should have received a copy of the GNU General Public License
%  along with this program.  If not, see <https://www.gnu.org/licenses/>.
%

%% creates asymmetrical filter windows
%% filter will always have negative weights
function y = awin(x,mu,L)
	l = max(mu-L,0);
	r = mu+L;
%	r = l+2*L;
	c = mu;
	b = [1;  
%             mu;
	     0;   % 1 value at left is zero
             0;   % 2 value at right is zero
             0;   % 3 d at left
%             0;   % 4 d at right
           %  0;   % 5 derivative at c
	   %  0;
           % 0;
             ];

ll = c-l;
rr = r-c;

variant = 1;
switch (variant)
case {-1} % piecewise linear kernel
	o = 1;
	A=[ vanderi_1d(l,c,o),vanderi_1d(c,r,o);
	    ll/(ll+rr)*murow(l,c,o),rr/(ll+rr)*murow(c,r,o);
	    vander_1d(l,o),zeros(1,o+1);
	    zeros(1,o+1),vander_1d(r,o)];
case {0}
	o = 3;
	A = [vanderi_1d(l,r,o)
             murow(l,r,o);
	     vander_1d(l,o);
	     vanderd_1d(c,o,1).;
             ];
case {1}
	% arbitrary derivative at left and right side, mean is mu
	o = 4;
	A = [vanderi_1d(l,r,o)
%             murow(l,r,o);
	     vanderd_1d(c,o,1).;
	     vander_1d(l,o);
	     vander_1d(r,o);
             ];
case {2}
	% zero d left, mean is mu
	o = 4;
	A = [
              vanderi_1d(l,r,o)
%	     murow(l,r,o); % mean
	     vander_1d(l,o); % f(l)
             vander_1d(r,o); % f(r)
             vanderd_1d(c,o,1).; % f'(c)
             vanderd_1d(l,o,1).; % f'(l)
	     ];  % int
case {2}
	% zero derivative left and right, mean is mu
	o = 7;
	A = [
              vanderi_1d(l,r,o)
	      murow(l,r,o)
	      vander_1d(l,o);
              vander_1d(r,o);
              vanderd_1d(c,o,1);
              vanderd_1d(l,o,1);
              vanderd_1d(r,o,1);
	];
case {3}
	% zero second derivative l and r and second derivative 0 at right
	o = 9;
	A = [
             vanderi_1d(l,r,o);
	     murow(l,r,o);
             vander_1d(l,o);
             vander_1d(r,o);
             vanderd_1d(c,o,1);
             vanderd_1d(l,o,1);
             vanderd_1d(r,o,1);
             vanderd_1d(l,o,2)
             vanderd_1d(r,o,2)
	];
case {4}
	A = [vander_1d(l,0);
             vander_1d(r,0);
             vanderd_1d(l,1,0);
             vanderd_1d(c,1,0);
             vanderd_1d(r,1,0);
             vanderi_1d(l,c,0)
             vanderi_1d(c,r,0)
             vanderd_1d(r,o,2).
             0, 0,  0,   6,    24*r,   90*r.^2, 168*r.^3,280*r^4,432*r^5];
case {5}
	A = [vander_1d(l,0);
             vander_1d(r,0);
             %vanderd_1d(l,1,0).;
             vanderd_1d(c,1,0);
             vanderd_1d(r,1,0);
             vanderi_1d(l,c,0)
             vanderi_1d(c,r,0)
             vanderd_1d(r,o,2).
             0, 0,  0,   6,    24*r,   90*r.^2, 168*r.^3,280*r^4,432*r^5
             0, 0,  0,   0,    24,   180*r,504*r.^2,1120*r^3,21602*r^4; ];
end

	c = A \ b;
%	A
%	A*c
%	c
	y = zeros(size(x));
	% regress coefficients
	if (variant >= 0)
		% expand
		fdx = (x>l & x<r);
		A = vander_1d(x(fdx),length(c)-1);
		y(fdx) = A*c;
	else
		fdx = x>l & x < mu;
		A = vander_1d(x(fdx),o);
		y(fdx) = A*c(1:o+1);
		fdx = x>=mu & x < r;
		A = vander_1d(x(fdx),o);
		y(fdx) = A*c(o+2:end);
	end
%	y = max(y,0);
%	y = y/sum(y);

	if (sum(y) <=0) %l+L>mu+sqrt(eps))
	[mu L]
	[sum(y)*(x(2)-x(1)),sum(y.*x)*(x(2)-x(1))./mu]
%	l+L
%	mu
	figure(1)
	clf
	plot(y)
	pause
	end
end

% 2016-03-03 18:15:02.635692885 +0100
% Karl Kastner, Berlin

function row = murow(l,r,order)
	row = zeros(1,order+1,class(l));
	for idx=1:order+1
		row(idx) = 1./(idx+1)*(r.^(idx+1)-l.^(idx+1));
	end
end

