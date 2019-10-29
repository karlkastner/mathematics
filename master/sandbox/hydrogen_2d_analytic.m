% Wed Jun 13 11:45:18 MSK 2012
% Karl Kästner, Berlin

% f = @(x,y) exp(-sqrt(x.^2 + y.^2)/2).^2
%  quad2d(f, -100, 100, -100, 100)/(2*pi) 

h=0.01
X =0:h:10;
% probability
y = (1/(2*pi))*exp(-X/2).^2
sum(2*pi*X.*y)
plot(X,y)


function v = hydrogen_analytic_2d(P, shell)
	% analytic solution of the ground state in 2D
% TODO the shell influence is not correct (this is the case for 3D)!!!
	v = 0.5*pi*exp(-2/level*sqrt(P(:,1).^2 + P(:,2).^2));
end

