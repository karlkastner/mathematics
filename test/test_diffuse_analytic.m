% Mon 30 May 10:15:38 CEST 2022
L = 10;
dx = 0.1;
nx = L/dx;
x = linspace(-L/2,L/2,nx)';
%x0 = L/2;
%x = x-x0;
x0 = -4;
x0 = 0;
% 1d

r = x;
y0 = zeros(nx,1);
[mv,mdx] = min(abs(x-x0));
x0 = x(mdx);
y0(mdx)=1;

t = 0.0001;
t = 0.001;
%t = 0.01;
%t = 0.1;
%t = 0.1;
e = 1;
clear y

% 1d
if (0)

D2     = derivative_matrix_2_1d(nx,L,2,'circular');
fun    = @(t,y) e*D2*y;
[tvec, y] = ode23(fun,[0,t],y0);

%clf
%imagesc(y)
%pause

y = y(end,:)';

y(:,2) = diffuse_analytic(t,y0,L,e,false);
y(:,3) = diffuse_analytic(t,y0,L,e,true);

dx = x(2)-x(1);
y(:,end+1) = 1/sqrt(4*pi*e*t)*exp(-abs(x-x0).^2/(4*e*t))*dx;
if (0)
y(:,end+1) = circshift(cconv(y0,y(:,end),nx),round(nx/2));
end

figure(1);
subplot(2,2,1)
plot(x,[y],'.-');
%^plot(x,[y0,y],'.-');
vline([x0])
legend('ode','da','daf','analytic')
%legend('y0','ode','da','daf','analytic')

subplot(2,2,2)
plot(x,[y-y(:,end)],'.-');

else
	y0=zeros(nx,nx);
	y0(round(nx)/2,round(nx)/2) = 1;
%	y0 = rand(nx,nx);

	% 2d
	[Dx,Dy,D2x,Dxy,D2y] = derivative_matrix_2d([nx,nx],[L,L],2,'circular');
	D2 = D2x+D2y;

	fun    = @(t,y) e*D2*y;
	[tvec, y] = ode23(fun,[0,t],y0(:));
	yy = reshape(y(end,:),nx,nx);	

	yy(:,:,2) = diffuse_analytic(t,y0,[L,L],e,true);
	yy(:,:,3) = diffuse_analytic(t,y0,[L,L],e,false);

	clf
	subplot(2,2,1)
	plot(squeeze(yy(:,round(end/2),:)));
	subplot(2,2,2)
	plot(squeeze(yy(round(end/2),:,:)));

	rms(flat(yy(:,:,3)-yy(:,:,2)))
	
end
