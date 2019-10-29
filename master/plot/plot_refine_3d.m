% Tue Aug  7 15:06:16 MSK 2012
% Karl Kästner, Berlin

function plot_refine_3d(pflag)

P = 0.5*[0 0 0;
         0 0 2;
         0 2 0;
         2 0 0;
	 0 0 1;
         0 1 0;
         1 0 0;
         0 1 1;
         1 0 1;
         1 1 0];

p1 = 1;
p2 = 2;
p3 = 3;
p4 = 4;
p12 = 5;
p13 = 6;
p14 = 7;
p23 = 8
p24 = 9;
p34 = 10; 

M = [    1.0000         0         0         0
         0    1.0000         0         0
   -0.3660    0.5000    0.8660         0
   -0.6830    0.5000    0.4330    0.7500];
P = P - 0.25;
P = [ones(10,1) P]*M;
P=P(:,2:4)

% 0 - split
T0 = [ p1 p2 p3 p4]

% 8 - split
T8 = [ p1 p12 p13 p14
       p2 p12 p23 p24
       p3 p13 p23 p34
       p4 p14 p24 p34
    p14, p23, p24, p34, 
    p13, p23, p14, p34,  
    p12, p14, p23, p24,  
    p12, p13, p23, p14 ];

% 4 - split
T4 = [p1 p12 p13 p4;
      p2 p12 p23 p4;
      p3 p13 p23 p4;
      p12 p13 p23 p4];

% 2 - split
T2 = [ p1 p13 p2 p4;
       p13 p2 p3 p4];

% 3-split
T3 = [ p3  p13  p23  p4;
       p13 p23  p2  p4;
       p13 p2  p1   p4];

clf
colormap([1 1 1])
alpha = 0;
subplot(1,5,1);
display_3d(P,T0,[],1,alpha)

subplot(1,5,2);
display_3d(P,T8,[],1,alpha)

subplot(1,5,3);
display_3d(P,T4,[],1,alpha)

subplot(1,5,4);
display_3d(P,T2,[],1,alpha)

subplot(1,5,5);
display_3d(P,T3,[],1,alpha)

for idx=1:5
	subplot(1,5,idx);
	set(gca,'visible','off')
	view([-72 22])
	zoom(1.3)
%	subplot(2,5,5+idx);
	h = title(idx-1)
	set(findall(gca, 'type', 'text'), 'visible', 'on')
	pos = get ( h, 'position' )
	set ( h, 'position', [0 +0.125 -0.5] )
end

if (nargin() > 0 && pflag > 0)
	preparePrint();
	print -deps ../img/refinement-3d.eps
	system('epstopdf ../img/refinement-3d.eps')
end

end

