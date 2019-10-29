P = [0 2;
     1 2;
     2 2;
     0 1;
     1 1;
     2 1;
     0 0;
     1 0;
     2 0;
     0.5 1.5;
     1   1.5;
     0.5 1];
     
T = [1 2 4;
     2 5 4;
     2 3 5;
     3 6 5;
     4 5 7;
     5 8 7;
     6 9 8];

for idx=1:size(T,1)
	line(	[ P(T(idx,1),1) P(T(idx,2),1) P(T(idx,3),1) P(T(idx,1),1) ], ...
		[ P(T(idx,1),2) P(T(idx,2),2) P(T(idx,3),2) P(T(idx,1),2) ], 'color', [0 0 0],'linewidth',2);
end

T = [10 11 12];
for idx=1:size(T,1)
	line(	[ P(T(idx,1),1) P(T(idx,2),1) P(T(idx,3),1) P(T(idx,1),1) ], ...
		[ P(T(idx,1),2) P(T(idx,2),2) P(T(idx,3),2) P(T(idx,1),2) ], 'color', [0 0 0],'linestyle','--','linewidth',2);
end

T = [1 2 10; % 10 4 2;
     2 3 11; % 3 5 11;
     5 7 12; % 12 7 4
     ];

for idx=1:size(T,1)
	line(	[ P(T(idx,1),1) P(T(idx,2),1) P(T(idx,3),1) P(T(idx,1),1) ], ...
		[ P(T(idx,1),2) P(T(idx,2),2) P(T(idx,3),2) P(T(idx,1),2) ], 'color', [0 0 0],'linestyle','--','linewidth',2);
end

set(gca,'xtick',[])
set(gca,'ytick',[])
axis([-0.005 2.005 -0.005 2.005])
axis square
axis off
box off
print -deps ../img/refine_2d.eps 
%axis tight

