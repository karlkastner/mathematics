% Thu 22 Nov 16:49:42 CET 2018

% in theory, one can also transform so that p1 = [-1,0] and p2 [+1,0]
% fit a quadratic and than select the point at x = 0

if (1)
syms x1 x2 y1 y2 x0 y0 l

l1 = sqrt((x1-x0)^2 + (y1-y0)^2)
l2 = sqrt((x2-x0)^2 + (y2-y0)^2)

eq = [l1^2 - l2^2, l1+l2 - l]
s  = solve(eq,[x0,y0])

s.x0 = simplify(s.x0);
s.y0 = simplify(s.y0)

clear fun
fun.x0 = matlabFunction(s.x0)
fun.y0 = matlabFunction(s.y0)

save('mat/optimum-chord-length.mat','fun');
end
% two solutions, one on the far side and one on the near side
% check which is nearer, choose nearer



