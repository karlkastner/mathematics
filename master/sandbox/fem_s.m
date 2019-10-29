% A-matrix
% Laplace-term
l : 0.5 * -1/(X(idx)-X(idx-1))
c : 0.5 * (1/(X(idx)-X(idx-1)) + 1/(X(idx+1)-X(idx)))
r : 0.5 * -1/(X(idx+1)-X(idx))
% potential term
l : X(idx-1)^2*log(X(idx)/X(idx-1))/(X(idx)-X(idx-1))^2 - 

% B-matrix
l: -1/3*(X(idx) - X(idx-1))
c:  1/3*(X(idx+1) - X(idx-1))
r: -1/3*(X(idx+1) - X(idx))

