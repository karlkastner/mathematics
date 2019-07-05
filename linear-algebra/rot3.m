% Mon Jun 24 12:59:43 UTC 2013
% Karl Kästner, Berlin
%
% 3D rotation matrix from angle
function R = rotM(phi,psi,theta)
	R =  [       1             0         0;
                     0      cos(phi) -sin(phi);
                     0      sin(phi)  cos(phi)] ...
           * [cos(psi)             0 -sin(psi);
                     0             1         0;
              sin(psi)             0 cos(psi)] ...
           * [cos(theta) -sin(theta)        0;
              sin(theta)  cos(theta)        0;
                       0           0        1];
end % rotM
