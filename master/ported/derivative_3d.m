% Thu Jul 12 16:06:08 MSK 2012
% Karl Kästner, Berlin

% returns the first derivative
% in : coefficients of 3D polynomial of degree p
% out : coefficients of 3D polynomial of degree p-1
% expects input C to be polynomial coefficients in rows
% 1 2 3 4 5  6  7  8  9  10
% 1 x y z x² y² z² xy xz yz ...
function [C_dx C_dy C_dz] = derivative_3d(C)
	switch (size(C,1))
		case (4)
			C_dx = C(2,:);
			C_dy = C(3,:);
			C_dz = C(4,:);
		case (10)
			C_dx = [  C( 2,:); 2*C( 5,:);   C( 8,:);   C( 9,:) ];
                        C_dy = [  C( 3,:);   C( 8,:); 2*C( 6,:);   C(10,:) ];
			C_dz = [  C( 4,:);   C( 9,:);   C(10,:); 2*C( 7,:) ];
		case (20)
			C_dx = [ C( 2,:); 2*C( 5,:);   C( 8,:);   C( 9,:); 3*C(11,:);   C(16,:);   C(18,:); 2*C(14,:); 2*C(15,:);   C(20,:) ];
                        C_dy = [ C( 3,:);   C( 8,:); 2*C( 6,:);   C(10,:);   C(14,:); 3*C(12,:);   C(19,:); 2*C(16,:);   C(20,:); 2*C(17,:) ];
			C_dz = [ C( 4,:);   C( 9,:);   C(10,:); 2*C( 7,:);   C(15,:);   C(17,:); 3*C(13,:);   C(20,:); 2*C(18,:); 2*C(19,:) ];
		case (35)
			C_dx = [    C( 2,:); 2*C( 5,:);   C( 8,:);   C( 9,:); 3*C(11,:);   C(16,:);   C(18,:); 2*C(14,:); 2*C(15,:);   C(20,:); ...
				  4*C(21,:);   C(26,:);   C(28,:); 3*C(24,:); 3*C(25,:); 2*C(30,:);   C(34,:); 2*C(31,:);   C(35,:); 2*C(33,:)];
                        C_dy = [    C( 3,:);   C( 8,:); 2*C( 6,:);   C(10,:);   C(14,:); 3*C(12,:);   C(19,:); 2*C(16,:);   C(20,:); 2*C(17,:); ...
				    C(24,:); 4*C(22,:);   C(29,:); 2*C(30,:);   C(33,:); 3*C(26,:); 3*C(27,:);   C(35,:); 2*C(32,:); 2*C(34,:)];
			C_dz = [    C( 4,:);   C( 9,:);   C(10,:); 2*C( 7,:);   C(15,:);   C(17,:); 3*C(13,:);   C(20,:); 2*C(18,:); 2*C(19,:);
				    C(25,:);   C(27,:); 4*C(23,:);   C(33,:); 2*C(31,:);   C(34,:); 2*C(32,:); 3*C(28,:); 3*C(29,:); 2*C(35,:);];
	otherwise
		error('derivative_3d','degree of the polynomial is not yet supported')
	end
end % derivative_3d

