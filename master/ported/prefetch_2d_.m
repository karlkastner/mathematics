% Sun Jun 24 19:58:01 MSK 2012
% Karl Kästner, Berlin

% fetches the element point coordinates
% constructs the vandermonde matrix
% inverts the vandermonde matrix (yields the test functions)

function [Phi area h_side s_angle C] = prefetch_2d(mesh)
	P = mesh.P;
	T = mesh.T;
	lt1 = size(T,1);
	lt2 = size(T,2);
	% allocate memory
	Phi = zeros(lt1,lt2,lt2);
	% determine number of the power
	nv = -1.5 + sqrt(2*lt2 + 0.25);
	for tdx=1:lt1
		%  prefetch triangle points
		A_ = P(T(tdx,:),:);
		
		% construct the Vandermonde matrix
		Va = vander_2d(A_,nv);
		% invert the Vandermonde matrix
		Phi(tdx,:,:) = inv(Va);

		% calculate element centre coordinate
		C(tdx,:) = 1/3*(P(T(tdx,1),:) + P(T(tdx,2),:) + P(T(tdx,3),:));

		% calculate element determinant
                A = [   1 P(T(tdx,1),:);
                        1 P(T(tdx,2),:);
                        1 P(T(tdx,3),:) ];
		determinant(tdx) = abs(det(A));

		% todo: warn for clockwise (area<0) and degenerated triangles (area/h_max < eps)

		% length of side opposit the points
		% todo, vectorise
		a = sqrt((P(T(tdx,2),1) - P(T(tdx,3),1))^2 + (P(T(tdx,2),2) - P(T(tdx,3),2))^2);
		b = sqrt((P(T(tdx,1),1) - P(T(tdx,3),1))^2 + (P(T(tdx,1),2) - P(T(tdx,3),2))^2);
		c = sqrt((P(T(tdx,1),1) - P(T(tdx,2),1))^2 + (P(T(tdx,1),2) - P(T(tdx,2),2))^2);
		h_side(tdx,1) = a;
		h_side(tdx,2) = b;
		h_side(tdx,3) = c;

		% sine of interior angles at the points
		% todo, vectorise
		sin_a = 0.5*determinant(tdx)/(0.5*b*c);
		sin_b = 0.5*determinant(tdx)/(0.5*a*c);
		sin_c = 0.5*determinant(tdx)/(0.5*a*b);

		s_angle(tdx,1) = sin_a;
		s_angle(tdx,2) = sin_b;
		s_angle(tdx,3) = sin_c;
	end
end

