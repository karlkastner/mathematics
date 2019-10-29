% Sat Jun 16 14:40:42 MSK 2012
% Karl Kästner, Berlin

% calculate partial derivatives inside the the elements
function dV = partial_derivative_2d(T, P, V)
	lt = size(T,1);
	dV = zeros(lt,2);
	for idx=1:lt
		A = [   1 P(T(idx,1),:);
                        1 P(T(idx,2),:);
                        1 P(T(idx,3),:) ];
		c = A \ [V(T(idx,1)); V(T(idx,2)); V(T(idx,3))];
		dV(idx,:) = c(2:3).';
	end
end

