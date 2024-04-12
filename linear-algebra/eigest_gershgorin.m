% Thu  4 Jan 13:33:33 CET 2024
function E = eigest_gershgorin(A)
	D = diag(A);
	B = abs(sum(A,2)-D);
	E = [D-B,D+B];
end

