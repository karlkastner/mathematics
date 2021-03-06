% Thu 11 Jan 10:52:01 CET 2018
%
%% inverse of stacked 4x4 matrices
function Ai = inv4x4(A)
	det = det4x4(A);
	Ai = [    A(2,2,:).*A(3,3,:).*A(4,4,:) ... % a11
		- A(2,2,:).*A(3,4,:).*A(4,3,:) ...
		- A(2,3,:).*A(3,2,:).*A(4,4,:) ...
		+ A(2,3,:).*A(3,4,:).*A(4,2,:) ...
		+ A(2,4,:).*A(3,2,:).*A(4,3,:) ...
		- A(2,4,:).*A(3,3,:).*A(4,2,:), ... % a12
		  A(1,2,:).*A(3,4,:).*A(4,3,:) ...
		- A(1,2,:).*A(3,3,:).*A(4,4,:) ...
		+ A(1,3,:).*A(3,2,:).*A(4,4,:) ...
		- A(1,3,:).*A(3,4,:).*A(4,2,:) ...
		- A(1,4,:).*A(3,2,:).*A(4,3,:) ...
		+ A(1,4,:).*A(3,3,:).*A(4,2,:), ... % a13
		  A(1,2,:).*A(2,3,:).*A(4,4,:) ...
		- A(1,2,:).*A(2,4,:).*A(4,3,:) ...
		- A(1,3,:).*A(2,2,:).*A(4,4,:) ...
		+ A(1,3,:).*A(2,4,:).*A(4,2,:) ...
		+ A(1,4,:).*A(2,2,:).*A(4,3,:) ...
		- A(1,4,:).*A(2,3,:).*A(4,2,:), ... % a14
		  A(1,2,:).*A(2,4,:).*A(3,3,:) ...
		- A(1,2,:).*A(2,3,:).*A(3,4,:) ...
		+ A(1,3,:).*A(2,2,:).*A(3,4,:) ...
		- A(1,3,:).*A(2,4,:).*A(3,2,:) ...
		- A(1,4,:).*A(2,2,:).*A(3,3,:) ...
		+ A(1,4,:).*A(2,3,:).*A(3,2,:); ...
	          A(2,1,:).*A(3,4,:).*A(4,3,:) ...  % a21
		- A(2,1,:).*A(3,3,:).*A(4,4,:) ...
		+ A(2,3,:).*A(3,1,:).*A(4,4,:) ...
		- A(2,3,:).*A(3,4,:).*A(4,1,:) ...
		- A(2,4,:).*A(3,1,:).*A(4,3,:) ...
		+ A(2,4,:).*A(3,3,:).*A(4,1,:), ... % a22
		  A(1,1,:).*A(3,3,:).*A(4,4,:) ...
		- A(1,1,:).*A(3,4,:).*A(4,3,:) ...
		- A(1,3,:).*A(3,1,:).*A(4,4,:) ...
		+ A(1,3,:).*A(3,4,:).*A(4,1,:) ...
		+ A(1,4,:).*A(3,1,:).*A(4,3,:) ...
		- A(1,4,:).*A(3,3,:).*A(4,1,:), ... % a23
		  A(1,1,:).*A(2,4,:).*A(4,3,:) ...
		- A(1,1,:).*A(2,3,:).*A(4,4,:) ...
		+ A(1,3,:).*A(2,1,:).*A(4,4,:) ...
		- A(1,3,:).*A(2,4,:).*A(4,1,:) ...
		- A(1,4,:).*A(2,1,:).*A(4,3,:) ...
		+ A(1,4,:).*A(2,3,:).*A(4,1,:), ... % a24
		  A(1,1,:).*A(2,3,:).*A(3,4,:) ...
		- A(1,1,:).*A(2,4,:).*A(3,3,:) ...
		- A(1,3,:).*A(2,1,:).*A(3,4,:) ...
		+ A(1,3,:).*A(2,4,:).*A(3,1,:) ...
		+ A(1,4,:).*A(2,1,:).*A(3,3,:) ...
		- A(1,4,:).*A(2,3,:).*A(3,1,:); ...
	          A(2,1,:).*A(3,2,:).*A(4,4,:) ... % a31
		- A(2,1,:).*A(3,4,:).*A(4,2,:) ...
		- A(2,2,:).*A(3,1,:).*A(4,4,:) ...
		+ A(2,2,:).*A(3,4,:).*A(4,1,:) ...
		+ A(2,4,:).*A(3,1,:).*A(4,2,:) ...
		- A(2,4,:).*A(3,2,:).*A(4,1,:), ... % a32
		  A(1,1,:).*A(3,4,:).*A(4,2,:) ...
		- A(1,1,:).*A(3,2,:).*A(4,4,:) ...
		+ A(1,2,:).*A(3,1,:).*A(4,4,:) ...
		- A(1,2,:).*A(3,4,:).*A(4,1,:) ...
		- A(1,4,:).*A(3,1,:).*A(4,2,:) ...
		+ A(1,4,:).*A(3,2,:).*A(4,1,:), ... % a33
		  A(1,1,:).*A(2,2,:).*A(4,4,:) ...
		- A(1,1,:).*A(2,4,:).*A(4,2,:) ...
		- A(1,2,:).*A(2,1,:).*A(4,4,:) ...
		+ A(1,2,:).*A(2,4,:).*A(4,1,:) ...
		+ A(1,4,:).*A(2,1,:).*A(4,2,:) ...
		- A(1,4,:).*A(2,2,:).*A(4,1,:), ... % a34
		  A(1,1,:).*A(2,4,:).*A(3,2,:) ...
		- A(1,1,:).*A(2,2,:).*A(3,4,:) ...
		+ A(1,2,:).*A(2,1,:).*A(3,4,:) ...
		- A(1,2,:).*A(2,4,:).*A(3,1,:) ...
		- A(1,4,:).*A(2,1,:).*A(3,2,:) ...
		+ A(1,4,:).*A(2,2,:).*A(3,1,:); ... % a41
		  A(2,1,:).*A(3,3,:).*A(4,2,:) ...
		- A(2,1,:).*A(3,2,:).*A(4,3,:) ...
		+ A(2,2,:).*A(3,1,:).*A(4,3,:) ...
		- A(2,2,:).*A(3,3,:).*A(4,1,:) ...
		- A(2,3,:).*A(3,1,:).*A(4,2,:) ...
		+ A(2,3,:).*A(3,2,:).*A(4,1,:), ... % a42
		  A(1,1,:).*A(3,2,:).*A(4,3,:) ...
		- A(1,1,:).*A(3,3,:).*A(4,2,:) ...
		- A(1,2,:).*A(3,1,:).*A(4,3,:) ...
		+ A(1,2,:).*A(3,3,:).*A(4,1,:) ...
		+ A(1,3,:).*A(3,1,:).*A(4,2,:) ...
		- A(1,3,:).*A(3,2,:).*A(4,1,:), ... % a43
		  A(1,1,:).*A(2,3,:).*A(4,2,:) ...
		- A(1,1,:).*A(2,2,:).*A(4,3,:) ...
		+ A(1,2,:).*A(2,1,:).*A(4,3,:) ...
		- A(1,2,:).*A(2,3,:).*A(4,1,:) ...
		- A(1,3,:).*A(2,1,:).*A(4,2,:) ...
		+ A(1,3,:).*A(2,2,:).*A(4,1,:), ... % a44
		  A(1,1,:).*A(2,2,:).*A(3,3,:) ...
		- A(1,1,:).*A(2,3,:).*A(3,2,:) ...
		- A(1,2,:).*A(2,1,:).*A(3,3,:) ...
		+ A(1,2,:).*A(2,3,:).*A(3,1,:) ...
		+ A(1,3,:).*A(2,1,:).*A(3,2,:) ...
		- A(1,3,:).*A(2,2,:).*A(3,1,:)];
	% normalise
	Ai = bsxfun(@times,1./det,Ai);
end % inv4x4

