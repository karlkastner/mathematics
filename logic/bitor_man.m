% Sun Jun  8 17:11:00 WIB 2014
% Karl Kastner, Berlin
%
%% bitwise OR of the numbers of the columns of A
%%
%% input:
%%	A (positive integer)
function B = bitor_man(A)
	B = A(1,:);
	for idx=2:size(A,1)
		B = bitor(B,A(idx,:));
	end
end

