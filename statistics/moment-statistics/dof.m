% Wed Jun 25 19:36:32 WIB 2014
% Karl Kastner, Berlin
%
%% mininum number of support points
%% for a polynomial of degree order in dim dimensions
%
function n = dof(dim,order)
	n = 1/factorial(dim)*prod((1:dim)+order);
end

