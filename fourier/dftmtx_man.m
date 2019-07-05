% Thu  4 May 12:35:41 CEST 2017
% Karl Kastner, Berlin
%
%% fourier matrix in matlab style with a limited number of rows,
%% columns of higher frequencies are omitted
%%
%% input :
%% n  : number of samples
%% nr : number of columns
%%
%% output :
%% F : fourier matrix
function F = dftmtx_man(n,nr)
	F = zeros(nr,n);
	F(1,1:n) = 1;
	id = (0:n-1);
	I  = (-1i).^(id*4/n);
	for idx=2:nr
		F(idx,:) = F(idx-1,:).*I;
	end
end
