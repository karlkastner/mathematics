% 2024-09-27 12:24:29.565923496 +0200
% Karl Kastner, Berlin
%% assemble a matrix, similar to sparse but works with symbolic variables
function A = assemble(buf1,buf2,buf3,n1,n2)
	for idx=1:length(buf1)
		A(buf1(idx),buf2(idx)) = buf3(idx);
	end

end

