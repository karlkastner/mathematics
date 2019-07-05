% 2016-01-28 14:51:38.868683343 +0100
%% estimation of the condition number
function c = condest_(A)
	t = [];  
	c = normest1(A,t);
	c = c*norm(A,1);
%	c = normest1(@condestf,t);
%   A_norm = norm(A,1);
%   c = Ainv_norm*A_norm;
%end
%v = v/norm(v,1);

    function f = condestf(flag, x)
        %CONDESTF   Function used by CONDEST.        
        if isequal(flag,'dim')
            f = max(size(A));
        elseif isequal(flag,'real')
            f = isreal(A);
        elseif isequal(flag,'notransp')
            f = A*x;
        elseif isequal(flag,'transp')
            f = A'*x;
        end
    end
end

