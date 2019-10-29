function A = cut_out(n, r0)
h= 1/(n+1);
x = 2*(((1:n)/(n+1))-0.5);
X = diag(sparse(x));
I = speye(n);

[XX YY] = meshgrid(x,x);
%XX = kron(I,X);
%YY = kron(X,I);
RR = sqrt(XX.^2 + YY.^2);
ID = reshape(cumsum(reshape(RR > r0,1,n^2)),n,n);

ik = [5 -10 10 -5 1];
ik = [4 -6 4 -1];
ik = [3 -3 1];
ik = [2 -1];

%XX, YY, RR, ID

m = ID(end,end);
A = spalloc(m, m, 5*ceil(sqrt(m)));
Rp = spalloc(m, m, m);
for idx=1:n
 for jdx=1:n
	if (RR(idx,jdx) > r0)
		id = ID(idx,jdx);
		Rp(id,id) = RR(idx,jdx)^2;
		A(id,id) = -4; % + 2/RR(idx);
		% left
		if (idx > 1)
			if (RR(idx-1,jdx) > r0)
				A(id,ID(idx-1,jdx)) =  1;
			else
				% interpolate
				for kdx=1:length(ik)
					jd = ID(idx+kdx-1,kdx);
					A(id,jd) = A(id,jd) + ik(kdx);
				end
			end
		end
		% right
		if (idx < n)
			if (RR(idx+1,jdx) > r0)
				A(id,ID(idx+1,jdx)) =  1;
			else
				% interpolate
				for kdx=1:length(ik)
					jd = ID(idx-kdx+1,kdx);
					%A(id,jd) = A(id,jd) + ik(length(ik)-kdx+1);
					A(id,jd) = A(id,jd) + ik(kdx);
				end
			end
		end
		% top
		if (jdx > 1)
			if (RR(idx,jdx-1) > r0)
				A(id,ID(idx,jdx-1)) =  1;
			else
				% interpolate
				for kdx=1:length(ik)
					jd = ID(idx,jdx+kdx-1);
					A(id,jd) = A(id,jd) + ik(kdx);
				end
			end
		end
		% bottom
		if (jdx < n)
			if (RR(idx,jdx+1) > r0)
				A(id,ID(idx,jdx+1)) =  1;
			else
				% interpolate
				for kdx=1:length(ik)
					jd = ID(idx,jdx-kdx+1);
					%A(id,jd) = A(id,jd) + ik(length(ik)-kdx+1);
					A(id,jd) = A(id,jd) + ik(kdx);
				end
			end
		end
	end % RR 
 end % jdx
end % idx

%full(A)
A = A/h^2 + Rp;

end % cut_out

