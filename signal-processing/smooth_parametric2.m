% 2018-02-25 14:09:11.339322390 
%% parametrically smooth the curve
function P2 = smooth__(P1,P2,P3,relax)
	p = relax;
	for idx=1:size(P1,2)
		P2(:,idx) = (1-p)*P2(:,idx) + p*smooth___(P1(:,idx),P2(:,idx),P3(:,idx));
	end
end

% smooth location of P2
function P2 = smooth___(P1,P2,P3)
	P0 = P1;  %+[1;0];
	% translate
	P1 = P1-P0;
	P2 = P2-P0;
	P3 = P3-P0;
	% scale
	s  = norm(P3);
	P3 = P3/s;
	P2 = P2/s;
	% rotate
	R = [ P3(1),P3(2),
             -P3(2),P3(1)];
	% test
%R*P3
%ause
	P2 = R*P2;
	%T  = ([-1,0; 1,0]*inv([P1,P3]))
	%P2 = T*P2;
	if (P2(1) < sqrt(eps) || P2(1) > 1-sqrt(eps))
		P2 = [0.5;0];
	else
		% set up the quadratic polynomial
		A = vander_1d([0,P2(1),1],2);
%if (all(flat(isfinite(A))) && cond(A) > 1e3)
%	P2
%pause
%end
		c = A \ [0;P2(2);0];
		% evaluate polynomial at 0.5
		A = vander_1d(0.5,2);
		P2 = [0.5;A*c]; %[0; c(3)];
	end
	% rotate  back
	P2 = R'*P2;
	% scale back
	P2 = s*P2;
	% translate back
	P2 = P2+P0;
end

