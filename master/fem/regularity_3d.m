% Wed Jun 13 20:45:42 MSK 2012
% Karl Kästner, Berlin

% volume check, degeneracy check
function [v_sum a_sum h_eff_max volume area h_eff h_max] = regularity_3d(P,T,Bc)

%%{	T = [1 2 3 4];
%	P = [rand(4,3)]
%	tdx = 1;
%		A = [  1 P(T(tdx,1),:)
%		       1 P(T(tdx,2),:)
%		       1 P(T(tdx,3),:)
%		       1 P(T(tdx,4),:)]';
%	E = [A(:,1) - A(:,2) A(:,1) - A(:,3) A(:,1) - A(:,4) A(:,2) - A(:,3) A(:,2) - A(:,4) A(:,3) - A(:,4)]
%	E = E(2:end,:)
%		svd(A)
%		svd(E)
%	pause
%	for idx=1:6, E_(:,idx) = E(:,idx)/norm(E(:,idx)); end
%	svd(E)
%	det(A)
%	A1= [  1 P(T(tdx,2),:)
%	       1 P(T(tdx,3),:)
%	       1 P(T(tdx,4),:)];
%	A2= [  1 P(T(tdx,1),:)
%	       1 P(T(tdx,3),:)
%	       1 P(T(tdx,4),:)];
%	A3= [  1 P(T(tdx,1),:)
%	       1 P(T(tdx,2),:)
%	       1 P(T(tdx,4),:)];
%	A4= [  1 P(T(tdx,1),:)
%	       1 P(T(tdx,2),:)
%	       1 P(T(tdx,3),:)]
%		det(A)/abs(det(A1)*det(A2)*det(A3)*det(A4))^0.25;
%	pause
%%}

	for tdx=1:size(T,1)
		A = [  1 P(T(tdx,1),:)
		       1 P(T(tdx,2),:)
		       1 P(T(tdx,3),:)
		       1 P(T(tdx,4),:)];
	        E = [P(T(tdx,1))-P(T(tdx,2)),...
                     P(T(tdx,1))-P(T(tdx,3)),...
                     P(T(tdx,1))-P(T(tdx,4)),...
                     P(T(tdx,2))-P(T(tdx,3)),...
                     P(T(tdx,2))-P(T(tdx,4)),...
                     P(T(tdx,3))-P(T(tdx,4))];
		%degen(idx)  = min(svd(E));
		volume(tdx) = 1/6*abs(det(A));
		h_max(tdx)  = sqrt(max(sum(E.^2)));
		h_eff(tdx)  = h_max(tdx)^4 / volume(tdx);
		%h_eff(tdx) = h_max^3 / volume(tdx);
		%degen(idx) = volume(idx)/abs(det(A1)*det(A2)*det(A3)*det(A4))^0.25;
	end
	% check boundary area
	area = zeros(size(Bc,1),1);
	a_sum = 0;
	for idx=1:size(Bc,1)
		%A = [  1 P(Bc(idx,1),:)
		%       1 P(Bc(idx,2),:)
		%       1 P(Bc(idx,3),:)];
		%area(idx,1) = 0.5*abs(det(A));
		%area_ = 0.5*norm( cross( P(Bc(idx,1),:)-P(Bc(idx,2),:), P(Bc(idx,1),:)-P(Bc(idx,3),:) ) )
		area(idx,1) = 0.5*sqrt( ...
			 det([ ones(3,1) P(Bc(idx,1:3),1) P(Bc(idx,1:3),2) ])^2 ...
			+ det([ ones(3,1) P(Bc(idx,1:3),1) P(Bc(idx,1:3),3) ])^2 ...
			+ det([ ones(3,1) P(Bc(idx,1:3),2) P(Bc(idx,1:3),3) ])^2 );
		if (area(idx,1) < 1e-12)
			Bc(idx,:)
			P(Bc(idx,1),:)
			P(Bc(idx,2),:)
			P(Bc(idx,3),:)
			[(1:size(P,1))' P]
			error('here');
		end
		a_sum = a_sum + area(idx,1);
	end

	v_sum = sum(volume);
	h_eff_max = max(h_eff);
end % function

