function Q = filter_twosided(Q,rho1,rho2)
Q_ = flipud(Q);
sum(Q)
for idx=3:size(Q,1)
	Q(idx,:) =   (1-rho1-rho2)*Q(idx,:) ...
		   + rho1*Q(idx-1,:) ...
		   + rho2*(2*Q(idx-1,:)-1*Q(idx-2,:));
	Q_(idx,:) = (1-rho1-rho2)*Q_(idx,:) ...
		   + rho1*Q_(idx-1,:) ...
		   + rho2*(2*Q_(idx-1,:)-1*Q_(idx-2,:));
%		   + rho*Qp;
%	Qp = 2*Q_(idx-1,:)-1*Q_(idx-2,:);
%	Q_(idx,:) = (1-rho)*Q_(idx,:) + rho*Qp;
end
Q_ = flipud(Q_);
Q = 0.5*(Q + Q_);
end
