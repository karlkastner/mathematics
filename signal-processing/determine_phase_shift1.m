% 2021-12-22 13:30:25.028490215 +0100
%
%% average phase and phase shift per time step of a train of waves
%
% rows    : time
% columns : space
function [phi, dphi_dt] = determine_phase(b,dt)
	nt = size(b,1);
	% 
	F = fft(b.').';

	% time differences
	% TODO use here the df_dt of the Rietkerk-model
	% df_dt  = diff(F,[],1)./dt;

	tdx   = (1:nt)';
	%dt  = (tdx(2)-tdx(1));

	% phase angle
	phi   = angle(F);

	% phase shift with time
	dphi_dt = diff(phi)./dt;

	dphi_dt = wrapTo2Pi(dphi_dt);

	% average over all time steps
	dphi_dt = mean(dphi_dt,1);

	% subtract shift over time
	phi = phi - dt*cvec(tdx)*rvec(dphi_dt);
	phi = wrapToPi(phi);

	% get the offset
	phi = mean(phi,1);
end
