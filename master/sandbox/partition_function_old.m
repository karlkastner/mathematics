% Tue Oct  4 20:47:45 MSD 2011
% Karl Kästner, Berlin

% todo degeneracy
% todo reduced mass of the electron
% todo exact constants incl. eps0 as 1/(my c^2)

function [ZZ EE Z__] = partition_function(n)
	% zero Kelvin in Celsius
	K0 = -273.15; % [K0] = °C
	% boltzmann factor
	kB = 1.3806488e-23; % [kB] = 1 J/K
	% room temperature
	T = 20.0 - K0; % [T] = Kelvin
	% boltzmann factor
	beta = 1.0/(kB * T) % [kB] = J/K
	% mass of electron
	me = 9.109e-31; % [me] = kg
	mr = me;
	m_star = me;
	% planck's constant
	h = 6.626e-34; % [h] = Js
	% charge of electron
	e = 1.602e-19; % [e] = As
	% permitivity in vacuum
	eps0 = 8.854e-12; % [eps0] = F/m

	Z = 0;
	ZZ = zeros(n,1);
	EE = zeros(n,1);
	for idx=1:n
		% get energy of current state
		% hydrogen atom
		E = -mr * e^4 / (8*(eps0*h*idx)^2);
		% potential well
		L = 1.2e-9; % 12 Angstroem
%		E = 0.5*h^2/m_star*(0.5*idx/L)^2;
		% canonical partition function
		Z_ = exp( beta*E); % why no - ?
		Z = Z + Z_;
		% store values
		EE(idx,1) = E;
		ZZ(idx,1) = Z;
		Z__(idx,1) = Z_;
	end
end % function partition function

