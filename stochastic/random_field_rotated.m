% Thu 20 Apr 20:04:29 CEST 2023
% Karl Kästner, Berlin
%
%  This program is free software: you can redistribute it and/or modify
%  it under the terms of the GNU General Public License as published by
%  the Free Software Foundation, either version 3 of the License, or
%  (at your option) any later version.
%
%  This program is distributed in the hope that it will be useful,
%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU General Public License for more details.
%
%  You should have received a copy of the GNU General Public License
%  along with this program.  If not, see <https://www.gnu.org/licenses/>.
%
%
% output : ex : random displacement of x
%          ey : random displamcenet of y
% input :
% n
% L 
% s  : scale of e
% fc : cutoff frequency of lowpass-filter
function [ex,ey,ea] = random_field_rotated(n,L,s,fc,order)
	% relative indices of 4 neighbours
	ix = [-0,-1,0,1];
	iy = [-1,0,1, 0];

	% local angle angle
	ea = 2*pi*rand(n);

	ea = lowpass2d_ideal(ea,L,fc,order);

	% displacements of grid points
	%id = 2:n-1;
	id = 2:n+1;

	ea_ = zeros(n+2);
	ea_(id,id)  = ea;
	ea_(id,1)   = ea(:,end);
	ea_(id,end) = ea(:,1);
	ea_(1,id)   = ea(end,:);
	ea_(end,id) = ea(1,:);
	sea = sin(ea_);
	cea = cos(ea_);
% if(0)
%	ex = zeros(n);
%	ey = zeros(n);
%	for idx=1:4
%		ex = ex + cos(ea_(id-ix(idx),id-iy(idx)))*ix(idx) + sin(ea_(id-ix(idx),id-iy(idx)))*iy(idx);
%		ey = ey - sin(ea_(id-ix(idx),id-iy(idx)))*ix(idx) + cos(ea_(id-ix(idx),id-iy(idx)))*iy(idx);
%	end
%	ex_ = ex;
%	ey_ = ey;
% end

	% TODO use circular cdiff
	ex =    - sea(id,id+1) + sea(id,id-1) - cea(id+1,id) + cea(id-1,id);
	ey =    - cea(id,id+1) + sea(id+1,id) + cea(id,id-1) - sea(id-1,id);

%	rms(ex_(:)-ex(:))
%	rms(ey_(:)-ey(:))

	if (0)
	ex__ = ex;
	ey__ = ey;

	ex = zeros(n+2);
	ey = zeros(n+2);
	% neighbour displacements for local rotation
	for idx=1:4
 		ex(id+ix(idx),id+iy(idx)) = ex(id+ix(idx),id+iy(idx)) + cos(ea)*ix(idx) + sin(ea)*iy(idx);
		ey(id+ix(idx),id+iy(idx)) = ey(id+ix(idx),id+iy(idx)) - sin(ea)*ix(idx) + cos(ea)*iy(idx);
		ex_ = ex_ + cos(ea_(id-ix(idx),id-iy(idx)))*ix(idx) + sin(ea_(id-ix(idx),id-iy(idx)))*iy(idx);
		ey_ = ey_ - sin(ea_(id-ix(idx),id-iy(idx)))*ix(idx) + cos(ea_(id-ix(idx),id-iy(idx)))*iy(idx);
	end

	% boundaries TODO cos y missing
	ex_        = ex(2:end-1,2:end-1);
	ex_(:,1)   = ex_(:,1) + ex(2:end-1,end);
	ex_(:,end) = ex_(:,end) + ex(2:end-1,1);
	ex_(1,:)   = ex_(1,:) + ex(end,2:end-1);
	ex_(end,:) = ex_(end,:) + ex(1,2:end-1);
	ex = ex_;
	ey_        = ey(2:end-1,2:end-1);
	ey_(:,1)   = ey_(:,1) + ey(2:end-1,end);
	ey_(:,end) = ey_(:,end) + ey(2:end-1,1);
	ey_(1,:)   = ey_(1,:) + ey(end,2:end-1);
	ey_(end,:) = ey_(end,:) + ey(1,2:end-1);
	ey = ey_;
	end

	% smooth the displacements
%	nf = 500;
%	for iex=1:4
	if (0)
	ex = lowpass2d_ideal(ex,L,fc,order);
	ey = lowpass2d_ideal(ey,L,fc,order);
	end
%	end
	ex = s*ex;
	ey = s*ey;
end

