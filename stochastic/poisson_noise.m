% Thu 19 Oct 09:13:10 CEST 2023
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
%% triggers events when the time step counter tdx reaches the value of next
%% sets the time of the next event according to the poisson random number with rate "rate"
%% there is the possibility that the random number is zero, in this case
%% several events happen simultaneously at the current step,
%% the counter "event" returns the number of events happening at the current time step
function [event,next] = poisson_noise(tdx,next,rate)
	event = (tdx >= next);
	fdx = find(event);
	event = double(event);
	% this iteration takes ~ log(n/rate) steps
	while (~isempty(fdx))
		%next(fdx) = next(fdx) + poissrnd(rate(fdx));
		next(fdx) = next(fdx) + exprnd(rate(fdx));
		fdx = fdx(tdx >= next(fdx));
		event(fdx) = event(fdx)+1;
	end
end

