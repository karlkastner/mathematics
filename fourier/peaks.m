% Mi 8. Apr 16:36:17 CEST 2015
% Karl Kastner, Berlin
%
%% peaks of the power spectrum of a disctrete fourier transform
%%
%% rule for peaks: there is no higher value left or right of the "peak"
%%                 until the signal drops to p*y_peak, p = 0.5
%%
%% works best, when spectrum has been smoothened
%%
%% input :
%% f : frequency
%% y : absolute value of fourier transform (power spectrum)
%% L : length in space or time of series
%%
%% output :
%%
%% a0 : amplitude
%% s0 : standard deviation (error?) of amplitude
%% w0 : width of peak
%% lambda = wave length (period?)
%% pdx : index of peak
%% f : frequency (if not given as input) 
%
function [f0, a0, s0, w0, lambda, pdx, f] = peaks(f,y,L)
	n = length(y);
	if (isempty(f))
		% TODO, times 2pi
		f = (0:n-1)'/L;
	end

	% TODO, no magic number
	p = 0.5;

%	pdx = [];
%	for idx=1:n
%		% find left threshold
%		fl  = find(y(1:idx-1) < p*y(idx),'last');
%		if (~isempty(fl) && max(y(fl:idx-1)) > y(idx))
%			continue;
%		end
%		% find right treshold
%		fr  = find(y(idx+1:end) < p*y(idx),'first')+idx;
%		% check that no point inbetween exceeds value
%		if (~isempty(fr) && max(y(idx+1:fr)) > y(idx))
%			continue;
%		end
%		pdx(end+1,1) = idx;
%	end

	% beam width
	w   = NaN(n,1);
	id1 = zeros(n,1);
	id2 = zeros(n,1);
	select = false(n,1);
	for idx=2:n-1
		if (y(idx) > y(idx-1) && y(idx) > y(idx+1))
			% do not remove square brackets
			l = max([1,find(y(1:idx-1) < p*y(idx),1,'last')]);
			r = min([n,idx + find(y(idx+1:end) < p*y(idx),1,'first')]);
			% only keep, if there are no higher peaks in-between
			if (max(y(idx+1:r)) <= y(idx) ...
			 && max(y(l:idx-1)) <= y(idx))
				select(idx) = true;
				% left half peak
				id1(idx) = l;
				% right half peak
				id2(idx) = r;
			end
		end % if
	end % for idx

	% peak width
	w = id2-id1-0.5;
	% amplitudesof peaks
	[a0, s0] = amplitude_from_peak(y,w);
	% normalise for data series length
	a0 = a0*1/n;

	id1 = id1(select);
	id2 = id2(select);
	a0 = a0(select);
	s0 = s0(select);
	y0 = y(select);
	f0 = f(select);
	w0 = w(select);

	% sort by magnitude
	% was y0
	[void, pdx] = sort(a0,1,'descend');
	a0       = a0(pdx);
	f0       = f0(pdx);
	s0	 = s0(pdx);
	id1      = id1(pdx);
	id2      = id2(pdx);
	w0       = w0(pdx);
	select   = find(select);
	pdx      = select(pdx);
	
if (0)
	% choose onl above detection threshold and with no large peaks near by
	select = false(n,1);
	block  = false(n,1);
	for idx=1:length(a0)
		if (~block(idx))
			select(idx) = true;
			% block adjacent cells
			block(id1(idx):id2(idx)) = true;
		end
	end
	a0  = a0(select);
	f0  = f0(select);
	s0  = s0(select);
	w0  = w0(select);
	pdx = pdx(select);

%	fdx = s0.^2 < f0.^2;
%	fdx =
%	a0  = a0(fdx);
%	p0  = p0(fdx);
%	s0  = s0(fdx);
%	f0  = f0(fdx);
end
	
	% bias reduced estimate of the wave length
	lambda = 1./f0;

	%lambda = 1./f0 - s0.^2./f0.^3;
end

