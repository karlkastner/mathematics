% 2013-07-22 09:37:48 UTC
% Karl Kästner, Berlin
%
%% convert wgs84 utm to latitute and longitude
% Kapuas is zone 49M
%
function [lat lon zone] = utm2latlon(east, north, zone)
	if(0 == length(east))
		lat  = [];
		lon  = [];
		zone = [];
	else
	mUtm=defaultm('utm');
	if (nargin() < 3 || isempty(zone))
		zone  = utmzone(east,north);
	end
	mUtm.zone  = zone;
	mUtm.geoid = almanac('earth','wgs84','meters');
	%mUtm.geoid = wgs84Ellipsoid('meters');
	mUtm       = defaultm(mUtm);
%        [eastutm2, northutm2] = mfwdtran(mUtm,44.3536,28.4981);
%default values taken from manual lat -> UTM conversion
%diffnorth2=4912239-northutm2
%diffeast2=619393.56-eastutm2
	% Perform inverse mapping to verify lat and lon
	%[lat, lon] = minvtran(mUtm, eastutm2,northutm2)
	[lat lon] = minvtran(mUtm, east, north);
	[a b] = minvtran(mUtm, lon, lat)
	end
end % utm2latlon

