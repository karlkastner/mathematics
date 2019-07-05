% 2013-08-29 19:38:43.51 +0700
% Karl Kastner, Berlin
%
%% transform latitude and longitude to WGS84 UTM
%
% function [easting northing zone] = latlon2utm(lat, lon, zone)
%
function [easting, northing, zone] = latlon2utm(lat, lon, zone)
	if(0 == length(lat))
		easting = [];
		northing =[];
		zone = [];
	else
		mUtm=defaultm('utm');
		if (nargin() < 3 || isempty(zone))
			% I missing, therefore no numeric conversion
			z = ['A','B','C','D','E','F','G','H','J','K','L',...
			     'M','N','P','Q','R','S','T','U','V','W','X'];
			id = find(isfinite(lon) & isfinite(lat),1,'first');
			zone_lon = floor((lon(id)+186)/6);
			zone_lat = ceil((lat(id)+98) / 8);
			zone = [num2str(zone_lon), z(zone_lat)];
		end
		mUtm.zone = zone;
		%mUtm.geoid = almanac('earth','wgs84','meters');
		mUtm.geoid = wgs84Ellipsoid('meters');
		mUtm       = defaultm(mUtm);
		[easting, northing] = mfwdtran(mUtm, lat, lon);
%		zone = mUtm.zone;
		%default values taken from manual lat -> UTM conversion
	
		% Perform inverse mapping to verify lat and lon
	%	[lat, lon] = minvtran(mUtm, easting, northing)
	end
end % latlon2utm


