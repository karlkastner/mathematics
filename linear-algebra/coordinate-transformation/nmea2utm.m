% Thu Jul 10 08:27:14 WIB 2014
% Karl Kastner, Berlin
%
%% convert nmea messages to utm coordinates
function [wgs84 utm] = nmea2utm(NMEAGGA, zone)
	% case insensitivity
	if (isfield(NMEAGGA,'Lat'))
		lat = NMEAGGA.Lat;
	else
		lat = NMEAGGA.lat;
	end
	% idiot proof, e.g. do not assume lat and long have same capitalisation
	if (isfield(NMEAGGA,'Long'))
		lon = NMEAGGA.Long;
	else
		lon = NMEAGGA.long;
	end

	lat         = double(lat);
	lon	    = double(lon);

	% set minus sign for southern latitudes
	if (isfield(NMEAGGA,'SN'))
		fdx = (NMEAGGA.SN == 'S');
	else
		% if no hemsiphere is given assume signed latitude
		% this assumption fails when the ship track crosses the eaquator
		% and the latitude is unsigned
		warning('Hemisphere information missing, unique lat-long conversion not feasible, assuming northern Hemispheror signed lat');
		fdx = false(size(lat));
	end
	lat(fdx)       = -lat(fdx);

	fdx = find(0 == lat);
	lat(fdx) = NaN;
	gdx = find(0 == lon);
	lon(fdx) = NaN;
	if (length(fdx) || length(gdx))
		warning('nmea2utm:zero positions ignored');
	end

	% take the median of the gps coordinates measured during the ADCP ping
	if (~isvector(lat))
		lat      = nanmedian(lat,2);
		lon      = nanmedian(lon,2);
	else
		% make sure this is a column vector
		lat = cvec(lat);
		lon = cvec(lon);
	end
	wgs84.lat = lat;
	wgs84.lon = lon;
	[utm.X, utm.Y] = latlon2utm(lat, lon, zone);
end % nmea2utm

