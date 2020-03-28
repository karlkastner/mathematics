% 2020-01-02 17:10:08.033102880 +0800

function [Z_egm96,dz] = ellipsoid2geoid(lat,lon,Z_wgs84)
	dz = cvec(geoidheight(lat, lon, 'EGM96'));
	Z_egm96 = Z_wgs84-dz;
end

