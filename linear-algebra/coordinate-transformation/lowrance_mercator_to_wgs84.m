% Thu Feb  6 16:50:59 WIB 2014
% Karl Kastner, berlin
%
%% convert lowrance coordinates to wgs84
%%
%% based on spreadsheet by D Whitney King and Patty B at Lowrance
% TODO, this should better be part of the lowrance class
function [ lat lon ] = lowrance_mercator_to_wgs84(x,y)
	%Function MerToGeoLong(xpos As Double)
	%'Merctor to Geographic in Degrees
	%'Input x position in meters (or whatever the standard unit), Output GeoLat in deg
	%'algrithm from Patty B at lowrance. WGS84 only
	%Dim RadToDeg As Double, DegToRad As Double, b As Double, PI As Double, HALF_PI As Double
	%RadToDeg = 57.2957795132
	%DegToRad = 0.0174532925199
	%PI = 3.141592654
	%HALF_PI = 1.570796327
	%MerToGeoLong = xpos * RadToDeg / b

	b = 6356752.3142;
	lat = rad2deg( 2 * atan( exp(y / b) ) - 0.5*pi );
	lon = rad2deg(x)/b;
end % function lowrance_mercator_to_wgs84()

