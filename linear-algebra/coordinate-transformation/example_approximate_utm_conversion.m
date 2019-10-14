% Thu 13 Sep 11:27:10 CEST 2018

function test

lat  = [0,0]
lon  = [109,110];
%l    = lon/dlon+1/2*nzone;
%x    = ((frac(l)-1/2)*c/nzone + x0)
[x,zone] = latlon2utm_simple(lat,lon) 
[x_,y,z]=latlon2utm(lat,lon)


% berlin
lat = 52.515891
lon = 13.408963
[x,zone] = latlon2utm_simple(lat,lon) 
[x_,y,z]=latlon2utm(lat,lon)

l=1:4:80; [x,y]=latlon2utm(l,zeros(1,length(l))); [x_,y_]=latlon2utm_simple(l,zeros(1,length(l))); [l;y/1e6;y_/1e6]

% two points of same latitude have similar y
% two points of same longitude have not similar x, excep they are in the centre of the zone(!)
end
