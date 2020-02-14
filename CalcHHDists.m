function [dHH,ia,ib]=CalcHHDists(x)

% Set radius of earth in m
r=6371.0088*1000;

% Find unique households (HHs)
[~,ia,ib]=unique(x.HHID);

% Convert HH latitudes and longitudes from degrees to radians
lat=x.latitude(ia)/180*pi;
lon=x.longitude(ia)/180*pi;

lat_diff=abs(bsxfun(@minus,lat,lat'));
lon_diff=abs(bsxfun(@minus,lon,lon'));

dHH=r*2*asin(sqrt((sin(lat_diff/2)).^2+cos(lat)*cos(lat').*(sin(lon_diff/2)).^2));

% d=dHH(ib,ib);