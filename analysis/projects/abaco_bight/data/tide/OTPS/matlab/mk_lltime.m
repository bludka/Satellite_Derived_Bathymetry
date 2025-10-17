%1) Lat: -52.664130°  Lon: -70.237260° (ADCP location)
%2) Lat: -52.690000°  Lon: -70.275000° (VMP, M5-2)
%3) Lat: -52.800000°  Lon: -70.610000°  (VMP, M6-2)
lat=[-52.664130 -52.690000 -52.800000];
lon=[-70.237260 -70.275000 -70.610000];
dt=1/24;
np=length(lat);
%
d1=datenum(2019,03,02);
d2=datenum(2019,03,08);
for k=1:np
fid=fopen(['lat_lon_time_IL_' int2str(k)],'w');
for t=d1:dt:d2
fprintf(fid,'%10.4f %10.4f %8d %4d %4d %4d %4d %4d\n',lat(k),lon(k),datevec(t));
end
fclose(fid);
end
