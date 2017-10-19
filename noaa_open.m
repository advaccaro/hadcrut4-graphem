%noaa_open.m
odir = '/home/scec-02/avaccaro/HadCRUT4.3/NOAA_GMT/data/';


%url = https://www.ncdc.noaa.gov/cag/time-series/global/globe/land_ocean/p12/12/1880-2015.csv

url = 'https://www.ncdc.noaa.gov/cag/time-series/global/globe/land_ocean/p12/12/1880-2015.csv';
localname = 'noaa_gmt.csv';
localpath = [odir localname];
urlwrite(url,localpath); %write data from url to local .csv



M = csvread('noaa_gmt.csv',3,0);
time = M(:,1); nt = length(time);
anom = M(:,2);
temp = nan(nt,1); tfrac = nan(nt,1); mon = nan(nt,1); tmdn = nan(nt,1);

%mean temperature estimates (1901-2000; Jan->Dec)
%url here: https://www.ncdc.noaa.gov/monitoring-references/faq/anomalies.php
climatology = [12,12.1,12.7,13.7,14.8,15.5,15.8,15.6,15.0,14,12.9,12.2];


for i = 1:nt
	mon(i) = mod(time(i),100); %extract month from time vector
	temp(i) = anom(i) + climatology(mon(i)); %add appropriate monthly value from climatology
	year(i) = (time(i) - mon(i))/100; %extract year
	tfrac(i) = year(i) + mon(i)/12 - 1/12; %rewrite time axis in fractional units
	tmdn(i) = datenum(year(i),mon(i),15);
end






%Assign to structure and save
NOAA.temp = temp;
NOAA.time = time;
NOAA.tfrac = tfrac;
NOAA.tmdn = tmdn;
NOAA.anom = anom;
NOAA.climatology = climatology;
NOAA.mon = mon;
matname = 'NOAA_GMT.mat';
matpath = [odir matname];

save(matpath, 'NOAA');
