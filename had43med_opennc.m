%had43med_opennc.m
%opens hadcrut4 ensemble median .nc file and saves as .mat

tic;
%addpath(genpath('/Users/adam/Desktop/HadCRUT4.2.0.0'))
addpath(genpath('/home/scec-02/avaccaro/HadCRUT4.3'))

%addpath(genpath('/Users/adam/Desktop/Treerings/matlib'))
%addpath(genpath('/home/geovault-02/avaccaro/matlib'))
%addpath('/home/scec-02/julieneg/matlib/RegEM_KCV')
%addpath(genpath('/home/scec-02/julieneg/matlib/'))
%indir = '/Users/adam/Desktop/HadCRUT4.2.0.0/data/rawdata/';
indir = '/home/scec-02/avaccaro/HadCRUT4.3/data/raw/';

%odir = '/Users/adam/Desktop/HadCRUT4.2.0.0/data/';


load('cw2013.mat') %satellite data

latname = 'latitude';
lonname = 'longitude';
timename = 'time';
tname = 'temperature_anomaly';
statusname = 'field_status';

filename = 'HadCRUT.4.3.0.0.median.nc';

%open file
ncidin = netcdf.open(filename);

%find correct indices for each variable
time_indx = netcdf.inqVarID(ncidin,timename);
lon_indx = netcdf.inqVarID(ncidin,lonname);
lat_indx = netcdf.inqVarID(ncidin,latname);
temp_indx = netcdf.inqVarID(ncidin,tname);
status_indx = netcdf.inqVarID(ncidin,statusname);

%now import
H43med.lat = netcdf.getVar(ncidin,lat_indx);
H43med.lon = netcdf.getVar(ncidin,lon_indx);
H43med.time = netcdf.getVar(ncidin,time_indx);
H43med.x = netcdf.getVar(ncidin,temp_indx);
H43med.stat = netcdf.getVar(ncidin, status_indx);

%encode NaN's
H43med.x(H43med.x < -99999) = NaN;

[nlon, nlat, ntime] = size(H43med.x);
nloc = nlat * nlon;

%reshape for RegEM
had43med = reshape(H43med.x, [nloc, ntime]);
had43med = had43med'; %time x space
had43med = double(had43med);

rawH = had43med;


[x,y] = meshgrid(H43med.lat,H43med.lon);
loc = [y(:),x(:)];

navail = sum(~isnan(had43med));
station = find(navail >= ntime/5);
nstation = find(navail < ntime/5);









cw = reshape(cw13.temp, [nloc,1973]);
cw = cw';
cw = double(cw);

%put on common time axis
starttime = datenum('Jan-1-1850');

cwtime = starttime + double(cw13.time);
cwtvec = datevec(cwtime);
cwtfrac = cwtvec(:,1) + cwtvec(:,2)/12 - 1/12;
cwendtime = cwtfrac(end);

htime = starttime + double(H43med.time);
htvec = datevec(htime);
htfrac = htvec(:,1) + htvec(:,2)/12 - 1/12;
H43med.tser = htime;
H43med.tvec = htvec;
H43med.tfrac = htfrac;

had43med((htfrac >= 1979 & htfrac <= cwendtime), nstation) = cw(cwtfrac >= 1979, nstation);
navail2 = sum(~isnan(had43med));
station2 = find(navail>=ntime/5);

odir = '/home/scec-02/avaccaro/HadCRUT4.3/data/';
save43tag = [odir 'had43med.mat'];
%save as .mat
save(save43tag, 'H43med' , 'had43med', 'loc', 'rawH')

toc;
