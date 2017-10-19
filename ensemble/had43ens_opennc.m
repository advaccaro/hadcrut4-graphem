%had43ens_opennc.m


%filenum = ; %filenum must be defined in pbs script!
addpath(genpath('/home/scec-02/avaccaro/HadCRUT4.3/'))
addpath('/home/scec-02/jianghaw/pseudoproxy/graphem_test/graphem/')


%%Initialize (0th step) (open raw data netCDf, update  w/ cw2013, format)
load('cw2013.mat')
indir = '/home/scec-02/avaccaro/HadCRUT4.3/ensemble/data/raw/';
infile = ['HadCRUT.4.3.0.0.anomalies.' num2str(filenum) '.nc'];
inpath = [indir infile];


%unpack netCDF file
ncidin = netcdf.open(inpath); %open netCDF file

%netCDF var names
latname = 'latitude'; 
lonname = 'longitude';
timename = 'time'; 
tname = 'temperature_anomaly';
statusname = 'field_status';

%find correct indices for each variable
time_indx = netcdf.inqVarID(ncidin,timename);
lon_indx = netcdf.inqVarID(ncidin,lonname);
lat_indx = netcdf.inqVarID(ncidin,latname);
temp_indx = netcdf.inqVarID(ncidin,tname);
status_indx = netcdf.inqVarID(ncidin,statusname);

%now import
H43ens.lat = netcdf.getVar(ncidin,lat_indx);
H43ens.lon = netcdf.getVar(ncidin,lon_indx);
H43ens.time = netcdf.getVar(ncidin,time_indx);
H43ens.x = netcdf.getVar(ncidin,temp_indx);
H43ens.stat = netcdf.getVar(ncidin, status_indx);

%encode NaN's
H43ens.x(H43ens.x < -99999) = NaN;

%Reshape for GraphEM
[nlon, nlat, ntime] = size(H43ens.x);
nloc = nlat * nlon;
had43ens = reshape(H43ens.x, [nloc, ntime]);
had43ens = had43ens';
had43ens = double(had43ens);

rawH43ens = had43ens;

[x,y] = meshgrid(H43ens.lat,H43ens.lon);
loc = [y(:),x(:)];

navail = sum(~isnan(had43ens));
station = find(navail >= ntime/5);
nstation = find(navail < ntime/5);

%load cw2013 data and update had4ens
%load('cw2013.mat') %interpolated data
cw = reshape(cw13.temp, [nloc,1973]);
cw = cw';
cw = double(cw);

%put on common time axis
starttime = datenum('Jan-1-1850');

cwtime = starttime + double(cw13.time);
cwtvec = datevec(cwtime);
cwtfrac = cwtvec(:,1) + cwtvec(:,2)/12 - 1/12;
cwendtime = cwtfrac(end);

htime = starttime + double(H43ens.time);
htvec = datevec(htime);
htfrac = htvec(:,1) + htvec(:,2)/12 - 1/12;
H43ens.tser = htime;
H43ens.tvec = htvec;
H43ens.tfrac = htfrac;

had43ens((htfrac >= 1979 & htfrac <= cwendtime), nstation) = cw(cwtfrac >= 1979, nstation);
navail2 = sum(~isnan(had43ens));
station2 = find(navail>=ntime/5);

odir = '/home/scec-02/avaccaro/HadCRUT4.3/ensemble/data/';
ofile = [infile(1:end-3) '.mat'];
opath = [odir ofile];

%save raw and updated hadcrut4 ensemble member as .mat
save(opath, 'H43ens' , 'had43ens', 'loc', 'rawH43ens')


