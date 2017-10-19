%cw13_ncopen.m
%opens the cowtan and way 2013 temperature series

%addpath(genpath('/Users/adam/Desktop/HadCRUT4.2.0.0'))
addpath(genpath('/home/scec-02/avaccaro/HadCRUT4.2'))

%addpath(genpath('/Users/adam/Desktop/Treerings/matlib'))
addpath(genpath('/home/geovault-02/avaccaro/matlib'))

%indir = '/Users/adam/Desktop/HadCRUT4.2.0.0/data/cw2013/';
indir = '/home/scec-02/avaccaro/HadCRUT4.2/data/cw2013/';

%odir = '/Users/adam/Desktop/HadCRUT4.2.0.0/data/';
odir = '/home/scec-02/avaccaro/HadCRUT4.2/data/';

infile = [indir 'had4_krig_v2_0_0.nc'];

latname = 'latitude';
lonname = 'longitude';
timename = 'time';
yearname = 'year';
monthname = 'month';
tempname = 'temperature_anomaly';


%open file
ncidin = netcdf.open(infile);

%find correct indices for each variable
lat_indx = netcdf.inqVarID(ncidin,latname);
lon_indx = netcdf.inqVarID(ncidin,lonname);
time_indx = netcdf.inqVarID(ncidin,timename);
year_indx = netcdf.inqVarID(ncidin,yearname);
month_indx = netcdf.inqVarID(ncidin,monthname);
temp_indx = netcdf.inqVarID(ncidin,tempname);

%now import
cw13.lat = netcdf.getVar(ncidin,lat_indx);
cw13.lon = netcdf.getVar(ncidin,lon_indx);
cw13.time = netcdf.getVar(ncidin, time_indx);
cw13.year = netcdf.getVar(ncidin, year_indx);
cw13.month = netcdf.getVar(ncidin, month_indx);
cw13.temp = netcdf.getVar(ncidin, temp_indx);

save([odir 'cw2013.mat'] 'cw13')
