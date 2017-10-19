%had43med_graphem_sp80_explore.m

addpath(genpath('/home/scec-02/avaccaro/HadCRUT4.3/'))

load('had43med.mat')

X = load('had43med_graphem_sp80_step2.mat')

[H43M80_EXP] = had43med_explore(X.Xf,H43med.lon,H43med.lat,X.tser,X.loc)


odir = '/home/scec-02/avaccaro/HadCRUT4.3/data/';
outname = 'H43med_SP80_EXP.mat';
outpath = [odir outname];
save(outpath,'H43M80_EXP');
clear all
%might as well do the others too!
%raw
load('had43med.mat')
odir = '/home/scec-02/avaccaro/HadCRUT4.3/data/';

[H43MRAW_EXP] = had43med_explore(rawH,H43med.lon,H43med.lat,H43med.tser,loc);
rawname = 'H43med_RAW_EXP.mat';
rawpath = [odir rawname];
save(rawpath,'H43MRAW_EXP')

clear all
%cowtan and way updated
load('had43med.mat')
odir = '/home/scec-02/avaccaro/HadCRUT4.3/data/';
[H43MUPD_EXP] = had43med_explore(had43med,H43med.lon,H43med.lat,H43med.tser,loc);
updname = ['H43med_UPD_EXP.mat'];
updpath = [odir updname];
save(updpath, 'H43MUPD_EXP');


%Cowtan and Way (kriging)
%load('cw2013.mat')
load('/home/scec-02/avaccaro/HadCRUT4.3fv/CW2013/data/cw2013.mat');
odir = '/home/scec-02/avaccaro/HadCRUT4.3/cw2013/data/';
cw2013_3d = double(cw13.temp); [nlon,nlat,nmon] = size(cw2013_3d);
cw2013_2d = reshape(cw2013_3d, [nlon*nlat,nmon])';
cwlon = double(cw13.lon); cwlat = double(cw13.lat);
cwtser = datenum('Jan-0-1850') + double(cw13.time);
[CW_EXP] = had43med_explore(cw2013_2d,cwlon,cwlat,cwtser,loc);
cwname = ['CW_EXP.mat'];
cwpath = [odir cwname];
save(cwpath, 'CW_EXP')




