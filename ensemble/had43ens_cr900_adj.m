%had43ens_cr900_adj.m

addpath(genpath('/home/scec-02/avaccaro/HadCRUT4.3/'))

load('had43med.mat')

adjR = neigh_radius_adj(loc,900);

odir = '/home/scec-02/avaccaro/HadCRUT4.3/ensemble/graphem_cr/data/';
adjtag = 'had43med_graphem_cr900_adj.mat';

adjpath = [odir adjtag];
save(adjpath, 'adjR')


