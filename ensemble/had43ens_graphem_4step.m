%had43ens_graphem_4step.m


%filenum = ; %filenum must be defined in pbs script!
addpath(genpath('/home/scec-02/avaccaro/HadCRUT4.3/'))
addpath('/home/scec-02/jianghaw/pseudoproxy/graphem_test/graphem/')


%%Initialize open .mat from (0th step)
indir = '/home/scec-02/avaccaro/HadCRUT4.3/ensemble/data/';
infile = ['HadCRUT.4.3.0.0.anomalies.' num2str(filenum) '.mat'];
inpath = [indir infile];
load(inpath)





%% GraphEM CR Part 1 (1st step)
target_cr = 900;
odir_cr = '/home/scec-02/avaccaro/HadCRUT4.3/ensemble/graphem_cr/data/';

cronename = ['had43ens_' num2str(filenum) '_cr' num2str(target_cr) '_step1.mat'];
cronepath = [odir_cr cronename];
crtwoname = ['had43ens_' num2str(filenum) '_cr' num2str(target_cr) '_step2.mat'];
crtwopath = [odir_cr crtwoname];



X = had43ens;

% Coordinates
np = length(loc);
lonlat = loc;
lon = lonlat(:,1);
lat = lonlat(:,2);

%time axes (split into recent/past)
tfrac = H43ens.tfrac;
trec = tfrac(tfrac >= 1925);
thist = tfrac(tfrac < 1925);
indrec = ismember(tfrac,trec);
indhist = ismember(tfrac,thist);
Xrec = X(indrec,:);
Xhist = X(indhist,:);

%GraphEM options
opt.regress = 'ols';
opt.stagtol = 5e-3;
opt.maxit = 40;
opt.lonlat = lonlat;
opt.useggm = 1;
ADJ = load('had43med_graphem_cr900_adj.mat');
opt.adj = ADJ.adjR;
if ~ismember(filenum,[1:8,10:20,23:25,27:29,31:39,41:42,44,46:49,51,53:62,64:72,74:75,76:77,80:82,84,86,100]) == 1
	[Xout, Mout, Cout] = graphem_JW(Xrec,opt);
	Xcon_cr = vertcat(Xhist,Xout);

	%save graphem_cr_step1
	save(cronepath, 'Xcon_cr','Xout', 'Mout', 'Cout', 'loc', 'H43ens')
else
	load(cronepath)
end

clear opt
%% GraphEM CR Part 2 (2nd step)
%GraphEM options
opt.regress = 'ols';
opt.stagtol = 5e-3;
opt.maxit = 40;
opt.lonlat = lonlat;
opt.useggm = 1;
opt.adj = ADJ.adjR;
opt.C0 = Cout;
if ~ismember(filenum,[1:8,10:20,23:25,27:29,31:39,41:42,44,46:49,51,53:62,64:70,72,74:75,76:77,80:82,84,86,100]) == 1
	[Xf_cr,Mf_cr,Cf_cr] = graphem_JW(Xcon_cr,opt);

	save(crtwopath, 'Xf_cr', 'Mf_cr', 'Cf_cr', 'loc', 'H43ens')
else
	load(crtwopath)
end


%save graphem_cr_step2

%% GraphEM SP Part 1 (3rd step)
target_spars = .8;
odir_sp = '/home/scec-02/avaccaro/HadCRUT4.3/ensemble/graphem_sp/data/';

spadjtag = [odir_sp 'had43ens_' num2str(filenum) '_sp' num2str(target_spars*100) '_adj.mat'];
sponetag = [odir_sp 'had43ens_' num2str(filenum) '_sp' num2str(target_spars*100) '_step1.mat'];
sptwotag = [odir_sp 'had43ens_' num2str(filenum) '_sp' num2str(target_spars*100) '_step2.mat'];

tcal = tfrac(tfrac >= 1960 & tfrac < 1991);
calib = ismember(tfrac,tcal);
Xcal = Xf_cr(calib,:);
S = corrcoef(Xcal); N = 150;
[spars,adj] = greedy_search_TT(S,target_spars, N);
spars_f = spars(end); adjM = adj{end};
save(spadjtag,'spars_f','adjM','target_spars')
clear adj spars

%GraphEM Options
opt.regress = 'ols';
opt.stagtol = 5e-3;
opt.maxit = 40;
opt.useggm = 1;
opt.err_export = 0;
opt.adj = adjM;

if ~ismember(filenum,[1:8,10:20,23:25,27:29,31:39,41:42,44,46:49,51,53:62,64:69,76:77,84,86,100]) == 1
	[Xo,Mo,Co] = graphem_JW(Xrec,opt);
	Xcon_sp = vertcat(Xhist,Xo);
	save(sponetag, 'Xo','Mo','Co','target_spars','loc','H43ens', 'Xrec', 'Xhist', 'Xcon_sp', 'spars_f')
else
	load(sponetag)
end
%save graphem_sp_step1

%% GraphEM SP Part 2 (4th step)
opt.regress = 'ols';
opt.stagtol = 5e-3;
opt.maxit = 40;
opt.useggm = 1;
opt.err_export = 0;
opt.adj = adjM;
opt.C0 = Co;
opt.lonlat = lonlat;

if ~ismember(filenum,[1:8,10:20,23:25,27:29,31:39,41:42,44,46:49,51,53:62,64:69,76:77,84,86,100]) == 1
	[Xf_sp,Mf_sp,Cf_sp] = graphem_JW(Xcon_sp,opt);


	%save graphem_sp_step2
	save(sptwotag, 'Xf_sp', 'Mf_sp', 'Cf_sp', 'target_spars', 'loc', 'H43ens', 'spars_f')
else
	display('already done!')
end




