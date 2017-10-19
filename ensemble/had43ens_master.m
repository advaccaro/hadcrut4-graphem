%had44ens_master.m
addpath(genpath('/home/scec-02/avaccaro/HadCRUT4.3/ensemble/'))

%% Initialize (define vars, preallocate, etc.)

%Load raw data (median)
raw = load('/home/scec-02/avaccaro/HadCRUT4.3/data/had43med.mat');
[nlon,nlat,nmon] = size(raw.H43med.x);
%nm = 1973; %nmons defined explicitly! (ONLY FOR INFILLED ENS MEMBERS, NOT MEDIAN)
lon = raw.H43med.lon;
lat = raw.H43med.lat;
[nt,ns] = size(raw.had43med); %nmonths x nstations
ny = floor(nt/12); %nyears
years = [1850:2014];
nm = ny*12; %number of months

%[x,y] = meshgrid(raw.H43med.lon,raw.H43med.lat);
%loc = [x(:),y(:)];
loc = raw.loc;
time = raw.H43med.tser;
%time = zeros(1973,1);
%time(1:1972) = 675700 + raw.H43med.time;
%time(end) = time(end-1)+30.5;

%weights = repmat(cosd(loc(:,2)), [1 ny]); %weights by lat for global mean
weights = repmat(cosd(loc(:,2)), [1 nmon]);


%% Unpack .mat files/Assemble 100-member ensemble
datadir = '/home/scec-02/avaccaro/HadCRUT4.3/ensemble/graphem_sp/data/';
datalist = dir([datadir '/*sp80_step2.mat']); %list of .mat files
ndata = length(datalist);

%ENS.globalmean = nan(ny,ndata);
%ENS.globalmean = nan(nmon,ndata);
%ENS.nino34 = nan(nmon,ndata);

ENS.globalmean = nan(nm,ndata);
ENS.nino34 = nan(nm,ndata);

for i = 1:ndata
filename = datalist(i).name;
E = had43ens_explore(filename); %function: had43ens_explore computes GMMT and NINO3.4 for each file
ENS.globalmean(:,i) = E.GMMT;
ENS.nino34(:,i) = E.NINO34m;
clear filename E
end

%E = load(filename); %Load ensemble member

%%% Global Annual mean

% %Annualize
% had43ens_ann = nan(ny,ns);
% for k = 1:ns
% for j = 1:ny
% had43ens_ann(j,k) = mean(E.Xf_sp(1+(j-1)*12:12*j,k));
% end
% end

%%Create global monthly mean time series
%temp = nan(nmon,1);
%for j = 1:nmon
%temp(j) = (E.Xf_sp(j,:)*weights(:,j))/sum(weights(:,j));
%end


% %Create global annual mean time series
% temp = nan(ny,1);
% for j = 1:ny
% temp(j) = (had43ens_ann(j,:)*weights(:,j))/sum(weights(:,j));
% end

%ENS.globalmean(:,i) = temp;

%% Nino3.4
%[nm,~] = size(E.Xf_sp);
%Xc = center(E.Xf_sp);
%Xc_3d = reshape(Xc', [nlon,nlat,nm]);
%[nino12,nino3,nino4,nino34] = nino_indices(Xc_3d, lon, lat);
%ENS.nino34(:,i) = nino34(1:end-1); %WARNING: ***EXTREMELY AD-HOC!!!***

%clear E filename temp had4ens_ann Xc Xc_3d
%end

save ENS.mat ENS

%% Calculate quantiles, median, mean
% %Global annual mean
% for i = 1:ny
% globalmean_sum(i,:) = quantile(ENS.globalmean(i,:),[.05 .25 .5 .75 .95]);
% end

%global monthly mean temperature
for i = 1:nm
globalmean_sum(i,:) = quantile(ENS.globalmean(i,:),[.05 .25 .5 .75 .95]);
end


%Nino34
for i = 1:nm
nino34_sum(i,:) = quantile(ENS.nino34(i,:),[.05 .25 .5 .75 .95]);
end



%% Plotting
%fig(1)
%hold on;
%plot(years,globalmean_sum)
%fancyplot_deco('12-Month Global Mean Temperature over Land & Ocean','Year','deg C')
%hepta_figprint('had4ens_globalmean')


%fig(2)
%hold on;
%plot(time,nino34_sum)
%datetick('x','yyyy')
%hleg = legend('')
%set(hleg, 'Location', 'SouthEast')
%legend('boxoff')
%fancyplot_deco('NINO3.4 SST Anomaly', 'Year', 'deg C')
%hepta_figprint('had4ens_nino34')




