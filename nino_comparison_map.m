%nino_comparison_map.m
%m_proj, m_grid, m_coast, m_line, m_pcolor

%% load raw HadCRUT4
R = load('/home/scec-02/avaccaro/HadCRUT4.3/data/had43med.mat')

%% load GraphEM-infilled HadCRUT4 (GLASSO)
H = load('/home/scec-02/avaccaro/HadCRUT4.3/graphem_sp/data/had43med_graphem_sp80_step2.mat')

%% load GraphEM-infilled HadCRUT4 (NEIGH)
N = load('/home/scec-02/avaccaro/HadCRUT4.3/graphem_cr/data/had43med_graphem_cr900_step2.mat')

%% load Cowtan and Way gridded dataset
load('/home/scec-02/avaccaro/HadCRUT4.3fv/CW2013/data/cw2013')
lat = double(cw13.lat); lon = double(cw13.lon); %C.lon(C.lon < 0) = C.lon(C.lon < 0) + 360;
[nlon,nlat,ntime] = size(cw13.temp);
C.tfrac = nan(ntime,1);
for i = 1:ntime
	C.tfrac(i) = 1850 + (i-1)/12;
end
C.grid = double(cw13.temp);


%% calculate 3 month means for 1877-1878 (DJF, MAM, JJA)
%DJF
cDJFt = C.tfrac(C.tfrac >= 1878-1/12 & C.tfrac <= 1878+1/12);
cDJFind = ismember(C.tfrac,cDJFt);
C.DJF = mean(C.grid(:,:,cDJFind),3);
hDJFt = H.tfrac(H.tfrac >= 1878-1/12 & H.tfrac <= 1878+1/12);
hDJFind = ismember(H.tfrac,hDJFt);
H.DJF = mean(H.Xf(hDJFind,:));
H.DJF = reshape(H.DJF, [nlon nlat]);



%MAM
cMAMt = C.tfrac(C.tfrac >= 1878+2/12 & C.tfrac <= 1878+4/12);
cMAMind = ismember(C.tfrac,cMAMt);
C.MAM = mean(C.grid(:,:,cMAMind),3);
hMAMt = H.tfrac(H.tfrac >= 1878+2/12 & H.tfrac <= 1878+4/12);
hMAMind = ismember(H.tfrac,hMAMt);
H.MAM = mean(H.Xf(hMAMind,:));
H.MAM = reshape(H.MAM, [nlon nlat]);

%JJA
cJJAt = C.tfrac(C.tfrac >= 1878+5/12 & C.tfrac <= 1878+7/12);
cJJAind = ismember(C.tfrac,cJJAt);
C.JJA = mean(C.grid(:,:,cJJAind),3);
hJJAt = H.tfrac(H.tfrac >= 1878+5/12 & H.tfrac <= 1878+7/12);
hJJAind = ismember(H.tfrac,hJJAt);
H.JJA = mean(H.Xf(hJJAind,:));
H.JJA = reshape(H.JJA, [nlon nlat]);





%% Plotting
fig('nino comparison maps')
hold on;
subplot(3,2,1); hold on;
m_proj('Robinson','clong',0);

%m_grid('xtick',[0:60:360],'tickdir','out','ytick',[-90:30:90], 'color','k', 'fontsize',10,'fontname','Times New Roman');
m_grid('box','off','xtick',9,'ytick',6,'xlabeldir','end','xticklabels',[],'yticklabels',[]);
m_coast('color','k');
h=m_pcolor(lon,lat,C.DJF'); caxis([-3 3]);
set(h,'EdgeAlpha',0)
title('CW - DJF')

subplot(3,2,2); hold on;
m_proj('Robinson','clong',0);
m_grid('box','off','xtick',9,'ytick',6,'xlabeldir','end','xticklabels',[],'yticklabels',[]);
m_coast('color','k');
h=m_pcolor(lon,lat,H.DJF'); caxis([-3 3]);
set(h,'EdgeAlpha',0)
title('GraphEM - DJF')

subplot(3,2,3); hold on;
m_proj('Robinson','clong',0);
m_grid('box','off','xtick',9,'ytick',6,'xlabeldir','end','xticklabels',[],'yticklabels',[]);
m_coast('color','k');
h=m_pcolor(lon,lat,C.MAM'); caxis([-3 3]);
set(h,'EdgeAlpha',0)
title('CW - MAM')

subplot(3,2,4); hold on;
m_proj('Robinson','clong',0);
m_grid('box','off','xtick',9,'ytick',6,'xlabeldir','end','xticklabels',[],'yticklabels',[]);
m_coast('color','k');
h=m_pcolor(lon,lat,H.MAM'); caxis([-3 3]);
set(h,'EdgeAlpha',0)
title('GraphEM - MAM')

subplot(3,2,5); hold on;
m_proj('Robinson','clong',0);
m_grid('box','off','xtick',9,'ytick',6,'xlabeldir','end','xticklabels',[],'yticklabels',[]);
m_coast('color','k');
h=m_pcolor(lon,lat,C.JJA'); caxis([-3 3]);
set(h,'EdgeAlpha',0)
title('CW - JJA')

subplot(3,2,6); hold on;
m_proj('Robinson','clong',0);
m_grid('box','off','xtick',9,'ytick',6,'xlabeldir','end','xticklabels',[],'yticklabels',[]);
m_coast('color','k');
h=m_pcolor(lon,lat,H.JJA'); caxis([-3 3]);
set(h,'EdgeAlpha',0)
title('GraphEM - JJA')



hepta_figprint('./figs/nino1878comparison.eps')






