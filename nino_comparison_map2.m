%nino_comparison_map2.m
%m_proj, m_grid, m_coast, m_line, m_pcolor

%% load raw HadCRUT4
R = load('/home/scec-02/avaccaro/HadCRUT4.3/data/had43med.mat')
[Rnx,Rny,Rntime] = size(R.H43med.x);
R.tfrac = nan(Rntime,1);
for i = 1:Rntime
	R.tfrac(i) = 1850 + (i-1)/12;
end
R.grid = double(R.H43med.x);

%% load GraphEM-infilled HadCRUT4 (GLASSO)
H = load('/home/scec-02/avaccaro/HadCRUT4.3/graphem_sp/data/had43med_graphem_sp80_step2.mat')

%% load GraphEM-infilled HadCRUT4 (NEIGH)
N = load('/home/scec-02/avaccaro/HadCRUT4.3/graphem_cr/data/had43med_graphem_cr900_step2.mat')
N.tfrac = H.tfrac;

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

rDJFt = R.tfrac(R.tfrac >= 1878-1/12 & R.tfrac <= 1878+1/12);
rDJFind = ismember(R.tfrac,rDJFt);
R.DJF = nmean(R.grid(:,:,rDJFind),3);

hDJFt = H.tfrac(H.tfrac >= 1878-1/12 & H.tfrac <= 1878+1/12);
hDJFind = ismember(H.tfrac,hDJFt);
H.DJF = mean(H.Xf(hDJFind,:));
H.DJF = reshape(H.DJF, [nlon nlat]);

nDJFt = N.tfrac(N.tfrac >= 1878-1/12 & N.tfrac <= 1878+1/12);
nDJFind = ismember(N.tfrac,nDJFt);
N.DJF = mean(N.Xf(nDJFind,:));
N.DJF = reshape(N.DJF,[nlon nlat]);

%MAM
cMAMt = C.tfrac(C.tfrac >= 1878+2/12 & C.tfrac <= 1878+4/12);
cMAMind = ismember(C.tfrac,cMAMt);
C.MAM = mean(C.grid(:,:,cMAMind),3);

rMAMt = R.tfrac(R.tfrac >= 1878+2/12 & R.tfrac <= 1878+4/12);
rMAMind = ismember(R.tfrac,rMAMt);
R.MAM = nmean(R.grid(:,:,rMAMind),3);

hMAMt = H.tfrac(H.tfrac >= 1878+2/12 & H.tfrac <= 1878+4/12);
hMAMind = ismember(H.tfrac,hMAMt);
H.MAM = mean(H.Xf(hMAMind,:));
H.MAM = reshape(H.MAM, [nlon nlat]);

nMAMt = N.tfrac(N.tfrac >= 1878+2/12 & N.tfrac <= 1878+4/12);
nMAMind = ismember(N.tfrac,nMAMt);
N.MAM = mean(N.Xf(nMAMind,:));
N.MAM = reshape(N.MAM, [nlon nlat]);

%JJA
cJJAt = C.tfrac(C.tfrac >= 1878+5/12 & C.tfrac <= 1878+7/12);
cJJAind = ismember(C.tfrac,cJJAt);
C.JJA = mean(C.grid(:,:,cJJAind),3);

rJJAt = R.tfrac(R.tfrac >= 1878+5/12 & R.tfrac <= 18978+7/12);
rJJAind = ismember(R.tfrac, rJJAt);
R.JJA = nmean(R.grid(:,:,rJJAind),3);

hJJAt = H.tfrac(H.tfrac >= 1878+5/12 & H.tfrac <= 1878+7/12);
hJJAind = ismember(H.tfrac,hJJAt);
H.JJA = mean(H.Xf(hJJAind,:));
H.JJA = reshape(H.JJA, [nlon nlat]);

nJJAt = N.tfrac(N.tfrac >= 1878+5/12 & N.tfrac <= 1878+7/12);
nJJAind = ismember(N.tfrac,nJJAt);
N.JJA = mean(N.Xf(nJJAind,:));
N.JJA = reshape(N.JJA, [nlon nlat]);

%Predefine arrays for easy plotting
ssta = cell(12,1);
ssta{1} = R.DJF';
ssta{2} = R.MAM';
ssta{3} = R.JJA';
ssta{4} = C.DJF';
ssta{5} = C.MAM';
ssta{6} = C.JJA';
ssta{7} = N.DJF';
ssta{8} = N.MAM';
ssta{9} = N.JJA';
ssta{10} = H.DJF';
ssta{11} = H.MAM';
ssta{12} = H.JJA';


titles = cell(12,1);
titles{1} = ['Raw - DJF'];
titles{2} = ['Raw - MAM'];
titles{3} = ['Raw - JJA'];
titles{4} = ['CW - DJF'];
titles{5} = ['CW - MAM'];
titles{6} = ['CW - JJA'];
titles{7} = ['GraphEM (Neigh) - DJF'];
titles{8} = ['GraphEM (Neigh) - MAM'];
titles{9} = ['GraphEM (Neigh) - JJA'];
titles{10} = ['GraphEM (GLASSO) - DJF'];
titles{11} = ['GraphEM (GLASSO) - MAM'];
titles{12} = ['GraphEM (GLASSO) - JJA'];








%% Plotting
fig('nino comparison maps'); clf;
%hold on;
%ha = tight_subplot(5,3,[.05 .00001], .05, .01);
ha = tight_subplot(5,3,[.01 .00001], .05, .01);
%ha = subplot_pos(5,3,[.1 .05],.05,.05)
%subplot(3,2,1); hold on;
for ii = 1:12
	axes(ha(ii));
	m_proj('Robinson', 'clong', 0);
	mpc = m_pcolor(lon,lat,ssta{ii}); caxis([-3 3]);
	set(mpc,'EdgeAlpha',0);
	coast = m_coast('color','k');
	m_grid('box','on','xtick',9,'ytick',6,'xlabeldir','end','xticklabels',[],'yticklabels',[]);
	%title(titles{ii});
	

end

%axes(ha(13)); title('DJF');

axes(ha(14)); %title('MAM');
caxis([-3 3]);
c = colorbar3('horiz',['^{\circ}C']);
%cpos2 = [.05 .05 .9 .03];
%cpos2 = [.15 .05 .7 .03]; 
cpos2 = [.15 .16 .7 .03];
c.Position = cpos2;

%axes(ha(15)); title('JJA');


hepta_figprint('./figs/nino1878comparison3.eps')





%mpc = m_pcolor(lon,lat,ssta{1}); caxis([-3 3]);

%cpos = c.Position;
%cpos(1) = cpos(1)/2;
%cpos(3) = cpos(3) * 3;
%cpos(4) = cpos(4) * 8;
%c.Position = cpos;

%set(mpc,'visible','off');





%ColorMap = get(gcf,'Colormap');
%m = size(ColorMap,1);
%CData = get(mpc,'CData');
%cmin = min(CData(:)); cmax = max(CData(:));
%idx = min(m,round((m-1)*(CData-cmin)/(cmax-cmin))+1);








%axes(ha(14)); caxis([-3 3]);
%set(ha(14),'CLim',[-3 3]);
%c = colorbar2('horiz'); caxis([-3 3]);
%set(c,'XTick',[-2 0 2]);



%set(ha,'XTickLabel',[]); set(ha,'YTickLabel',[]);
%set(ha,'Box','off');
%set(ha,'MinorGridAlpha',0)
%set(ha,'xtick',[]); set(ha,'ytick',[]);
%set(ha,'color','none')
%set(ha,'visible','off')
%set(ha,'layer','bottom')












