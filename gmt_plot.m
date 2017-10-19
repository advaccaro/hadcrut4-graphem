%gmt_plot.m


%Load GMT datasets
load_gmt_datasets

%Load ensemble
load('/home/scec-02/avaccaro/HadCRUT4.3/ensemble/data/ENS_SP80.mat')

%Abbreviate dataset names
HR = H43MRAW_EXP;
HU = H43MUPD_EXP;
HG = H43M80_EXP;
G = GISTEMP_GMT;
CW = CW_EXP;
N = NOAA_GMT;
E = ENS;

%Flip datenum time axes (because it's a pain in the ass)
HR.tmdn = flipud(HR.tmdn);
HU.tmdn = flipud(HU.tmdn);
HG.tmdn = flipud(HG.tmdn);
%G.tmdn = flipud(G.tmdn); %already oriented in the proper direction!  Hooray!
CW.tmdn = flipud(CW.tmdn);



%% Calculate differences (to HR)
%HadCRUT family
HU.GMTm_diff = HU.GMTm - HR.GMTm;
HG.GMTm_diff = HG.GMTm - HR.GMTm;

HR.tmin = min(HR.tmdn); HR.tmax = max(HR.tmdn);

%Cowtan and Way
CW.tmin = min(CW.tmdn); CW.tmax = max(CW.tmdn);
CW.diffind1 = find(CW.tmdn >= HR.tmin & CW.tmdn <= HR.tmax);
CW.diffind2 = find(HR.tmdn >= CW.tmin & HR.tmdn <= CW.tmax);
CW.tdiff = intersect(CW.tmdn,HR.tmdn);
CW.GMTm_diff = CW.GMTm(CW.diffind1) - HR.GMTm(CW.diffind2);


% GISTEMP
G.tmin = min(G.tmdn); G.tmax = max(G.tmdn);
G.diffind1 = find(G.tmdn >= HR.tmin & G.tmdn <= HR.tmax);
G.diffind2 = find(HR.tmdn >= G.tmin & HR.tmdn <= G.tmax);
G.tdiff = intersect(G.tmdn,HR.tmdn);
G.GMTm_diff = G.GMTm(G.diffind1) - HR.GMTm(G.diffind2);

%% Get indices of recent (1985-present)
HG.trec_ind = find(HG.tmfrac >= 1985); 
G.trec_ind = find(G.tmfrac >= 1985);
N.trec_ind = find(N.tmfrac >= 1985);


%% Plotting
load JEG_graphics

%GMT Figure 1 (HadCRUT family w/ differences)
fig('HadCRUT family GMT w/ differences'); clf;
subplot(2,1,1);
hold on; %xlim = ([1850,2020]);
plot(HR.tmdn,HR.GMTm,'color','g');
plot(CW.tmdn,CW.GMTm,'color','r');
plot(HG.tmdn,HG.GMTm,'color',skyblue );
datetick('x','yyyy')

[hleg,objh] = legend('HadCRUT4.3 median raw', 'HadCRUT4.3 Cowtan & Way', 'HadCRUT4.3 + GraphEM (0.8% sparsity)');
set(hleg,'Location','SouthEast')
set(objh,'linewidth',2)
legend('boxoff')
ylabel_str = ['Temperature (' char(176) 'C)'];
fancyplot_deco('Monthly Global Mean Temperature', 'Time (Year)', ylabel_str,14,'Helvetica')
%axis([1849 2020 -1.5 1])
subplot(2,1,2);
hold on; %xlim = ([1850,2020]);
plot(CW.tdiff,CW.GMTm_diff,'color','r')
plot(HG.tmdn,HG.GMTm_diff,'color',skyblue)
datetick('x','yyyy')
%xlim = ([1850,2020])
[hleg,objh] = legend('HadCRUT4.3 Cowtan & Way', 'HadCRUT4.3 + GraphEM (0.8% sparsity)');
set(hleg,'Location','SouthEast')
set(objh,'linewidth',2)
legend('boxoff')
fancyplot_deco('Same, difference to raw HadCRUT4.3 median', 'Time (Year)', '{\Delta}T ({\circ}C)',14,'Helvetica')
%axis([1849 2020 -1.5 0.5])
odir = '/home/scec-02/avaccaro/HadCRUT4.3/figs/';
fig1name = 'had43med_gmt_plot1.eps';
fig1path = [odir fig1name];
hepta_figprint(fig1path)


%GMT Figure 2 (GISTEMP w/ recent)
ind1 = HG.trec_ind; ind2 = G.trec_ind; ind3 = N.trec_ind;
fig('The Gs (GraphEM and GISS)'); clf;
subplot(2,1,1);
hold on;
plot(HG.tmdn,HG.GMTm,'color',skyblue); 
plot(G.tmdn,G.GMTm,'color',ornj);
plot(N.tmdn,N.GMTm,'color',maroon);
datetick('x','yyyy')
[hleg,objh] = legend('HadCRUT4.3 + GraphEM (0.8% sparsity)', 'GISTEMP', 'NOAA');
set(hleg,'Location','SouthEast')
set(objh,'linewidth',2)
legend('boxoff')
fancyplot_deco('Monthly Global Mean Temperature', 'Time (Year)', ylabel_str, 14, 'Helvetica')
subplot(2,1,2)
hold on;
plot(HG.tmdn(ind1), HG.GMTm(ind1), 'color',skyblue);
plot(G.tmdn(ind2), G.GMTm(ind2), 'color', ornj);
plot(N.tmdn(ind3),N.GMTm(ind3),'color',maroon);
datetick('x','yyyy')
[hleg,objh] = legend('HadCRUT4.3+ GraphEM (0.8% sparsity)', 'GISTEMP','NOAA');
set(hleg,'Location','SouthEast')
set(objh,'linewidth',2)
pause;
legend('boxoff')
fancyplot_deco('Same, 1985 - present', 'Time (Year)', ylabel_str, 14, 'Helvetica')

fig2name = 'had43med_gmt_plot2.eps';
fig2path = [odir fig2name];
hepta_figprint(fig2path)
%plot(




