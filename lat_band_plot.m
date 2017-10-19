%lat_band_plot.m


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

%Flip datenum time axes
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
group{1} = HG.LAT_60N_90N;
group2{1} =  HR.LAT_60N_90N;
group3{1} = CW.LAT_60N_90N;
leg_str{1} =  ['60' char(176) 'N - 90' char(176) 'N'];
group{2} = HG.LAT_30N_60N;
group2{2} = HR.LAT_30N_60N;
group3{2} = CW.LAT_30N_60N;
leg_str{2} = ['30' char(176) 'N - 60' char(176) 'N'];
group{3} = HG.LAT_0_30N;
group2{3} = HR.LAT_0_30N;
group3{3} = CW.LAT_0_30N;
leg_str{3} = ['0' char(176) 'N - 30' char(176) 'N'];
group{4} = HG.LAT_30S_0;
group2{4} = HR.LAT_30S_0;
group3{4} = CW.LAT_30S_0;
leg_str{4} = ['30' char(176) 'S - 0' char(176) 'S'];
group{5} = HG.LAT_60S_30S;
group2{5} = HR.LAT_60S_30S;
group3{5} = CW.LAT_60S_30S;
leg_str{5} = ['60' char(176) 'S - 30' char(176) 'S'];
group{6} = HG.LAT_90S_60S;
group2{6} = HR.LAT_90S_60S;
group3{6} = CW.LAT_90S_60S;
leg_str{6} = ['90' char(176) 'S - 60' char(176) 'S'];
title_str{1} = ['Zonal Means of 30' char(176) ' Latitude Bands'];
for i = 2:6
title_str{i} = '';
end
%GMT Figure 1 (HadCRUT family w/ differences)
fig('HadCRUT latitude bands'); clf;
for i = 1:6
hAx(i) = subplot(7,1,i); %grid on;
hold on; xlim([1850,2017]); %ylim([-12,12]);
set(gca,'YGrid','on');
if ismember(i,[1,6])
ylim([-12,12]);
else
ylim([-3,3]);
end

plot(HR.tmdn,group2{i},'color','g'); %grid on;
plot(CW.tmdn,group3{i},'color','r'); %grid on;
plot(HG.tmdn,group{i},'color',skyblue ); %grid on;

datetick('x','yyyy')

%[hleg,objh] = legend('HadCRUT4.3 median raw', 'HadCRUT4.3 Cowtan & Way', 'HadCRUT4.3 + GraphEM (0.8% sparsity)');
%leg_str = ['60' char(176) 'N - 90' char(176) 'N'];
%[hleg,objh] = legend(leg_str{i});
%set(hleg,'Location','NorthWest')
%set(objh,'linewidth',2)
%legend('boxoff')
t=title(leg_str{i});
set(gca,'Fontsize',6.5);
set(t,'FontSize',7);
set(t,'Position',[t.Position(1) t.Position(2)-1.5,t.Position(3)]);

%ylabel_str = ['Temperature Anomaly (' char(176) 'C)'];
%fancyplot_deco(title_str{i}, 'Time (Year)', ylabel_str,8,'Helvetica')
if i == 6
[hL,iL,~,~] = legend('HadCRUT4.3 median raw','HadCRUT4.3 Cowtan & way', 'HadCRUT4.3 + GraphEM (0.8% sparsity)', 'Location', 'Best');
set(hL,'FontSize',8);%,'Orientation','horizontal');
end

end
hAx(7) = subplot(7,1,7);

p1=get(hAx(1),'position');  % position of UL-most axes
p2=get(hAx(7),'position'); % ditto for LR-most
hAxOuter=axes('position',[p1(1) p2(2) p2(1)+p2(3)-p1(1)  p1(2)+p1(4)-p2(2)], ...
              'color','none','visible','off'); % an outer axis for global p
delete(hAx(7))            % remove the unwanted last subplot axes
y_str = ['Temperature anomaly (' char(176) 'C)'];
x_str = ['Time (Year)'];
hX=text(-0.1,0.5,y_str,'rotation',90,'fontsize',11, ...
        'horizontalalignment','center','verticalalignment','bottom');
hY=text(0.5,0.1,x_str,'rotation',0,'fontsize',11, ...
        'horizontalalignment','center','verticalalignment','top');
suptitle(['Zonal Means for 30' char(176) 'Latitude Bands']);



set(hL,'Position',[0.6, 0.08, 0.85, 0.05]);
set(iL, 'LineWidth',2);
set(hL,'Box','Off');
odir = '/home/scec-02/avaccaro/HadCRUT4.3/figs/';
fig1name = 'had43med_lat_bands.eps';
fig1path = [odir fig1name];
hepta_figprint(fig1path)






