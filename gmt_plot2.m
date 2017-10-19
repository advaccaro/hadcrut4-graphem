%gmt_plot2.m


%Load GMT datasets
load_gmt_datasets

%Load ensemble
load('/home/scec-02/avaccaro/HadCRUT4.3/ensemble/data/ENS_SP80.mat')

%Abbreviate dataset names
HR = H43MRAW_EXP;
%HU = H43MUPD_EXP;
HG = H43M80_EXP;
G = GISTEMP_GMT;
CW = CW_EXP;
N = NOAA_GMT;
%E = ENS; %screen broken records out here (for now)
%ind_in = setdiff(1:97,19); %#19, (file 26) screened out
%E.GMTm = ENS.GMTm(:,ind_in);
%E.NINO34m = ENS.GMTm(:,ind_in);
E.GMTm = ENS.GMTm;
%E.NINO34m = ENS.G


%Flip datenum time axes (because it's a pain in the ass)
HR.tmdn = flipud(HR.tmdn);
%HU.tmdn = flipud(HU.tmdn);
HG.tmdn = flipud(HG.tmdn);
%G.tmdn = flipud(G.tmdn); %already oriented in the proper direction!  Hooray!
CW.tmdn = flipud(CW.tmdn);



%% Calculate differences (to HR)
%HadCRUT family
%HU.GMTm_diff = HU.GMTm - HR.GMTm;
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




% ENS
[nt,nf] = size(E.GMTm);
E.tmfrac = nan(nt,1); %set up time axis
for i = 1:nt
	E.tmfrac(i) = 1850 + (i-1)/12;
end
E.tmdn = HR.tmdn;
E.tmdn(end+1) = E.tmdn(end) + 30;
%calculate differences
[dind1,dind2] = find_overlap(E.tmfrac,HR.tmfrac);
E.GMTm_diff = E.GMTm(dind1,:) - repmat(HR.GMTm(dind2),[1 nf]);
%calculate quantiles, median
for i = 1:nt
	GMsum(i,:) = quantile(E.GMTm(i,:), [.025 .25 .5 .75 .975]);
end

for i = 1:nt-1
	GMsum_diff(i,:) = quantile(E.GMTm_diff(i,:), [.025 .25 .5 .75 .975]);
end

%% Get indices of recent (1985-present)
HG.trec_ind = find(HG.tmfrac >= 1985); ind1 = HG.trec_ind;
G.trec_ind = find(G.tmfrac >= 1985); ind2 = G.trec_ind;
N.trec_ind = find(N.tmfrac >= 1985); ind3 = N.trec_ind;
E.trec_ind = find(E.tmfrac >= 1985); ind4 = E.trec_ind;


E.GMTm_rec = E.GMTm(ind4,:); nrec = length(ind4);
for i = 1:nrec
	GMsum_rec(i,:) = quantile(E.GMTm_rec(i,:), [.025 .25 .5 .75 .975]);
end




% Calculate envelopes for ensemble
%(Outer envelope = 5-95% ||| inner envelope = 25-75%)
[vertxsqO,y2O,ZO] = fill3_prep(E.tmdn',GMsum(:,5)',GMsum(:,1)');
[vertxsqI,y2I,ZI] = fill3_prep(E.tmdn',GMsum(:,4)',GMsum(:,2)');
[vertxsqO_diff,y2O_diff,ZO_diff] = fill3_prep(HR.tmdn',GMsum_diff(:,5)',GMsum_diff(:,1)');
[vertxsqI_diff,y2I_diff,ZI_diff] = fill3_prep(HR.tmdn',GMsum_diff(:,4)',GMsum_diff(:,2)');
[vertxsqO_rec,y2O_rec,ZO_rec] = fill3_prep(E.tmdn(ind4)',GMsum_rec(:,5)',GMsum_rec(:,1)');
[vertxsqI_rec,y2I_rec,ZI_rec] = fill3_prep(E.tmdn(ind4)',GMsum_rec(:,4)',GMsum_rec(:,2)');




%% Plotting
load JEG_graphics

%GMT Panel 1 (HadCRUT family w/ differences)
fig('GMT plot 1'); clf;
%subplot(4,1,1);
subplot(2,1,1);
hold on; %xlim = ([1850,2020]);
hp1 = fill3(vertxsqO,y2O,ZO,ligr); %plot outer envelope
hp2 = fill3(vertxsqI,y2I,ZI,ligr); %plot inner envolope
%alpha(hp1,.3); alpha(hp2,.3); %reduce opacity of envelopes
set(hp1,'EdgeColor',ligr); set(hp2,'EdgeColor',ligr);
%set(hp1,'EdgeAlpha',.3); set(hp2,'EdgeAlpha',.3); %remove edge lines
p1=plot(HG.tmdn,HG.GMTm,'color','k' );
p2=plot(HR.tmdn,HR.GMTm,'color','g');
p3=plot(CW.tmdn,CW.GMTm,'color','r');
datetick('x','yyyy')
[hleg,objh] = legend([p1,p2,p3],{'GraphEM ensemble median (0.8% sparsity)', 'HadCRUT4.3 median raw', 'HadCRUT4.3 Cowtan & Way'});
%[hleg,objh] = legend('HadCRUT4.3 median raw', 'HadCRUT4.3 Cowtan & Way', 'HadCRUT4.3 + GraphEM (0.8% sparsity)');
set(hleg,'Location','SouthEast')
set(objh,'linewidth',2)
legend('boxoff')
ylabel_str = ['Temperature (' char(176) 'C)'];
fancyplot_deco('Monthly Global Mean Temperature', 'Time (Year)', ylabel_str,14,'Helvetica')




%Panel 2 (HadCRUT family differences)
%subplot(4,1,2);
subplot(2,1,2);
hold on; %xlim = ([1850,2020]);
hp3 = fill3(vertxsqO_diff,y2O_diff,ZO_diff,ligr);
hp4 = fill3(vertxsqI_diff,y2I_diff,ZI_diff,ligr);
set(hp3, 'EdgeColor',ligr); set(hp4,'EdgeColor',ligr);
p4=plot(HG.tmdn,HG.GMTm_diff,'color','k');
p5=plot(CW.tdiff,CW.GMTm_diff,'color','r');

datetick('x','yyyy')

[hleg,objh] = legend([p4,p5],{'GraphEM ensemble median (0.8% sparsity)','HadCRUT4.3 Cowtan & Way'});
set(hleg,'Location','SouthEast')
set(objh,'linewidth',2)
legend('boxoff')
fancyplot_deco('Same, difference to raw HadCRUT4.3 median', 'Time (Year)', '{\Delta}T ({\circ}C)',14,'Helvetica')

%save/print
odir = '/home/scec-02/avaccaro/HadCRUT4.3/figs/';
figname = 'gmt_plot1.eps';
figpath = [odir figname];
hepta_figprint(figpath)




%Panel 3 (GISTEMP w/ recent)
%ind1 = HG.trec_ind; ind2 = G.trec_ind; ind3 = N.trec_ind;
%fig('The Gs (GraphEM and GISS)'); clf;
fig('GMT plot 2');
%subplot(4,1,3);
subplot(2,1,1);
hold on;
hp5 = fill3(vertxsqO,y2O,ZO,ligr);
hp6 = fill3(vertxsqI,y2I,ZI,ligr);
set(hp5,'EdgeColor',ligr); set(hp6,'EdgeColor',ligr);
p6=plot(HG.tmdn,HG.GMTm,'color','k'); 
p7=plot(G.tmdn,G.GMTm,'color',ornj);
p8=plot(N.tmdn,N.GMTm,'color',maroon);
datetick('x','yyyy')
[hleg,objh] = legend([p6,p7,p8], {'GraphEM ensemble median (0.8% sparsity)', 'GISTEMP', 'NOAA'});
set(hleg,'Location','SouthEast')
set(objh,'linewidth',2)
legend('boxoff')
fancyplot_deco('Monthly Global Mean Temperature', 'Time (Year)', ylabel_str, 14, 'Helvetica')

%Panel 4 (recent)
%subplot(4,1,4)
subplot(2,1,2)
hold on;
hp7 = fill3(vertxsqO_rec,y2O_rec,ZO_rec,ligr);
hp8 = fill3(vertxsqI_rec,y2I_rec,ZI_rec,ligr);
set(hp7,'EdgeColor',ligr); set(hp8,'EdgeColor',ligr);
p9=plot(HG.tmdn(ind1), HG.GMTm(ind1), 'color','k');
p10=plot(G.tmdn(ind2), G.GMTm(ind2), 'color', ornj);
p11=plot(N.tmdn(ind3),N.GMTm(ind3),'color',maroon);
datetick('x','yyyy')
[hleg,objh] = legend([p9,p10,p11],{'GraphEM ensemble median (0.8% sparsity)', 'GISTEMP','NOAA'});
set(hleg,'Location','SouthEast')
set(objh,'linewidth',2)
legend('boxoff')
fancyplot_deco('Same, 1985 - present', 'Time (Year)', ylabel_str, 14, 'Helvetica')



%save/print
odir = '/home/scec-02/avaccaro/HadCRUT4.3/figs/';
figname = 'gmt_plot2.eps';
figpath = [odir figname];
hepta_figprint(figpath)







