%had43ens_gmmtplot.m
%server version
%makes GMMT plot for 100-member GraphEM-infilled HadCRUT4.3 ensemble
%plots against GISTEMP, raw HadCRUT4.3, and updated HadCRUT4.3

%% LOAD DATA
%addpath(genpath('/Users/adam/Desktop/HadCRUT4.3/'))
addpath(genpath('/home/scec-02/avaccaro/HadCRUT4.3/'))

%raw had43med
load Had43med_RAW_NINO

%had43med updated
load Had43med_UPD_NINO

%GISTEMP
load GISTEMP

%GraphEM infilled had43ens
load ENS

%addpath(genpath('/Users/adam/Desktop/Treerings'))
load JEG_graphics

%% Useful metadata
years = [1850:2014]; ny = length(years);
mons = flipud(datenum(Had43med_RAW.tm)); %monthly fractional time axis GOING FORWARD IN TIME
[nm,nlat,nlon] = size(Had43med_RAW.SSTmonthly);
nloc = nlat * nlon;

%% Prepare ENS
%Screen out corrupt records
nrecs = size(ENS.globalmean,2); %number of records in ensemble
nums = [1:nrecs];
exclude = [2,7,20,28,38,42,65,67,72]; %indices of corrupt datafiles
include = setdiff(nums,exclude); %records kept in
ENSf.globalmean = ENS.globalmean(:,include);

%Calculate quantiles, median
for i =1:nm
    GMsum(i,:) = quantile(ENSf.globalmean(i,:),[.05 .25 .5 .75 .95]);
end
ENSf.GMsum = GMsum;

%% Had43med differences
Had43med_RAW.tmf = flipud(datenum(Had43med_RAW.tm));
Had43med_RAW.tvec = datevec(Had43med_RAW.tmf);
Had43med_RAW.tfrac = Had43med_RAW.tvec(:,1) + Had43med_RAW.tvec(:,2)/12 - 1/12;
Had43med_RAW.tmin = min(Had43med_RAW.tfrac); Had43med_RAW.tmax = max(Had43med_RAW.tfrac);
lowerlim = Had43med_RAW.tmin; upperlim = Had43med_RAW.tmax;
Had43med_UPD.diff = Had43med_UPD.GMMT - Had43med_RAW.GMMT;

for i = 1:5
	ENSf.diff(:,i) = GMsum(:,i) - Had43med_RAW.GMMT;
end

%% GISTEMP
GISTEMP.tvec = datevec(GISTEMP.tm);
GISTEMP.tfrac = GISTEMP.tvec(:,1) + GISTEMP.tvec(:,2)/12 - 1/12;
GISTEMP.tmin = min(GISTEMP.tfrac); GISTEMP.tmax = max(GISTEMP.tfrac);
GISTEMP.diffind1 = find(GISTEMP.tfrac >= lowerlim & GISTEMP.tfrac <= upperlim);
GISTEMP.diffind2 = find(Had43med_RAW.tfrac >= GISTEMP.tmin & Had43med_RAW.tfrac <= GISTEMP.tmax);
GISTEMP.tdiff = intersect(GISTEMP.tfrac,Had43med_RAW.tfrac);
GISTEMP.tind = ismember(GISTEMP.tfrac, GISTEMP.tdiff);
GISTEMP.diff = GISTEMP.GMMT(GISTEMP.diffind1)*0.01 - Had43med_RAW.GMMT(GISTEMP.diffind2);


%% Cowtan and Way
load('CW_EXP.mat')


%% Annualize (No point?!)
% raw_ann = nan(ny,1); %raw_ann_diff = nan(ny,1);
% upd_ann = nan(ny,1); upd_ann_diff = nan(ny,1);
% GISTEMP_ann = nan(ny,1); GISTEMP_ann_diff = nan(ny,1);
% for i = 1:5
%     ENSf.ann(:,i) = nan(ny,1); ENSf.ann_diff(:,i) = nan(ny,1);
% end
% for i = 1:ny
%     raw_ann(i) = nmean(Had43med_RAW.GMMT(1+(i-1)*12:12*i));
%     raw_ann_diff(i) = nmean(Had43med_RAW.diff
%     upd_ann(i) = nmean(Had43med_UPD.GMMT(1+(i-1)*12:12*i));
%     upd_ann_diff(i) = nmean(Had43med_UPD.diff(1+(i-1)*12:12*i));
%     GISTEMP_ann(i) = 
%     GISTEMP_ann_diff(i) = 
%     for j = 1:5
%        ENSf.ann(i,j) = 



%% Decadal smoothing (going to skip this for now)
rawmean_smooth = hepta_smooth(Had43med_RAW.GMMT, 1/120);
updmean_smooth = hepta_smooth(Had43med_UPD.GMMT, 1/120);
updmean_diff_smooth = hepta_smooth(Had43med_UPD.diff, 1/120);
for i = 1:5
GMsum_smooth(:,i) = hepta_smooth(GMsum(:,i), 1/120);
ENSf.diff_smooth(:,i) = hepta_smooth(ENSf.diff(:,i),  1/120);
end
gistemp_smooth = hepta_smooth(GISTEMP.GMMT, 1/120);
gistemp_diff_smooth = hepta_smooth(GISTEMP.diff, 1/120);
%% Trend lines (going to skip this for now)


%% Calculate envelopes for ensemble
%(Outer envelope = 5-95% ||| inner envelope = 25-75%)
XcasesO = mons'; %row vector 
NcasesO = length(XcasesO); SqcasesO = [XcasesO XcasesO];
XcasesI = XcasesO; %will need to change this if there is ever a difference
NcasesI = NcasesO; %see above
SqcasesI = SqcasesO; %see above

topO = GMsum(:,5); topO = topO'; %row vector
botO = GMsum(:,1); botO = botO'; %row vector
topI = GMsum(:,4); topI = topI'; %row vector
botI = GMsum(:,2); botI = botI'; %row vector

vO = zeros(2*NcasesO,1); vI = zeros(2*NcasesI,1);

SvertsO = zeros(2*NcasesO,1); SvertsI = zeros(2*NcasesI,1);

for j = 1:NcasesO
SvertsO(2*j) = XcasesO(j);
SvertsO(1+(2*(j-1))) = XcasesO(j);
vO(2*j) = topO(j);
vO(1+2*(j-1)) = botO(j);
end

for j = 1:NcasesI
SvertsI(2*j) = XcasesI(j);
SvertsI(1+(2*(j-1))) = XcasesI(j);
vI(2*j) = topI(j);
vI(1+2*(j-1)) = botI(j);
end

yO = [botO, fliplr(topO)]; y2O = [yO, yO(1)];
yI = [botI, fliplr(topI)]; y2I = [yI, yI(1)];

vertxO = [XcasesO, fliplr(XcasesO)];
vertxsqO = [vertxO, vertxO(1)];
ZO = zeros(size(vertxsqO));

vertxI = [XcasesI, fliplr(XcasesI)];
vertxsqI = [vertxI, vertxI(1)];
ZI = zeros(size(vertxsqI));



%% Plotting

%% Test plots below:

%HadCRUT4.3 GMMT plot (w/ envelopes)
fig('had43ens gmmt plot w/ envelopes'); clf;
hold on;
hp1 = fill3(vertxsqO,y2O,ZO,'k'); %plot outer envelope
hp2 = fill3(vertxsqI,y2I,ZI,'b'); %plot inner envelope
alpha(hp1,.3); alpha(hp2,.3); %reduce opacity of envelopes
set(hp1,'EdgeAlpha',0.3); set(hp2,'EdgeAlpha',0.3); %remove edge lines
set(hp2,'EdgeColor','b');
plot(mons,GMsum(:,3),'color',skyblue)
plot(mons,Had43med_RAW.GMMT,'g')
plot(mons,Had43med_UPD.GMMT,'r')
plot(GISTEMP.tm,GISTEMP.GMMT*0.01,'color',ornj)
plot(flipud(CW_EXP.tmdn),CW_EXP.GMTm,'color',dark_violet);
datetick('x','yyyy')
[hleg,objh] = legend('GraphEM ensemble 5-95% C.I.', 'GraphEM ensemble 25-75% C.I.', 'GraphEM ensemble median', 'Raw HadCRUT4.3 median', 'Updated HadCRUT4.3 median', 'GISTEMP', 'Cowtan & Way');
set(hleg,'Location','SouthEast')
set(objh, 'linewidth', 2)
legend('boxoff')
fancyplot_deco('Global Mean Monthly Temperature','Time (Year)','Temperature ({\circ}C)',14,'Helvetica')
%hepta_figprint('./figs/had43ens_gmmtplot_all_errorbars.eps')
print -painters -dpdf -cmyk -r1000 './figs/had43ens_gmmtplot_all_errorbars.pdf'


%test plot (only hadcrut family + ensemble w/ errorbars)
fig('had43ens gmmt plot had fam w/ envelopes'); clf;
hold on;
hp1 = fill3(vertxsqO,y2O,ZO,'k'); %plot outer envelope
hp2 = fill3(vertxsqI,y2I,ZI,'b'); %plot inner envelope
alpha(hp1,.3); alpha(hp2,.3); %reduce opacity of envelopes
%set(hp1,'EdgeAlpha',0); set(hp2,'EdgeAlpha',0); %remove edge lines
set(hp2,'EdgeColor',b);
plot(mons,GMsum(:,3),'color',skyblue)
plot(mons,Had43med_RAW.GMMT,'g')
plot(mons,Had43med_UPD.GMMT,'r')
%plot(GISTEMP.tm,GISTEMP.GMMT*0.01,'color',ornj)
datetick('x','yyyy')
[hleg,objh] = legend('GraphEM ensemble 5-95% C.I.', 'GraphEM ensemble 25-75% C.I.', 'GraphEM ensemble median', 'Raw HadCRUT4.3 median', 'Updated HadCRUT4.3 median');
set(hleg,'Location','SouthEast')
set(objh, 'linewidth', 2)
legend('boxoff')
fancyplot_deco('Global Mean Monthly Temperature','Time (Year)','Temperature ({\circ}C)',14,'Helvetica')
hepta_figprint('./figs/had43ens_gmmtplot_hadfam_errorbars.eps')

%test plot (only hadcrut ensemble w/ errorbars)
fig('had43ens gmmt plot only ens w/ envelopes'); clf;
hold on;
hp1 = fill3(vertxsqO,y2O,ZO,'k'); %plot outer envelope
hp2 = fill3(vertxsqI,y2I,ZI,'b'); %plot inner envelope
alpha(hp1,.3); alpha(hp2,.3); %reduce opacity of envelopes
set(hp1,'EdgeAlpha',0.3); set(hp2,'EdgeAlpha',0.3); %remove edge lines
set(hp2,'EdgeColor','b');
plot(mons,GMsum(:,3),'color',skyblue)
%plot(mons,Had43med_RAW.GMMT,'g')
%plot(mons,Had43med_UPD.GMMT,'r')
%plot(GISTEMP.tm,GISTEMP.GMMT*0.01,'color',ornj)
datetick('x','yyyy')
[hleg,objh] = legend('GraphEM ensemble 5-95% C.I.', 'GraphEM ensemble 25-75% C.I.', 'GraphEM ensemble median');
set(hleg,'Location','SouthEast')
set(objh, 'linewidth', 2)
legend('boxoff')
fancyplot_deco('Global Mean Monthly Temperature','Time (Year)','Temperature ({\circ}C)',14,'Helvetica')
%hepta_figprint('./figs/had43ens_gmmtplot_ens_errorbars.eps')
print -painters -dpdf -cmyk -r1000 './figs/had43ens_gmmtplot_ens_errorbars.pdf'

%test plot (only errorbars)
fig('had43ens gmmt plot only ens w/ envelopes'); clf;
hold on;
hp1 = fill3(vertxsqO,y2O,ZO,'k'); %plot outer envelope
hp2 = fill3(vertxsqI,y2I,ZI,'b'); %plot inner envelope
alpha(hp1,.3); alpha(hp2,.3); %reduce opacity of envelopes
set(hp1,'EdgeAlpha',0.3); set(hp2,'EdgeAlpha',0.3); %remove edge lines
set(hp2,'EdgeColor','b');
%plot(mons,GMsum(:,3),'color',skyblue)
%plot(mons,Had43med_RAW.GMMT,'g')
%plot(mons,Had43med_UPD.GMMT,'r')
%plot(GISTEMP.tm,GISTEMP.GMMT*0.01,'color',ornj)
datetick('x','yyyy')
[hleg,objh] = legend('GraphEM ensemble 5-95% C.I.', 'GraphEM ensemble 25-75% C.I.');
set(hleg,'Location','SouthEast')
set(objh, 'linewidth', 2)
legend('boxoff')
fancyplot_deco('Global Mean Monthly Temperature','Time (Year)','Temperature ({\circ}C)',14,'Helvetica')
%hepta_figprint('./figs/had43ens_gmmtplot_only_errorbars.eps')
print -painters -dpdf -cmyk -r1000 './figs/had43ens_gmmtplot_only_errorbars.pdf'





%% Real plots below: (commented out for now)

%fig('Had43ens sp70 gmmt plot all')
%subplot(2,1,1)
%ciplot_maud(GMsum(:,1),GMsum(:,5),mons,'k')
%hold on;
%ciplot_maud(GMsum(:,2),GMsum(:,4),mons,skyblue)
%plot(mons,GMsum(:,3),'b')
%plot(mons, Had43med_RAW.GMMT, 'color', 'g')
%plot(mons,Had43med_UPD.GMMT,'r')
%plot(GISTEMP.tm,GISTEMP.GMMT.*0.01, 'Color', ornj)

%datetick('x','yyyy')
%%legtext = ['GraphEM ensemble 5-95% C.I.', 'GraphEM ensemble 25-75% C.I.', 'GraphEM ensemble median', 'HadCRUT4.3 ensemble median', 'HadCRUT4.3 median w/ interpoalted data update', 'GISTEMP'];
%%[hleg, objh] = legend(legtext);
%[hleg,objh] = legend('GraphEM ensemble 5-95% C.I.', 'GraphEM ensemble 25-75% C.I.', 'GraphEM ensemble median', 'HadCRUT4.3 ensemble median', 'HadCRUT4.3 median w/ interpoalted data update', 'GISTEMP');
%set(hleg,'Location','SouthEast')
%set(objh, 'linewidth', 2)
%legend('boxoff')
%fancyplot_deco('Global Mean Monthly Temperature','Time (Year)','Temperature ({\circ}C)',14,'Helvetica')
%subplot(2,1,2)
%ciplot_maud(ENSf.diff(:,1),ENSf.diff(:,5),mons,'k');
%hold on;
%ciplot(ENSf.diff(:,2),ENSf.diff(:,4),mons,skyblue);
%plot(mons,ENSf.diff(:,3),'b');
%plot(mons, Had43med_UPD.diff,'color','r')
%plot(GISTEMP.tm(GISTEMP.tind), GISTEMP.diff,'color',ornj)
%datetick('x','yyyy')
%%legtext2 = ['GraphEM ensemble 5-95% C.I.', 'GraphEM ensemble 25-75% C.I.', 'GraphEM ensemble median', 'HadCRUT4.3 median w/ interpolated data update', 'GISTEMP'];
%%[hleg,objh] = legend(legtext2);
%[hleg2,objh2] = legend('GraphEM ensemble 5-95% C.I.', 'GraphEM ensemble 25-75% C.I.', 'GraphEM ensemble median', 'HadCRUT4.3 median w/ interpolated data update', 'GISTEMP');
%set(hleg2,'Location','SouthEast');
%set(objh2,'linewidth',2);
%legend('boxoff')
%fancyplot_deco('Same, difference from HadCRUT4.3 median','Time (Year)', '{\Delta}T ({\circ}C)',14,'Helvetica')
%hepta_figprint('had43ens_gmmtplot_all.eps')
%print -painters -dpdf -cmyk -r1000 had43ens_gmmtplot_diff_all.pdf
