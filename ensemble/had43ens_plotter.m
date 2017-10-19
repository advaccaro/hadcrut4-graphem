%had43ens_plotter.m
addpath(genpath('/home/scec-02/avaccaro/HadCRUT4.3/'))
load('ENS.mat'); %contains globalmean/nino34 of 100-member ensemble
load('had43med.mat'); %contains raw had43median and raw+sat update temperature fields
addpath(genpath('/home/scec-02/avaccaro/HadCRUT4.2/'))
load('NINO34_indices.mat'); %contains ERSSTv3, Kaplan, HadSST2 NINO34
load JEG_graphics
years = [1850:2014]; mons = [1:1980]; %THIS IS INCORRECT, NEED TO FIX THIS TO BE MONTHLY-RESOLVED
ny = length(years);
nm = ny*12;

%% TEMPORARY MEASURE: EXCLUDE CORRUPT ENTRIES
nums = [1:71];
exclude = [4,10,15,29,36,46,71,73]; %corrupt datafiles
include = setdiff(nums,exclude);
ENSf.globalmean = ENS.globalmean(:,include);
ENSf.nino34 = ENS.globalmean(:,include);


%% Calculate quantiles, median, mean of 100-member ensemble
%Global annual mean
for i = 1:nm
globalmean_sum(i,:) = quantile(ENSf.globalmean(i,:),[.05 .25 .5 .75 .95]);
end


%Nino34
for i = 1:nm
nino34_sum(i,:) = quantile(ENSf.nino34(i,:),[.05 .25 .5 .75 .95]);
end

%% Raw HadCRUT4 Median
%Global annual mean
[nlon, nlat, ntime] = size(H43med.x);
nloc = nlat*nlon;

%Hrawa = nan(ny,nloc); %WANT THIS IN MONTHLY RESOLUTION
%for i = 1:nloc
%for j = 1:ny
%Hrawa(j,i)=nmean(rawH(1+(j-1)*12:12*j,i));
%end
%end
%temp_1 = nan(ny,1); %WANT THIS IN MONTHLY RESOLUTION
%weights1 = repmat(cosd(loc(:,2)),[1,ny]);
%for i = 1:ny %WANT THIS IN MONTHLY RESOLUTION 
%weights1(isnan(Hrawa(i,:)),i) = 0;
%Hrawa(i,(isnan(Hrawa(i,:)))) = 0;
%temp_1(i) = (Hrawa(i,:)*weights1(:,i))/sum(weights1(:,i));
%end

%Nino34
Xin = had43med;
Xin = center(Xin);
Xc_3d = reshape(Xin',[nlon,nlat,ntime]);
[nt, np] = size(Xin);
lon = H43med.lon;
lat = H43med.lat;
Raw_3d = reshape(rawH', [nlon,nlat,ntime]);
time = 675700 + H43med.time;
[satnino12,satnino3,satnino4,satnino34] = nino_indices(Xc_3d, lon, lat);
[rawnino12,rawnino3,rawnino4,rawnino34] = nino_indices(Raw_3d, lon, lat);

%% Raw Had4 Median + Sattelite update
%Global annual mean
Hsat = had43med; %2d already

%Annualize
Hsata = nan(ny,nloc); %WANT THIS IN MONTHLY RESOLUTION

for i = 1:nloc
for j = 1:ny %WANT THIS IN MONTHLY RESOLUTION
Hsata(j,i) = nmean(Hsat(1+(j-1)*12:12*j,i));
end
end

%create global average time series
temp_2 = nan(ny,1); %WANT THIS IN MONTHLY RESOLUTION
weights2 = repmat(cosd(loc(:,2)),[1,ny]); %WANT THIS IN MONTHLY RESOLUTION

for i = 1:ny %WANT THIS IN MONTHLY RESOLUTION
weights2(isnan(Hsata(i,:)),i) = 0;
Hsata(i,(isnan(Hsata(i,:)))) = 0;
temp_2(i) = (Hsata(i,:)*weights2(:,i))/sum(weights2(:,i));
end

%Nino34

%% Plotting  (NOTHING REALLY WORKS BELOW THIS POINT!  FIX EVERYTHING!!!)
fig('1')
hold on;
%plot(years,globalmean_sum)
ciplot_maud(globalmean_sum(:,1),globalmean_sum(:,5),mons,ligr)
ciplot_maud(globalmean_sum(:,2),globalmean_sum(:,4),mons,skyblue)
plot(mons,globalmean_sum(:,3),'Color',blue)
plot(mons,temp_1,'Color',bright_green) %RAW
plot(mons,temp_2,'Color',red) %RAW + SAT UPDATE
hleg = legend('GraphEM ensemble 5-95% C.I.','GraphEM ensemble 25-75% C.I.', 'GraphEM ensemble median','HadCRUT4.3 ensemble median','HadCRUT4.3 median w/ Satellite' )
set(hleg, 'Location', 'NorthWest')
legend('boxoff')
fancyplot_deco('12-Month Global Mean Temperature over Land & Ocean','Year','deg C')
hepta_figprint('had43ens_globalmean_v1')

%decadal smoothing
rawmean_smooth = hepta_smooth(temp_1,1/10);
satmean_smooth = hepta_smooth(temp_2,1/10);
globalmean_sum_smooth = nan(ny,5);
for i = 1:5
globalmean_sum_smooth(:,i) = hepta_smooth(globalmean_sum(:,i),1/10);
end

%trend lines
rawmean_coeffs = polyfit(mons',rawmean_smooth,1);
satmean_coeffs = polyfit(mons',satmean_smooth,1);
globalmean_sum_coeffs = nan(5,2);

for i = 1:5
globalmean_sum_coeffs(i,:) = polyfit(mons',globalmean_sum_smooth(:,i),1);
end

rawmean_trend = polyval(rawmean_coeffs,years);
satmean_trend = polyval(satmean_coeffs,years);
globalmean_sum_trend = nan(5,nm);
for i = 1:5
globalmean_sum_trend(i,:) = polyval(globalmean_sum_coeffs(i,:),mons);
end

fig('2')
hold on;
subplot(2,1,1)
hold on;
ciplot_maud(globalmean_sum(:,1),globalmean_sum(:,5),years,ligr)
ciplot_maud(globalmean_sum(:,2),globalmean_sum(:,4),years,skyblue)
plot(years,globalmean_sum(:,3),'Color',blue)
plot(years,temp_1,'Color',bright_green) %RAW
plot(years,temp_2,'Color',red) %RAW + SAT UPDATE
hleg = legend('GraphEM ensemble 5-95% C.I.','GraphEM ensemble 25-75% C.I.', 'GraphEM ensemble median','HadCRUT4.3 ensemble median','HadCRUT4.3 median w/ interpoalted data update' )
set(hleg, 'Location', 'NorthWest')
legend('boxoff')
fancyplot_deco('Global Mean Monthly Temperature over Land & Ocean','Year','deg C')
subplot(2,1,2)
hold on;
ciplot_maud(globalmean_sum_smooth(:,1),globalmean_sum_smooth(:,5),years,ligr)
ciplot_maud(globalmean_sum_smooth(:,2),globalmean_sum_smooth(:,4),years,skyblue)
plot(years,globalmean_sum_smooth(:,3),'Color',blue)
plot(years,rawmean_smooth,'Color',bright_green)
plot(years,satmean_smooth,'Color',red)

%add trend lines
plot(years,globalmean_sum_trend(3,:),'Color',blue,'LineStyle','--')
plot(years,rawmean_trend,'Color',bright_green,'LineStyle','--')
plot(years,satmean_trend,'Color',red,'LineStyle','--')
fancyplot_deco('idem, 10-year lowpass filter', 'Year', 'deg C')

hepta_figprint('had43ens_globalmean_v2')


fig('3')
hold on;
ciplot_maud(nino34_sum(1:end-1,1),nino34_sum(1:end-1,5),time,ligr)
ciplot_maud(nino34_sum(1:end-1,2),nino34_sum(1:end-1,4),time,skyblue)
plot(time,nino34_sum(1:end-1,3),'Color',blue)
plot(time,rawnino34,'Color',bright_green)
plot(time,satnino34,'Color',red)
datetick('x','yyyy')
hleg2 = legend('GraphEM ensemble 5-95% C.I.','GraphEM ensemble 25-75% C.I.', 'GraphEM ensemble median','HadCRUT4.3 ensemble median','HadCRUT4.3 median w/ Satellite' )
legend('boxoff')
fancyplot_deco('NINO3.4 SST Anomaly', 'Year', 'deg C')
hepta_figprint('had43ens_nino34')

%annualize
rawnino34_ann = nan(ny,1);
satnino34_ann = nan(ny,1);
nino34_sum_ann = nan(ny,5);

for i = 1:ny
rawnino34_ann(i) = nmean(rawnino34(1+(i-1)*12:12*i));
satnino34_ann(i) = nmean(satnino34(1+(i-1)*12:12*i));

for j = 1:5
nino34_sum_ann(i,j) = nmean(nino34_sum(1+(i-1)*12:12*i,j));
end
end


%get rid of nan's
rawnino34_ann_nmean1 = nmean(rawnino34_ann);
satnino34_ann_nmean1 = nmean(satnino34_ann);

rawnino34_nanind = isnan(rawnino34_ann);
satnino34_nanind = isnan(satnino34_ann);

rawnino34_ann(rawnino34_nanind) = rawnino34_ann_nmean1;
satnino34_ann(satnino34_nanind) = satnino34_ann_nmean1;

%10-year smoothing
nino34_sum_annsmooth = nan(ny,5);

rawnino34_annsmooth = hepta_smooth(rawnino34_ann,1/10);
satnino34_annsmooth = hepta_smooth(satnino34_ann,1/10);

for i = 1:5
nino34_sum_annsmooth(:,i) = hepta_smooth(nino34_sum_ann(:,i),1/10);
end


%Trend lines
rawnino34_coeffs = polyfit(years',rawnino34_annsmooth,1);
satnino34_coeffs = polyfit(years',satnino34_annsmooth,1);
nino34_sum_coeffs = nan(5,2);

for i = 1:5
nino34_sum_coeffs(i,:) = polyfit(years',nino34_sum_annsmooth(:,i),1);
end

rawnino34_trend = polyval(rawnino34_coeffs,years);
satnino34_trend = polyval(satnino34_coeffs,years);
nino34_sum_trend = nan(5,ny);


for i = 1:5
nino34_sum_trend(i,:) = polyval(nino34_sum_coeffs(i,:),years);
end

fig('4')
hold on;
ciplot_maud(nino34_sum_ann(:,1),nino34_sum_ann(:,5),years,ligr)
ciplot_maud(nino34_sum_ann(:,2),nino34_sum_ann(:,4),years,skyblue)
plot(years,nino34_sum_ann(:,3),'Color',blue)
plot(years,rawnino34_ann,'Color',bright_green)
plot(years,satnino34_ann,'Color',red)
%datetick('x','yyyy')
hleg4 = legend('GraphEM ensemble 5-95% C.I.','GraphEM ensemble 25-75% C.I.', 'GraphEM ensemble median','HadCRUT4.3 ensemble median','HadCRUT4.3 median w/ Satellite' )
legend('boxoff')
fancyplot_deco('NINO3.4 SST Anomaly', 'Year', 'deg C')
hepta_figprint('had43ens_nino34_v2')


fig('5')
hold on;
subplot(2,1,1)
hold on;
ciplot_maud(nino34_sum_ann(:,1),nino34_sum_ann(:,5),years,ligr)
ciplot_maud(nino34_sum_ann(:,2),nino34_sum_ann(:,4),years,skyblue)
plot(years,nino34_sum_ann(:,3),'Color',blue)
plot(years,rawnino34_ann,'Color',bright_green)
plot(years,satnino34_ann,'Color',red)
%datetick('x','yyyy')
hleg5 = legend('GraphEM ensemble 5-95% C.I.','GraphEM ensemble 25-75% C.I.', 'GraphEM ensemble median','HadCRUT4.3 ensemble median','HadCRUT4.3 median w/ Satellite' )
legend('boxoff')
fancyplot_deco('NINO3.4 SST Anomaly', 'Year', 'deg C')

subplot(2,1,2)
hold on;
ciplot_maud(nino34_sum_annsmooth(:,1),nino34_sum_annsmooth(:,5),years,ligr)
ciplot_maud(nino34_sum_annsmooth(:,2),nino34_sum_annsmooth(:,4),years,skyblue)
plot(years,nino34_sum_annsmooth(:,3),'Color',blue)
plot(years,rawnino34_annsmooth,'Color',bright_green)
plot(years,satnino34_annsmooth,'Color',red)

%add trend lines
plot(years,nino34_sum_trend(3,:),'Color',blue,'LineStyle','--')
plot(years,rawnino34_trend,'Color',bright_green,'LineStyle','--')
plot(years,satnino34_trend,'Color',red,'LineStyle','--')
fancyplot_deco('idem, 10-year lowpass filter', 'Year', 'deg C')

hepta_figprint('had43ens_nino34_v3')


fig('6')
hold on;
subplot(2,1,1)
hold on;
%ciplot_maud(nino34_sum_ann(:,1),nino34_sum_ann(:,5),years,ligr)
%ciplot_maud(nino34_sum_ann(:,2),nino34_sum_ann(:,4),years,skyblue)
plot(years,nino34_sum_ann(:,3),'Color',blue)
plot(years,rawnino34_ann,'Color',bright_green)
plot(years,satnino34_ann,'Color',red)
plot(ERSSTv3.ta,ERSSTv3.NINO34,'Color',firebrick)
plot(Kaplan.ta,Kaplan.NINO34,'Color',pink)
plot(HadSST2i.ta,HadSST2i.NINO34,'Color',maroon)
%datetick('x','yyyy')
hleg6 = legend('GraphEM','HadCRUT4.3','HadCRUT4.3 w/ Satellite', 'ERSSTv3', 'Kaplan', 'HadSST2' )
legend('boxoff')
fancyplot_deco('NINO3.4 SST Anomaly', 'Year', 'deg C')

subplot(2,1,2)
hold on;
%ciplot_maud(nino34_sum_annsmooth(:,1),nino34_sum_annsmooth(:,5),years,ligr)
%ciplot_maud(nino34_sum_annsmooth(:,2),nino34_sum_annsmooth(:,4),years,skyblue)
plot(years,nino34_sum_annsmooth(:,3),'Color',blue)
plot(years,rawnino34_annsmooth,'Color',bright_green)
plot(years,satnino34_annsmooth,'Color',red)
plot(ERSSTv3.ta,ERSSTv3.NINO34_lw,'Color', firebrick)
plot(Kaplan.ta,Kaplan.NINO34_lw,'Color',pink)
plot(HadSST2i.ta,HadSST2i.NINO34_lw,'Color',maroon)

%add trend lines
plot(years,nino34_sum_trend(3,:),'Color',blue,'LineStyle','--')
plot(years,rawnino34_trend,'Color',bright_green,'LineStyle','--')
plot(years,satnino34_trend,'Color',red,'LineStyle','--')
plot(ERSSTv3.ta,ERSSTv3.NINO34_trend,'Color', firebrick,'LineStyle','--')
plot(Kaplan.ta,Kaplan.NINO34_trend,'Color',pink,'LineStyle','--')
plot(HadSST2i.ta,HadSST2i.NINO34_trend,'Color',maroon,'LineStyle','--')
fancyplot_deco('idem, 10-year lowpass filter', 'Year', 'deg C')

hepta_figprint('had43ens_nino34_v5', 400)
