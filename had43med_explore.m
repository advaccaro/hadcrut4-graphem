function [NAME_EXP] = had43med_explore(grid_2d,lon_in,lat_in,time_in,loc)
%% nino_explore.m
% INPUTS: sst_in - 3d sst data (lonxlatxtime)
% 		  lon_in - vector of longitudes
%		  lat_in - vector of latitudes
%		  time_in - time axis in standard serial time (days since Jan-0-0000)



[nmon,nloc] = size(grid_2d);
%sst_3d = reshape(had4med',[nlon, nlat, nmon]);
%sst_3d = sst_in;

% fix grids
lon  = double(lon_in); lon(lon<0) = lon(lon<0) + 360;
lat=double(lat_in); 
nlat = length(lat); nlon = length(lon);

grid_3d = reshape(grid_2d', [nlon nlat nmon]);


%tr = double(cobe.time); % raw time, in units of days since 1891-1-1 00:00:00
%tn = tr + datenum(1891,1,15);  % convert to standard serial date (i.e. add 1800 years)
tn = time_in; %time in standard serial date (days since jan-0-0000)
tv = datevec(tn); 
t = tv(:,1)+tv(:,2)/12-1/12;
%sstr  = double(sst_3d); % may need to increase Java Heap memory space to load (Preferences panel)
%keep everything out to 2013
%sst = sstr(:,:,t<2014);
%sst = sstr;

%t = t(t<2014);
% Dimensions
nx=length(lon); % 5 by 5 deg
ny=length(lat);
ns=nx*ny;
nt=length(t);  nyrs = floor(nt/12);


% GMMT
GRID = grid_2d;
temp = nan(nt,1);
weights = repmat(cosd(loc(:,2)),[1,nt]);
for i = 1:nt
	weights(isnan(GRID(i,:))) = 0;
	GRID(i,(isnan(GRID(i,:)))) = 0;
	temp(i) = (GRID(i,:)*weights(:,i))/sum(weights(:,i));
end
GMMT = temp;

for i = 1:nyrs
	GMMTa_yr(i) = nmean(GMMT(12*(i-1)+1:12*i));
end


% compute climatoly
%ind = find(t >= 1961 & t < 1991); 
%sstp = permute(sst, [3 2 1]);
%sstp12 = sstp(1:nyrs*12,:,:);
%sst_clim = climatology(sstp(ind,:,:));  
%ssta = sstp12 - repmat(sst_clim,[nyrs 1 1]);
%ssta = sst; %already anomalies
grid_3dp = permute(grid_3d, [3 2 1]);

%  Find land points and replace by NaN's
%land = isnan(ssta);

% NINO3/4 computation
NINO3_x=find(lon >=210 & lon <=270);
NINO3_y=find(lat>=-5 & lat <=5);
NINO3=nmean(nmean(grid_3dp(:,NINO3_y,NINO3_x),2),3);
NINO34_x=find(lon >=190 & lon <=240);  % (170W to 120 W)
NINO34=nmean(nmean(grid_3dp(:,NINO3_y,NINO34_x),2),3);
NINO34m = NINO34;


%12 month mean
years12 = [t(1): t(1) + nyrs - 1];
for i = 1:nyrs
	NINO34a_yr(i) = nmean(NINO34(12*(i-1)+1:12*i));
end

% Latitude Band computations
LAT_90S_60S_y = find(lat>= -90 & lat <= -60);
LAT_90S_60S=nmean(nmean(grid_3dp(:,LAT_90S_60S_y,:),2),3);
LAT_60S_30S_y = find(lat>=-60 & lat <= -30);
LAT_60S_30S=nmean(nmean(grid_3dp(:,LAT_60S_30S_y,:),2),3);
LAT_30S_0_y = find(lat >= -30 & lat <= 0);
LAT_30S_0=nmean(nmean(grid_3dp(:,LAT_30S_0_y,:),2),3);
LAT_0_30N_y = find(lat >= 0 & lat <= 30);
LAT_0_30N=nmean(nmean(grid_3dp(:,LAT_0_30N_y,:),2),3);
LAT_30N_60N_y = find(lat >= 30 & lat <= 60);
LAT_30N_60N = nmean(nmean(grid_3dp(:,LAT_30N_60N_y,:),2),3);
LAT_60N_90N_y = find(lat >= 60 & lat <= 90);
LAT_60N_90N = nmean(nmean(grid_3dp(:,LAT_60N_90N_y,:),2),3);

% Lat Bands averaged overtime
lats = [-90:10:90]; nlats = length(lats);
%time_groups = ;

% flip it to most recent first
GRID_3DK=flipdim(grid_3dp,1);
tk=flipud(t);
NINO34f = flipud(NINO34); GMMTf = flipud(GMMT);

% convert time axis using datenum technology
year=floor(tk); 
month=round((tk-year)*12)+1;
tv=[year , month  , repmat(15,nt,1)];
tn=datenum(tv);
tm=datestr(tn);
%  declare arrays
ta=[max(year):-1:min(year)+1]'; na=length(ta);


% Average over DJF period (NH Winter, W)
for k=1:na
	feb=find(tn==datenum([ta(k),2,15,0,0,0]));
	dec=find(tn==datenum([ta(k)-1,12,15,0,0,0]));
	GRID_3DKw(k,:,:)= nmean(GRID_3DK(feb:dec,:,:),1); 
	%SST_Kw(k,:,:)	
	NINO34ka(k)=nmean(NINO34f(feb:dec));
	GMMTa_djf(k)=nmean(GMMTf(feb:dec));
end


NINO34a_djf = NINO34ka;
		
% Time convention  DJF[n]=Dec[n-1] + Jan[n] + Feb[n] , n from 1857 to 2000.

%%%%%%%%%%%%%%%%%%%%%%%%%
% select ocean points and reshape
%tmp=reshape(SST_Kw,na,ns);
%ocean=~isnan(tmp(1,:));
%SSTo=tmp(:,ocean); 
%[nyears,nr]=size(SSTo);
%s=[1:nr];


% NINO3 check
%NINO3w=nmean(nmean(SST_Kw(:,NINO3_y,NINO3_x),2),3);


% Decadal smoothing
NINO34d_djf_in = hepta_smooth(NINO34a_djf(find(~isnan(NINO34a_djf)==1)),1/10);
NINO34d_yr_in = hepta_smooth(NINO34a_yr(find(~isnan(NINO34a_yr)==1)),1/10);
NINO34d_full_in = hepta_smooth(NINO34m(find(~isnan(NINO34m)==1)),1/120);
GMMTd_yr = hepta_smooth(GMMTa_yr,1/10);
GMMTd_full = hepta_smooth(GMMT,1/120);
GMMTd_djf = hepta_smooth(GMMTa_djf,1/10);

NINO34d_djf = nan(length(NINO34a_djf),1);
NINO34d_djf(find(~isnan(NINO34a_djf)==1)) = NINO34d_djf_in;
NINO34d_yr = nan(length(NINO34a_yr),1);
NINO34d_yr(find(~isnan(NINO34a_yr)==1)) = NINO34d_yr_in;
NINO34d_full = nan(length(NINO34m),1);
NINO34d_full(find(~isnan(NINO34m)==1)) = NINO34d_full_in;



%Trend lines

trendcoeffs = polyfit(t,GMMT,1);
trend = polyval(trendcoeffs,t);
trendcoeffsd = polyfit(t,GMMTd_full,1);
trendd = polyval(trendcoeffsd,t);

%Trend from 1998-end
trind = find(t >= 1998);
tmr=t(trind);
GMMTr = GMMT(trind,:);
trendcoeffsr = polyfit(tmr,GMMTr,1);
trendr = polyval(trendcoeffsr,tmr);
GMTrd = GMMTd_full(trind,:);
trendcoeffsrd = polyfit(tmr,GMTrd,1);
trendrd = polyval(trendcoeffsrd,tmr);


%20th century trend
t20ind = find(t >= 1900 & t<= 1999);
tm20 = t(t20ind);
GMT20=GMMT(t20ind,:);
trendcoeffs20 = polyfit(tm20,GMT20,1);
trend20 = polyval(trendcoeffs20,tm20);

%1951-2012
t1951 = find(t >= 1951 & t <= 2012);
tm1951 = t(t1951);
GMT1951 = GMMT(t1951,:);
trendcoeffs1951 = polyfit(tm1951,GMT1951,1);
trend1951 = polyval(trendcoeffs1951,tm1951);
% Graphic check
%fig('Kaplan sanity check')
%subplot(3,1,1:2)
%pcolor(s,ta,SSTo)
%subplot(3,1,3)
%plot(t,NINO3), hold on, axis tight
%plot(ta,NINO3w,'g-',ta,NINO3a,'r-','linewidth',[2]),hold off





%  Assign to structure and save
%
DATASET.SSTm=flipud(GRID_3DK); % monthly values
DATASET.SSTa_djf=flipud(GRID_3DKw);  %  DJF avg
DATASET.lon=lon;
DATASET.lat=lat;
DATASET.ta_yr = years12;
DATASET.ta_djf=ta;
DATASET.tm=flipud(tm(1:length(NINO34m),:));
DATASET.tmdn = flipud(datenum(DATASET.tm));
DATASET.tmfrac = t(1:length(NINO34m));
DATASET.tmr = tmr;
DATASET.tm20 = tm20;
DATASET.tm1951 = tm1951;
DATASET.td_djf = ta;
DATASET.td_yr = years12;
DATASET.td_full = t(1:length(NINO34d_full));
DATASET.NINO34m=NINO34m;
DATASET.NINO34a_djf=NINO34a_djf;
DATASET.NINO34a_yr = NINO34a_yr;
DATASET.NINO34d_djf = NINO34d_djf;
DATASET.NINO34d_yr = NINO34d_yr;
DATASET.NINO34d_full = NINO34d_full;
DATASET.GMTm = GMMT;
DATASET.GMTa_djf = GMMTa_djf;
DATASET.GMTa_yr = GMMTa_yr;
DATASET.GMTd_djf = GMMTd_djf;
DATASET.GMTd_yr =GMMTd_yr;
DATASET.GMTd_full = GMMTd_full;
DATASET.LAT_90S_60S = LAT_90S_60S;
DATASET.LAT_60S_30S = LAT_60S_30S;
DATASET.LAT_30S_0 = LAT_30S_0;
DATASET.LAT_0_30N = LAT_0_30N;
DATASET.LAT_30N_60N = LAT_30N_60N;
DATASET.LAT_60N_90N = LAT_60N_90N;
DATASET.trendcoeffs = trendcoeffs;
DATASET.trendcoeffsd = trendcoeffsd;
DATASET.trendcoeffsr = trendcoeffsr;
DATASET.trendcoeffsrd = trendcoeffsrd;
DATASET.trendcoeffs20 = trendcoeffs20;
DATASET.trendcoeffs1951 = trendcoeffs1951;
DATASET.trend = trend;
DATASET.trendd = trendd;
DATASET.trendr = trendr;
DATASET.trendrd = trendrd;
DATASET.trend20 = trend20;
DATASET.trend1951 = trend1951;


%COBE.GMMT=GMMT;
NAME_EXP = DATASET;
