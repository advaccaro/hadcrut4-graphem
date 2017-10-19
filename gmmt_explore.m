%gmmt_explore.m

function [NAME_NINO] = gmmt_explore(sst_in,lon_in,lat_in,time_in)

% INPUTS: sst_in - 3d sst data (lonxlatxtime)
% 		  lon_in - vector of longitudes
%		  lat_in - vector of latitudes
%		  time_in - time axis in standard serial time (days since Jan-0-0000)



[nlon,nlat,nmon] = size(sst_in);
%sst_3d = reshape(had4med',[nlon, nlat, nmon]);
sst_3d = sst_in;

% fix grids
lon  = double(lon_in); lon(lon<0) = lon(lon<0) + 360;
lat=double(lat_in); 
%tr = double(cobe.time); % raw time, in units of days since 1891-1-1 00:00:00
%tn = tr + datenum(1891,1,15);  % convert to standard serial date (i.e. add 1800 years)
tn = time_in; %time in standard serial date (days since jan-0-0000)
tv = datevec(tn); 
t = tv(:,1)+tv(:,2)/12-1/12;
sstr  = double(sst_3d); % may need to increase Java Heap memory space to load (Preferences panel)
%keep everything out to 2013
%sst = sstr(:,:,t<2014);
sst = sstr;

%t = t(t<2014);
% Dimensions
nx=length(lon); % 5 by 5 deg
ny=length(lat);
ns=nx*ny;
nt=length(t);  nyrs = floor(nt/12);


% compute climatology
ind = find(t >= 1961 & t < 1991); 
sstp = permute(sst, [3 2 1]);
sstp12 = sstp(1:nyrs*12,:,:);
sst_clim = climatology(sstp(ind,:,:));  
ssta = sstp12 - repmat(sst_clim,[nyrs 1 1]);

%  Find land points and replace by NaN's
%land = isnan(ssta);

% NINO3/4 computation
NINO3_x=find(lon >=210 & lon <=270);
NINO3_y=find(lat>=-5 & lat <=5);
NINO3=nmean(nmean(ssta(:,NINO3_y,NINO3_x),2),3);
NINO34_x=find(lon >=190 & lon <=240);  % (170W to 120 W)
NINO34=nmean(nmean(ssta(:,NINO3_y,NINO34_x),2),3);
NINO34m = NINO34;

%12 month mean
years12 = [t(1): t(1) + nyrs - 1];
for i = 1:nyrs
	NINO34a_yr(i) = nmean(NINO34(12*(i-1)+1:12*i));
end


% flip it to most recent first
SST_K=flipdim(ssta,1);
tk=flipud(t);
NINO3f = flipud(NINO3); NINO34f = flipud(NINO34);

% convert time axis using datenum technology
year=floor(tk); 
month=round((tk-year)*12)+1;
tv=[year , month  , repmat(15,nt,1)];
tn=datenum(tv);
tm=datestr(tn);
%  declare arrays
ta=[max(year):-1:min(year)+1]'; na=length(ta);
SST_Kw = zeros(na,ny,nx);
NINO3a = zeros(na,1);
% Average over DJF period (NH Winter, W)
for k=1:na
	feb=find(tn==datenum([ta(k),2,15,0,0,0]));
	dec=find(tn==datenum([ta(k)-1,12,15,0,0,0]));
	SST_Kw(k,:,:)= nmean(SST_K(feb:dec,:,:),1); 
	%SST_Kw(k,:,:)	
	NINO3ka(k)=nmean(NINO3f(feb:dec));
	NINO34ka(k)=nmean(NINO34f(feb:dec));
end

NINO3a_djf = NINO3ka;
NINO34a_djf = NINO34ka;
		
% Time convention  DJF[n]=Dec[n-1] + Jan[n] + Feb[n] , n from 1857 to 2000.

%%%%%%%%%%%%%%%%%%%%%%%%%
% select ocean points and reshape
tmp=reshape(SST_Kw,na,ns);
ocean=~isnan(tmp(1,:));
SSTo=tmp(:,ocean); 
[nyears,nr]=size(SSTo);
s=[1:nr];


% NINO3 check
NINO3w=nmean(nmean(SST_Kw(:,NINO3_y,NINO3_x),2),3);


% Decadal smoothing
NINO34d_djf = hepta_smooth(NINO34a_djf,1/10);
NINO34d_yr = hepta_smooth(NINO34a_yr,1/10);
NINO34d_full = hepta_smooth(NINO34m,1/120);


% Graphic check
%fig('Kaplan sanity check')
%subplot(3,1,1:2)
%pcolor(s,ta,SSTo)
%subplot(3,1,3)
%plot(t,NINO3), hold on, axis tight
%plot(ta,NINO3w,'g-',ta,NINO3a,'r-','linewidth',[2]),hold off





%  Assign to structure and save
%
DATASET.SSTm=flipud(SST_K); % monthly values
DATASET.SSTa_djf=flipud(SST_Kw);  %  DJF avg
DATASET.SSTo_djf=flipud(SSTo);    %   DJF avg on ocean pts only (no NaNs)
DATASET.lon=lon;
DATASET.lat=lat;
DATASET.ta_yr = years12;
DATASET.ta_djf=ta;
DATASET.tm=flipud(tm(1:length(NINO34m),:));
DATASET.tmdn = flipud(datenum(DATASET.tm));
DATASET.tmfrac = t(1:length(NINO34m));
DATASET.td_djf = ta;
DATASET.td_yr = years12;
DATASET.td_full = t(1:length(NINO34d_full));
DATASET.NINO3=NINO3a;
DATASET.NINO34m=NINO34m;
DATASET.NINO34a_djf=NINO34a_djf;
DATASET.NINO34a_yr = NINO34a_yr;
DATASET.NINO34d_djf = NINO34d_djf;
DATASET.NINO34d_yr = NINO34d_yr;
DATASET.NINO34d_full = NINO34d_full;
DATASET.ocean_pts=ocean;
%COBE.GMMT=GMMT;
NAME_NINO = DATASET;
