%lat_band_contourf.m
clear all; close all;
addpath(genpath('/home/scec-02/avaccaro/HadCRUT4.3/cbrewer/'));

load('/home/scec-02/avaccaro/HadCRUT4.3/graphem_sp/data/had43med_graphem_sp80_step2.mat');

lat = loc(:,2); nlat = length(lat); %get latitudes and length
years = tvec(:,1); %get years
minYear = min(years); maxYear = max(years); %find min and max year
yearBounds = [minYear:5:maxYear]; %define boundaries for years
nInts = length(yearBounds) - 1; %number of 5-year intervals
latStep = 20;
latBounds = fliplr([-90:latStep:90]); %define boundaries for latitude bands
nBands = length(latBounds) - 1; %number of latitude bands

zonalMeans = zeros(nBands,nInts);
pLats = zeros(nBands,1); pTime = zeros(nInts,1);
for k = 1:nBands
	for j = 1:nInts
		clear timeInd spaceInd
		%get indices:
		timeInd = find(years >= yearBounds(j) & years <= yearBounds(j+1));
		spaceInd = find(lat <= latBounds(k) & lat >= latBounds(k+1));
		zonalMeans(k,j) = nmean(nmean(Xf(timeInd,spaceInd)));
		pTime(j) = mean([yearBounds(j:j+1)]);
	end
	pLat(k) = mean([latBounds(k:k+1)]);
end


%plotting

pTime = repmat(pTime,1,nBands)';
pLat = repmat(pLat',1,nInts);

title_str = ['5-Year Zonal Means for ' int2str(latStep) char(176) ' Latitude Bands']; 
xlbl_str = 'Time (year)';
ylbl_str = 'Latitude (degrees)';
cb_str = ['Temperature anomaly (' char(176) 'C)'];

%fig 1
fig('lat bands pcolor'); clf;
pco = pcolor(pTime,pLat,zonalMeans);
title(title_str);%title(['Zonal Averages']);
caxis([-2.5,2.5]);
cbmap = cbrewer('div','RdBu',21);
colormap(flipud(cbmap)); 
hcb = colorbar;
xlabel(hcb, cb_str);
xlabel(xlbl_str);  ylabel(ylbl_str);
outFile = ['./figs/had43med_lat_band_pcolor_' int2str(latStep) 'd.eps'];
%hepta_figprint(outFile);


%fig 2
fig('lat bands pcolor interp'); clf;
pco2=pcolor(pTime,pLat,zonalMeans); 
shading('interp');
title(title_str);
caxis([-2.5,2.5]); colormap(flipud(cbmap));
hcb = colorbar;
xlabel(hcb,cb_str);
xax = get(gca,'xaxis');yax = get(gca,'yaxis'); 
xlabel(xlbl_str); ylabel(ylbl_str);
outFile2 = ['./figs/had43med_lat_band_pcolor_' int2str(latStep) 'd_interp.eps'];
%hepta_figprint(outFile2);

%fig 3 
[ntime,ncells] = size(Xf);
lats = unique(lat); nlats = length(lats);
%Xf_sm = hepta_smooth(Xf,1/(12*5));
latMeans = zeros(ntime,nlats);
cLat = zeros(ntime,nlats);
cTime = zeros(ntime,nlats);
for i = 1:ntime
	for j = 1:nlats
		latInd = find(lat == lats(j));
		latMeans(i,j) = nmean(Xf(i,latInd));
		cLat(i,j) = lats(j);
		cTime(i,j) = tser(i);
	end
end



cLims = [2.5]; Ns = [16]; YRs = [10];
for k = 1:length(YRs);
	group{k}.latMeans_sm = zeros(ntime,nlats);
	yrsm = YRs(k);
	for j = 1:nlats
		group{k}.latMeans_sm(:,j) = hepta_smooth(latMeans(:,j),1/(12*yrsm));
	end
end
for i = 1:length(cLims);
for j = 1:length(Ns);
for k = 1:length(YRs);
	cLim = cLims(i); n = Ns(j); yrsm = YRs(k);
	fig('lat bands contourf'); clf;
	ttl_str = ['Zonal Means with ' int2str(yrsm) '-year smoothing']; 
	[ctf,ln] = contourf(cTime,cLat,group{k}.latMeans_sm,n);
	datetick('x','yyyy');
	shading('interp'); title(ttl_str);
	caxis([-cLim,cLim]); colormap(flipud(cbmap));
	hcb = colorbar;
	xlabel(hcb,cb_str); xlabel(xlbl_str); ylabel(ylbl_str);
	xlim([tser(1) tser(end)]);
	figPath = ['./figs/zonal/had43med_zonal_' int2str(yrsm) 'yr_' int2str(cLim) 'c_' int2str(n) 'n.eps'];
	hepta_figprint(figPath);
end
end
end


for i = 1:length(cLims);
for j = 1:length(Ns);
for k = 1:length(YRs);
	cLim = cLims(i); n = Ns(j); yrsm = YRs(k);
	fig('lat bands pcolor2'); clf;
	ttl_str = ['Zonal Means with ' int2str(yrsm) '-year smoothing']; 
	%[ctf,ln] = contourf(cTime,cLat,group{k}.latMeans_sm,n);
	pco2 = pcolor(cTime,cLat,group{k}.latMeans_sm);	
	datetick('x','yyyy');
	shading('interp'); title(ttl_str);
	caxis([-cLim,cLim]); colormap(flipud(cbmap));
	hcb = colorbar;
	xlabel(hcb,cb_str); xlabel(xlbl_str); ylabel(ylbl_str);
	xlim([tser(1) tser(end)]);
	figPath = ['./figs/zonal/had43med_zonal_' int2str(yrsm) 'yr_' int2str(cLim) 'c_' int2str(n) 'n.eps'];
	%hepta_figprint(figPath);
end
end
end


