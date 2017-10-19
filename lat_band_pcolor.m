%lat_band_pcolor.m
clear all; close all;
%rmpath('/home/s%cec-02/avaccaro/HadCRUT4.3/opendap/loaddap-3.7.3/testsuite/matlab/grid.m')l

load('/home/scec-02/avaccaro/HadCRUT4.3/graphem_sp/data/had43med_graphem_sp80_step2.mat');

lat = loc(:,2); nlat = length(lat); %get latitudes and length
years = tvec(:,1); %get years
minYear = min(years); maxYear = max(years); %find min and max year
yearBounds = [minYear:5:maxYear]; %define boundaries for years
nInts = length(yearBounds) - 1; %number of 5-year intervals
latStep = 20;
latBounds = fliplr([-90:latStep:90]); %define boundaries for latitude bands
nBands = length(latBounds) - 1; %number of latitude bands
%cLim = 3;
cLims = [2, 2.5, 3];
%ngroups = 11; %number of groups for colormap
NGroups = [9,13,17,21];
latSteps = [10,20];
%zonalMeans = zeros(nInts,nBands);
%pLats = zeros(nBands,1); pTime = zeros(nInts,1);
%for j = 1:nInts
%	for k = 1:nBands
%		clear timeInd spaceInd
%		%get indices:
%		timeInd = find(years >= yearBounds(j) & years <= yearBounds(j+1));
%		spaceInd = find(lat <= latBounds(k) & lat >= latBounds(k+1));
%		zonalMeans(j,k) = nmean(nmean(Xf(timeInd,spaceInd)));
%		pLat(k) = mean([latBounds(k:k+1)]);
%	end
%	pTime(j) = mean([yearBounds(j:j+1)]);
%end
%Xf_


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


%zonalMeans = flipud(zonalMeans);
%plotting

pTime = repmat(pTime,1,nBands)';
pLat = repmat(pLat',1,nInts);

xlbl_str = 'Time (year)';
ylbl_str = 'Latitude (degrees)';
cb_str = ['Temperature anomaly (' char(176) 'C)'];


for i = 1:length(latSteps)
latStep = latSteps(i);
for j = 1:length(NGroups)
ngroups = NGroups(j);
for k = 1:length(cLims)
cLim = cLims(k)
cbmap = cbrewer('div','RdBu',ngroups);
title_str = ['5-Year Zonal Means for ' int2str(latStep) char(176) ' Latitude Bands']; 

%pcolor(pTime,pLats,zonalMeans);
%fig('lat bands pcolor'); clf;
%pco = pcolor(pTime,pLat,zonalMeans);
%title(title_str);%title(['Zonal Averages']);
%caxis([-cLim,cLim]);
%cbmap = cbrewer('div','RdBu',ngroups);
%colormap(flipud(cbmap)); 
%hcb = colorbar;
%%title(hcb,['Temperature anomaly (' char(176) 'C)']);
%xlabel(hcb, cb_str);
%xax = get(gca,'xaxis'); %set(xax,'TickValues',[]);
%yax = get(gca,'yaxis'); %set(yax,'TickValues',[]);
%xlabel(xlbl_str); %ylabel(['Latitude (' int2str(latStep) char(176) ' Latitude Bands)']);
%ylabel(ylbl_str);
%outFile = ['./figs/had43med_lat_band_pcolor_' int2str(latStep) 'd.eps'];
%hepta_figprint(outFile);



%zonalMeans = flipud(zonalMeans);
fig('lat bands pcolor interp'); clf;
%pcolor(pLat,pTime,zonalMeans);
pco2=pcolor(pTime,pLat,zonalMeans); 
shading('interp');
%grid(gca,'on');
title(title_str);%title(['Zonal Averages']);
caxis([-cLim,cLim]); colormap(flipud(cbmap));
hcb = colorbar;
xlabel(hcb,cb_str);
xax = get(gca,'xaxis'); %set(xax,'TickValues',[]);
%set(xax,'TickLabels',[1850:25:2000, 2015]);
yax = get(gca,'yaxis'); %set(yax,'TickValues',[]);
%set(yax,'TickLabels',[-90,-70,-50,-30,0,30,50,70,90]);
xlabel(xlbl_str); %ylabel(['Latitude (' int2str(latStep) char(176) ' Latitude Bands)']);
ylabel(ylbl_str);
%outFile2 = ['./figs/had43med_lat_band_pcolor_' int2str(latStep) 'd_interp.eps'];
%outFile3 = ['./figs/had43med_lat_band_pcolor_' int2str(latStep) 'd_interp.jpeg'];
outFile4 = ['./figs/had43med_lat_band_pcolor_interp_' int2str(latStep) 'd_' int2str(cLim*10) 'c_' int2str(ngroups) 'n.jpeg'];
print(outFile4,'-djpeg','-r500');
end
end
end
%yrsm = 10;
%[ntime,nspace] = size(Xf);
%Xf_sm = zeros(ntime,nspace);
%for i = 1:ntime
%for j = 1:nspace
%Xf_sm(:,j) = hepta_smooth(Xf(:,j),1/(12*yrsm));
%end
%zonalMeans_sm = zeros(nBands,nInts);
%pLats = zeros(nBands,1); pTime = zeros(nInts,1);
%for k = 1:nBands
%	for j = 1:nInts
%		clear timeInd spaceInd
%		%get indices:
%		timeInd = find(years >= yearBounds(j) & years <= yearBounds(j+1));
%		spaceInd = find(lat <= latBounds(k) & lat >= latBounds(k+1));
%		zonalMeans_sm(k,j) = nmean(nmean(Xf_sm(timeInd,spaceInd)));
%		%pTime(j) = mean([yearBounds(j:j+1)]);
%	end
%	%pLat(k) = mean([latBounds(k:k+1)]);
%end


%fig('10yr smoothed 5year bins 20d latitude');
%pco3=pcolor(pTime,pLat,zonalMeans_sm); 
%shading('interp');
%grid(gca,'on');
%title(title_str);%title(['Zonal Averages']);
%caxis([-cLim,cLim]); colormap(flipud(cbmap));
%hcb = colorbar;
%xlabel(hcb,cb_str);
%xax = get(gca,'xaxis'); %set(xax,'TickValues',[]);
%%set(xax,'TickLabels',[1850:25:2000, 2015]);
%yax = get(gca,'yaxis'); %set(yax,'TickValues',[]);
%%set(yax,'TickLabels',[-90,-70,-50,-30,0,30,50,70,90]);
%xlabel(xlbl_str); %ylabel(['Latitude (' int2str(latStep) char(176) ' Latitude Bands)']);
%ylabel(ylbl_str);
%outFile3 = ['./figs/had43med_lat_band_pcolor_' int2str(latStep) 'd_interp.jpeg'];
%print(outFile3,'-djpeg','-r500');












%hepta_figprint(outFile2);


%zonalMeans2 = nan(nBands);
%pLat2 = nan(10,34);
%pTime2 = nan(10,34);
%for i = 1:9
%	for j = 1:33
%		zonalMeans2(i,j) = zonalMeans(i,j);
%		pLat2(i,j) = pLat(i,j);
%		pTime2(i,j) = pTime(i,j);
%	end
%end

%fig('3');clf;
%pcolor(pTime2,pLat2,zonalMeans2);
%fig('3');clf;
%subplot(211);pcolor(zonalMeans);
%subplot(212);pcolor(zonalMeans);shading('interp');

%fig('4');clf;
%subplot(211);pcolor(flipud(zonalMeans));
%subplot(212);pcolor(flipud(zonalMeans));shading('interp');

%fig('5');clf;
%pcolor(flipud(zonalMeans)); shading('interp');

%fig('6');clf;
%pcolor(flipud(pTime),flipud(pLat),zonalMeans);shading('interp');
