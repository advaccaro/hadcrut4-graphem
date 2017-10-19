%had43med_trend_map.m
%m_proj, m_grid, m_coast, m_line, m_pcolor

startyear = 2010;

%% load raw HadCRUT4
R = load('/home/scec-02/avaccaro/HadCRUT4.3/data/had43med.mat')
[Rnx,Rny,Rntime] = size(R.H43med.x);
R.tfrac = nan(Rntime,1);
for i = 1:Rntime
	R.tfrac(i) = 1850 + (i-1)/12;
end
R.grid = double(R.H43med.x);
Rnyears = floor(Rntime/12); R.years = [1850:1:1849+Rnyears]';
R.grida = nan(Rnx,Rny,Rnyears);
R.trec = R.tfrac(R.tfrac >= startyear);
R.recind = ismember(R.tfrac,R.trec);
R.fullco = cell(Rnx,Rny); R.recco = cell(Rnx,Rny);
R.full = nan(Rnx,Rny); R.rec = nan(Rnx,Rny);



%% load GraphEM-infilled HadCRUT4 (GLASSO)
H = load('/home/scec-02/avaccaro/HadCRUT4.3/graphem_sp/data/had43med_graphem_sp80_step2.mat')
Ht = length(H.tfrac);
H.trec = H.tfrac(H.tfrac >= startyear);
H.recind = ismember(H.tfrac,H.trec);
H.grid = reshape(H.Xf', [Rnx Rny Ht]);
Hnyears = floor(Ht/12); H.years = [1850:1:1849+Hnyears]';
H.grida = nan(Rnx,Rny,Hnyears);
H.fullco = cell(Rnx,Rny); H.recco = cell(Rnx,Rny);
H.full = nan(Rnx,Rny); H.rec = nan(Rnx,Rny);

%% load GraphEM-infilled HadCRUT4 (NEIGH)
N = load('/home/scec-02/avaccaro/HadCRUT4.3/graphem_cr/data/had43med_graphem_cr900_step2.mat')
N.tfrac = H.tfrac; Nt = length(N.tfrac);
N.trec = N.tfrac(N.tfrac >= startyear);
N.recind = ismember(N.tfrac,N.trec);
N.grid = reshape(N.Xf', [Rnx Rny Nt]);
Nnyears = floor(Nt/12); N.years = [1850:1:1849+Nnyears]';
N.grida = nan(Rnx,Rny,Nnyears);
N.fullco = cell(Rnx,Rny); N.recco = cell(Rnx,Rny);
N.full = nan(Rnx,Rny); N.rec = nan(Rnx,Rny);

%% load Cowtan and Way gridded dataset
load('/home/scec-02/avaccaro/HadCRUT4.3fv/CW2013/data/cw2013')
lat = double(cw13.lat); lon = double(cw13.lon); %C.lon(C.lon < 0) = C.lon(C.lon < 0) + 360;
[nlon,nlat,ntime] = size(cw13.temp);
C.tfrac = nan(ntime,1);
for i = 1:ntime
	C.tfrac(i) = 1850 + (i-1)/12;
end
C.trec = C.tfrac(C.tfrac >= startyear);
C.recind = ismember(C.tfrac,C.trec);
C.grid = double(cw13.temp);
Cnyears = floor(ntime/12); C.years = [1850:1:1849+Cnyears]';
C.grida = nan(Rnx,Rny,Cnyears);
C.fullco = cell(Rnx,Rny); C.recco = cell(Rnx,Rny);
C.full = nan(Rnx,Rny); C.rec = nan(Rnx,Rny);
for i = 1:Rnx
	for j = 1:Rny
		for k =1:Rnyears
			R.grida(i,j,k) = nmean(R.grid(i,j,1+(k-1)*12:k*12));
		end
		for k = 1:Hnyears
			H.grida(i,j,k) = nmean(H.grid(i,j,1+(k-1)*12:k*12));
		end
		for k = 1:Cnyears
			C.grida(i,j,k) = nmean(C.grid(i,j,1+(k-1)*12:k*12));
		end
		for k = 1:Nnyears
			N.grida(i,j,k) = nmean(N.grid(i,j,1+(k-1)*12:k*12));
		end
	end
end

R.fullaco = cell(Rnx,Rny); R.fulla = nan(Rnx,Rny);
H.fullaco = cell(Rnx,Rny); H.fulla = nan(Rnx,Rny);
C.fullaco = cell(Rnx,Rny); C.fulla = nan(Rnx,Rny);
N.fullaco = cell(Rnx,Rny); N.fulla = nan(Rnx,Rny);


for i = 1:Rnx
	for j = 1:Rny
		Rfulla = reshape(R.grida(i,j,:), [numel(R.grida(i,j,:)),1]);
		rnan = ~isnan(Rfulla);
		R.fullaco{i,j} = polyfit(R.years(rnan),Rfulla(rnan),1);
		R.fulla(i,j) = R.fullaco{i,j}(1);

		Nfulla = reshape(N.grida(i,j,:), [numel(N.grida(i,j,:)),1]);
		N.fullaco{i,j} = polyfit(N.years,Nfulla,1);
		N.fulla(i,j) = N.fullaco{i,j}(1);

		Hfulla = reshape(H.grida(i,j,:), [numel(H.grida(i,j,:)),1]);
		H.fullaco{i,j} = polyfit(H.years,Hfulla,1);
		H.fulla(i,j) = H.fullaco{i,j}(1);

		Cfulla = reshape(C.grida(i,j,:), [numel(C.grida(i,j,:)),1]);
		C.fullaco{i,j} = polyfit(C.years,Cfulla,1);
		C.fulla(i,j) = C.fullaco{i,j}(1);
	end
end




for i = 1:Rnx
	for j = 1:Rny
		Cfull = reshape(C.grid(i,j,:), [numel(C.grid(i,j,:)),1]);
		C.fullco{i,j} = polyfit(C.tfrac,Cfull,1);		
		C.full(i,j) = C.fullco{i,j}(1);
		Crec = reshape(C.grid(i,j,C.recind), [numel(C.grid(i,j,C.recind)),1]);
		C.recco{i,j} = polyfit(C.trec,Crec,1);
		C.rec(i,j) = C.recco{i,j}(1);

		Nfull = reshape(N.grid(i,j,:), [numel(N.grid(i,j,:)),1]);		
		N.fullco{i,j} = polyfit(N.tfrac,Nfull,1);
		N.full(i,j) = N.fullco{i,j}(1);
		Nrec = reshape(N.grid(i,j,N.recind), [numel(N.grid(i,j,N.recind)),1]);
		N.recco{i,j} = polyfit(N.trec,Nrec,1);
		N.rec(i,j) = N.recco{i,j}(1);

		Rfull = reshape(R.grid(i,j,:), [numel(R.grid(i,j,:)),1]);
		nanind = ~isnan(Rfull);
		R.fullco{i,j} = polyfit(R.tfrac(nanind),Rfull(nanind),1);
		R.full(i,j) = R.fullco{i,j}(1);
		Rrec = reshape(R.grid(i,j,R.recind), [numel(R.grid(i,j,R.recind)),1]);
		nanind2 = ~isnan(Rrec);
		R.recco{i,j} = polyfit(R.trec(nanind2),Rrec(nanind2),1);
		R.rec(i,j) = R.recco{i,j}(1);

		Hfull = reshape(H.grid(i,j,:), [numel(H.grid(i,j,:)),1]);
		H.fullco{i,j} = polyfit(H.tfrac,Hfull,1);
		H.full(i,j) = H.fullco{i,j}(1);
		Hrec = reshape(H.grid(i,j,H.recind), [numel(H.grid(i,j,H.recind)),1]);
		H.recco{i,j} = polyfit(H.trec,Hrec,1);
		H.rec(i,j) = H.recco{i,j}(1);

		clear Cfull Nfull Rfull Hfull Crec Nrec Rrec Hrec
	end
end


%Predefine arrays for easy plotting
sfull = cell(4,1);
sfull{1} = R.full'*100;
sfull{2} = C.full'*100;
sfull{3} = N.full'*100;
sfull{4} = H.full'*100;

srec = cell(4,1);
srec{1} = R.rec';
srec{2} = C.rec';
srec{3} = N.rec';
srec{4} = H.rec';




tfull = cell(4,1);
tfull{1} = 'HadCRUT4.3 median raw';
tfull{2} = 'HadCRUT4.3 median Cowtan & Way';
tfull{3} = 'GraphEM (Neighborhood graph)';
tfull{4} = 'GraphEM (GLASSO)';

%% Full time period



%% Plotting
fig('entire period linear trend map'); clf;

ha = tight_subplot(3,2,[.05 .00001], .05, .01);

for ii = 1:4
	axes(ha(ii));
	m_proj('Robinson', 'clong', 0);
	mpc = m_pcolor(lon,lat,sfull{ii}); %caxis([-.25 .25]);
	set(mpc,'EdgeAlpha',0);
	coast = m_coast('color','k');
	m_grid('box','on','xtick',9,'ytick',6,'xlabeldir','end','xticklabels',[],'yticklabels',[]);
	%title(titles{ii});
	title(tfull{ii});
	

end

axes(ha(5));
c = colorbar3('horiz',['^{\circ}C']);
cpos2 = [.15 .16 .7 .03];
c.Position = cpos2;

%figname1 = ['./figs/linear_trend_full.eps'];
figname1 = ['./figs/linear_trend_1850_to_2015.eps'];
hepta_figprint(figname1)


%fig('raw linear trend map 4 entire period'); clf;
%m_proj('Robinson','clong',0);
%mpc=m_pcolor(lon,lat,sfull{4});
%set(mpc,'EdgeAlpha',0);
%coast=m_coast('color','k');
%m_grid('box','on','xtick',9,'ytick',6,'xlabeldir','end','xticklabels',[],'yticklabels',[]);
%title(tfull{4})



%% Recent period


%% Plotting
fig('recent period linear trend map'); clf;

ha = tight_subplot(3,2,[.05 .00001], .05, .01);

for ii = 1:4
	axes(ha(ii));
	m_proj('Robinson', 'clong', 0);
	mpc = m_pcolor(lon,lat,srec{ii}); %caxis([-2 3])
	set(mpc,'EdgeAlpha',0);
	coast = m_coast('color','k');
	m_grid('box','on','xtick',9,'ytick',6,'xlabeldir','end','xticklabels',[],'yticklabels',[]);
	title(tfull{ii});

end

axes(ha(5));
c = colorbar3('horiz',['^{\circ}C']);
cpos2 = [.15 .16 .7 .03];
c.Position = cpos2;

%figname2 = ['./figs/linear_trend_' num2str(startyear) '.eps'];
figname2 = ['./figs/linear_trend_' num2str(startyear) '_to_2015.eps'];
hepta_figprint(figname2);


%% testing
R.ttest = R.tfrac( R.tfrac >= 1900 & R.tfrac <= 2009);
R.testind = ismember(R.tfrac,R.ttest);
R.testco = cell(Rnx,Rny);
R.test = nan(Rnx,Rny);

C.ttest = C.tfrac(C.tfrac>=1900 & C.tfrac <= 2009);
C.testind = ismember(C.tfrac,C.ttest);
C.testco = cell(Rnx,Rny);
C.test = nan(Rnx,Rny);

N.ttest = N.tfrac(N.tfrac>=1900 & N.tfrac <= 2009);
N.testind = ismember(N.tfrac,N.ttest);
N.testco = cell(Rnx,Rny);
N.test = nan(Rnx,Rny);

H.ttest = H.tfrac(H.tfrac>=1900 & H.tfrac <= 2009);
H.testind = ismember(H.tfrac,H.ttest);
H.testco = cell(Rnx,Rny);
H.test = nan(Rnx,Rny);


for i = 1:Rnx
	for j = 1:Rny
		Rtest = reshape(R.grid(i,j,R.testind), [numel(R.grid(i,j,R.testind)),1]);
		nanind3 = ~isnan(Rtest);
		if sum(nanind3) == 0
			R.test(i,j) = NaN;
		else
			Rmdl = fitlm(R.ttest(nanind3),Rtest(nanind3));
			R.test(i,j) = Rmdl.Coefficients.Estimate(2);	
		end	
		%R.testco{i,j} = polyfit(R.ttest(nanind3),Rtest(nanind3),1);
		%R.test(i,j) = R.testco{i,j}(1);
		
		Ctest = reshape(C.grid(i,j,C.testind), [numel(C.grid(i,j,C.testind)),1]);
		%nanind3 = ~isnan(Ctest);
		Cmdl = fitlm(C.ttest,Ctest);
		C.test(i,j) = Cmdl.Coefficients.Estimate(2);		
		%C.testco{i,j} = polyfit(C.ttest,Ctest,1);
		%C.test(i,j) = C.testco{i,j}(1);
		
		Ntest = reshape(N.grid(i,j,N.testind), [numel(N.grid(i,j,N.testind)),1]);
		%N.testco{i,j} = polyfit(N.ttest,Ntest,1);
		%N.test(i,j) = N.testco{i,j}(1);
		Nmdl = fitlm(N.ttest,Ntest);
		N.test(i,j) = Nmdl.Coefficients.Estimate(2);		


		Htest = reshape(H.grid(i,j,H.testind), [numel(H.grid(i,j,H.testind)),1]);
		%H.testco{i,j} = polyfit(H.ttest,Htest,1);
		%H.test(i,j) = H.testco{i,j}(1);
		Hmdl = fitlm(H.ttest,Htest);
		H.test(i,j) = Hmdl.Coefficients.Estimate(2);
		
	end
end

stest = cell(4,1);
stest{1} = R.test';
stest{2} = C.test';
stest{3} = N.test';
stest{4} = H.test';
		
fig('test period linear trend map'); clf;

ha = tight_subplot(3,2,[.05 .00001], .05, .01);

for ii = 1:4
	axes(ha(ii));
	m_proj('Robinson', 'clong', 0);
	mpc = m_pcolor(lon,lat,stest{ii}); %caxis([-2 3]);
	set(mpc,'EdgeAlpha',0);
	coast = m_coast('color','k');
	m_grid('box','on','xtick',9,'ytick',6,'xlabeldir','end','xticklabels',[],'yticklabels',[]);
	title(tfull{ii});

end


axes(ha(5));
c = colorbar3('horiz',['^{\circ}C']);
cpos2 = [.15 .16 .7 .03];
c.Position = cpos2;

figname3 = ['./figs/linear_trend_1900_to_2009.eps'];
hepta_figprint(figname3);





%% Plotting
fig('entire period annualized linear trend map'); clf;

ha = tight_subplot(3,2,[.05 .00001], .05, .01);
sfulla = cell(4,1);
sfulla{1} = R.fulla';
sfulla{2} = C.fulla';
sfulla{3} = N.fulla';
sfulla{4} = H.fulla';
for ii = 1:4
	axes(ha(ii));
	m_proj('Robinson', 'clong', 0);
	mpc = m_pcolor(lon,lat,sfulla{ii}); %caxis([-.25 .25]);
	set(mpc,'EdgeAlpha',0);
	coast = m_coast('color','k');
	m_grid('box','on','xtick',9,'ytick',6,'xlabeldir','end','xticklabels',[],'yticklabels',[]);
	%title(titles{ii});
	title(tfull{ii});
	

end

axes(ha(5));
c = colorbar3('horiz',['^{\circ}C']);
cpos2 = [.15 .16 .7 .03];
c.Position = cpos2;

%figname1 = ['./figs/linear_trend_fulla.eps'];
figname4 = ['./figs/linear_trend_1850_to_2015_annualized.eps'];
hepta_figprint(figname4)

