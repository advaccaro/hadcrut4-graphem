%make_gmt_table.m
%% Initialize
addpath(genpath('/home/scec-02/avaccaro/HadCRUT4.3/'))

load_gmt_datasets

%shorten dataset names
HR = H43MRAW_EXP;
HU = H43MUPD_EXP;
HG = H43M80_EXP;
G = GISTEMP_GMT;
CW = CW_EXP;
N = NOAA_GMT;

nSets = 4;

%preallocate arrays (columns)
name = cell(nSets,1); period = cell(nSets,1); 
rho_m = NaN(nSets,1); rho_msig = NaN(nSets,1); rho_m_cell = cell(nSets,1);
rho_a = NaN(nSets,1); rho_asig = NaN(nSets,1); rho_a_cell = cell(nSets,1);
rho_d = NaN(nSets,1); rho_dsig = NaN(nSets,1); rho_d_cell = cell(nSets,1);
trend = NaN(nSets,1); trend_rec = NaN(nSets,1); 
trend_cell = cell(nSets,1);
trend_rec_cell = cell(nSets,1);

%% HadCRUT4.3 raw
n=1; 
X = HR;
name{n} = '$HadCRUT4.3_{raw}$';
time1 = datevec(min(X.tmdn));
mon1 = time1(2); year1 = time1(1);
timef = datevec(max(X.tmdn));
monf = timef(2); yearf = timef(1);
period{n} = [num2str(mon1) '/' num2str(year1) '--' num2str(monf) '/' num2str(yearf)];
[rho_m(n), rho_msig(n)] = series_corr(X.tmfrac,X.GMTm,HG.tmfrac,HG.GMTm);
[rho_a(n), rho_asig(n)] = series_corr(X.ta_yr,X.GMTa_yr,HG.ta_yr,HG.GMTa_yr);
[rho_d(n), rho_dsig(n)] = series_corr(X.tmfrac,X.GMTd_full,HG.tmfrac,HG.GMTd_full);
trend(n) = X.trendcoeffs(1)*10; trend_rec(n) = X.trendcoeffsr(1)*10;
clear time1 mon1 year1 timef monf yearf X


%% HadCRUT4.3 CW update
n=n+1;
X = CW;
name{n} = '$HadCRUT4.3_{CW}$';
time1 = datevec(min(X.tmdn));
mon1 = time1(2); year1 = time1(1);
timef = datevec(max(X.tmdn));
monf = timef(2); yearf = timef(1);
period{n} = [num2str(mon1) '/' num2str(year1) '--' num2str(monf) '/' num2str(yearf)];
[rho_m(n), rho_msig(n)] = series_corr(X.tmfrac,X.GMTm,HG.tmfrac,HG.GMTm);
[rho_a(n), rho_asig(n)] = series_corr(X.ta_yr,X.GMTa_yr,HG.ta_yr,HG.GMTa_yr);
[rho_d(n), rho_dsig(n)] = series_corr(X.tmfrac,X.GMTd_full,HG.tmfrac,HG.GMTd_full);
trend(n) = X.trendcoeffs(1)*10; trend_rec(n) = X.trendcoeffsr(1)*10;
clear time1 mon1 year1 timef monf yearf X


%% GISTEMP
n=n+1;
X = G;
name{n} = 'GISTEMP';
time1 = datevec(min(X.tmdn));
mon1 = time1(2); year1 = time1(1);
timef = datevec(max(X.tmdn));
monf = timef(2); yearf = timef(1);
period{n} = [num2str(mon1) '/' num2str(year1) '--' num2str(monf) '/' num2str(yearf)];
[rho_m(n), rho_msig(n)] = series_corr(X.tmfrac,X.GMTm,HG.tmfrac,HG.GMTm);
[rho_a(n), rho_asig(n)] = series_corr(X.ta_yr,X.GMTa_yr,HG.ta_yr,HG.GMTa_yr);
[rho_d(n), rho_dsig(n)] = series_corr(X.tmfrac,X.GMTd_full,HG.tmfrac,HG.GMTd_full);
trend(n) = X.trendcoeffs(1)*10; trend_rec(n) = X.trendcoeffsr(1)*10;
clear time1 mon1 year1 timef monf yearf X


%% NOAA
n=n+1;
X = N;
name{n} = 'NOAAGlobalTemp';
time1 = datevec(min(X.tmdn));
mon1 = time1(2); year1 = time1(1);
timef = datevec(max(X.tmdn));
monf = timef(2); yearf = timef(1);
period{n} = [num2str(mon1) '/' num2str(year1) '--' num2str(monf) '/' num2str(yearf)];
[rho_m(n), rho_msig(n)] = series_corr(X.tmfrac,X.GMTm,HG.tmfrac,HG.GMTm);
[rho_a(n), rho_asig(n)] = series_corr(X.ta_yr,X.GMTa_yr,HG.ta_yr,HG.GMTa_yr);
[rho_d(n), rho_dsig(n)] = series_corr(X.tmfrac,X.GMTd_full,HG.tmfrac,HG.GMTd_full);
trend(n) = X.trendcoeffs(1)*10; trend_rec(n) = X.trendcoeffsr(1)*10;
clear time1 mon1 year1 timef monf yearf X


%rho_m = mat2cell(rho_m, [1 1 1]); 
%rho_a = mat2cell(rho_a, [1 1 1]); rho_a_cell
%rho_d = mat2cell(rho_d, [1 1 1]); rho_d_cell
%trend = mat2cell(trend, [1 1 1]); trend_cell
%trend_rec = mat2cell(trend_rec, [1 1 1]); trend_rec_cell

for i = 1:n
	rho_m_cell{i} = ['\textbf{',num2str(rho_m(i),'%+3.2f'),'}'];
	rho_a_cell{i} = ['\textbf{',num2str(rho_a(i),'%+3.2f'),'}'];
	rho_d_cell{i} = ['\textbf{',num2str(rho_d(i),'%+3.2f'),'}'];
	trend_cell{i} = num2str(trend(i),'%3.4f');
	trend_rec_cell{i} = num2str(trend_rec(i),'%3.4f');
end

M=[name period rho_m_cell rho_a_cell rho_d_cell];

%  still a problem with Row Labels
cLabels = {'Dataset', 'Period', '$\rho(GMTm)$', '$\rho(GMTa)$', '$\rho(GMTd)$'};
latextable(M,'Horiz',cLabels,'name','gmt_db_table_input.tex','Hline',[1]);
