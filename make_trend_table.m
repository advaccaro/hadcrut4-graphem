%make_trend_table.m

addpath(genpath('/home/scec-02/avaccaro/HadCRUT4.3/'))

load_gmt_datasets

%shorten dataset names
HR = H43MRAW_EXP;
HU = H43MUPD_EXP;
HG = H43M80_EXP;
G = GISTEMP_GMT;
CW = CW_EXP;
N = NOAA_GMT;

nSets = 5;

%preallocate arrays (columns)
name = cell(nSets,1); period = cell(nSets,1); 
trend = NaN(nSets,1); trend_rec = NaN(nSets,1);
trend20 = NaN(nSets,1); trend1951 = NaN(nSets,1);
trend_cell = cell(nSets,1);
trend_rec_cell = cell(nSets,1);
trend20_cell = cell(nSets,1);
trend1951_cell = cell(nSets,1);


n = 0;

%% HadCRUT4.3 GraphEM
n=n+1;
X = HG;
name{n} = '$HadCRUT4.3_{GraphEM}$';
time1 = datevec(min(X.tmdn));
mon1 = time1(2); year1 = time1(1);
timef = datevec(max(X.tmdn));
monf = timef(2); yearf = timef(1);
period{n} = [num2str(mon1) '/' num2str(year1) '--' num2str(monf) '/' num2str(yearf)];
trend(n) = X.trendcoeffs(1)*10; %warming per decade
trend_rec(n) = X.trendcoeffsr(1)*10;
trend20(n) = X.trendcoeffs20(1)*10;
trend1951(n) = X.trendcoeffs1951(1)*10;
clear X time1 timef mon1 monf



%% HadCRUT4.3 raw
n=n+1;
X = HR;
name{n} = '$HadCRUT4.3_{raw}$';
time1 = datevec(min(X.tmdn));
mon1 = time1(2); year1 = time1(1);
timef = datevec(max(X.tmdn));
monf = timef(2); yearf = timef(1);
period{n} = [num2str(mon1) '/' num2str(year1) '--' num2str(monf) '/' num2str(yearf)];
trend(n) = X.trendcoeffs(1)*10; %warming per decade
trend_rec(n) = X.trendcoeffsr(1)*10;
trend20(n) = X.trendcoeffs20(1)*10;
trend1951(n) = X.trendcoeffs1951(1)*10;
clear X time1 timef mon1 monf


%% HadCRUT4.3 Cowtan & Way
n=n+1;
X = CW;
name{n} = '$HadCRUT4.3_{CW}$';
time1 = datevec(min(X.tmdn));
mon1 = time1(2); year1 = time1(1);
timef = datevec(max(X.tmdn));
monf = timef(2); yearf = timef(1);
period{n} = [num2str(mon1) '/' num2str(year1) '--' num2str(monf) '/' num2str(yearf)];
trend(n) = X.trendcoeffs(1)*10; %warming per decade
trend_rec(n) = X.trendcoeffsr(1)*10;
trend20(n) = X.trendcoeffs20(1)*10;
trend1951(n) = X.trendcoeffs1951(1)*10;
clear X time1 timef mon1 monf

%% GISTEMP
n=n+1;
X = G;
name{n} = 'GISTEMP';
time1 = datevec(min(X.tmdn));
mon1 = time1(2); year1 = time1(1);
timef = datevec(max(X.tmdn));
monf = timef(2); yearf = timef(1);
period{n} = [num2str(mon1) '/' num2str(year1) '--' num2str(monf) '/' num2str(yearf)];
trend(n) = X.trendcoeffs(1)*10; %warming per decade
trend_rec(n) = X.trendcoeffsr(1)*10;
trend20(n) = X.trendcoeffs20(1)*10;
trend1951(n) = X.trendcoeffs1951(1)*10;
clear X time1 timef mon1 monf


%% NOAA
n=n+1;
X = N;
name{n} = 'NOAA';
time1 = datevec(min(X.tmdn));
mon1 = time1(2); year1 = time1(1);
timef = datevec(max(X.tmdn));
monf = timef(2); yearf = timef(1);
period{n} = [num2str(mon1) '/' num2str(year1) '--' num2str(monf) '/' num2str(yearf)];
trend(n) = X.trendcoeffs(1)*10; %warming per decade
trend_rec(n) = X.trendcoeffsr(1)*10;
trend20(n) = X.trendcoeffs20(1)*10;
trend1951(n) = X.trendcoeffs1951(1)*10;
clear X time1 timef mon1 monf


for i = 1:n
	trend_cell{i} = num2str(trend(i), '%+3.3f');
	trend_rec_cell{i} = num2str(trend_rec(i), '%+3.3f');
	trend20_cell{i} = num2str(trend20(i), '%+3.3f');
	trend1951_cell{i} = num2str(trend1951(i), '%+3.3f');
end

M = [name period trend_cell trend_rec_cell trend20_cell trend1951_cell];
cLabels = {'Dataset', 'Period', 'm(Period)', 'm(1998-end)', 'm(20th century)', 'm(1951-2012)'};
latextable(M,'Horiz',cLabels,'name','trend_db_table_input.tex','Hline',[1]);





