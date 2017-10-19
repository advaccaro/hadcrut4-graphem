%had43med_avail_map.m
addpath(genpath('/home/scec-02/avaccaro/HadCRUT4.3/'))


%load data set
load('had43med.mat')

%custom colors
dkgr = [ 0 .3922 0];
grey = [.3922 .3922 .3922];


%bar graph w/ calculations
[ntime,nstations] = size(rawH);

for i = 1:ntime
	navail(i) = sum(~isnan(rawH(i,:)));
end


subplot(3,2,5:6); hold on;
hb = bar(H43med.tfrac,navail);
hb.BarWidth = 1;
hb.EdgeColor = dkgr;
axis([1850 2015 0 2100])
title('HadCRUT4 data availability over time');
xlabel('Time (Year)'); ylabel('# observations')


%%%%%%%%%%%%%%%%%%%%%
ind1850 = find(H43med.tfrac == 1850);
ind1900 = find(H43med.tfrac == 1900);
ind1950 = find(H43med.tfrac == 1950);
ind2000 = find(H43med.tfrac == 2000);

non_nan1850 = ~isnan(rawH(ind1850,:));
non_nan1900 = ~isnan(rawH(ind1900,:));
non_nan1950 = ~isnan(rawH(ind1950,:));
non_nan2000 = ~isnan(rawH(ind2000,:));

nan1850 = isnan(rawH(ind1850,:));
nan1900 = isnan(rawH(ind1900,:));
nan1950 = isnan(rawH(ind1950,:));
nan2000 = isnan(rawH(ind2000,:));

lons = loc(:,1);
lats = loc(:,2);



plons{1} = lons(non_nan1850); plats{1} = lats(non_nan1850);
plons{2} = lons(non_nan1900); plats{2} = lats(non_nan1900);
plons{3} = lons(non_nan1950); plats{3} = lats(non_nan1950);
plons{4} = lons(non_nan2000); plats{4} = lats(non_nan2000);

nlons{1} = lons(nan1850); nlats{1} = lats(nan1850);
nlons{2} = lons(nan1900); nlats{2} = lats(nan1900);
nlons{3} = lons(nan1950); nlats{3} = lats(nan1950);
nlons{4} = lons(nan2000); nlats{4} = lats(nan2000);

titles{1} = 'January 1850'; titles{2} = 'January 1900'; titles{3} = 'January 1950'; titles{4} = 'January 2000';


for i = 1:4
	subplot(3,2,i); hold on;
	m_proj('Robinson','clong',0);
	%m_grid('xtick',[-180:60:180],'tickdir','out','ytick',[-90:30:90], 'color','k', 'fontsize',10,'fontname','Times New Roman');
	m_grid('xtick',[-180:60:180],'xticklabels',[],'tickdir','out','ytick',[-90:30:90], 'yticklabels',[],'color','k', 'fontsize',10);
	%m_grid('xtick','no','ytick','no', 'color','k')% 'fontsize',10,'fontname','Times New Roman');
	%m_grid('color','k');	
	m_coast('color','k');
	%hl=m_line(plons{i},plats{i},'color',dkgr,'marker','^','MarkerFaceColor',dkgr,'MarkerSize',8.5,'LineStyle','none');
	hl1=m_line(plons{i},plats{i},'color',dkgr,'marker','*','MarkerFaceColor',dkgr,'MarkerSize',1.5,'LineStyle','none');
	hl2=m_line(nlons{i},nlats{i},'color',grey,'marker','*','MarkerFaceColor',dkgr,'MarkerSize',1,'LineStyle','none');
	title(titles{i}); ;
end


odir = '/home/scec-02/avaccaro/HadCRUT4.3/figs/';
oname = 'had43med_avail_map.pdf'; opath = [odir oname];
print -painters -dpdf -cmyk -r1000 '/home/scec-02/avaccaro/HadCRUT4.3/figs/had43med_avail_map.pdf'
print -painters -dpdf -cmyk -r1000 '/home/scec-02/avaccaro/HadCRUT4.3/figs/had43med_avail_map.eps'

%for j = 1:3
%    np = length(group(j).nproxy);
%    for i = 1:np
%        nproxy(j,i) = group(j).nproxy(i);
%    end
    %group(j).lon(group(j).lon>180) = group(j).lon(group(j).lon>180) - 360;
%    tm(j,:) = group(1).tm;
%end

%for j = 1:3
%    proxycolor(j,:) = icon{j,1};
%    hl(j)=m_line(group(j).lon,group(j).lat,'color',icon{j,1},'marker',icon{j,2},'MarkerFaceColor',icon{j,1},'MarkerSize',8.5,'LineStyle','none');
%end



%[LEGH,OBJH]=legend(hl(:),icon{:,3});%pause; 
%set(LEGH,'FontName','Times','FontSize',14,'Location','SouthEast');
%set( findobj(OBJH,'Type','line'), 'Markersize', 12)
%legend('boxon')
