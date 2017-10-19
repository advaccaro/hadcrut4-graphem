%nino_plot.m

%Load NINO3.4 datasets
load_nino_datasets

% Abbreviate dataset names
B = BUNGE_NINO;
C = COBE_NINO;
CW = CW_EXP;
E = ERSSTv4_NINO;
HG = H43M80_EXP;
HR = H43MRAW_EXP;
HU = H43MUPD_EXP;
K = Kext_NINO;

%flip datenum time axes (because it's a pain in the rear end)
C.tmdn = flipud(C.tmdn);
CW.tmdn = flipud(CW.tmdn);
E.tmdn = flipud(E.tmdn);
HR.tmdn = flipud(HR.tmdn);
HU.tmdn = flipud(HU.tmdn);
HG.tmdn = flipud(HG.tmdn);

%% Calculate differences to Cowtan & Way
%HadCRUT family
[CWind1,CWind2] = find_overlap(HR.tmfrac,CW.tmfrac);
HG.NINO34m_Cdiff = HG.NINO34m(CWind1) - CW.NINO34m(CWind2);
HR.NINO34m_Cdiff = HR.NINO34m(CWind1) - CW.NINO34m(CWind2);





%% Calculate differences to Kaplan


%Bunge
%clear X
[B.Kind1,B.Kind2] = find_overlap(B.tmfrac,K.tmfrac);
B.NINO34m_Kdiff = B.NINO34m(B.Kind1) - K.NINO34m(B.Kind2);



% ERSSTv4
[E.Kind1,E.Kind2] = find_overlap(E.tmfrac,K.tmfrac);
E.NINO34m_Kdiff = E.NINO34m(E.Kind1) - K.NINO34m(E.Kind2);

%COBE
[C.Kind1,C.Kind2] = find_overlap(C.tmfrac,K.tmfrac);
C.NINO34m_Kdiff = C.NINO34m(C.Kind1) - K.NINO34m(C.Kind2);

%HadCRUT4.3 GraphEM
[HG.Kind1,HG.Kind2] = find_overlap(HG.tmfrac,K.tmfrac);
HG.NINO34m_Kdiff = HG.NINO34m(HG.Kind1) - K.NINO34m(HG.Kind2);





%% Plotting
load JEG_graphics

%NINO3.4 Figure 1 (HadCRUT family w/ differences)
fig('HadCRUT family NINO3.4 w/ differences'); clf;
subplot(2,1,1);
hold on;
plot(HR.tmdn,HR.NINO34m,'color','g');
plot(CW.tmdn,CW.NINO34m,'color','r');
plot(HG.tmdn,HG.NINO34m,'color',skyblue);
datetick('x','yyyy')
[hleg,objh] = legend('HadCRUT4.3 median raw', 'HadCRUT4.3 Cowtan & Way', 'HadCRUT4.3 + GraphEM (0.8% sparsity)');
set(hleg,'Location','SouthEast')
set(objh,'linewidth',2)
legend('boxoff')
ylabel_str = ['Temperature (' char(176) 'C)'];
fancyplot_deco('Monthly NINO3.4 index for various temperature products', 'Time (Year)', ylabel_str,14,'Helvetica')
subplot(2,1,2); hold on;
plot(HR.tmdn(CWind1),HR.NINO34m_Cdiff,'color','g');
plot(HG.tmdn(CWind1),HG.NINO34m_Cdiff,'color',skyblue);
datetick('x','yyyy')
[hleg,objh] = legend('HadCRUT4.3 median raw', 'HadCRUT4.3 + GraphEM (0.8% sparsity)');
set(hleg,'Location','SouthEast')
set(objh,'linewidth',2)
legend('boxoff')
fancyplot_deco('Same, difference from HadCRUT 4.3 Cowtan & Way', 'Time (Year)', '{\Delta}T ({\circ}C)',14,'Helvetica')

odir = '/home/scec-02/avaccaro/HadCRUT4.3/figs/';
ninofig1 = 'had43med_nino_plot1.eps';
nino1path = [odir ninofig1];
hepta_figprint(nino1path)


%NINO3.4 Figure 2 (all the pals)
fig('NINO3.4 many datasets')
subplot(2,1,1); hold on;
plot(K.tmdn,K.NINO34m,'color',dark_green);
plot(B.tmdn,B.NINO34m,'color',dark_violet);
plot(E.tmdn,E.NINO34m,'color',firebrick);
plot(C.tmdn,C.NINO34m,'color',pink);
plot(HG.tmdn,HG.NINO34m,'color',skyblue);
datetick('x','yyyy')
[hleg,objh] = legend('Kaplan','Bunge & Clarke','ERSSTv4','COBE','HadCRUT4.3 + GraphEM (0.8% sparsity)');
set(hleg,'Location','SouthEast')
set(objh,'linewidth',2)
pause;
legend('boxoff')
fancyplot_deco('Monthly NINO3.4 index for various temperature products', 'Time (Year)', ylabel_str,14,'Helvetica')
subplot(2,1,2); hold on;
plot(B.tmdn(B.Kind1),B.NINO34m_Kdiff,'color',dark_violet);
plot(E.tmdn(E.Kind1),E.NINO34m_Kdiff,'color',firebrick);
plot(C.tmdn(C.Kind1),C.NINO34m_Kdiff,'color',pink);
plot(HG.tmdn(HG.Kind1),HG.NINO34m_Kdiff,'color',skyblue);
datetick('x','yyyy');
[hleg,objh] = legend('Bunge & Clarke', 'ERSSTv4', 'COBE', 'HadCRUT4.3 + GraphEM (0.8% sparsity)');
set(hleg,'Location','SouthEast')
set(objh,'linewidth',2)
legend('boxoff')
fancyplot_deco('Same, difference to Kaplan', 'Time (Year)', '{\Delta}T ({\circ}C)',14,'Helvetica')

ninofig2 = 'had43med_nino_plot2.eps';
nino2path = [odir ninofig2];
hepta_figprint(nino2path)






