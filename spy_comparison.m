%spy_comparison.m
C = load('/home/scec-02/avaccaro/HadCRUT4.3/graphem_cr/data/had43med_grapem_cr900_adj.mat');

S = load('/home/scec-02/avaccaro/HadCRUT4.3/graphem_sp/data/had43med_graphem_sp80_adj.mat');

fig('CR900'); clf;
spy(C.adjR);
titlestr = ['Neighborhood graph cutoff radius: 900km'];	
	title(titlestr, 'FontSize', 9)
hepta_figprint('./figs/had43med_graphem_cr900_spy.eps')
clear titlestr

fig('SP80');clf;
spy(S.adjM);
titlestr = ['GLASSO sparsity: 0.8%'];
title(titlestr,'FontSize',9)
hepta_figprint('./figs/had43med_graphem_sp80_spy.eps')
clear titlestr

adjO = C.adjR.*S.adjM;
fig('Overlap');clf;
spy(adjO);
titlestr = ['Overlapping points'];
title(titlestr,'FontSize',9)
hepta_figprint('./figs/had43med_graphem_overlap_spy.eps')
clear titlestr

adjCS = C.adjR - S.adjM;
adjCS(adjCS == -1) = 0;
fig('CR900 - SP80');clf;
spy(adjCS)
titlestr = ['Points unique to neighborhood graph'];
title(titlestr,'FontSize',9)
hepta_figprint('./figs/had43med_graphem_cr900-sp80.eps')
clear titlestr

adjSC = S.adjM - C.adjR;
adjSC(adjSC == -1) = 0;
fig('SP80 - CR900');clf;
spy(adjSC)
titlestr = ['Points unique to GLASSO graph'];
title(titlestr,'FontSize',9)
hepta_figprint('./figs/had43med_graphem_sp80-cr900.eps')
clear titlestr
