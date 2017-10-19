%ghcn_open.m

B = dlmread('grid-mntp-1880-current-v3.3.0.dat' ,'\t');
B(find(B == -9999)) = NaN;
