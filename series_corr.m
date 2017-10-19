%series_corr.m
function [r,signif, p] = series_corr(t1,x1,t2,x2)
% find overlapping time interval and calculates correlations 
% calls on corr_sig.m
% INPUTS: t1/t2 - time axes using same dating scheme
	%     x1/x2 - data series to be compared



min1 = min(t1); max1 = max(t1); min2 = min(t2); max2 = max(t2);

if min1 > min2
	tmin = min1;
elseif min1 < min2
	tmin = min2;
elseif min1 == min2
	tmin = min1;
else
	display('error: no minimum time step')
end

if max1 > max2
	tmax = max2;
elseif max2 > max1
	tmax = max1;
elseif max1 == max2
	tmax = max1;
else
	display('error: no maximum time step')
end

tind1 = find( t1 >= tmin & t1 <= tmax);
tind2 = find( t2 >= tmin & t2 <= tmax);
xt1 = x1(tind1); xt2 = x2(tind2);

nxt1 = length(xt1); nxt2 = length(xt2);
if ~isequal(nxt1,nxt2) == 1
	display('error: inequal lengths')
end

nindf = nan;
n = 1;
for i = 1:nxt1
	if (~isnan(xt1(i)) == 1 & ~isnan(xt2(i)) == 1) == 1
		nindf(n) = i;
		n = n + 1;
	end
end	




[r,signif,p] = corr_sig(xt1(nindf),xt2(nindf));
