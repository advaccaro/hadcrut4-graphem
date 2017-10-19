%find_overlap.m
function [ind1 ind2] = find_overlap(t1,t2);

t1max = max(t1); t1min = min(t1);
t2max = max(t2); t2min = min(t2);


if t1max > t2max
	tmax = t2max;
elseif t2max > t1max
	tmax = t1max;
elseif t1max == t2max
	tmax = t1max;
else
	display('error: no max time')
end

if t1min > t2min
	tmin = t1min;
elseif t2min > t1min
	tmin = t2min;
elseif t1min == t2min
	tmin = t1min;
else
	display('error: no min time');
end

ind1 = find(t1 <= tmax & t1 >= tmin);
ind2 = find(t2 <= tmax & t2 >= tmin);


