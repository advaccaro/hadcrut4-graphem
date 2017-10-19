%fill3_prep.m
%Adam Vaccaro
%prepares quantiles for envelope plotting
%inputs: x axis (time axis/used for for plotting) as row vector,
%	upper and lower boundaries as row vectors
%outputs: vertxsq,y2,Z (used with fill3 e.g.: fill3(vertxsq,y2,Z,'k') )

function [vertxsq,y2,Z] = fill3_prep(xaxis,upperbound,lowerbound)


Xcases = xaxis; %row vector
Ncases = length(Xcases); Sqcases = [Xcases Xcases];
top = upperbound; bot = lowerbound; %row vector
v = zeros(2*Ncases,1);
Sverts = zeros(2*Ncases,1);

for j = 1:Ncases
Sverts(2*j) = Xcases(j);
Sverts(1+(2*(j-1))) = Xcases(j);
v(2*j) = top(j);
v(1+2*(j-1)) = bot(j);
end


y = [bot, fliplr(top)]; y2 = [y, y(1)];


vertx = [Xcases, fliplr(Xcases)];
vertxsq = [vertx, vertx(1)];
Z = zeros(size(vertxsq));






