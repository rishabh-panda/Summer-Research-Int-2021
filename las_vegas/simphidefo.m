function [phi_defo,lindefo] = simphidefo(btemp,NPS,stdDEFO)
error(nargchk(1,3,nargin))
if (nargin<2)
  NPS = 1;% [m] default
end
if (nargin<3)
  stdDEFO = 10.0;% [mm/y] default
end
if (max(stdDEFO)<=0)
  phi_defo = 0;
  lindefo  = zeros(NPS,1);
  return;
end
if (length(stdDEFO)~=1 & length(stdDEFO)~=NPS)
  error('stdDEFO should specify DEM for each PS');
end
% --- Get displacement rate per point ---------------
if (length(stdDEFO)==NPS)
  disp('Using given input displacement rates');
  lindefo    = reshape(stdDEFO,1,NPS);% force lying
else
  randn('state',sum(100*clock));% change state of random number generator
  lindefo    = stdDEFO.*randn(1,NPS);% per PS point
end