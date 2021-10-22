function [phi_topo,DEMerror] = simphitopo(bperp,NPS,maxDEM)
error(nargchk(1,3,nargin))
if (nargin<2)
  NPS = 1;% [m] default
end
if (nargin<3)
  maxDEM = 40.0;% [m] default
end
if (max(maxDEM)<=0)
  phi_topo = 0;
  DEMerror  = zeros(NPS,1);
  return;
end
if (length(maxDEM)~=1 & length(maxDEM)~=NPS)
  error('maxDEM should specify DEM for each PS');
end
% --- Get height per point --------------------------------
if (length(maxDEM)==NPS)
  disp('Using given input DEM errors');
  DEMerror   = reshape(maxDEM,1,NPS);% force lying
else
  rand('state',sum(100*clock));% change state of random number generator
  DEMerror   = maxDEM*rand(1,NPS)-maxDEM./2;
end
% --- Initialize variables --------------------------------
wavelength = 0.056;
slantrange = 850000;
inc_angle  = 23.0*pi/180;% [rad]
KK         = -4*pi/wavelength;
% --- Use same height conversion factor for all points ----
h2p        = KK.*bperp./(slantrange*sin(inc_angle));
phi_topo   = h2p*DEMerror;% interferometric phase
%EOF