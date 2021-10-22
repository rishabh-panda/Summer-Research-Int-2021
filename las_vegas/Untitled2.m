% --- Set up functional model for ILS ------------------------------
% ------------------------------------------------------------------
disp('Setting up design matrix for DEM error and lin.defo');
% --- Use same height conversion factor for all points ----
wavelength = 0.056;% [m]
slantrange = 850000;% [m]
inc_angle  = 23.0*pi/180;% [rad]
KK         = -4*pi/wavelength;
h2p        = KK.*Bperp./(slantrange*sin(inc_angle));% [1/m]
v2p        = KK*Btemp*1e-3;% [y/mm]
B          = [h2p, v2p, repmat(KK, [NIFG,1])];% design matrix for float parameters
NPM        = size(B,2);% number of float parameters (and pseudo-obs)