clc;
clear;
w = warning ('off','all');
% NIFG=19;
% NPS=7043;
load('Bperp.mat');
load('Btemp.mat');
load('acq_times.mat');
% load('X_rg.mat');
% load('Y_az.mat');
ps_ij = csvread('./PSC.csv');
X = ps_ij(:,2)+1;
Y = ps_ij(:,1)+1;
% load ../lv/xy.mat;
% X = xy(:,1);
% Y = xy(:,2);
Btemp=Btemp./365.0;
clear ps_ij;
load('PHASE_updated.mat');
load('F_updated.mat');
NIFG = numel(Bperp);
NPS = numel(X);

rand('state',30000); % To produce reproducible result
% --- Sparsify the points to obtain reference network points --
NO_SPARSIFICATION=0;
if (NO_SPARSIFICATION==1)
  warning('Not performing sparsification by user request')
  IDX_REF = 1:length(X);
else
  F       = rand(size(X));% e.g., inverse amplitude dispersion index
  F=1./F;
  grid_x = 25*2; % 20 m x 4 m for ERS-1 (slant res. x azimuth res.)
  grid_y = 125*2;
  IDX1    = sparsify(X,Y,[grid_x,grid_y], F,0);%
  IDX2    = sparsify(X(IDX1),Y(IDX1),[grid_x,grid_y], F(IDX1),1);%
  IDX3    = sparsify(X(IDX1(IDX2)),Y(IDX1(IDX2)),[grid_x,grid_y], F(IDX1(IDX2)),2);%
  IDX_REF = IDX1(IDX2(IDX3));% index in X,Y
end
X_REF    = X(IDX_REF);
Y_REF    = Y(IDX_REF);

% --- Create a network using Delaunay --------------------------
% --- Create arcs index vectors that specify this network ------
TRI      = delaunay(X_REF,Y_REF);
all_from = [TRI(:,1); TRI(:,2); TRI(:,3)].';
all_to   = [TRI(:,2); TRI(:,3); TRI(:,1)].';
% --- Make a list from-to for all arcs, sort it ---
% --- Make sure from is always less than to -------
qfrom    = min(all_from,all_to);% index in X_REF
qto      = max(all_from,all_to);% index in X_REF
ft       = [qfrom.', qto.'];
ft       = unique(ft,'rows');% remove double arcs
IDX_from = ft(:,1);% index in ?_REF (sorted)
IDX_to   = ft(:,2);% index in ?_REF
clear ft qfrom qto TRI;

% --- Remove long arcs -----------------------------
all_arclengths     = sqrt((X_REF(IDX_from)-X_REF(IDX_to)).^2 + ...
                          (Y_REF(IDX_from)-Y_REF(IDX_to)).^2);
max_length         = 5000;
too_long           = find(all_arclengths > max_length);
IDX_from(too_long) = [];% remove these arcs
IDX_to(too_long)   = [];%  +assume all points remain
all_arclengths     = sqrt((X_REF(IDX_from)-X_REF(IDX_to)).^2 + ...
                          (Y_REF(IDX_from)-Y_REF(IDX_to)).^2);

% % --- Store the "true" parameters at the arcs ------
% % --- of the reference network ---------------------
% allDEM_arc_true   = simtopo(IDX_REF(IDX_to))-simtopo(IDX_REF(IDX_from));
% allDEFO_arc_true  = simdefo(IDX_REF(IDX_to))-simdefo(IDX_REF(IDX_from));

% --- Report stats ---------------------------------
NARC = length(IDX_from);
NREF = length(unique([IDX_from;IDX_to]));%
disp(' ');
disp('Overview of reference network:');
disp('------------------------------');
disp(['  Number of points in network:            ', num2str(NREF)]);
disp(['  Number of arcs:                         ', num2str(NARC)]);
disp(['  Average number of arcs per point:       ', num2str(2*NARC/NREF)]);% 1 arc is 2ps
disp(['  Minimum arc length:                     ', num2str(min(all_arclengths))]);
disp(['  Maximum arc length:                     ', num2str(max(all_arclengths))]);
disp(['  Mean arc length:                        ', num2str(mean(all_arclengths))]);
disp(['  Standard deviation arc length:          ', num2str(std(all_arclengths))]);
disp(' ');
disp('Please <Press a key> to continue');

% --- Ask user if VCE is desired --------------------------------------
NO_VCE = input('Enter "1" to *not* perform VCE [0]:     ');
if (isempty(NO_VCE)) NO_VCE=0; end
if (NO_VCE==1)
  warning('only for quicker testing')
end

% --- Create a design matrix for "integration" ------
% --- do not reduce it for reference point yet ------
LS_DESIGN_C = zeros(NARC, NREF);% chapter 4.4.2
for arc=1:NARC
  LS_DESIGN_C(arc,IDX_from(arc)) = -1;
  LS_DESIGN_C(arc,IDX_to(arc))   =  1;
end
NARC_per_point = sum(abs(LS_DESIGN_C));

% --- For VCE, find a set of arcs such that each point is used only once ---
disp(' ');
disp('Obtaining a set of arcs for VCE that do not use common points');
tmp = LS_DESIGN_C;
for p=1:NREF
  used = find(tmp(:,p));% e.g., in rows [1,2]
  tmp(used(2:length(used)),:)=[];% remove all rows that also have this point
end
NARC_vce     = size(tmp,1);
IDX_from_vce = zeros(1,NARC_vce);
IDX_to_vce   = zeros(1,NARC_vce);
for arc=1:NARC_vce
  IDX_from_vce(arc) = find(tmp(arc,:)==-1);
  IDX_to_vce(arc)   = find(tmp(arc,:)==1);
end
clear tmp;% not needed anymore

% --- Compute distances of each arc for VCE ------------------------
arclengths = sqrt((X_REF(IDX_from_vce)-X_REF(IDX_to_vce)).^2 + ...
                  (Y_REF(IDX_from_vce)-Y_REF(IDX_to_vce)).^2);
disp(['  Number of arcs for VCE:                 ', num2str(NARC_vce)]);
disp(['  Mean arc length for VCE:                ', num2str(mean(arclengths))]);
disp(['  Standard deviation arc length for VCE:  ', num2str(std(arclengths))]);

% --- Plot acquisitions/network ------------------------------------
figure(1)
  %%% Panel 1: baseline distribution
  subplot(1,3,1)
    plot([Bperp],acq_times,'r+');% add master
    title('Baseline distribution');
    xlabel('perpendicular baseline');
    ylabel('acquisition time');
  subplot(1,3,2)
  %%% Panel 3: spatial distribution of points, network
    plotarc(X_REF, Y_REF, IDX_from, IDX_to);
    hold on
    plot(X, Y, 'k.', 'MarkerSize',2);
    hold off
    axis('tight')
    title('Reference Network');
    xlabel('Range');
    ylabel('Azimuth');
  %%% Panel 4: arcs used for VCE
  subplot(1,3,3)
    plotarc(X_REF, Y_REF, IDX_from_vce, IDX_to_vce);
    title('Arcs for Variance Component Estimation');
    xlabel('Range');
    ylabel('Azimuth');
disp(' ');
disp('Figure 1 show the data');

% ------------------------------------------------------------------
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
B          = [h2p, v2p];% design matrix for float parameters
NPM        = size(B,2);% number of float parameters (and pseudo-obs)

% --- Set up stochastic model for ILS ------------------------------
% --- A priori vc-matrix of the dd observations --------------------
disp('Qy: set to a priori vc-matrix of the DD-observations');
disp('See section 4.3 of the book');
[Qy, Fk_init] = psivcmtx(NIFG);%

% --- Regularization of the system of equations --------------------
% --- Add pseudo-observations and their variances ------------------
disp(' ');
disp('Setting up the model of observation equations:');
disp('   y1 = A1*a + B1*b + e1;  Qy1=D{e1},              (Eq. 1)');
disp('where:');
disp('  y1 [N x 1] vector with N wrapped phase observations.');
disp('  A1 [N x N] design matrix for integer parameters.');
disp('  a  [N x 1] vector of N unknown integer parameters (the ambiguities).');
disp('  B1 [N x 1] design matrix for float parameters.');
disp('  b  [2 x 1] vector of unknown float parameters.');
disp('  e1 [N x 1] vector of measurement noise.');
disp(' Qy1 [N x N] vc-matrix of the noise (D{e} is dispersion of e).');
disp(' ');
disp('This system of equations (Eq. 1) is under-determined,');
disp('i.e., there are more unknown parameters than observations.');
disp('Therefor, (Eq. 1) is regularized using 2 pseudo-observations');
disp('for the float parameters, with value 0 and certain standard deviation.');
disp('  [y1] = [A1] * a + [B1] * b + [e1];  Qy=[D{e1}  0 ],  (Eq. 2)');
disp('  [y2]   [A2]       [B2]       [e2]      [ 0  D{e2}]');
disp(' ');

% --- Add pseudo-observation -------------------------------
%y1     = y;
%y2     = 0;% zero pseudo-observation
y2      = zeros(NPM,1);% zero pseudo-observation
%y      = [y1; y2];

% --- Set up design matrix for integer parameters ----------
A1      = -2*pi*eye(NIFG);
A2      = zeros(NPM,NIFG);
A       = [A1; A2];

% --- Set up design matrix for float parameters ------------
B1      = B;%
B2      = eye(NPM);
B       = [B1; B2];

maxDEM = 25;
stdDEFO = 10;
% --- Set up vc-matrix -------------------------------------
var_pseudo_obs = [0.5*maxDEM, 2*stdDEFO].^2;% set a priori "soft-bounds"
Qy1     = Qy;
Qy2     = diag(var_pseudo_obs);
Qy      = [Qy1,     zeros(NIFG,NPM); ...
           zeros(NPM,NIFG),    Qy2];
       
       % --- Perform the estimation using the LAMBDA method -------
% --- Decorrelation of system of equations for faster and --
% --- more robust estimation -------------------------------
% --- See chapter 3 ----------------------------------------
P_B     = eye(NIFG+NPM) - B * inv(B.'*inv(Qy)*B) * B.'*inv(Qy);
A_      = P_B*A;% reduced design matrix
Qahat   = inv(A_.'*inv(Qy)*A_);%

% --- Decorrelation using Z-transform, function "zt" -------
% --- Speed-up, compute Z-transform outside loop (independent of observations)
[Z,L,D] = zt(Qahat);
Linv    = inv(L);% pre-computed outside of loop
Dinv    = 1./D;% pre-computed outside of loop
invZt   = inv(Z.');% for fixed solution a_check
invLtDL = inv(L.'*diag(D)*L);% pre-computed outside of loop for "ebs"

% --- 7: Estimate float parameter using unwrapped data -----
% --- Least-squares solution of E{y}=Ax; D{y}=Qy -----------
% --- is given by xhat=inv(A.'*inv(Qy)*A)*A.'*inv(Qy)*y ----
PROJ_LS     = inv(B1.'*inv(Qy1)*B1)*B1.'*inv(Qy1);% for speed


% --- Wait for user ----------------------------------------
% --- Then perform a VCE -------------------------------------------
% --- Use these arcs to estimate using a priori stochastic model ---
disp('Set up of matrices completed');
disp(' ');
disp(' ');
disp('Start of estimate the variance components.');
disp(['  First estimate parameters at ', num2str(NARC_vce), ' independent arcs']);
disp('  Using the Extended Bootstrap using the a priori vc-matrix (wrapped data).');


% --- Form the double-differences at the arcs ------------------------
allEST_init = zeros(NPM, NARC_vce);
allsqnorm   = zeros(1,NARC_vce);
allRES      = zeros(NIFG,NARC_vce);
for arc=1:NARC_vce
  % --- 0: Double-difference phase observations at arc ---------------
  y1   = wrap(PHASE(:,IDX_to_vce(arc))-PHASE(:,IDX_from_vce(arc)));
  %y    = [y1;y2];% add pseudo-observations
  % --- 1: float solution for ambiguities ----------------------------
  afloat = y1./(-2.*pi);% float solution (special case for PSI)
  zfloat = Z.' * afloat;% decorrelate float solution for faster ILS search
  % --- 2: Obtain bound for ILS search, function "ebs" ---------------
  [zfixed,sqnorm] = ebs(zfloat,L,D,1,invLtDL);% 1 candidate, reuse invLtDL
  % --- 4: Fixed solution a_check, Inverse Z-transform ---------------
  afixed = invZt*zfixed;
  % --- 5: Unwrap data using the estimated ambiguities ---------------
  y_uw   = y1 + 2.*pi.*afixed;
  % --- 6: Ordinary weighted least-squares using unwrapped data ------
  bhat   = PROJ_LS*y_uw;% estimated DEM error, defo differences
  yhat   = B1*bhat;
  ehat   = wrap(y_uw-yhat);% residuals are always wrapped, closest solution
  allRES(:,arc)      = ehat;
  allsqnorm(arc)     = sqnorm;
  allEST_init(:,arc) = bhat;% store to get feeling for parameter values
end
disp(' ...done');
disp(' ');

% --- After the initial estimation at these arcs, estimate the -------
% --- variance components and use them to estimate the parameters ----
% --- at all arcs of the reference network ---------------------------
% --- Estimate the variance components -------------------------------
disp('Estimate the variance components (can take some time)');
% disp('Please <Press a key>');
% if (do_profile~=1) pause; end;
IDX_ok      = find(allsqnorm<mean(allsqnorm)+std(allsqnorm));
disp(['Using ', num2str(length(IDX_ok)), ' arcs of ', num2str(NARC_vce), ' for estimation of Fk']);
if (NO_VCE==1)
  warning('Not performing VCE by user request');
  Fk = Fk_true;
else
  [Fk, varFk] = psivce(Qy1,B1,allRES(:,IDX_ok));
end
disp(' ...done');



% --- ILS Estimate parameters with this stochastic model --------------
disp(' ');
disp('Estimate the parameters at arcs of network using the a posteriori vc-matrix.');

% --- Update vc-matrix -------------------------------------------------
% --- Set bounds on parameters using preliminary estimates -------------
% --- Qy1: a posteriori vc-matrix for double-diff observations ---------
if (NO_VCE~=1)
  Qy1     = psivcmtx(Fk);%<-- evaluate estimated components!
end

% --- Qy2: vc-matrix for pseudo-observations ---------------------------
var_pseudo_obs = (1.5*std(allEST_init(:,IDX_ok).')).^2;%
disp(['Std. of pseudo-observations for DEM error, displ. rate: ', num2str(sqrt(var_pseudo_obs))]);
Qy2     = diag(var_pseudo_obs);%<-- derived from preliminary estimates
Qy      = [Qy1,     zeros(NIFG,NPM); ...
           zeros(NPM,NIFG),    Qy2];



% --- Update with new model --------------------------------------------
Qahat   = inv(A_.'*inv(Qy)*A_);% functional model same, stochastic model updated

% --- Decorrelation using Z-transform, function "zt" ---------------------
% --- speed-up: computations outside loop (independent of observations) --
[Z,L,D] = zt(Qahat);
Linv    = inv(L);% pre-computed outside of loop
Dinv    = 1./D;% pre-computed outside of loop
invZt   = inv(Z.');% for fixed solution a_check
invLtDL = inv(L.'*diag(D)*L);% pre-computed outside of loop for "ebs"



% --- Projector to estimate float parameter using unwrapped data -------
Qb_hat     = inv(B1.'*inv(Qy1)*B1);% vc-matrix of estimated parameters
PROJ_LS    = Qb_hat*B1.'*inv(Qy1);% for speed outside loop
disp(['Estimation at ', num2str(NARC), ' arcs of the network']);
disp('Qb_hat, vc-matrix of estimated parameters follows: ');
Qb_hat = Qb_hat
disp(' ');
allESTarc = zeros(NPM, NARC);
allsqnorm = zeros(1,NARC);
allRESarc = zeros(NIFG,NARC);
allVARarc = zeros(1,NARC);
invQy1 = inv(Qy1);
for arc=1:NARC
if (mod(arc-1,50)==0)
disp([' ILS: arc ', num2str(arc), ' of ', num2str(NARC)]);
end
% --- 0: Double-difference phase observations at arc ---------------
y1 = wrap(PHASE(:,IDX_REF(IDX_to(arc))) - ...
PHASE(:,IDX_REF(IDX_from(arc))));
%y = [y1;y2];% add pseudo-observations
% --- 1: float solution for ambiguities ----------------------------
afloat = y1./(-2.*pi);% float solution (special case for PSI)
zfloat = Z.' * afloat;% decorrelate float solution (for faster search)
% --- 2: Obtain bound for ILS search, function "ebs" ---------------
[zfixed_bs, sqnorm_bs] = ebs(zfloat,L,D,1,invLtDL);% 1 candidate
% --- 3: Fixed solution for decorrelated ambiguities using "ils"----
[zfixed_ils, sqnorm_ils, ierr] = ils(zfloat,Linv,Dinv,sqnorm_bs,1);
% --- Check ierr status flag ---------------------------------------
if (bitand(ierr,1))
warning('not enough candiates found in "ils"');
end
if (bitand(ierr,2))
warning('maxloop reached in "ils", using bootstrapped estimate');
sqnorm_ils = sqnorm_bs;
zfixed_ils = zfixed_bs;
end
% --- 4: Fixed solution a_check, Inverse Z-transform ---------------
afixed = invZt*zfixed_ils;
% --- 5: Unwrap data using the estimated ambiguities ---------------
y_uw = y1 + 2.*pi.*afixed;
% --- 6: Ordinary weighted least-squares using unwrapped data ------
bhat = PROJ_LS*y_uw;% estimated DEM error, defo differences
yhat = B1*bhat;
ehat = wrap(y_uw-yhat);% residuals are always wrapped, closest solution
allRESarc(:,arc) = ehat; % Residuals per arc
allsqnorm(arc) = sqnorm_ils;
allESTarc(:,arc) = bhat;% store estimated parameters at arcs
allVARarc(arc) = (ehat.'*invQy1*ehat)./(NIFG-NPM);
end
disp(' ...done');
disp(' ');

IDX_refpnt      = 223;
C               = LS_DESIGN_C;% integration of estimates at arcs
C(:,IDX_refpnt) = [];% remove singularity, DEM_refpnt==0; DEFO_refpnt==0
disp(['Selected reference point with IDX = ', num2str(IDX_refpnt)]);
% --- remove estimated with large sqnorm -------------------------------
% --- This is not a good way of hypothesis testing ---------------------
% --- and no test for points -------------------------------------------
IDX_ok    = find(allsqnorm<mean(allsqnorm)+2*std(allsqnorm)+eps);
y         = allESTarc(:,IDX_ok).';% remove bad arcs
C         = C(IDX_ok,:);% remove bad arcs
IDX_from  = IDX_from(IDX_ok);% reduce for removed arcs
IDX_to    = IDX_to(IDX_ok);% reduce for removed arcs
disp(['Removing ', num2str(NARC-length(IDX_ok)), ' arcs of ', ...
  num2str(NARC),' that had a bad fit before inversion']);
disp('  if many arcs are removed like this, it may lead');
disp('  to points that are not connected with any arc (which should be taken out),')
disp('  or to isolated networks that require additional reference points')
disp('  Here, this is not tested, i.e., the following system of equations');
disp('  may be singular, which will then lead to a crash.');
points_with_no_arcs = find(sum(abs(C))==0);
disp(['  Number of unconnected points:         ', ...
  num2str(length(points_with_no_arcs))]);

mean_VAR= zeros(1,size(C,2));
for i = 1:size(C,2)
    con_arc_ix=find(abs(C(:,i))==1);
    mean_VAR(i)= mean(allVARarc(con_arc_ix));
end
coh_points=find(mean_VAR<=3);
C=C(:,coh_points);


m = NARC; % no of observations/eqns
n = NREF - 1; % no of unknowns
q = m - n; % degrees of freedom (dof)
alpha = 0.05; % 95% probablity
df = q;
critical_val = chi2inv(1 - alpha,df)

OMT = 0;
for i=1:NIFG
tmp = allRESarc(i,:)*inv(diag(allVARarc))*allRESarc(i,:).';
OMT = OMT + tmp;
if OMT > critical_val
fprintf('Null hypothesis is rejected before seeing all interferograms\n');
disp(i)
break;
end
end
% OMT1 =ehat.'*inv(Qy1)*ehat;

if (OMT > critical_val)
fprintf('Null hypothesis is rejected\n');
end

%good_ix = find(OMT <= critical_val);
%all_bad_arc= allRESarc;%bad arcs along good arcs
%all_bad_arc(good_ix,:)=[];%bad arcs
%[num_bad_arc,col]=size(all_bad_arc);%calc num of good arcs
%rejecting H_0 greater than k_alpha
%for i=1:num_bad_arc
%T= (all_bad_arc(:,i).^2)./(allVARarc);%calc T_q
%idx(i)= length(find(T>critical_val));  
%end
critical_val = chi2inv(1 - alpha,1);
%all_bad_arc= allRESarc;
idx=zeros(m,1);
for i=1:m
T= (allRESarc(:,i).^2)./(var(allRESarc(:,i)));
idx(i,1)=length(find(T>critical_val));
end