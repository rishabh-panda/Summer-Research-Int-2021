%% Detect
m =  NIFG; % no of observations/eqns
n = NPM; % no of unknowns
q = m - n; % degrees of freedom (dof)
OMT = diag(allRESarc.'*inv(Qy1)*allRESarc);
% OMT1 =ehat.'*inv(Qy1)*ehat;

alpha = 0.05; % 95% probablity
df = q;
critical_val = chi2inv(1 - alpha,df)
good_ix = find(OMT <= critical_val);
bad_ix = NARC - good_ix;