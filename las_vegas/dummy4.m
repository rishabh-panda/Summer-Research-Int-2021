AllBadArcs = allRESarc;
AllBadArcs(:,good_ix) = [];
[NoGoodArcs,colnum]=size(good_ix);
NbadArcs = NARC - NoGoodArcs;
Qy1_diag = diag(Qy1);
idx = zeros(NbadArcs,1);
alpha = 0.05; % 95% probablity
df = 1;
critical_val = chi2inv(1 - alpha,df);
for j = 1:NbadArcs
    T_ = AllBadArcs(:,j)./Qy1_diag;
    idx(j) = length(find(T_ <= critical_val));
end