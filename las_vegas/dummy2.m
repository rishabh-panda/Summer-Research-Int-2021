%m_NIFG = 43;
%m_NPM = 2
m_Qbhat= C.'*inv(Qy_est)*C;
m_Qyhat = C*m_Qbhat*C.';
m_Qehat = Qy_est- m_Qyhat;
m_bhat= m_Qbhat*C.'*inv(Qy_est)*y;
m_yhat= C*m_bhat;
m_ehat= y-m_yhat;
%done lis
m_omt = m_ehat.'*inv(Qy_est)*m_ehat;
%m_q = m_NIFG - m_NPM;

m_T1 = m_ehat(:,1).'*inv(Qy_est)*C;
m_T2 = inv(C.'*inv(Qy_est)*m_Qehat*inv(Qy_est)*C);
m_T3 = C.'*inv(Qy_est)*m_ehat(:,1);
m_T = m_T1*m_T2*m_T3;