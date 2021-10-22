function pdf = pdf_ifg_phase(gamma,phi)

    rho = gamma/(1+gamma);
    beta = abs(rho)*cosd(phi);
    pdf = (1 - rho*rho)*(1 + beta.* acos(-beta)/sqrt(1 - beta.*beta))./(2*pi*(1 - beta.*beta)); 