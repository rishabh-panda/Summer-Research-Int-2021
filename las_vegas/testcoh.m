coh = zeros(5000*200,1);
for i = 1:5000*200
    coh(i) = double(enscoh(phi(i,:)));
    ml(i) = gamma_ml(phi(i,:));
end