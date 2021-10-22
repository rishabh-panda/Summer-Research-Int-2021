width = (dlmread('width.txt'));
height = (dlmread('len.txt'));
% corr = zeros(height,width);
% for i=1:NPS
%     tcoh = enscoh(residualph(:,i));
%     corr(Y(i),X(i)) = tcoh;
% end
% 
% fid=fopen('tcorr.bin','w');
% fwrite(fid,corr.','float');
% fclose(fid);
% 
% costscale=100;
% nshortcycle=200;
% maxshort=32000;

fprintf('Unwrapping in 2D space...\n')
ph_uw = zeros(NIFG,NPS);
for k=1:NIFG
    fprintf('Unwrapping interferogram: %d\n',k)
    
    fid=fopen('snaphu.conf','w');
    fprintf(fid,'INFILE  snaphu.in\n');
    fprintf(fid,'OUTFILE snaphu.out\n');
    %fprintf(fid,'COSTINFILE snaphu.costinfile\n'); % TODO: add cost
    %function based on distance between pts
    fprintf(fid,'STATCOSTMODE  DEFO\n');
    fprintf(fid,'INFILEFORMAT  FLOAT_DATA\n');
    fprintf(fid,'OUTFILEFORMAT FLOAT_DATA\n');
    fclose(fid);

    % fid=fopen('snaphu.costinfile','w');
    % fwrite(fid,rowcost','int16');
    % fwrite(fid,colcost','int16');
    % fclose(fid);

    ifgw = zeros(height,width);
    for i=1:NPS
        ifgw(Y(i),X(i)) = residualph(k,i);
    end
    fid=fopen('snaphu.in','w');
    fwrite(fid,ifgw.','float');
    fclose(fid);

    cmdstr=['snaphu -d -f snaphu.conf ',num2str(width),' >& snaphu.log'];
    
    [a,b] =system(cmdstr); % Running SNAPHU
    fid=fopen('snaphu.out');
    ifguw=fread(fid,[width,inf],'float').';
    fclose(fid);
    for i=1:NPS
        ph_uw(k,i) = ifguw(Y(i),X(i)) + modelled_phase_uw(k,i); % Adding back unwrapped residual to model phase 
        %(remember we have not corrected atmosphere, or any other noise yet)
    end
end

sm_cov=eye(NIFG);
master_day = 19970613;
acq_day = sort([acq_times;19970613]);
master_ix = find(acq_day == master_day);
bt = [Btemp([1:master_ix-1]);0;Btemp([master_ix+1:end])];
G = [ones(size(acq_times)),Btemp];
lambda = wavelength;
m = lscov(G,double(ph_uw'),sm_cov); % May give error as we have not considered master - master ifg
ph_all = -m(2,:)'*365.25/4/pi*lambda*1000;