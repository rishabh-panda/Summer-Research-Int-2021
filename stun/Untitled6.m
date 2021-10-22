figure(2)
  subplot(1,2,1)
plotps(X(IDX_ok), Y(IDX_ok), allESTPS(1,IDX_ok), [-40,40]);
hold on
plot(X_REF(IDX_refpnt), Y_REF(IDX_refpnt), 'r+', 'MarkerSize',8);
hold off
title('Estimated DEM error [m] at PS points')
subplot(1,2,2)
plotps(X(IDX_ok), Y(IDX_ok), allESTPS(2,IDX_ok), [-40,40]);
hold on
plot(X_REF(IDX_refpnt), Y_REF(IDX_refpnt), 'r+', 'MarkerSize',8);
hold off
title('Estimated Displ. rates [mm/y] at PS points')
min(allESTPS(2,IDX_ok))
max(allESTPS(2,IDX_ok))
mean(allESTPS(2,IDX_ok))
std(allESTPS(2,IDX_ok))