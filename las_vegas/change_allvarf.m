IDX_ok1 = find(allvarf<1.0); %% vary for four consecutive values : 1,2,3,4
IDX_ok2 = find(allvarf<2.0); %% vary for four consecutive values : 1,2,3,4
IDX_ok3 = find(allvarf<3.0); %% vary for four consecutive values : 1,2,3,4
IDX_ok4 = find(allvarf<4.0); %% vary for four consecutive values : 1,2,3,4
% --- Plot the estimated parameters at the points ----------------------

disp("DEF data for IDX_ok1: ")
% Deformation
min(allESTPS(2,IDX_ok1))
max(allESTPS(2,IDX_ok1))
mean(allESTPS(2,IDX_ok1))
std(allESTPS(2,IDX_ok1))
disp("DEM data for IDX_ok1: ")
% DEM
min(allESTPS(1,IDX_ok1))
max(allESTPS(1,IDX_ok1))
mean(allESTPS(1,IDX_ok1))
std(allESTPS(1,IDX_ok1))
figure(1)
  subplot(1,2,1)
plotps(X(IDX_ok1), Y(IDX_ok1), allESTPS(1,IDX_ok1), [-40,40]);
hold on
plot(X_REF(IDX_refpnt), Y_REF(IDX_refpnt), 'r+', 'MarkerSize',8);
hold off
title('Estimated DEM error [m] at PS points')
subplot(1,2,2)
plotps(X(IDX_ok1), Y(IDX_ok1), allESTPS(2,IDX_ok1), [-40,40]);
hold on
plot(X_REF(IDX_refpnt), Y_REF(IDX_refpnt), 'r+', 'MarkerSize',8);
hold off
title('Estimated Displ. rates [mm/y] at PS points')


disp("DEF data for IDX_ok2: ")
% Deformation
min(allESTPS(2,IDX_ok2))
max(allESTPS(2,IDX_ok2))
mean(allESTPS(2,IDX_ok2))
std(allESTPS(2,IDX_ok2))
disp("DEM data for IDX_ok2: ")
% DEM
min(allESTPS(1,IDX_ok2))
max(allESTPS(1,IDX_ok2))
mean(allESTPS(1,IDX_ok2))
std(allESTPS(1,IDX_ok2))
figure(2)
  subplot(1,2,1)
plotps(X(IDX_ok2), Y(IDX_ok2), allESTPS(1,IDX_ok2), [-40,40]);
hold on
plot(X_REF(IDX_refpnt), Y_REF(IDX_refpnt), 'r+', 'MarkerSize',8);
hold off
title('Estimated DEM error [m] at PS points')
subplot(1,2,2)
plotps(X(IDX_ok2), Y(IDX_ok2), allESTPS(2,IDX_ok2), [-40,40]);
hold on
plot(X_REF(IDX_refpnt), Y_REF(IDX_refpnt), 'r+', 'MarkerSize',8);
hold off
title('Estimated Displ. rates [mm/y] at PS points')

disp("DEF data for IDX_ok3: ")
% Deformation
min(allESTPS(2,IDX_ok3))
max(allESTPS(2,IDX_ok3))
mean(allESTPS(2,IDX_ok3))
std(allESTPS(2,IDX_ok3))
disp("DEM data for IDX_ok3: ")
% DEM
min(allESTPS(1,IDX_ok3))
max(allESTPS(1,IDX_ok3))
mean(allESTPS(1,IDX_ok3))
std(allESTPS(1,IDX_ok3))
figure(3)
  subplot(1,2,1)
plotps(X(IDX_ok3), Y(IDX_ok3), allESTPS(1,IDX_ok3), [-40,40]);
hold on
plot(X_REF(IDX_refpnt), Y_REF(IDX_refpnt), 'r+', 'MarkerSize',8);
hold off
title('Estimated DEM error [m] at PS points')
subplot(1,2,2)
plotps(X(IDX_ok3), Y(IDX_ok3), allESTPS(2,IDX_ok3), [-40,40]);
hold on
plot(X_REF(IDX_refpnt), Y_REF(IDX_refpnt), 'r+', 'MarkerSize',8);
hold off
title('Estimated Displ. rates [mm/y] at PS points')


disp("DEF data for IDX_ok4: ")
% Deformation
min(allESTPS(2,IDX_ok4))
max(allESTPS(2,IDX_ok4))
mean(allESTPS(2,IDX_ok4))
std(allESTPS(2,IDX_ok4))
disp("DEM data for IDX_ok4: ")
% DEM
min(allESTPS(1,IDX_ok4))
max(allESTPS(1,IDX_ok4))
mean(allESTPS(1,IDX_ok4))
std(allESTPS(1,IDX_ok4))
figure(4)
  subplot(1,2,1)
plotps(X(IDX_ok4), Y(IDX_ok4), allESTPS(1,IDX_ok4), [-40,40]);
hold on
plot(X_REF(IDX_refpnt), Y_REF(IDX_refpnt), 'r+', 'MarkerSize',8);
hold off
title('Estimated DEM error [m] at PS points')
subplot(1,2,2)
plotps(X(IDX_ok4), Y(IDX_ok4), allESTPS(2,IDX_ok4), [-40,40]);
hold on
plot(X_REF(IDX_refpnt), Y_REF(IDX_refpnt), 'r+', 'MarkerSize',8);
hold off
title('Estimated Displ. rates [mm/y] at PS points')