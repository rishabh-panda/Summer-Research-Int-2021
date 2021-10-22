disp(['  max DEM error:  ', num2str(max(allESTPS(1,IDX_ok)))]);
disp(['  min DEM error:  ', num2str(min(allESTPS(1,IDX_ok)))]);
disp(['  mean DEM error:  ', num2str(mean(allESTPS(1,IDX_ok)))]);
disp(['  std DEM error:  ', num2str(std(allESTPS(1,IDX_ok)))]);
disp(['  max DEF:  ', num2str(max(allESTPS(2,IDX_ok)))]);
disp(['  min DEF:  ', num2str(min(allESTPS(2,IDX_ok)))]);
disp(['  mean DEF:  ', num2str(mean(allESTPS(2,IDX_ok)))]);
disp(['  std DEF:  ', num2str(std(allESTPS(2,IDX_ok)))]);
figure(1)
 subplot(1,2,1)
  plot(allESTPS(1,IDX_ok),'.');
  title(' DEM error [m]');
  xlabel('arcs');
  ylabel('[m]');
 subplot(1,2,2) 
  plot(allESTPS(2,IDX_ok),'.');
  title('DEF [mm/y]');
  xlabel('arcs');
  ylabel('[mm/y]');  
