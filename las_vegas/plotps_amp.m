function h = plotps_amp(X,Y,Z, amp)
%PLOTPS   plot scattered data triplets.
%   PLOTPS(X,Y,Z) plots the points using the current colormap
%   to the current figure.
%
%   PLOTPS(X,Y,Z, ZLIM) where ZLIM is a two element vectors
%   clips the Z data to the given range.
%
%   H=PLOTPS(...) returns handle H to the figure.
%
%   Example:
%     x = round(rand(100,1)*10000);
%     y = round(rand(100,1)*10000);
%     z = round(randn(100,1)*10);
%     plotps(x,y,z,[-20,20]);
%
%   See also STUN.

% This function is part of the STUN toolbox which
% accompanies the book:
% Kampes, B.M., "Radar Interferometry -- The Persistent
% Scatterer Technique", published by Springer, 2006.

% You are allowed to use and/or modify this code for your
% own purposes.

% The STUN toolbox is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

% Written by (c)Bert Kampes, 21-Oct-2005
% Tested with Matlab version: 5.3.0.10183 (R11)
% $Revision: 1.3 $  $Date: 2005/12/07 08:23:06 $



% % --- Check input -------------------------------------------
% error(nargchk(3,4,nargin))
% if (length(X)~=length(Y))
%   error('length of X and Y must be equal');
% end
% if (length(X)~=length(Z))
%   error('length of X and Z must be equal');
% end
% if (nargin<4)
%   ZLIM = [floor(min(Z)-eps), ceil(max(Z)+eps)];% default
% end
% if (length(ZLIM)~=2)
%   error('ZLIM should be a length 2 vector');
% end
% Z = [Z(:);ZLIM(1);ZLIM(2)]; % add dummies for scaling
% Z(find(Z<=ZLIM(1))) = ZLIM(1);
% Z(find(Z>=ZLIM(2))) = ZLIM(2);
% 
% 
% 
% % --- Create the plot ---------------------------------------
% cmap = colormap(flipud(jet(64)));% current map
% NCOL = size(cmap,1);
% % --- Scale Z to range 1:NCOL ---
% clipz_min = min(Z);% plotting range
% clipz_max = max(Z);% plotting range
% if (clipz_min==clipz_max)
%   ZC   = ones(size(Z));
% else
%   ZC   = round(1+(NCOL-1-eps).*(Z-clipz_min)./(clipz_max-clipz_min));% [1:64]
% end
% 
% 
% % --- Plot points with different colors per point -----------
% hold off
% for i=1:length(X)
%   plot(X(i),Y(i),'s', ...
%     'Color',cmap(ZC(i),:), ...
%     'MarkerFaceColor',cmap(ZC(i),:), ...
%     'MarkerSize',4);
%   hold on
% end
% hold off
% % --- Add colorbar with corrected annotation ---------------- AxesH =
% % axes('CLim', [1, NCOL]);
% 
% caxis([1,NCOL]);
% c = colorbar;
% % colormap(flipud(jet(64)));
% zero = round(1+(NCOL-eps)*(-1)*clipz_min./(clipz_max-clipz_min));
% c.YTick = [1,zero,NCOL];
% c.YTickLabel = [clipz_min,0,clipz_max];
% % set(c,'YTick',[1,NCOL]); set(c,'YTIckLabel',[clipz_min,clipz_max]);
% 
% axis tight;
% 
% 
% % --- Return figure handle if requested ---------------------
% if (nargout==1)
%   h = gcf;
% end

%%% EOF
clipz_min = min(Z);% plotting range
clipz_max = max(Z);% plotting range
c=flipud(jet(64));
amp_mean=uint16(amp);
[n,m]=size(amp_mean);
c=[gray(256);c];

est_range = clipz_max - clipz_min;
col_ix=round(((Z-clipz_min)*63/est_range)+1);

col_ix(col_ix>64)=64;
col_ix(col_ix<1)=1;

pixel_aspect_ratio = 4.858995206750195*2;
az_ix=pixel_aspect_ratio:pixel_aspect_ratio:m*pixel_aspect_ratio;
plot_pix_m=40;
mean_x_post = 19.833000000000002;
plot_pixel_size=round(abs(plot_pix_m/mean_x_post));

pixel_margin1=floor((plot_pixel_size-1)/2);
pixel_margin2=ceil((plot_pixel_size-1)/2);
    
for i=1 : length(X)
    ix1 = Y(i)-floor((pixel_aspect_ratio*plot_pixel_size/2)-1):Y(i)+ceil(pixel_aspect_ratio*plot_pixel_size/2);
    ix1 = ix1(ix1>0&ix1<=n);
    ix2=X(i)-pixel_margin1+1:X(i)+pixel_margin2+1;
    ix2=ix2(ix2>0&ix2<=m);
    if ~(isnan(col_ix(i)))
        amp_mean(ix1,ix2)=256+col_ix(i);
    end
end

amp_mean=(fliplr(amp_mean));

image(az_ix,[1:size(amp_mean,1)],amp_mean)
set(gca,'xtick',[],'ytick',[])
% axis equal
axis tight
hold off
colormap(c);

plotlims =[clipz_min, clipz_max];
textcolor = [0 0 0.004];
textsize = 8;
units = 'mm/y';
h = colorbar('peer',gca);
ylim=get(h,'ylim');
set(h,'ylim',[ylim(2)-64,ylim(2)])
set(h,'ytick',[ylim(2)-64,ylim(2)],'yticklabel',plotlims,'xcolor',textcolor,'ycolor',textcolor,'fontweight','bold','color',textcolor,'FontSize',abs(textsize))
set(get(h,'ylabel'),'String',units,'FontSize',abs(textsize))  
    