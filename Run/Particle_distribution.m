File = 'Euler.nc';
Z_r = ncread(File, 'Z_r');
Days = ncread(File, 'Day');
N_ind=ncread(File, 'N_ind');
N_cell=ncread(File, 'N_cell');

t = tiledlayout(1,2,'TileSpacing','Compact');

%plot number of super-individuals
nexttile;
h = pcolor(Days, Z_r, N_ind);
xticks([0 180 365 540 730])
axHdl = get(h, 'Parent');
z = get(axHdl, 'CLim');
drawnow
set(axHdl, 'CLim', [5, 12]);
set(axHdl,'TickDir','out'); 
colorbar
shading flat;
hold on
title('Number of super-individuals m^{-3}')

%Plot number of cells
nexttile;
N_cell = N_cell./10^10;
h = pcolor(Days, Z_r, N_cell);
xticks([0 180 365 540 730])
axHdl = get(h, 'Parent');
z = get(axHdl, 'CLim');
drawnow
set(axHdl, 'CLim', [0, 1]);
set(axHdl,'TickDir','out'); 
colorbar
shading flat;
hold on
title('10^{10} cells m^{-3}')

xlabel(t, 'Time')
ylabel(t, 'Depth (m)')

% Print to a pdf file
set(gcf,'PaperUnits','normalized');
set(gcf,'PaperPositionMode', 'auto');
set(gcf,'PaperOrientation','landscape');
set(gcf,'PaperType','a4');

print('-dpdf','Particle_distribution.pdf');
close all;
%% 

