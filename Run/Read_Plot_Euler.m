File = 'Euler.nc';
Z_r = ncread(File, 'Z_r');
Days = ncread(File, 'Day');
NO3 = ncread(File, 'NO3');
PN = ncread(File, 'PN');
PC = ncread(File, 'PC');
CHL = ncread(File, 'CHL');
DET = ncread(File, 'DET');

tiledlayout(2,2)
nexttile;
h = pcolor(Days, Z_r, NO3);
axHdl = get(h, 'Parent');
z = get(axHdl, 'CLim');
drawnow
set(axHdl, 'CLim', [0, .8]);
colorbar
shading flat;
hold on
title('Nitrate (mmol m^-3)')

nexttile;
h = pcolor(Days, Z_r, PC);
axHdl = get(h, 'Parent');
z = get(axHdl, 'CLim');
drawnow
set(axHdl, 'CLim', [0, 4]);
colorbar
shading flat;
title('Phyto C')

nexttile;
h = pcolor(Days, Z_r, PN);
axHdl = get(h, 'Parent');
z = get(axHdl, 'CLim');
drawnow
set(axHdl, 'CLim', [0, .6]);
colorbar
shading flat;
title('Phyto N')

nexttile;

h = pcolor(Days, Z_r, CHL);
axHdl = get(h, 'Parent');
z = get(axHdl, 'CLim');
drawnow
set(axHdl, 'CLim', [0, 2]);
colorbar
shading interp
title('Chl')

% Print to a pdf file
set(gcf,'PaperUnits','normalized');
set(gcf,'PaperPositionMode', 'auto');
set(gcf,'PaperOrientation','landscape');
set(gcf,'PaperType','a4');

print('-dpdf','Euler_output.pdf');
close all;
%% 

