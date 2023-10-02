File = 'Euler.nc';
Z_r = ncread(File, 'Z_r');
Z_w=ncread(File, 'Z_w');
Days = ncread(File, 'Day');
Temp=ncread(File, 'Temp');
Kv  = ncread(File, 'Kv');
NO3 = ncread(File, 'NO3');
PN = ncread(File, 'PN');
PC = ncread(File, 'PC');
CHL = ncread(File, 'CHL');
DET = ncread(File, 'DET');
ZOO = ncread(File, 'ZOO');
TZOO = squeeze(sum(ZOO, 1));
NPP = ncread(File, 'NPP');

%Check how many years of the model run
NYear = Days(end)/365;
Years = double(Days)/365;

t = tiledlayout(3,3,'TileSpacing','Compact');

%plot temperature
nexttile;
h = pcolor(Years, Z_r, Temp);
xticks([0:1:NYear])
axHdl = get(h, 'Parent');
z = get(axHdl, 'CLim');
drawnow
set(axHdl, 'CLim', [16, 28]);
set(axHdl,'TickDir','out'); 
colorbar
shading flat;
hold on
title('Temperature (ºC)')

%Plot Kv
nexttile;
h = pcolor(Years, Z_w, Kv);
xticks([0:1:NYear])
axHdl = get(h, 'Parent');
z = get(axHdl, 'CLim');
drawnow
set(axHdl, 'CLim', [0, .1]);
set(axHdl,'TickDir','out'); 
colorbar
shading flat;
hold on
title('Diffusivity (m^2 s^{-1})')

%Plot nitrate
nexttile;
h = pcolor(Years, Z_r, NO3);
xticks([0:1:NYear])
axHdl = get(h, 'Parent');
z = get(axHdl, 'CLim');
drawnow
set(axHdl, 'CLim', [0, .8]);
set(axHdl,'TickDir','out'); 
colorbar
shading flat;
hold on
title('Nitrate (mmol m^-3)')

nexttile;
h = pcolor(Years, Z_r, PC);
xticks([0:1:NYear])
axHdl = get(h, 'Parent');
z = get(axHdl, 'CLim');
drawnow
set(axHdl, 'CLim', [0, 4]);
set(axHdl,'TickDir','out'); 
colorbar
shading flat;
title('Phyto C')

nexttile;
h = pcolor(Years, Z_r, PN);
xticks([0:1:NYear])
axHdl = get(h, 'Parent');
z = get(axHdl, 'CLim');
drawnow
set(axHdl, 'CLim', [0, .6]);
set(axHdl,'TickDir','out'); 
colorbar
shading flat;
title('Phyto N')

nexttile;

h = pcolor(Years, Z_r, CHL);
xticks([0:1:NYear])
axHdl = get(h, 'Parent');
z = get(axHdl, 'CLim');
drawnow
set(axHdl, 'CLim', [0, 2]);
set(axHdl,'TickDir','out'); 
colorbar
shading interp
title('Chl')

nexttile;

h = pcolor(Years, Z_r, TZOO);
xticks([0:1:NYear])
axHdl = get(h, 'Parent');
z = get(axHdl, 'CLim');
drawnow
set(axHdl, 'CLim', [0, .5]);
set(axHdl,'TickDir','out'); 
colorbar
shading interp
title('Total Zoo')

%Plot detritus
nexttile;

h = pcolor(Years, Z_r, DET);
xticks([0:1:NYear])
axHdl = get(h, 'Parent');
z = get(axHdl, 'CLim');
drawnow
set(axHdl, 'CLim', [0, .1]);
set(axHdl,'TickDir','out'); 
colorbar
shading interp
title('Detritus')

%Plot NPP
nexttile;

h = pcolor(Years, Z_r, NPP);
xticks([0:1:NYear])
axHdl = get(h, 'Parent');
z = get(axHdl, 'CLim');
drawnow
set(axHdl, 'CLim', [0, 50]);
set(axHdl,'TickDir','out'); 
colorbar
shading interp
title('NPP')

xlabel(t, 'Year')
ylabel(t, 'Depth (m)')


% Print to a pdf file
set(gcf,'PaperUnits','normalized');
set(gcf,'PaperPositionMode', 'manual');
set(gcf,'PaperPosition', [0 0 1 1]);

set(gcf,'PaperOrientation','landscape');
set(gcf,'PaperType','a4');

print('-dpdf','Euler_output.pdf');
close all;
%% 

