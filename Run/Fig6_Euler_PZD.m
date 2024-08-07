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

%Only plot final year
finalyear = 1;

if (finalyear == 1)
    Temp  = Temp(:,  (Days(end)-365):Days(end));
    Kv = Kv(:, (Days(end)-365):Days(end));
    NO3 = NO3(:, (Days(end)-365):Days(end));
    PN = PN(:, (Days(end)-365):Days(end));
    PC = PC(:, (Days(end)-365):Days(end));
    DET = DET(:, (Days(end)-365):Days(end));
    CHL = CHL(:, (Days(end)-365):Days(end));
    TZOO = TZOO(:, (Days(end)-365):Days(end));
    NPP = NPP(:, (Days(end)-365):Days(end));
    xvar = 0:365;
    xlab = [0:180:365];

else
    xvar = Years;
    xlab = [0:1:NYear];


end;

%Check how many years of the model run
NYear = Days(end)/365;
Years = double(Days)/365;

t = tiledlayout(2,2,'TileSpacing','Compact');

%plot phytoplankton carbon
nexttile;
h = pcolor(xvar, Z_r, PC);
xticks(xlab)
axHdl = get(h, 'Parent');
z = get(axHdl, 'CLim');
drawnow
set(axHdl, 'CLim', [0, 3]);
set(axHdl,'TickDir','out'); 
colorbar
shading flat;
title('(a) Phyto carbon (mmol C m^{-3})')
ax = gca;
ax.TitleHorizontalAlignment = 'left';
ax.TitleFontSizeMultiplier = 1;
ax.FontSize = 8; 
nexttile;
h = pcolor(xvar, Z_r, PN);
xticks(xlab)
axHdl = get(h, 'Parent');
z = get(axHdl, 'CLim');
drawnow
set(axHdl, 'CLim', [0, .5]);
set(axHdl,'TickDir','out'); 
colorbar
shading flat;
title('(b) Phyto nitrogen (mmol N m^{-3})')
ax = gca;
ax.TitleHorizontalAlignment = 'left';
ax.TitleFontSizeMultiplier = 1;
ax.FontSize = 8; 

nexttile;

h = pcolor(xvar, Z_r, TZOO);
xticks(xlab)
axHdl = get(h, 'Parent');
z = get(axHdl, 'CLim');
drawnow
set(axHdl, 'CLim', [0, .8]);
set(axHdl,'TickDir','out'); 
colorbar
shading interp
title('(c) Total Zooplankton (mmol N m^{-3})')
ax = gca;
ax.TitleHorizontalAlignment = 'left';
ax.TitleFontSizeMultiplier = 1;
ax.FontSize = 8; 

%Plot detritus
nexttile;

h = pcolor(xvar, Z_r, DET);
xticks(xlab)
axHdl = get(h, 'Parent');
z = get(axHdl, 'CLim');
drawnow
set(axHdl, 'CLim', [0, .2]);
set(axHdl,'TickDir','out'); 
colorbar
shading interp
title('(d) Detritus (mmol N m^{-3})')
ax = gca;
ax.TitleHorizontalAlignment = 'left';
ax.TitleFontSizeMultiplier = 1;
ax.FontSize = 8; 
 xlabel(t, 'Day','FontSize',9);
 ylabel(t, 'Depth (m)','FontSize',9);

% Print to a pdf file
set(gcf,'PaperUnits','normalized');
set(gcf,'PaperPositionMode', 'manual');
set(gcf,'PaperPosition', [0 0 1 1]);

set(gcf,'PaperOrientation','landscape');
set(gcf,'PaperType','a4');
exportgraphics(gcf,'Fig6_PZD.pdf','ContentType','vector');
% c=datestr(datetime('today'));
% if (finalyear == 1)
%   print('-dpdf',['PZD_output_finalyear', c, '.pdf']);
% else
%    print('-dpdf',['PZD_output_allyear', c, '.pdf']);
% end
close all;
%% 

