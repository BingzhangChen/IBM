File = 'Euler.nc';
Z_r = ncread(File, 'Z_r');
Days = ncread(File, 'Day');
CDiv_avg = ncread(File, 'CDiv_avg');
CDiv_var = ncread(File, 'CDiv_var');
Topt_avg = ncread(File, 'Topt_avg');
Topt_var = ncread(File, 'Topt_var');
Lnalpha_avg = ncread(File, 'Lnalpha_avg');
Lnalpha_var = ncread(File, 'Lnalpha_var');
TLnV_cov = ncread(File, 'TLnV_cov');
Talp_cov = ncread(File, 'Talp_cov');
ALnV_cov = ncread(File, 'ALnV_cov');

%Check how many years of the model run
NYear = Days(end)/365;
Years = double(Days)/365;

%Convert CDiv_avg to Volume
Vavg = (12 .* exp(CDiv_avg) ./10^(-.69) ).^(1/.88);

%Convert Volume to ESD
ESD_avg = (6/pi .* Vavg).^(1/3);

%Only plot final year
finalyear = 1;

if (finalyear == 1)
    ESD_avg  = ESD_avg(:,  (Days(end)-365):Days(end));
    CDiv_var = CDiv_var(:, (Days(end)-365):Days(end));
    Topt_var = Topt_var(:, (Days(end)-365):Days(end));
    Topt_avg = Topt_avg(:, (Days(end)-365):Days(end));
    Lnalpha_avg = Lnalpha_avg(:, (Days(end)-365):Days(end));
    Lnalpha_var = Lnalpha_var(:, (Days(end)-365):Days(end));
    TLnV_cov = TLnV_cov(:, (Days(end)-365):Days(end));
    ALnV_cov = ALnV_cov(:, (Days(end)-365):Days(end));
    Talp_cov = Talp_cov(:, (Days(end)-365):Days(end));
    xvar = 0:365;
    xlab = [0:180:365];
else
    xvar = Years;
    xlab = [0:1:NYear];
end;

t = tiledlayout(3,3,'TileSpacing','Compact');

%plot phytoplankton mean size
nexttile;
h = pcolor(xvar, Z_r, ESD_avg);
xticks(xlab)
axHdl = get(h, 'Parent');
z = get(axHdl, 'CLim');
drawnow
set(axHdl, 'CLim', [1, 120]);
set(axHdl, 'ColorScale', 'log')
set(axHdl,'TickDir','out'); 

colorbar
shading flat;
hold on
title('(a) Mean ESD (\mum)')
ax = gca;
ax.TitleHorizontalAlignment = 'left';
ax.FontSize = 7; 
ax.TitleFontSizeMultiplier = 1;

%Plot phyto. size diversity
nexttile;
h = pcolor(xvar, Z_r, CDiv_var);
xticks(xlab)
axHdl = get(h, 'Parent');
z = get(axHdl, 'CLim');
drawnow
set(axHdl, 'CLim', [0.5, 15]);
set(axHdl,'TickDir','out'); 
colorbar
shading flat;
hold on
title('(b) Var(CDiv) (log pmol C cell^{-1})^2')
ax = gca;
ax.TitleHorizontalAlignment = 'left';
ax.FontSize = 7; 
ax.TitleFontSizeMultiplier = 1;

%Plot mean Topt
nexttile;
h = pcolor(xvar, Z_r, Topt_avg);
xticks(xlab)
axHdl = get(h, 'Parent');
z = get(axHdl, 'CLim');
drawnow
set(axHdl, 'CLim', [20, 26]);
set(axHdl,'TickDir','out'); 
colorbar
shading flat;
hold on;
title('(c) Mean T_{opt} (ºC)')
ax = gca;
ax.TitleHorizontalAlignment = 'left';
ax.FontSize = 7; 
ax.TitleFontSizeMultiplier = 1;

%Plot  Topt  variance
nexttile;
h = pcolor(xvar, Z_r, Topt_var);
xticks(xlab)
axHdl = get(h, 'Parent');
z = get(axHdl, 'CLim');
drawnow
set(axHdl, 'CLim', [0.5, 100]);
set(axHdl,'TickDir','out'); 
set(axHdl, 'ColorScale', 'log')
colorbar
shading flat;
hold on;
title('(d) Var(T_{opt}) (ºC^2)');
ax = gca;
ax.TitleHorizontalAlignment = 'left';
ax.FontSize = 7; 
ax.TitleFontSizeMultiplier = 1;

%Plot mean alphaChl
nexttile;
h = pcolor(xvar, Z_r, Lnalpha_avg);
xticks(xlab)
axHdl = get(h, 'Parent');
z = get(axHdl, 'CLim');
drawnow
set(axHdl, 'CLim', [-2, -1]);
set(axHdl,'TickDir','out'); 
colorbar
shading flat;
hold on;
title('(e) Mean Ln \alpha^{Chl}');
ax = gca;
ax.TitleHorizontalAlignment = 'left';
ax.FontSize = 7; 
ax.TitleFontSizeMultiplier = 1;

%Plot  Lnalpha  variance
nexttile;
h = pcolor(xvar, Z_r, Lnalpha_var);
xticks(xlab)
axHdl = get(h, 'Parent');
z = get(axHdl, 'CLim');
drawnow
set(axHdl, 'CLim', [0.01, .1]);
set(axHdl,'TickDir','out'); 
colorbar
shading flat;
hold on;
title('(f) Var(Ln \alpha^{Chl})')
ax = gca;
ax.TitleHorizontalAlignment = 'left';
ax.FontSize = 7; 
ax.TitleFontSizeMultiplier = 1;

%Plot  covariance between Topt and size (lnV)
nexttile;
h = pcolor(xvar, Z_r, TLnV_cov);
xticks(xlab)
axHdl = get(h, 'Parent');
z = get(axHdl, 'CLim');
drawnow
set(axHdl, 'CLim', [-1, 1]);
set(axHdl,'TickDir','out'); 
colorbar
shading flat;
hold on
title('(g) Cov(T_{opt}, size)')
ax = gca;
ax.TitleHorizontalAlignment = 'left';
ax.FontSize = 7; 
ax.TitleFontSizeMultiplier = 1;

%Plot  covariance between Topt and Lnalpha
nexttile;
h = pcolor(xvar, Z_r, Talp_cov);
xticks(xlab)
axHdl = get(h, 'Parent');
z = get(axHdl, 'CLim');
drawnow
set(axHdl, 'CLim', [-.2, .1]);
set(axHdl,'TickDir','out'); 
colorbar
shading flat;
hold on
title('(h) Cov(T_{opt}, Ln \alpha^{Chl})')
ax = gca;
ax.TitleHorizontalAlignment = 'left';
ax.FontSize = 7; 
ax.TitleFontSizeMultiplier = 1;

%Plot  covariance between size and Lnalpha
nexttile;
h = pcolor(xvar, Z_r, ALnV_cov);
xticks(xlab)
axHdl = get(h, 'Parent');
z = get(axHdl, 'CLim');
drawnow
set(axHdl, 'CLim', [-.2, 0.6]);
set(axHdl,'TickDir','out'); 
colorbar
shading flat;
hold on
title('(i) Cov(size, Ln \alpha^{Chl})')
ax = gca;
ax.TitleHorizontalAlignment = 'left';
ax.FontSize = 7; 
ax.TitleFontSizeMultiplier = 1;

if (finalyear == 1)
  xlabel(t, 'Day', 'FontSize',9)
else
  xlabel(t, 'Year', 'FontSize',9)
end
ylabel(t, 'Depth (m)', 'FontSize',9)

% Print to a pdf file
set(gcf,'PaperUnits','normalized');
set(gcf,'PaperPositionMode', 'manual');
set(gcf,'PaperPosition', [0 0 1 1]);

set(gcf,'PaperOrientation','landscape');
set(gcf,'PaperType','a4');
%c=datestr(datetime('today'));
exportgraphics(gcf,'Fig8_PHY_trait.pdf','ContentType','vector');
close all;
%% 

