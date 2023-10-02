File = 'Euler.nc';
Z_r = ncread(File, 'Z_r');
Days = ncread(File, 'Day');
ZOO = ncread(File, 'ZOO');
TZOO = squeeze(sum(ZOO, 1));
NZOO = size(ZOO, 1);
%Check how many years of the model run
NYear = Days(end)/365;
Years = double(Days)/365;


t = tiledlayout(4,5,'TileSpacing','Compact');

%plot temperature
for i = 1:NZOO,
  cff = squeeze(ZOO(i,:,:));
  nexttile;
  h = pcolor(Years, Z_r, cff);
  xticks([0:1:NYear])
  axHdl = get(h, 'Parent');
  z = get(axHdl, 'CLim');
  drawnow
  set(axHdl, 'CLim', [0, .1]);
  set(axHdl,'TickDir','out'); 
  shading flat;
  hold on
  title(['ZOO', num2str(i)])
end

colorbar
xlabel(t, 'Year')
ylabel(t, 'Depth (m)')

% Print to a pdf file
set(gcf,'PaperUnits','normalized');
set(gcf,'PaperPositionMode', 'manual');
set(gcf,'PaperPosition', [0 0 1 1]);
set(gcf,'PaperOrientation','landscape');
set(gcf,'PaperType','a4');

print('-dpdf','ZOO_output.pdf');
close all;
%% 

