close all
clear all
%==========================================================================

%==========================================================================
%% ------------------- Data obtained from R-script --------------------- %%
%==========================================================================

%==========================================================================
%% Data abundance by season and group:
%==========================================================================
data   = xlsread('PlanktonAbundance_WIN.csv');
%--------------------------------------------------------------------------
VOLUME = data(:,1); % Volume size classes
PHYWIN = data(:,2); % Phytoplantkon abundance in winter
ZOOWIN = data(:,3); % Zooplantkon abundance in winter
%==========================================================================
data   = xlsread('PlanktonAbundance_SUM.csv');
%--------------------------------------------------------------------------
PHYSUM = data(:,2); % Phytoplantkon abundance in summer
ZOOSUM = data(:,3); % Zooplantkon abundance in summer
%==========================================================================

%==========================================================================
%% Regression data:
%==========================================================================
% Phytoplankton in winter
data = xlsread('Regression_PHYWIN.csv');
%--------------------------------------------------------------------------
VFIT_PHYWIN = data(:,1); % Volume
FIT_PHYWIN  = data(:,2); % Abundance predicted
CIL_PHYWIN  = data(:,3); % Confidence interval low limit
CIU_PHYWIN  = data(:,4); % Confidence interval upper limit
PTE_PHYWIN  = data(1,5); % Slope
%==========================================================================
% Phytoplankton in summer
data = xlsread('Regression_PHYSUM.csv');
%--------------------------------------------------------------------------
VFIT_PHYSUM = data(:,1);
FIT_PHYSUM  = data(:,2);
CIL_PHYSUM  = data(:,3);
CIU_PHYSUM  = data(:,4);
PTE_PHYSUM  = data(1,5);
%==========================================================================
% Zooplankton in winter
data = xlsread('Regression_ZOOWIN.csv');
%--------------------------------------------------------------------------
VFIT_ZOOWIN = data(:,1);
FIT_ZOOWIN  = data(:,2);
CIL_ZOOWIN  = data(:,3);
CIU_ZOOWIN  = data(:,4);
PTE_ZOOWIN  = data(1,5);
%==========================================================================
% Zooplankton in summer
data = xlsread('Regression_ZOOSUM.csv');
%--------------------------------------------------------------------------
VFIT_ZOOSUM = data(:,1);
FIT_ZOOSUM  = data(:,2);
CIL_ZOOSUM  = data(:,3);
CIU_ZOOSUM  = data(:,4);
PTE_ZOOSUM  = data(1,5);
%==========================================================================

%==========================================================================
%% FIGURE:
%==========================================================================
make_it_tight = true;
% m = distance between subplots
% n = distance to bottom and top
% p = distance to left and right
subplot = @(m,n,p) subtightplot (m, n, p, [0.01 0.03], [0.15 0.10], [0.08 0.04]);
if ~make_it_tight,  clear subplot;  end
%..........................................................................
Fig = figure;
set(gcf, 'Color','white');
set(Fig, 'Position',[400,200,1000,400]);
%==========================================================================
subplot(1,2,1)
hold on
box on
grid on
%..........................................................................
% Size spectra
plot(log10(VOLUME),log10(PHYWIN),'o','Color',[227  26  28]/255)
plot(log10(VOLUME(3:end)),log10(ZOOWIN(3:end)),'o','Color',[ 34  94 168]/255)
%..........................................................................
% Rregression line
plot(VFIT_PHYWIN,FIT_PHYWIN,'-', 'Color',[227  26  28]/255);
plot(VFIT_ZOOWIN,FIT_ZOOWIN,'--','Color',[ 34  94 168]/255);
%..........................................................................
% Confidence intervals
fill([VFIT_PHYWIN; flipud(VFIT_PHYWIN)],[CIL_PHYWIN; flipud(CIU_PHYWIN)], ...
     [227 26 28]/255, 'FaceAlpha', 0.2, 'EdgeColor', 'none','HandleVisibility','off');  % Phyto winter
fill([VFIT_ZOOWIN; flipud(VFIT_ZOOWIN)],[CIL_ZOOWIN; flipud(CIU_ZOOWIN)], ...
     [34 94 168]/255, 'FaceAlpha', 0.2, 'EdgeColor', 'none','HandleVisibility','off'); % Zoo Winter
%..........................................................................
xlim([log10(0.1) log10(1e11)]) 
xlabel('Volume (\mum^{3})')
set(gca,'XTick',[log10(1) log10(10e1) log10(10e3) log10(10e5) log10(10e7) log10(10e9)],...
        'XTickLabel',{'0' '10^2' '10^3' '10^5' '10^7' '10^9'})
%..........................................................................
ylim([log10(0.4) log10(10e12)])
ylabel('Abundance (ind m^{-3})')
set(gca,'YTick',[log10(1) log10(10e2) log10(10e5) log10(10e8) log10(10e11)],...
        'YTickLabel',{'0' '10^2' '10^5' '10^8' '10^{11}'})
%..........................................................................
title('(a) Size Spectra for February 20^{th}')
%..........................................................................
str = ['Slope = ',sprintf('%.2f',PTE_PHYWIN)];
text(1.70,9.90,str,'LineStyle','none','Color',[227  26  28]/255,'HorizontalAlignment','Left');
%..........................................................................
str = ['Slope = ',sprintf('%.2f',PTE_ZOOWIN)];
text(6.89,4.76,str,'LineStyle','none','Color',[ 34  94 168]/255,'HorizontalAlignment','Left');
%..........................................................................
legend(' Phy abundance',...
       ' Zoo abundance',...
       ' Phy regression line',...
       ' Zoo regression line',...
       ' Phy confidence interval',...
       ' Zoo confidence interval',...
       ' Location','best')
legend boxoff
%..........................................................................
set(gca,'TickLength',[0.005, 0.005])
set(gca,'FontSize',11)
%==========================================================================
subplot(1,2,2)
hold on
box on
grid on
%..........................................................................
plot(log10(VOLUME),log10(PHYSUM),'o','Color',[227  26  28]/255)
plot(log10(VOLUME(3:end)),log10(ZOOSUM(3:end)),'o','Color',[ 34  94 168]/255)
%..........................................................................
plot(VFIT_PHYSUM,FIT_PHYSUM,'-', 'Color',[227  26  28]/255);
plot(VFIT_ZOOSUM,FIT_ZOOSUM,'--','Color',[ 34  94 168]/255);
%..........................................................................
fill([VFIT_PHYSUM; flipud(VFIT_PHYSUM)],[CIL_PHYSUM; flipud(CIU_PHYSUM)], ...
     [227 26 28]/255, 'FaceAlpha', 0.2, 'EdgeColor', 'none');  % Phyto summer
fill([VFIT_ZOOSUM; flipud(VFIT_ZOOSUM)],[CIL_ZOOSUM; flipud(CIU_ZOOSUM)], ...
     [34 94 168]/255, 'FaceAlpha', 0.2, 'EdgeColor', 'none'); % Zoo summer
%..........................................................................
xlim([log10(0.1) log10(1e11)]) 
xlabel('Volume (\mum^{3})')
set(gca,'XTick',[log10(1) log10(10e1) log10(10e3) log10(10e5) log10(10e7) log10(10e9)],...
        'XTickLabel',{'0' '10^2' '10^3' '10^5' '10^7' '10^9'})
%..........................................................................
ylim([log10(0.4) log10(10e12)])
% ylabel('Abundance (ind m^{-3})')
set(gca,'YTick',[log10(1) log10(10e2) log10(10e5) log10(10e8) log10(10e11)],...
        'YTickLabel',[])
%..........................................................................
title('(b) Size Spectra for August 20^{th}')
%..........................................................................
str = ['Slope = ',sprintf('%.2f',PTE_PHYSUM)];
text(2.00,10.14,str,'LineStyle','none','Color',[227  26  28]/255,'HorizontalAlignment','Left');
%..........................................................................
str = ['Slope = ',sprintf('%.2f',PTE_ZOOSUM)];
text(6.87,4.88,str,'LineStyle','none','Color',[ 34  94 168]/255,'HorizontalAlignment','Left');
%..........................................................................
set(gca,'TickLength',[0.005, 0.005])
set(gca,'FontSize',11)
%==========================================================================
%% Save figure:
%..........................................................................
set(gcf,'PaperPositionMode','auto');
print(gcf, '-dpng', 'FIG10_SIZESPECTRA.png');
%--------------------------------------------------------------------------
print('FIG10_SIZESPECTRA', '-depsc'); % Save as an EPS file
%--------------------------------------------------------------------------
% Get the current figure position in inches
fig = gcf;
fig.Units = 'inches';
figPosition = fig.Position;
%..........................................................................
% Set PaperSize to match the figure dimensions
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [figPosition(3) figPosition(4)]);
%..........................................................................
% Save the figure as a PDF using the print command
print(gcf, 'FIG10_SIZESPECTRA.pdf', '-dpdf', '-bestfit');
%==========================================================================
return