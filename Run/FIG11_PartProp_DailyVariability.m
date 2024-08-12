close all
clear all
format long g
warning off
%==========================================================================

%==========================================================================
%% Script to plot daily variability of cellular content and environmental conditions
% during winter and summer conditions of a selected super-individual
%--------------------------------------------------------------------------
% Writen by Iria Sala
% Created on 20.02.2023
%==========================================================================

%==========================================================================
%% Load Lagrangian fields:
%==========================================================================
%..........................................................................
LAG  = 'ParY6.nc';
%..........................................................................
HOUR = 1:8737;
TIME = repmat(HOUR',1,20000);           % Time (hours)
ZDEP = double(nc_varget(LAG,'Z'));      % Depth (m)
%..........................................................................
CELL = double(nc_varget(LAG,'ID'));     % Particle ID
CNUM = double(nc_varget(LAG,'N_cell')); % Cellular abundance
%..........................................................................
CPAR = double(nc_varget(LAG,'PAR'));    % Environmental PAR (W m-2)
TEMP = double(nc_varget(LAG,'Temp'));   % Environmental temperature (degree C)
CNO3 = double(nc_varget(LAG,'NO3'));    % Environmental NO3 (mmol N m-3)
%..........................................................................
PHYN = double(nc_varget(LAG,'PN'));     % Nutrogen content (pmol N cell-1)
PHYC = double(nc_varget(LAG,'PC'));     % Carbon content (pmol C cell-1)     
PCHL = double(nc_varget(LAG,'CHL'));    % Chlorophyll content (pg Chl cell-1)
%==========================================================================

%==========================================================================
%% Locate super-individual #30477:
%==========================================================================
NUMSI = 149306;
idx   = find(CELL == NUMSI);
TIMEi = TIME(idx);
%--------------------------------------------------------------------------
ZDEPi = ZDEP(idx);
CPARi = CPAR(idx);
TEMPi = TEMP(idx);
CNO3i = CNO3(idx);
%--------------------------------------------------------------------------
CNUMi = CNUM(idx);
PHYNi = PHYN(idx);
PHYCi = PHYC(idx);
PCHLi = PCHL(idx);
%==========================================================================

%==========================================================================
%% Winter period (1:2160):
%==========================================================================
TiW = 1;
TeW = 169;
%--------------------------------------------------------------------------
TIMEw = TIMEi(TiW:TeW);
TIMEw = 1:1/24:8;
TVECw = 1:1:8; % Time vector for plotting
TOUTw = 1:1:8;
%--------------------------------------------------------------------------
ZDEPw = ZDEPi(TiW:TeW);
CPARw = CPARi(TiW:TeW);
TEMPw = TEMPi(TiW:TeW);
CNO3w = CNO3i(TiW:TeW);
%--------------------------------------------------------------------------
CNUMw = CNUMi(TiW:TeW);
PHYNw = PHYNi(TiW:TeW);
PHYCw = PHYCi(TiW:TeW);
PCHLw = PCHLi(TiW:TeW);
%--------------------------------------------------------------------------
RAQNw = PHYNw ./ PHYCw;
CHLCw = PCHLw ./ PHYCw;
%==========================================================================

%==========================================================================
%% Summer period (4321:6552):
%==========================================================================
TiS = 4321;
TeS = 4489;
%--------------------------------------------------------------------------
TIMEs = TIMEi(TiS:TeS);
TIMEs = 1:1/24:8;
TVECs = 1:1:8; % Time vector for plotting
TOUTs = 1:1:8;
%--------------------------------------------------------------------------
ZDEPs = ZDEPi(TiS:TeS);
CPARs = CPARi(TiS:TeS);
TEMPs = TEMPi(TiS:TeS);
CNO3s = CNO3i(TiS:TeS);
%--------------------------------------------------------------------------
CNUMs = CNUMi(TiS:TeS);
PHYNs = PHYNi(TiS:TeS);
PHYCs = PHYCi(TiS:TeS);
PCHLs = PCHLi(TiS:TeS);
%--------------------------------------------------------------------------
RAQNs = PHYNs ./ PHYCs;
CHLCs = PCHLs ./ PHYCs;
%==========================================================================


%==========================================================================
%% FIGURE:
%==========================================================================
make_it_tight = true;
% m = distance between subplots
% n = distance to bottom and top
% p = distance to left and right
subplot = @(m,n,p) subtightplot (m, n, p, [0.04 0.07], [0.11 0.03], [0.07 0.02]);
if ~make_it_tight,  clear subplot;  end
%..........................................................................
Fig = figure;
set(gcf, 'Color','white');
set(Fig, 'Position',[220,100,1200,500]);
%==========================================================================
subplot(2,4,1)
box on
grid on
hold on
%..........................................................................
plot(TIMEw,ZDEPw,'Color',[ 34  94 168]/255,'LineStyle','--','LineWidth',1.5);
plot(TIMEs,ZDEPs,'Color',[227  26  28]/255,'LineStyle','-','LineWidth',1.5);
%..........................................................................
xlim([TIMEw(1) TIMEw(end)]);
% xlabel('Time (months)');
set(gca,'XTick',TVECw,'XTickLabel',[],'FontSize',11)
%..........................................................................
ylim([-50 0]);
ylabel('Depth (m)')
% set(gca,'YTick',[-175,-140,-105,-70,-35,0],'YTickLabel',[-175,-140,-105,-70,-35,0],'FontSize',11)
%..........................................................................
str = '(a)';
text(7.10,-3,str,'LineStyle','none','Color','k','HorizontalAlignment','Left');
%..........................................................................
set(gca,'TickLength',[0.005, 0.005])
set(gca,'FontSize',11,'FontName','Arial');
%..........................................................................
hold off
%==========================================================================
subplot(2,4,2)
box on
grid on
hold on
%..........................................................................
plot(TIMEw,TEMPw,'Color',[ 34  94 168]/255,'LineStyle','--','LineWidth',1.5);
plot(TIMEs,TEMPs,'Color',[227  26  28]/255,'LineStyle','-','LineWidth',1.5);
%..........................................................................
xlim([TIMEw(1) TIMEw(end)]);
% xlabel('Time (months)');
set(gca,'XTick',TVECw,'XTickLabel',[],'FontSize',11)
%..........................................................................
ylim([19.00 28.00]);
ylabel('Temperature (ºC)')
% set(gca,'YTick',[18.85:0.20:19.85],'YTickLabel',[18.85:0.20:19.85],'FontSize',11)
%..........................................................................
legend('Winter period','Summer period','Location','best','FontSize',11)
legend box off
%..........................................................................
str = '(b)';
text(7.10,27.4,str,'LineStyle','none','Color','k','HorizontalAlignment','Left');
%..........................................................................
set(gca,'TickLength',[0.005, 0.005])
set(gca,'FontSize',11,'FontName','Arial');
%..........................................................................
hold off
%==========================================================================
subplot(2,4,3)
box on
grid on
hold on
%..........................................................................
plot(TIMEw,CPARw,'Color',[ 34  94 168]/255,'LineStyle','--','LineWidth',1.5);
plot(TIMEs,CPARs,'Color',[227  26  28]/255,'LineStyle','-','LineWidth',1.5);
%..........................................................................
xlim([TIMEw(1) TIMEw(end)]);
% xlabel('Time (months)');
set(gca,'XTick',TVECw,'XTickLabel',[],'FontSize',11)
%..........................................................................
ylim([0 200]);
ylabel('PAR (W m^{-2})')
%..........................................................................
str = '(c)';
text(7.10,185,str,'LineStyle','none','Color','k','HorizontalAlignment','Left');
%..........................................................................
set(gca,'TickLength',[0.005, 0.005])
set(gca,'FontSize',11,'FontName','Arial');
%..........................................................................
hold off
%==========================================================================
subplot(2,4,4)
box on
grid on
hold on
%..........................................................................
plot(TIMEw,CNO3w,'Color',[ 34  94 168]/255,'LineStyle','--','LineWidth',1.5);
plot(TIMEs,CNO3s,'Color',[227  26  28]/255,'LineStyle','-','LineWidth',1.5);
%..........................................................................
xlim([TIMEw(1) TIMEw(end)]);
% xlabel('Time (months)');
set(gca,'XTick',TVECw,'XTickLabel',[],'FontSize',11)
%..........................................................................
ylim([0.50 1.30]);
ylabel('DIN (mmol m^{-3})')
set(gca,'yticklabel',num2str(get(gca,'ytick')','%.2f'))
%..........................................................................
str = '(d)';
text(7.10,1.24,str,'LineStyle','none','Color','k','HorizontalAlignment','Left');
%..........................................................................
set(gca,'TickLength',[0.005, 0.005])
set(gca,'FontSize',11,'FontName','Arial');
%..........................................................................
hold off
%==========================================================================
subplot(2,4,5)
box on
grid on
hold on
%..........................................................................
plot(TIMEw,CNUMw,'Color',[ 34  94 168]/255,'LineStyle','--','LineWidth',1.5);
plot(TIMEs,CNUMs,'Color',[227  26  28]/255,'LineStyle','-','LineWidth',1.5);
%..........................................................................
xlim([TIMEw(1) TIMEw(end)]);
xlabel('Time (days)');
set(gca,'XTick',TVECw,'XTickLabel',TOUTw,'FontSize',11)
%..........................................................................
ylim([0.1e8 4.5e8]);
ylabel('Cell abundance');
set(gca,'ytick',1e8:1e8:4e8,'yticklabel',{'1e8' '2e8' '3e8' '4e8'})
%..........................................................................
str = '(e)';
text(7.10,4.2e8,str,'LineStyle','none','Color','k','HorizontalAlignment','Left');
%..........................................................................
set(gca,'TickLength',[0.005, 0.005])
set(gca,'FontSize',11,'FontName','Arial');
%..........................................................................
hold off
%==========================================================================
subplot(2,4,6)
box on
grid on
hold on
%..........................................................................
plot(TIMEw,PHYCw,'Color',[ 34  94 168]/255,'LineStyle','--','LineWidth',1.5);
plot(TIMEs,PHYCs,'Color',[227  26  28]/255,'LineStyle','-','LineWidth',1.5);
%..........................................................................
xlim([TIMEw(1) TIMEw(end)]);
xlabel('Time (days)');
set(gca,'XTick',TVECw,'XTickLabel',TOUTw,'FontSize',11)
%..........................................................................
ylim([0.05 0.13]);
ylabel('P_C (pmol C cell^{-1})')
set(gca,'yticklabel',num2str(get(gca,'ytick')','%.2f'))
%..........................................................................
str = '(f)';
text(7.10,0.123,str,'LineStyle','none','Color','k','HorizontalAlignment','Left');
%..........................................................................
set(gca,'TickLength',[0.005, 0.005])
set(gca,'FontSize',11,'FontName','Arial');
%..........................................................................
hold off
%==========================================================================
subplot(2,4,7)
box on
grid on
hold on
%..........................................................................
plot(TIMEw,RAQNw,'Color',[ 34  94 168]/255,'LineStyle','--','LineWidth',1.5);
plot(TIMEs,RAQNs,'Color',[227  26  28]/255,'LineStyle','-','LineWidth',1.5);
%..........................................................................
xlim([TIMEw(1) TIMEw(end)]);
xlabel('Time (days)');
set(gca,'XTick',TVECw,'XTickLabel',TOUTw,'FontSize',11)
%..........................................................................
ylim([0.15 0.22]);
ylabel('Q^N (pmol N pmol C^{-1})')
set(gca,'yticklabel',num2str(get(gca,'ytick')','%.2f'))
%..........................................................................
str = '(g)';
text(7.10,0.216,str,'LineStyle','none','Color','k','HorizontalAlignment','Left');
%..........................................................................
set(gca,'TickLength',[0.005, 0.005])
set(gca,'FontSize',11,'FontName','Arial');
%..........................................................................
hold off
%==========================================================================
subplot(2,4,8)
box on
grid on
hold on
%..........................................................................
plot(TIMEw,CHLCw,'Color',[ 34  94 168]/255,'LineStyle','--','LineWidth',1.5);
plot(TIMEs,CHLCs,'Color',[227  26  28]/255,'LineStyle','-','LineWidth',1.5);
%..........................................................................
xlim([TIMEw(1) TIMEw(end)]);
xlabel('Time (days)');
set(gca,'XTick',TVECw,'XTickLabel',TOUTw,'FontSize',11)
%..........................................................................
ylim([0.28 0.55]);
ylabel('Chl:C (pg Chl pmol C^{-1})')
set(gca,'yticklabel',num2str(get(gca,'ytick')','%.2f'))
%..........................................................................
str = '(h)';
text(7.10,0.528,str,'LineStyle','none','Color','k','HorizontalAlignment','Left');
%..........................................................................
set(gca,'TickLength',[0.005, 0.005])
set(gca,'FontSize',11,'FontName','Arial');
%..........................................................................
hold off
%==========================================================================
%% Save figure:
%..........................................................................
set(gcf,'PaperPositionMode','auto');
print(gcf, '-dpng','FIG11_PartProp_DailyVariability.png');
%--------------------------------------------------------------------------
print('FIG11_PartProp_DailyVariability', '-depsc'); % Save as an EPS file
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
print(gcf, 'FIG11_PartProp_DailyVariability.pdf', '-dpdf', '-bestfit');
%==========================================================================
return
