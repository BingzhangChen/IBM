close all
clear all
format long g
%==========================================================================

%==========================================================================
%% Load data from the 5000 super-individuals simulation:
%==========================================================================
path = 'G:\My Drive\RESEARCH\PROJECT LEVERHULME\CODE\1D-Model\OUTPUTS\v033_RUN049_05K\';
%..........................................................................
EUL  = fullfile(path,'Euler.nc');
%..........................................................................
DAYS = double(nc_varget(EUL,'Day'));
tfin = length(DAYS) - 1;
tini = tfin - 364;
%..........................................................................
ZDEP = double(nc_varget(EUL,'Z_w'));
ZDIF = diff(ZDEP)';
ZDIF = repmat(ZDIF,365,1);
%..........................................................................
DINN = double(nc_varget(EUL,'NO3'));
%..........................................................................
PHYN = double(nc_varget(EUL,'PN'));
PHYC = double(nc_varget(EUL,'PC'));
PCHL = double(nc_varget(EUL,'CHL'));
%..........................................................................
NPP = double(nc_varget(EUL,'NPP'));
%==========================================================================
% Compute integrated values
DINN_05K = sum((DINN(tini:tfin,:) .* ZDIF),2);
PHYC_05K = sum((PHYC(tini:tfin,:) .* ZDIF),2);
PHYN_05K = sum((PHYN(tini:tfin,:) .* ZDIF),2);
PCHL_05K = sum((PCHL(tini:tfin,:) .* ZDIF),2);
NPP_05K  = sum((NPP (tini:tfin,:) .* ZDIF),2);
%==========================================================================
% Load Rao index
RAO_05K  = xlsread('Rao05K.csv');
%==========================================================================
return
%==========================================================================
%% Load data from the 10000 super-individuals simulation:
%==========================================================================
path  = 'G:\My Drive\RESEARCH\PROJECT LEVERHULME\CODE\1D-Model\OUTPUTS\v033_RUN050_10K\';
%..........................................................................
EUL  = fullfile(path,'Euler.nc');
%..........................................................................
DINN = double(nc_varget(EUL,'NO3'));
%..........................................................................
PHYN = double(nc_varget(EUL,'PN'));
PHYC = double(nc_varget(EUL,'PC'));
PCHL = double(nc_varget(EUL,'CHL'));
%..........................................................................
NPP  = double(nc_varget(EUL,'NPP'));
%==========================================================================
% Compute integrated values
DINN_10K = sum((DINN(tini:tfin,:) .* ZDIF),2);
PHYC_10K = sum((PHYC(tini:tfin,:) .* ZDIF),2);
PHYN_10K = sum((PHYN(tini:tfin,:) .* ZDIF),2);
PCHL_10K = sum((PCHL(tini:tfin,:) .* ZDIF),2);
NPP_10K  = sum((NPP (tini:tfin,:) .* ZDIF),2);
%==========================================================================
% Load Rao index
RAO_10K  = xlsread('Rao10K.csv');
%==========================================================================

%==========================================================================
%% Load data from the 20000 super-individuals simulation:
%==========================================================================
path  = 'G:\My Drive\RESEARCH\PROJECT LEVERHULME\CODE\1D-Model\OUTPUTS\v033_RUN048_20K\';
%..........................................................................
EUL  = fullfile(path,'Euler.nc');
%..........................................................................
DINN = double(nc_varget(EUL,'NO3'));
%..........................................................................
PHYN = double(nc_varget(EUL,'PN'));
PHYC = double(nc_varget(EUL,'PC'));
PCHL = double(nc_varget(EUL,'CHL'));
%..........................................................................
NPP  = double(nc_varget(EUL,'NPP'));
%==========================================================================
% Compute integrated values
DINN_20K = sum((DINN(tini:tfin,:) .* ZDIF),2);
PHYC_20K = sum((PHYC(tini:tfin,:) .* ZDIF),2);
PHYN_20K = sum((PHYN(tini:tfin,:) .* ZDIF),2);
PCHL_20K = sum((PCHL(tini:tfin,:) .* ZDIF),2);
NPP_20K  = sum((NPP (tini:tfin,:) .* ZDIF),2);
%==========================================================================
% Load Rao index
RAO_20K  = xlsread('Rao20K.csv');
%==========================================================================

%==========================================================================
%% Load data from the 20000 super-individuals simulation:
%==========================================================================
path  = 'G:\My Drive\RESEARCH\PROJECT LEVERHULME\CODE\1D-Model\OUTPUTS\v033_RUN051_50K\';
%..........................................................................
EUL  = fullfile(path,'Euler.nc');
%..........................................................................
DINN = double(nc_varget(EUL,'NO3'));
%..........................................................................
PHYN = double(nc_varget(EUL,'PN'));
PHYC = double(nc_varget(EUL,'PC'));
PCHL = double(nc_varget(EUL,'CHL'));
%..........................................................................
NPP  = double(nc_varget(EUL,'NPP'));
%==========================================================================
% Compute integrated values
DINN_50K = sum((DINN(tini:tfin,:) .* ZDIF),2);
PHYC_50K = sum((PHYC(tini:tfin,:) .* ZDIF),2);
PHYN_50K = sum((PHYN(tini:tfin,:) .* ZDIF),2);
PCHL_50K = sum((PCHL(tini:tfin,:) .* ZDIF),2);
NPP_50K  = sum((NPP (tini:tfin,:) .* ZDIF),2);
%==========================================================================
% Load Rao index
RAO_50K  = xlsread('Rao50K.csv');
%==========================================================================

%==========================================================================
%% Define time vector:
%==========================================================================
DAYS = 1:365;
Tvec = [1,32,60,91,121,152,182,213,244,274,305,335];
Tout = {'J' 'F' 'M' 'A' 'M' 'J' 'J' 'A' 'S' 'O' 'N' 'D'};
%==========================================================================


%==========================================================================
%% FIGURE: Total
%==========================================================================
make_it_tight = true;
% m = distance between subplots
% n = distance to bottom and top
% p = distance to left and right
subplot = @(m,n,p) subtightplot (m, n, p, [0.05 0.08], [0.10 0.04], [0.07 0.02]);
if ~make_it_tight,  clear subplot;  end
%..........................................................................
Fig = figure;
set(gcf, 'Color','white');
set(Fig, 'Position',[80, 50,1000,500]);
%==========================================================================
subplot(2,3,1);
%..........................................................................
hold on
box on
%..........................................................................
plot(DAYS,DINN_05K,'Color',[215 025 028]/255,'LineStyle','-','LineWidth',2);
plot(DAYS,DINN_10K,'Color',[253 174 097]/255,'LineStyle','--','LineWidth',2);
plot(DAYS,DINN_20K,'Color',[171 221 164]/255,'LineStyle','-.','LineWidth',2);
plot(DAYS,DINN_50K,'Color',[ 43 131 186]/255,'LineStyle',':','LineWidth',2);
%..........................................................................
xlim([1 365]);
% xlabel('Time (days)')
set(gca,'XTick',Tvec,'XTickLabel',[],'FontSize',10)
%..........................................................................
ylim([325 360]);
ylabel('DIN (mmol m^{-2})');
% set(gca,'YTick',[-250:50:0],'YTickLabel',[-250:50:0],'FontSize',10)
%..........................................................................
str = '(a)';
text(15,358,str,'LineStyle','none','Color','k','HorizontalAlignment','Left');
%..........................................................................
set(gca,'FontSize',11,'FontName','Times'); 
%==========================================================================
subplot(2,3,2);
%..........................................................................
hold on
box on
%..........................................................................
plot(DAYS,PHYC_05K,'Color',[215 025 028]/255,'LineStyle','-','LineWidth',2);
plot(DAYS,PHYC_10K,'Color',[253 174 097]/255,'LineStyle','--','LineWidth',2);
plot(DAYS,PHYC_20K,'Color',[171 221 164]/255,'LineStyle','-.','LineWidth',2);
plot(DAYS,PHYC_50K,'Color',[ 43 131 186]/255,'LineStyle',':','LineWidth',2);
%..........................................................................
xlim([1 365]);
% xlabel('Time (months)')
set(gca,'XTick',Tvec,'XTickLabel',[],'FontSize',10)
%..........................................................................
ylim([90 220]);
ylabel('P_C (mmol C m^{-2})');
% set(gca,'YTick',[-250:50:0],'YTickLabel',[-250:50:0],'FontSize',10)
%..........................................................................
str = '(b)';
text(15,213,str,'LineStyle','none','Color','k','HorizontalAlignment','Left');
%..........................................................................
% legend(' 1000 super-individuals',' 5000 super-individuals','10000 super-individuals')
% legend box off
%..........................................................................
set(gca,'FontSize',11,'FontName','Times'); 
%==========================================================================
subplot(2,3,3);
%..........................................................................
hold on
box on
%..........................................................................
plot(DAYS,PHYN_05K,'Color',[215 025 028]/255,'LineStyle','-','LineWidth',2);
plot(DAYS,PHYN_10K,'Color',[253 174 097]/255,'LineStyle','--','LineWidth',2);
plot(DAYS,PHYN_20K,'Color',[171 221 164]/255,'LineStyle','-.','LineWidth',2);
plot(DAYS,PHYN_50K,'Color',[ 43 131 186]/255,'LineStyle',':','LineWidth',2);
%..........................................................................
xlim([1 365]);
% xlabel('Time (months)')
set(gca,'XTick',Tvec,'XTickLabel',[],'FontSize',10)
%..........................................................................
ylim([15 35]);
ylabel('P_N (mmol N m^{-2})');
% set(gca,'YTick',[-250:50:0],'YTickLabel',[-250:50:0],'FontSize',10)
%..........................................................................
str = '(c)';
text(15,34,str,'LineStyle','none','Color','k','HorizontalAlignment','Left');
%..........................................................................
legend('  5 K',' 10 K',' 20 K',' 50 K','Location','best',...
       'FontSize',12);
legend box off
%..........................................................................
set(gca,'FontSize',11,'FontName','Times'); 
%==========================================================================
subplot(2,3,4);
%..........................................................................
hold on
box on
%..........................................................................
plot(DAYS,PCHL_05K,'Color',[215 025 028]/255,'LineStyle','-','LineWidth',2);
plot(DAYS,PCHL_10K,'Color',[253 174 097]/255,'LineStyle','--','LineWidth',2);
plot(DAYS,PCHL_20K,'Color',[171 221 164]/255,'LineStyle','-.','LineWidth',2);
plot(DAYS,PCHL_50K,'Color',[ 43 131 186]/255,'LineStyle',':','LineWidth',2);
%..........................................................................
xlim([1 365]);
xlabel('Time (months)')
set(gca,'XTick',Tvec,'XTickLabel',Tout,'FontSize',10)
%..........................................................................
ylim([35 80]);
ylabel('CHL (mg Chl m^{-2})');
% set(gca,'YTick',[-250:50:0],'YTickLabel',[-250:50:0],'FontSize',10)
%..........................................................................
str = '(d)';
text(15,77,str,'LineStyle','none','Color','k','HorizontalAlignment','Left');
%..........................................................................
set(gca,'FontSize',11,'FontName','Times'); 
%==========================================================================
subplot(2,3,5);
%..........................................................................
hold on
box on
%..........................................................................
plot(DAYS,NPP_05K,'Color',[215 025 028]/255,'LineStyle','-','LineWidth',2);
plot(DAYS,NPP_10K,'Color',[253 174 097]/255,'LineStyle','--','LineWidth',2);
plot(DAYS,NPP_20K,'Color',[171 221 164]/255,'LineStyle','-.','LineWidth',2);
plot(DAYS,NPP_50K,'Color',[ 43 131 186]/255,'LineStyle',':','LineWidth',2);
%..........................................................................
xlim([1 365]);
xlabel('Time (months)')
set(gca,'XTick',Tvec,'XTickLabel',Tout,'FontSize',10)
%..........................................................................
ylim([50 350]);
ylabel('NPP (mg C m^{-2} d^{-1})');
% set(gca,'YTick',[-250:50:0],'YTickLabel',[-250:50:0],'FontSize',10)
%..........................................................................
str = '(e)';
text(15,333,str,'LineStyle','none','Color','k','HorizontalAlignment','Left');
%..........................................................................
set(gca,'FontSize',11,'FontName','Times'); 
%==========================================================================
subplot(2,3,6);
%..........................................................................
hold on
box on
%..........................................................................
plot(DAYS,RAO_05K,'Color',[215 025 028]/255,'LineStyle','-','LineWidth',2);
plot(DAYS,RAO_10K,'Color',[253 174 097]/255,'LineStyle','--','LineWidth',2);
plot(DAYS,RAO_20K,'Color',[171 221 164]/255,'LineStyle','-.','LineWidth',2);
plot(DAYS,RAO_50K,'Color',[ 43 131 186]/255,'LineStyle',':','LineWidth',2);
%..........................................................................
xlim([1 365]);
xlabel('Time (months)')
set(gca,'XTick',Tvec,'XTickLabel',Tout,'FontSize',10)
%..........................................................................
ylim([0.15 0.25]);
ylabel('Rao index');
% set(gca,'YTick',[-250:50:0],'YTickLabel',[-250:50:0],'FontSize',10)
%..........................................................................
str = '(f)';
text(15,0.245,str,'LineStyle','none','Color','k','HorizontalAlignment','Left');
%..........................................................................
set(gca,'FontSize',11,'FontName','Times'); 
%==========================================================================
%% Save figure:
%..........................................................................
set(gcf,'PaperPositionMode','auto');
print(gcf, '-dpng', 'FIG12_SIMCOMPARISON.png');
%--------------------------------------------------------------------------
print('FIG12_SIMCOMPARISON', '-depsc'); % Save as an EPS file
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
print(gcf, 'FIG12_SIMCOMPARISON.pdf', '-dpdf', '-bestfit');
%==========================================================================
return