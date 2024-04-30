clear all;


master_path = ['sirs_model.SIRS_JSF_particle=2000_intMag=1.0/'];
JSF_colour = 1/255*[131, 96, 150];

%%
tW = 3;
f1 = figure;

subplot(1,3,1)
beta_dist = readmatrix([master_path,'betaCoef_distribution.csv']);
violinplot(beta_dist,'beta', 'HalfViolin','right','Showdata',false,...
    'ShowBox', false,'ShowMean',false,'ShowMedian', false,...
    'ViolinColor',JSF_colour,'ViolinAlpha',0.9);
view(90,-90)
ylabel('$\beta$','FontSize',18,'Interpreter','latex')
% xlabel(,'FontSize',22,'Interpreter','latex')
plot([1 2],[0 0], 'r--',LineWidth=2)
plot([1 2],[4 4], 'r--',LineWidth=2)
plot([1 2],[2 2], 'k:',LineWidth=tW)

axis([1 1.32 -0.25 4.25])

subplot(1,3,2)
rho_dist = readmatrix([master_path,'gammaCoef_distribution.csv']);
violinplot(rho_dist,'beta', 'HalfViolin','right','Showdata',false, ...
    'ShowBox', false,'ShowMean',false,'ShowMedian', false,...
    'ViolinColor',JSF_colour,'ViolinAlpha',0.9);
view(90,-90)
ylabel('$\gamma$','FontSize',18,'Interpreter','latex')
plot([1 2],[0 0], 'r--',LineWidth=2)
plot([1 2],[3 3], 'r--',LineWidth=2)
plot([1 2],[1 1], 'k:',LineWidth=tW)
axis([1 1.32 -0.25 3.25])


subplot(1,3,3)
pi_dist = readmatrix([master_path,'omegaCoef_distribution.csv']);
violinplot(pi_dist,'beta', 'HalfViolin','right','Showdata',false, ...
    'ShowBox', false,'ShowMean',false,'ShowMedian', false,...
    'ViolinColor',JSF_colour,'ViolinAlpha',0.9);
view(90,-90)
ylabel('$\omega$','FontSize',18,'Interpreter','latex')
plot([1 2],[0 0], 'r--',LineWidth=2)
plot([1 2],[2 2], 'r--',LineWidth=2)
plot([1 2],[1 1], 'k:',LineWidth= tW)

axis([1 1.32 -0.25 2.75])


%%
