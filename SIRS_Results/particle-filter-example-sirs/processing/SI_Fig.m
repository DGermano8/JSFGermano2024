master_path = '2000_100/sirs_model.SIRS_JSF_particle=2000_intMag=1/';
% master_path = 'sirs_model.SIRS_JSF_particle=6000_intMag=0.875/';
% master_path = 'sirs_model.SIRS_JSF_particle=6000_intMag=0.75/';

% raw_data = readmatrix(['sim_df.csv']);
omega_ci_data =  readmatrix([master_path,'omegaEst.csv']);
beta_ci_data =  readmatrix([master_path,'betaEst.csv']);
gamma_ci_data =  readmatrix([master_path,'gammaEst.csv']);

%%

% mean
omega_est(1) = omega_ci_data(end-8,5);
% 95% lower
omega_est(2) = omega_ci_data(end-8,5)-omega_ci_data(end,4);
% 95% upper
omega_est(3) = omega_ci_data(end,5)-omega_ci_data(end-8,5);
% 50% lower
omega_est(4) = omega_ci_data(end-8,5)-omega_ci_data(end-4,4);
% 50% upper
omega_est(5) = omega_ci_data(end-4,5)-omega_ci_data(end-8,5);

beta_est(1) = beta_ci_data(end-8,5);
beta_est(2) = beta_ci_data(end-8,5)-beta_ci_data(end,4);
beta_est(3) = beta_ci_data(end,5)-beta_ci_data(end-8,5);
beta_est(4) = beta_ci_data(end-8,5)-beta_ci_data(end-4,4);
beta_est(5) = beta_ci_data(end-4,5)-beta_ci_data(end-8,5);

gamma_est(1) = gamma_ci_data(end-8,5);
gamma_est(2) = gamma_ci_data(end-8,5)-gamma_ci_data(end,4);
gamma_est(3) = gamma_ci_data(end,5)-gamma_ci_data(end-8,5);
gamma_est(4) = gamma_ci_data(end-8,5)-gamma_ci_data(end-4,4);
gamma_est(5) = gamma_ci_data(end-4,5)-gamma_ci_data(end-8,5);
    
%%
% subplot(2,1,1)
JSF_colour = 1/255*[131, 96, 150];
Tau_colour = 1/255*[237, 123, 123];
Gil_colour = 1/255*[240, 184, 110];

f=figure;
hold on;
errorbar(1,omega_est(1),omega_est(2),omega_est(3),'o','linewidth',2.5,'color',Gil_colour)
errorbar(1,omega_est(1),omega_est(4),omega_est(5),'o','linewidth',2.5,'color',Tau_colour)
scatter(1,omega_est(1), 100 ,JSF_colour,'filled')
plot([0 2],[0 0],':','color',[0.5 0.5 0.5],'linewidth',2.5)
plot([0 2],[2 2],':','color',[0.5 0.5 0.5],'linewidth',2.5)
axis([0.5 1.5 -0.1 2.1])

f=figure;
hold on;
errorbar(1,beta_est(1),beta_est(2),beta_est(3),'o','linewidth',2.5,'color',Gil_colour)
errorbar(1,beta_est(1),beta_est(4),beta_est(5),'o','linewidth',2.5,'color',Tau_colour)
scatter(1,beta_est(1), 100 ,JSF_colour,'filled')
plot([0 2],[0 0],':','color',[0.5 0.5 0.5],'linewidth',2.5)
plot([0 2],[4 4],':','color',[0.5 0.5 0.5],'linewidth',2.5)
axis([0.5 1.5 -0.1 4.1])

f=figure;
hold on;
errorbar(1,gamma_est(1),gamma_est(2),gamma_est(3),'o','linewidth',2.5,'color',Gil_colour)
errorbar(1,gamma_est(1),gamma_est(4),gamma_est(5),'o','linewidth',2.5,'color',Tau_colour)
scatter(1,gamma_est(1), 100 ,JSF_colour,'filled')
plot([0 2],[0 0],':','color',[0.5 0.5 0.5],'linewidth',2.5)
plot([0 2],[3 3],':','color',[0.5 0.5 0.5],'linewidth',2.5)
axis([0.5 1.5 -0.1 3.1])



