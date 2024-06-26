
patient_list = ['432192'; '443108'; '444332';'444391';'445602';'451152'];

patient = patient_list(5,:);
master_path = ['Data/6000/',patient,'/src.tiv.RefractoryCellModel_JSF_6000/'];

beta_ci_data =  readmatrix([master_path,'betaEst.csv']);
delta_ci_data =  readmatrix([master_path,'deltaEst.csv']);
lnV0_ci_data =  readmatrix([master_path,'lnV0Est.csv']);
phi_ci_data =  readmatrix([master_path,'phiEst.csv']);
pi_ci_data =  readmatrix([master_path,'piEst.csv']);
rho_ci_data =  readmatrix([master_path,'rhoEst.csv']);

%%
CI_val = [0, 12, 25, 37, 50, 62, 75, 87, 95];

[beta_CI_min_ii,beta_CI_max_ii, beta_CI_time] = GetCIs(beta_ci_data);
[delta_CI_min_ii,delta_CI_max_ii, delta_CI_time] = GetCIs(delta_ci_data);
[lnV0_CI_min_ii,lnV0_CI_max_ii, lnV0_CI_time] = GetCIs(log10(exp(lnV0_ci_data)));
[phi_CI_min_ii,phi_CI_max_ii, phi_CI_time] = GetCIs(phi_ci_data);
[pi_CI_min_ii,pi_CI_max_ii, pi_CI_time] = GetCIs(pi_ci_data);
[rho_CI_min_ii,rho_CI_max_ii, rho_CI_time] = GetCIs(rho_ci_data);

%%
f = figure;
sgtitle(patient)
% subplot(2,1,1)
JSF_colour = 1/255*[131, 96, 150];
Tau_colour = 1/255*[237, 123, 123];
Gil_colour = 1/255*[240, 184, 110];


subplot(2,3,1)
hold on;
for jj = 1:2:length(CI_val)
    ciplot(lnV0_CI_min_ii(jj,:),lnV0_CI_max_ii(jj,:),lnV0_CI_time(jj,:),JSF_colour,(length(CI_val)-jj+1)/length(CI_val))
end
plot([0 100],[0 0],'--','color',Tau_colour,'linewidth',2.5)
plot([0 100],[log10(exp(5)) log10(exp(5))],'--','color',Tau_colour,'linewidth',2.5)
title('Prediction of $$\log_{10}(V_0)$$','fontsize',14,'Interpreter','latex')
% xlabel('Number of data used','fontsize',14,'Interpreter','latex')
% legend({'0\% CI', '25\% CI','50\% CI', '75\% CI', '95\% CI',...
%     'Prior'},...
%     'Location','northwest','NumColumns',2,'fontsize',12,'Interpreter','latex')
axis([0 7 -0.5 1.05*log10(exp(5))])



subplot(2,3,2)
hold on;
for jj = 1:2:length(CI_val)
    ciplot(beta_CI_min_ii(jj,:),beta_CI_max_ii(jj,:),beta_CI_time(jj,:),JSF_colour,(length(CI_val)-jj+1)/length(CI_val))
end
plot([0 100],[0 0],'--','color',Tau_colour,'linewidth',2.5)
plot([0 100],[10 10],'--','color',Tau_colour,'linewidth',2.5)
title('Prediction of $$\beta \,(\times 10^{-9})$$','fontsize',14,'Interpreter','latex')
% xlabel('Number of data used','fontsize',14,'Interpreter','latex')
% legend({'0\% CI', '25\% CI','50\% CI', '75\% CI', '95\% CI',...
%     'Prior'},...
%     'Location','northwest','NumColumns',2,'fontsize',12,'Interpreter','latex')
axis([0 14 -1 11])
ax = gca;
% ax.YAxis(1).Scale ="log";


subplot(2,3,3)
hold on;
for jj = 1:2:length(CI_val)
    ciplot(delta_CI_min_ii(jj,:),delta_CI_max_ii(jj,:),delta_CI_time(jj,:),JSF_colour,(length(CI_val)-jj+1)/length(CI_val))
end
plot([0 100],[1 1],'--','color',Tau_colour,'linewidth',2.5)
plot([0 100],[4 4],'--','color',Tau_colour,'linewidth',2.5)
title('Prediction of $$\delta$$','fontsize',14,'Interpreter','latex')
% xlabel('Number of data used','fontsize',14,'Interpreter','latex')
% legend({'0\% CI', '25\% CI','50\% CI', '75\% CI', '95\% CI',...
%     'Prior'},...
%     'Location','northwest','NumColumns',2,'fontsize',12,'Interpreter','latex')
axis([0 14 0.5 4.5])


subplot(2,3,4)
hold on;
for jj = 1:2:length(CI_val)
    ciplot(pi_CI_min_ii(jj,:),pi_CI_max_ii(jj,:),pi_CI_time(jj,:),JSF_colour,(length(CI_val)-jj+1)/length(CI_val))
end
plot([0 100],[350 350],'--','color',Tau_colour,'linewidth',2.5)
plot([0 100],[410 410],'--','color',Tau_colour,'linewidth',2.5)
title('Prediction of $$\pi$$','fontsize',14,'Interpreter','latex')
% xlabel('Number of data used','fontsize',14,'Interpreter','latex')
% legend({'0\% CI', '25\% CI','50\% CI', '75\% CI', '95\% CI',...
%     'Prior'},...
%     'Location','northwest','NumColumns',2,'fontsize',12,'Interpreter','latex')
axis([0 14 345 415])


subplot(2,3,5)
hold on;
for jj = 1:2:length(CI_val)
    ciplot(phi_CI_min_ii(jj,:),phi_CI_max_ii(jj,:),phi_CI_time(jj,:),JSF_colour,(length(CI_val)-jj+1)/length(CI_val))
end
plot([0 100],[0 0],'--','color',Tau_colour,'linewidth',2.5)
plot([0 100],[10 10],'--','color',Tau_colour,'linewidth',2.5)
title('Prediction of $$\phi \,(\times 10^{-5})$$','fontsize',14,'Interpreter','latex')
% xlabel('Number of data used','fontsize',14,'Interpreter','latex')
% legend({'0\% CI', '25\% CI','50\% CI', '75\% CI', '95\% CI',...
%     'Prior'},...
%     'Location','northwest','NumColumns',2,'fontsize',12,'Interpreter','latex')
axis([0 14 -0.5 10.5])


subplot(2,3,6)
hold on;
for jj = 1:2:length(CI_val)
    ciplot(rho_CI_min_ii(jj,:),rho_CI_max_ii(jj,:),rho_CI_time(jj,:),JSF_colour,(length(CI_val)-jj+1)/length(CI_val))
end
plot([0 100],[0 0],'--','color',Tau_colour,'linewidth',2.5)
plot([0 100],[0.024 0.024],'--','color',Tau_colour,'linewidth',2.5)
title('Prediction of $$\rho$$','fontsize',14,'Interpreter','latex')
% xlabel('Number of data used','fontsize',14,'Interpreter','latex')
% legend({'0\% CI', '25\% CI','50\% CI', '75\% CI', '95\% CI',...
%     'Prior'},...
%     'Location','northwest','NumColumns',2,'fontsize',12,'Interpreter','latex')
axis([0 14 -0.005 0.029])

%%
function [CI_min_ii,CI_max_ii, CI_time] = GetCIs(ci_data)
    ci_time = ci_data(:,2);
    ci_val =  ci_data(:,3);
    ci_min = ci_data(:,4);
    ci_max = ci_data(:,5);

    % CIs: [ 0, 12.5, 25, 37.5, 50, 62.5, 75, 87.5, 95 ]
    CI_val = [0, 12, 25, 37, 50, 62, 75, 87, 95];

    CI_min_ii = zeros(length(CI_val),size(ci_time,1)/length(CI_val));
    CI_max_ii = zeros(length(CI_val),size(ci_time,1)/length(CI_val));
    CI_time = zeros(length(CI_val),size(ci_time,1)/length(CI_val));

    for jj = 1:length(CI_val)
        CI_min_ii(jj,:) = ci_min(jj:9:end);
        CI_max_ii(jj,:) = ci_max(jj:9:end);
        CI_time(jj,:) = ci_time(jj:9:end);
    end
end