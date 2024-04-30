clear all;

patient_list = ['432192'; '443108'; '444332';'444391';'445602';'451152'];

raw_data_mat = zeros(length(patient_list),14,2);
ci_data_mat = zeros(length(patient_list),153,5);
ext_data_mat = zeros(length(patient_list),16,2);

fc_time = 15;
int_time = 15;

beta_ci_mat =  zeros(length(patient_list),135,5);
delta_ci_mat =  zeros(length(patient_list),135,5);
lnV0_ci_mat =  zeros(length(patient_list),135,5);
phi_ci_mat =  zeros(length(patient_list),135,5);
pi_ci_mat =  zeros(length(patient_list),135,5);
rho_ci_mat =  zeros(length(patient_list),135,5);

for p_id=1:length(patient_list)
    
    patient = patient_list(p_id,:);

    master_path = [patient,'/src.tiv.RefractoryCellModel_JSF_2000/'];
    raw_data = readmatrix(['PatientData/CSVs/',patient,'.csv']);
    raw_data_mat(p_id,:,:) = raw_data;
    ci_data =  readmatrix([master_path,'plt_df.csv']);
    ci_data_mat(p_id,:,:) = ci_data;

    ext_data = readmatrix([master_path,'ext_prdf.csv']);
    ext_data(fc_time,:) = [];
    ext_data_mat(p_id,:,:) = ext_data;
    
    beta_ci_data =  readmatrix([master_path,'betaEst.csv']);
    delta_ci_data =  readmatrix([master_path,'deltaEst.csv']);
    lnV0_ci_data =  readmatrix([master_path,'lnV0Est.csv']);
    phi_ci_data =  readmatrix([master_path,'phiEst.csv']);
    pi_ci_data =  readmatrix([master_path,'piEst.csv']);
    rho_ci_data =  readmatrix([master_path,'rhoEst.csv']);
    
    beta_ci_mat(p_id,:,:) = beta_ci_data(:,1:5);
    delta_ci_mat(p_id,:,:) = delta_ci_data(:,1:5);
    lnV0_ci_mat(p_id,:,:) = lnV0_ci_data(:,1:5);
    phi_ci_mat(p_id,:,:) = phi_ci_data(:,1:5);
    pi_ci_mat(p_id,:,:) = pi_ci_data(:,1:5);
    rho_ci_mat(p_id,:,:) = rho_ci_data(:,1:5);
end
%%

beta_est = zeros(length(patient_list),5);
delta_est = zeros(length(patient_list),5);
lnV0_est = zeros(length(patient_list),5);
phi_est = zeros(length(patient_list),5);
pi_est = zeros(length(patient_list),5);
rho_est = zeros(length(patient_list),5);

for p_id=1:length(patient_list)
    % mean
    beta_est(p_id,1) = beta_ci_mat(p_id,end-8,5);
    % 95% lower
    beta_est(p_id,2) = beta_ci_mat(p_id,end-8,5)-beta_ci_mat(p_id,end,4);
    % 95% upper
    beta_est(p_id,3) = beta_ci_mat(p_id,end,5)-beta_ci_mat(p_id,end-8,5);
    % 50% lower
    beta_est(p_id,4) = beta_ci_mat(p_id,end-8,5)-beta_ci_mat(p_id,end-4,4);
    % 50% upper
    beta_est(p_id,5) = beta_ci_mat(p_id,end-4,5)-beta_ci_mat(p_id,end-8,5);
    
    delta_est(p_id,1) = delta_ci_mat(p_id,end-8,5);
    delta_est(p_id,2) = delta_ci_mat(p_id,end-8,5)-delta_ci_mat(p_id,end,4);
    delta_est(p_id,3) = delta_ci_mat(p_id,end,5)-delta_ci_mat(p_id,end-8,5);
    delta_est(p_id,4) = delta_ci_mat(p_id,end-8,5)-delta_ci_mat(p_id,end-4,4);
    delta_est(p_id,5) = delta_ci_mat(p_id,end-4,5)-delta_ci_mat(p_id,end-8,5);
    
    
    lnV0_est(p_id,1) = log10(exp(lnV0_ci_mat(p_id,end-8,5)));
    lnV0_est(p_id,2) = log10(exp(lnV0_ci_mat(p_id,end-8,5)-lnV0_ci_mat(p_id,end,4)));
    lnV0_est(p_id,3) = log10(exp(lnV0_ci_mat(p_id,end,5)-lnV0_ci_mat(p_id,end-8,5)));
    lnV0_est(p_id,4) = log10(exp(lnV0_ci_mat(p_id,end-8,5)-lnV0_ci_mat(p_id,end-4,4)));
    lnV0_est(p_id,5) = log10(exp(lnV0_ci_mat(p_id,end-4,5)-lnV0_ci_mat(p_id,end-8,5)));
    
    phi_est(p_id,1) = phi_ci_mat(p_id,end-8,5);
    phi_est(p_id,2) = phi_ci_mat(p_id,end-8,5)-phi_ci_mat(p_id,end,4);
    phi_est(p_id,3) = phi_ci_mat(p_id,end,5)-phi_ci_mat(p_id,end-8,5);
    phi_est(p_id,4) = phi_ci_mat(p_id,end-8,5)-phi_ci_mat(p_id,end-4,4);
    phi_est(p_id,5) = phi_ci_mat(p_id,end-4,5)-phi_ci_mat(p_id,end-8,5);
    
    pi_est(p_id,1) = pi_ci_mat(p_id,end-8,5);
    pi_est(p_id,2) = pi_ci_mat(p_id,end-8,5)-pi_ci_mat(p_id,end,4);
    pi_est(p_id,3) = pi_ci_mat(p_id,end,5)-pi_ci_mat(p_id,end-8,5);
    pi_est(p_id,4) = pi_ci_mat(p_id,end-8,5)-pi_ci_mat(p_id,end-4,4);
    pi_est(p_id,5) = pi_ci_mat(p_id,end-4,5)-pi_ci_mat(p_id,end-8,5);
    
    rho_est(p_id,1) = rho_ci_mat(p_id,end-8,5);
    rho_est(p_id,2) = rho_ci_mat(p_id,end-8,5)-rho_ci_mat(p_id,end,4);
    rho_est(p_id,3) = rho_ci_mat(p_id,end,5)-rho_ci_mat(p_id,end-8,5);
    rho_est(p_id,4) = rho_ci_mat(p_id,end-8,5)-rho_ci_mat(p_id,end-4,4);
    rho_est(p_id,5) = rho_ci_mat(p_id,end-4,5)-rho_ci_mat(p_id,end-8,5);
    
end

%%
JSF_colour = 1/255*[131, 96, 150];
Tau_colour = 1/255*[237, 123, 123];
Gil_colour = 1/255*[240, 184, 110];
myblue_colour = 1/255*[95, 130, 237];

f=figure;
hold on;
errorbar(1:6,lnV0_est(:,1),lnV0_est(:,2),lnV0_est(:,3),'o','linewidth',2.5,'color',Gil_colour)
errorbar(1:6,lnV0_est(:,1),lnV0_est(:,4),lnV0_est(:,5),'o','linewidth',2.5,'color',Tau_colour)
scatter(1:6,lnV0_est(:,1), 100 ,JSF_colour,'filled')

ylabel('$$\log_{10}(V_0)$$','fontsize',16,'Interpreter','latex')
xlabel('Patient','fontsize',16,'Interpreter','latex')

plot([0 7],[0 0],':','color',[0.5 0.5 0.5],'linewidth',2.5)
plot([0 7],[log10(exp(5)) log10(exp(5))],':','color',[0.5 0.5 0.5],'linewidth',2.5)
axis([0 7 -0.5 1.05*log10(exp(6))])
xticks([1:6])
xticklabels({patient_list})
xtickangle(45)
hold off;
f.Position = [391   324   303   217];
exportgraphics(gcf,'Figures/covid_fit_lnV0.png','Resolution',1000)



f=figure;
hold on;
errorbar(1:6,beta_est(:,1),beta_est(:,2),beta_est(:,3),'o','linewidth',2.5,'color',Gil_colour)
errorbar(1:6,beta_est(:,1),beta_est(:,4),beta_est(:,5),'o','linewidth',2.5,'color',Tau_colour)
scatter(1:6,beta_est(:,1), 100 ,JSF_colour,'filled')

ylabel('$$\beta \,(\times 10^{-9})$$','fontsize',16,'Interpreter','latex')
xlabel('Patient','fontsize',16,'Interpreter','latex')

plot([0 7],[0 0],':','color',[0.5 0.5 0.5],'linewidth',2.5)
plot([0 7],[10 10],':','color',[0.5 0.5 0.5],'linewidth',2.5)
axis([0 7 -1.25 11.5])
xticks([1:6])
xticklabels({patient_list})
xtickangle(45)
hold off;
f.Position = [391   324   303   217];
exportgraphics(gcf,'Figures/covid_fit_beta.png','Resolution',1000)


f=figure;
hold on;
errorbar(1:6,delta_est(:,1),delta_est(:,2),delta_est(:,3),'o','linewidth',2.5,'color',Gil_colour)
errorbar(1:6,delta_est(:,1),delta_est(:,4),delta_est(:,5),'o','linewidth',2.5,'color',Tau_colour)
scatter(1:6,delta_est(:,1), 100 ,JSF_colour,'filled')

ylabel('$$\delta$$','fontsize',16,'Interpreter','latex')
xlabel('Patient','fontsize',16,'Interpreter','latex')

plot([0 7],[1 1],':','color',[0.5 0.5 0.5],'linewidth',2.5)
plot([0 7],[4 4],':','color',[0.5 0.5 0.5],'linewidth',2.5)
axis([0 7 0.5 4.5])
xticks([1:6])
xticklabels({patient_list})
xtickangle(45)
hold off;
f.Position = [391   324   303   217];
exportgraphics(gcf,'Figures/covid_fit_delta.png','Resolution',1000)


f=figure;
hold on;
errorbar(1:6,pi_est(:,1),pi_est(:,2),pi_est(:,3),'o','linewidth',2.5,'color',Gil_colour)
errorbar(1:6,pi_est(:,1),pi_est(:,4),pi_est(:,5),'o','linewidth',2.5,'color',Tau_colour)
scatter(1:6,pi_est(:,1), 100 ,JSF_colour,'filled')

ylabel('$$\pi$$','fontsize',16,'Interpreter','latex')
xlabel('Patient','fontsize',16,'Interpreter','latex')

plot([0 7],[350 350],':','color',[0.5 0.5 0.5],'linewidth',2.5)
plot([0 100],[410 410],':','color',[0.5 0.5 0.5],'linewidth',2.5)
axis([0 7 342.5 417.5])
xticks([1:6])
xticklabels({patient_list})
xtickangle(45)
hold off;
f.Position = [391   324   303   217];
exportgraphics(gcf,'Figures/covid_fit_pi.png','Resolution',1000)


f=figure;
hold on;
errorbar(1:6,phi_est(:,1),phi_est(:,2),phi_est(:,3),'o','linewidth',2.5,'color',Gil_colour)
errorbar(1:6,phi_est(:,1),phi_est(:,4),phi_est(:,5),'o','linewidth',2.5,'color',Tau_colour)
scatter(1:6,phi_est(:,1), 100 ,JSF_colour,'filled')

ylabel('$$\phi \,(\times 10^{-5})$$','fontsize',16,'Interpreter','latex')
xlabel('Patient','fontsize',16,'Interpreter','latex')

plot([0 7],[0 0],':','color',[0.5 0.5 0.5],'linewidth',2.5)
plot([0 7],[10 10],':','color',[0.5 0.5 0.5],'linewidth',2.5)
axis([0 7 -1.0 11])
xticks([1:6])
xticklabels({patient_list})
xtickangle(45)
hold off;
f.Position = [391   324   303   217];
exportgraphics(gcf,'Figures/covid_fit_phi.png','Resolution',1000)


f=figure;
hold on;
errorbar(1:6,rho_est(:,1),rho_est(:,2),rho_est(:,3),'o','linewidth',2.5,'color',Gil_colour)
errorbar(1:6,rho_est(:,1),rho_est(:,4),rho_est(:,5),'o','linewidth',2.5,'color',Tau_colour)
scatter(1:6,rho_est(:,1), 100 ,JSF_colour,'filled')

ylabel('$$\rho$$','fontsize',16,'Interpreter','latex')
xlabel('Patient','fontsize',16,'Interpreter','latex')
plot([0 7],[0 0],':','color',[0.5 0.5 0.5],'linewidth',2.5)
plot([0 7],[0.024 0.024],':','color',[0.5 0.5 0.5],'linewidth',2.5)
hold off;
xticks([1:6])
xticklabels({patient_list})
xtickangle(45)
axis([0 7 -0.005 0.029])
f.Position = [391   324   303   217];
exportgraphics(gcf,'Figures/covid_fit_rho.png','Resolution',1000)


legend('95\% CI', '50\% CI', 'Mean', 'Prior','NumColumns',2,'fontsize',15,'Interpreter','latex','location','north')
axis([-1 0 -1 0.5])
title('')
axis off;
exportgraphics(gcf,'Figures/covid_fit_legend.png','Resolution',1000)


%%
