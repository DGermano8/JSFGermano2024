clear all;

patient_list = ['432192'; '443108'; '444332';'444391';'445602';'451152'];

raw_data_mat = zeros(length(patient_list),14,2);
ci_data_mat = zeros(length(patient_list),207,5);
ext_data_mat = zeros(length(patient_list),22,2);

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

    master_path = ['Data/',patient,'/src.tiv.RefractoryCellModel_JSF_2000/'];
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

figure;
subplot(2,6,7)
hold on;
errorbar(1:6,lnV0_est(:,1),lnV0_est(:,2),lnV0_est(:,3),'o','linewidth',2.5)
errorbar(1:6,lnV0_est(:,1),lnV0_est(:,4),lnV0_est(:,5),'o','linewidth',2.5)
title('$$\log_{10}(V_0)$$','fontsize',14,'Interpreter','latex')
plot([0 7],[0 0],':','color',Tau_colour,'linewidth',2.5)
plot([0 7],[log10(exp(5)) log10(exp(5))],':','color',Tau_colour,'linewidth',2.5)
axis([0 7 -0.5 1.05*log10(exp(6))])
xticks([1:6])
xticklabels({patient_list})
xtickangle(45)
hold off;


subplot(2,6,8)
hold on;
errorbar(1:6,beta_est(:,1),beta_est(:,2),beta_est(:,3),'o','linewidth',2.5)
errorbar(1:6,beta_est(:,1),beta_est(:,4),beta_est(:,5),'o','linewidth',2.5)
title('$$\beta \,(\times 10^{-9})$$','fontsize',14,'Interpreter','latex')
plot([0 7],[0 0],':','color',Tau_colour,'linewidth',2.5)
plot([0 7],[10 10],':','color',Tau_colour,'linewidth',2.5)
axis([0 7 -1.25 11.5])
xticks([1:6])
xticklabels({patient_list})
xtickangle(45)
hold off;

subplot(2,6,9)
hold on;
errorbar(1:6,delta_est(:,1),delta_est(:,2),delta_est(:,3),'o','linewidth',2.5)
errorbar(1:6,delta_est(:,1),delta_est(:,4),delta_est(:,5),'o','linewidth',2.5)
title('$$\delta$$','fontsize',14,'Interpreter','latex')
plot([0 7],[1 1],':','color',Tau_colour,'linewidth',2.5)
plot([0 7],[4 4],':','color',Tau_colour,'linewidth',2.5)
axis([0 7 0.5 4.5])
xticks([1:6])
xticklabels({patient_list})
xtickangle(45)
hold off;

subplot(2,6,10)
hold on;
errorbar(1:6,pi_est(:,1),pi_est(:,2),pi_est(:,3),'o','linewidth',2.5)
errorbar(1:6,pi_est(:,1),pi_est(:,4),pi_est(:,5),'o','linewidth',2.5)
title('$$\pi$$','fontsize',14,'Interpreter','latex')
plot([0 7],[350 350],':','color',Tau_colour,'linewidth',2.5)
plot([0 100],[410 410],':','color',Tau_colour,'linewidth',2.5)
axis([0 7 342.5 417.5])
xticks([1:6])
xticklabels({patient_list})
xtickangle(45)
hold off;

subplot(2,6,11)
hold on;
errorbar(1:6,phi_est(:,1),phi_est(:,2),phi_est(:,3),'o','linewidth',2.5)
errorbar(1:6,phi_est(:,1),phi_est(:,4),phi_est(:,5),'o','linewidth',2.5)
title('$$\phi \,(\times 10^{-5})$$','fontsize',14,'Interpreter','latex')
plot([0 7],[0 0],':','color',Tau_colour,'linewidth',2.5)
plot([0 7],[10 10],':','color',Tau_colour,'linewidth',2.5)
axis([0 7 -1.0 11])
xticks([1:6])
xticklabels({patient_list})
xtickangle(45)
hold off;

subplot(2,6,12)
hold on;
errorbar(1:6,rho_est(:,1),rho_est(:,2),rho_est(:,3),'o','linewidth',2.5)
errorbar(1:6,rho_est(:,1),rho_est(:,4),rho_est(:,5),'o','linewidth',2.5)
title('$$\rho$$','fontsize',14,'Interpreter','latex')
plot([0 7],[0 0],':','color',Tau_colour,'linewidth',2.5)
plot([0 7],[0.024 0.024],':','color',Tau_colour,'linewidth',2.5)
hold off;
xticks([1:6])
xticklabels({patient_list})
xtickangle(45)
axis([0 7 -0.005 0.029])

%%

for p_id=1:length(patient_list)
    subplot(2,6,p_id)

    patient = patient_list(p_id,:);
    master_path = ['Data/',patient,'/src.tiv.RefractoryCellModel_JSF_2000/'];
    raw_data = readmatrix(['PatientData/CSVs/',patient,'.csv']);
    ci_data =  readmatrix([master_path,'plt_df.csv']);
    ext_data = readmatrix([master_path,'ext_prdf.csv']);
    fc_time = 15;
    int_time = 15;
    ext_data(fc_time,:) = [];
    time = raw_data(:,1)+1;
    V_data = raw_data(:,2);
    
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

    BC_CI_min_ii = CI_min_ii(:,1:(fc_time+1));
    BC_CI_max_ii = CI_max_ii(:,1:(fc_time+1));
    BC_CI_time = CI_time(:,1:(fc_time+1));

    FC_CI_min_ii = CI_min_ii(:,(fc_time+2):end);
    FC_CI_max_ii = CI_max_ii(:,(fc_time+2):end);
    FC_CI_time = CI_time(:,(fc_time+2):end);
    
    colororder([0 0 0; myblue_colour])


    hold on;
    yyaxis left

    for jj = [1:2:length(CI_val)]
        ciplot(BC_CI_min_ii(jj,:),BC_CI_max_ii(jj,:),BC_CI_time(jj,:),Tau_colour,1.0*jj/length(CI_val))
    end
    for jj = [1:2:length(CI_val)]
        ciplot(FC_CI_min_ii(jj,:),FC_CI_max_ii(jj,:),FC_CI_time(jj,:),JSF_colour,1.0*jj/length(CI_val))
    end

    % xlabel('Time (days)','fontsize',14,'Interpreter','latex')
    % ylabel('Viral Load','fontsize',14,'Interpreter','latex')
    title(patient,'fontsize',14,'Interpreter','latex')


    axis([0 14 10^(-1) 10^(9)])

    ax = gca;
    ax.YAxis(1).Scale ="log";

    scatter(time(time<=fc_time),10.^(V_data(time<=fc_time)),50,'ko','filled')
    scatter(time(time>fc_time),10.^(V_data(time>fc_time)),50,'ko')
    plot([-1 15],[10^(-0.65) 10^(-0.65)],'--','linewidth',2.5,'color',[0.5 0.5 0.5])

    hold off

    yyaxis right
    hold on
    plot(ext_data((ext_data(:,1) >= 0),1),ext_data((ext_data(:,1) >= 0),2),'-','linewidth',2.5,'color',myblue_colour)
    axis([0 21 0 1])
    % ylabel('Probability of Viral Clearance','fontsize',14,'Interpreter','latex')

    plot([int_time int_time],[-1 10^5],':','linewidth',2.5,'color','black')
    hold off


end

