clear all;

patient_list = ['432192'; '443108'; '444332';'444391';'445602';'451152'];

raw_data_mat = zeros(length(patient_list),14,2);
ci_data_mat = zeros(length(patient_list),207,5);
ext_data_mat = zeros(length(patient_list),22,2);

fc_time = 14;
int_time = 14;
plot_time = 20;

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
    
    
end
%%
JSF_colour = 1/255*[131, 96, 150];
Tau_colour = 1/255*[237, 123, 123];
Gil_colour = 1/255*[240, 184, 110];
myblue_colour = 1/255*[95, 130, 237];


%%

for p_id=1:length(patient_list)
    f=figure;

    patient = patient_list(p_id,:);
    master_path = ['Data/',patient,'/src.tiv.RefractoryCellModel_JSF_2000/'];
    raw_data = readmatrix(['PatientData/CSVs/',patient,'.csv']);
    ci_data =  readmatrix([master_path,'plt_df.csv']);
    ext_data = readmatrix([master_path,'ext_prdf.csv']);
    fc_time = 14;
    int_time = 14;
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

    FC_CI_min_ii = CI_min_ii(:,(fc_time+1):end);
    FC_CI_max_ii = CI_max_ii(:,(fc_time+1):end);
    FC_CI_time = CI_time(:,(fc_time+1):end);
    
    colororder([0 0 0; myblue_colour])


    hold on;
    yyaxis left

    for jj = [1:2:length(CI_val)-1]
        ciplot(BC_CI_min_ii(jj,:),BC_CI_max_ii(jj,:),BC_CI_time(jj,:),Tau_colour,1.0*jj/(length(CI_val)-0))
    end
    plot(BC_CI_time(length(CI_val),:),BC_CI_min_ii(length(CI_val),:),'-','color',Tau_colour,'linewidth',1.5,'MarkerFaceColor',JSF_colour)
    for jj = [1:2:length(CI_val)-1]
        ciplot(FC_CI_min_ii(jj,:),FC_CI_max_ii(jj,:),FC_CI_time(jj,:),JSF_colour,1.0*jj/(length(CI_val)-0))
    end
    plot(FC_CI_time(length(CI_val),:),FC_CI_min_ii(length(CI_val),:),'-','color',JSF_colour,'linewidth',1.5,'MarkerFaceColor',JSF_colour)

    xlabel('Time (days)','fontsize',16,'Interpreter','latex')
    ylabel('Viral Load','fontsize',16,'Interpreter','latex')
    title(patient,'fontsize',16,'Interpreter','latex')
    
    yticks([10^0 10^3 10^6 10^9])
    xticks([0 5 10 15 20])
%     xticklabels({patient_list})
    


    axis([0 plot_time 10^(-1) 10^(9)])

    ax = gca;
    ax.YAxis(1).Scale ="log";

    scatter(time(time<=fc_time),10.^(V_data(time<=fc_time)),50,'ko','filled')
%     scatter(time(time>fc_time),10.^(V_data(time>fc_time)),50,'ko')
    plot([-1 plot_time],[10^(-0.65) 10^(-0.65)],'--','linewidth',2.5,'color',[0.5 0.5 0.5])
    

    hold off

    yyaxis right
    hold on
    plot(ext_data(1:fc_time,1),ext_data(1:fc_time,2),'-','linewidth',2.5,'color',myblue_colour)
    plot(ext_data(fc_time:end,1),ext_data(fc_time:end,2),'--','linewidth',2.5,'color',myblue_colour)
    axis([0 plot_time -0.01 1.01])
    
%     plot([int_time int_time],[10^(-2) 10^9],':','linewidth',2.5,'color','black')

    % ylabel('Probability of Viral Clearance','fontsize',14,'Interpreter','latex')

%     plot([int_time int_time],[-1 10^5],':','linewidth',2.5,'color','black')
    hold off

    f.Position = [391   324   303   217];
    exportgraphics(gcf,['Figures2/covid_p',num2str(p_id),'.png'],'Resolution',1000)
end

%%
% title('')
% legend({'95\% CI','75\% CI','50\% CI','25\% CI','0\% CI','95\% CI','75\% CI','50\% CI','25\% CI','0\% CI','Nasal viral load data','Limit of detection','Prob. of viral clearance','Estimated Prob. of viral clearance'},'NumColumns',3,'fontsize',15,'Interpreter','latex')
% axis([15 16 10^(9) 10^(10)])
% axis off;
% exportgraphics(gcf,['Figures2/covid_legend.png'],'Resolution',1000)