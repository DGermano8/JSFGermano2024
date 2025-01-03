patient_list = ['432192'; '443108'; '444332';'444391';'445602';'451152'];

patient = patient_list(5,:);

master_path = ['Data/5000/',patient,'/src.tiv.RefractoryCellModel_JSF_5000/'];
% master_path = 'sirs_model.SIRS_JSF_particle=6000_intMag=0.875/';
% master_path = 'sirs_model.SIRS_JSF_particle=6000_intMag=0.75/';

raw_data = readmatrix(['PatientData/CSVs/',patient,'.csv']);
ci_data =  readmatrix([master_path,'V_plt_df.csv']);

ext_data = readmatrix([master_path,'ext_prdf.csv']);

fc_time = 14;
int_time = 14;
plot_time = 20;

ext_data(fc_time,:) = [];

%%

time = raw_data(:,1)+1;
V_data = raw_data(:,2);

% plot(time,I_data,'o')

%%
% CIs: [ 0, 12.5, 25, 37.5, 50, 62.5, 75, 87.5, 95 ]

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
%%
f=figure;
% subplot(2,1,1)

title(['Patient ID ',patient],'fontsize',20,'Interpreter','latex')

JSF_colour = 1/255*[131, 96, 150];
Tau_colour = 1/255*[237, 123, 123];
Gil_colour = 1/255*[240, 184, 110];
myblue_colour = 1/255*[95, 130, 237];

colororder([0 0 0; myblue_colour])


hold on;
yyaxis left

for jj = [1:2:length(CI_val)]
    ciplot(BC_CI_min_ii(jj,1:(fc_time+1)),BC_CI_max_ii(jj,1:(fc_time+1)),BC_CI_time(jj,1:(fc_time+1)),Tau_colour,1.0*jj/length(CI_val))
end
for jj = [1:2:length(CI_val)]
    ciplot(BC_CI_min_ii(jj,(fc_time+1):end),BC_CI_max_ii(jj,(fc_time+1):end),BC_CI_time(jj,(fc_time+1):end),JSF_colour,1.0*jj/length(CI_val))
end
for jj = [1:2:length(CI_val)]
    ciplot(FC_CI_min_ii(jj,:),FC_CI_max_ii(jj,:),FC_CI_time(jj,:),JSF_colour,1.0*jj/length(CI_val))
end

xlabel('Time (days)','fontsize',18,'Interpreter','latex')
ylabel('Viral Load','fontsize',18,'Interpreter','latex')
% axis([0 14 10^(-1) 1.05*max([max([FC_CI_max_ii BC_CI_max_ii]) 10.^V_data'])])
axis([0 plot_time 10^(-1) 10^(9)])

yticks([10^0 10^3 10^6 10^9])
xticks([0 5 10 15 20])

ax = gca;
ax.YAxis(1).Scale ="log";

scatter(time(time<=fc_time),10.^(V_data(time<=fc_time)),50,'Wo','filled')
scatter(time(time>=fc_time),10.^(V_data(time>=fc_time)),50,'wo')
plot([-1 plot_time],[10^(-0.65) 10^(-0.65)],'--','linewidth',2.5,'color',[0.5 0.5 0.5])

plot_darkmode
hold off

yyaxis right
hold on
% plot(ext_data((ext_data(1:fc_time,1) >= 0),1),ext_data((ext_data(1:fc_time,1) >= 0),2),'-','linewidth',2.5,'color',myblue_colour)
plot(ext_data(1:fc_time,1),ext_data(1:fc_time,2),'-','linewidth',2.5,'color',myblue_colour)
plot(ext_data(fc_time:end,1),ext_data(fc_time:end,2),'--','linewidth',2.5,'color',myblue_colour)
axis([0 plot_time -0.02 1.02])
ylabel('Probability of Viral Clearance','fontsize',18,'Interpreter','latex')
% xlabel('time (days)')
yticks([0 0.25 0.5 0.75 1])

plot([int_time int_time],[-1 10^5],':','linewidth',2.5,'color','black')
plot_darkmode
hold off
% legend({'Inference 95\% CI', 'Inference 75\% CI','Inference 50\% CI', 'Inference 25\% CI', 'Inference 0\% CI',...
%     'Prediction 95\% CI','Prediction 75\% CI','Prediction 50\% CI', 'Prediction 25\% CI', 'Prediction 0\% CI',...
%     'Inference data','Future data', 'Probability of fade-out','Prediction start'},...
%     'Location','northeast','NumColumns',3,'fontsize',12,'Interpreter','latex')

% f.Position = [420   324   303   217];
f.Position = [420   324   420   320];
exportgraphics(gca,['figs/',patient,'.pdf'],'Resolution',200)