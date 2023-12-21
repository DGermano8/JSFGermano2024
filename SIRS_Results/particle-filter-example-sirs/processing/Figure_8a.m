master_path = '2000_100/sirs_model.SIRS_JSF_particle=2000_intMag=1/';
% master_path = 'sirs_model.SIRS_JSF_particle=6000_intMag=0.875/';
% master_path = 'sirs_model.SIRS_JSF_particle=6000_intMag=0.75/';

raw_data = readmatrix(['2000_100/sim_df.csv']);
ci_data =  readmatrix([master_path,'plt_df.csv']);

ext_data = readmatrix([master_path,'ext_prdf.csv']);

fc_time = 100;
int_time = 100;

%%

time = raw_data(:,1);
I_data = raw_data(:,2);

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

FC_CI_min_ii = CI_min_ii(:,(fc_time+3):end);
FC_CI_max_ii = CI_max_ii(:,(fc_time+3):end);
FC_CI_time = CI_time(:,(fc_time+3):end);
%%
figure;
% subplot(2,1,1)
JSF_colour = 1/255*[131, 96, 150];
Tau_colour = 1/255*[237, 123, 123];
Gil_colour = 1/255*[240, 184, 110];
myblue_colour = 1/255*[95, 130, 237];

colororder([0 0 0; myblue_colour])


hold on;
yyaxis left

for jj = [1:2:length(CI_val)]
    ciplot(BC_CI_min_ii(jj,:),BC_CI_max_ii(jj,:),BC_CI_time(jj,:),Tau_colour,1.0*jj/length(CI_val))
end
for jj = [1:2:length(CI_val)]
    ciplot(FC_CI_min_ii(jj,:),FC_CI_max_ii(jj,:),FC_CI_time(jj,:),JSF_colour,1.0*jj/length(CI_val))
end

xlabel('Time (days)','fontsize',14,'Interpreter','latex')
ylabel('Number of infectious','fontsize',14,'Interpreter','latex')
axis([-2 max(time)+1 -100 1.05*max(max([FC_CI_max_ii BC_CI_max_ii]))])

scatter(time(time<=fc_time),I_data(time<=fc_time),10,'ko','filled')
scatter(time(time>fc_time),I_data(time>fc_time),10,'ko')

yyaxis right
plot(ext_data((ext_data(:,1) >= int_time),1),ext_data((ext_data(:,1) >= int_time),2),'-','linewidth',2.5,'color',myblue_colour)
axis([0 max(ext_data(:,1)) 0 1])
ylabel('Probability of fade-out','fontsize',14,'Interpreter','latex')
% xlabel('time (days)')

plot([int_time int_time],[-1 10^5],':','linewidth',1.5,'color','black')


legend({'Inference 95\% CI', 'Inference 75\% CI','Inference 50\% CI', 'Inference 25\% CI', 'Inference 0\% CI',...
    'Prediction 95\% CI','Prediction 75\% CI','Prediction 50\% CI', 'Prediction 25\% CI', 'Prediction 0\% CI',...
    'Inference data','Future data', 'Probability of fade-out','Prediction start'},...
    'Location','northeast','NumColumns',3,'fontsize',12,'Interpreter','latex')

