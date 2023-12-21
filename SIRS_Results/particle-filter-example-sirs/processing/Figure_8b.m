

clear all;

ext_data_all = zeros(301,2,7);
int_vals = [0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0];
for ii=1:7
    ext_data_ii = readmatrix(['2000_100/sirs_model.SIRS_JSF_particle=2000_intMag=',num2str(int_vals(ii)),'/ext_prdf.csv']);

    ext_data_all(:,:,ii) = ext_data_ii;
end
% ext_data_70 = readmatrix(['2000_100/sirs_model.SIRS_JSF_particle=2000_intMag=0.7/ext_prdf.csv']);
% 
% ext_data_75 = readmatrix(['2000_100/sirs_model.SIRS_JSF_particle=2000_intMag=0.75/ext_prdf.csv']);
% 
% ext_data_80 = readmatrix(['2000_100/sirs_model.SIRS_JSF_particle=2000_intMag=0.8/ext_prdf.csv']);
% 
% ext_data_85 = readmatrix(['2000_100/sirs_model.SIRS_JSF_particle=2000_intMag=0.85/ext_prdf.csv']);
% 
% ext_data_90 = readmatrix(['2000_100/sirs_model.SIRS_JSF_particle=2000_intMag=0.9/ext_prdf.csv']);
% 
% ext_data_95 = readmatrix(['2000_100/sirs_model.SIRS_JSF_particle=2000_intMag=0.95/ext_prdf.csv']);
% 
% ext_data_no = readmatrix(['2000_100/sirs_model.SIRS_JSF_particle=2000_intMag=1.0/ext_prdf.csv']);

int_time = 100;
%%
JSF_colour = 1/255*[131, 96, 150];
Tau_colour = 1/255*[237, 123, 123];

myblue_colour = 1/255*[95, 130, 237];

c_map = myColour2Gradient(7,Tau_colour,myblue_colour);
figure;

hold on;
for ii=1:7
    plot(ext_data_all((ext_data_all(:,1,ii) >= int_time),1,ii),ext_data_all((ext_data_all(:,1,ii) >= int_time),2,ii),'-','linewidth',2.5,'color',c_map(ii,:))
end
% ii=7;
% plot(ext_data_all((ext_data_all(:,1,ii) >= int_time),1,ii),ext_data_all((ext_data_all(:,1,ii) >= int_time),2,ii),'-','linewidth',2.5,'color',myblue_colour)


axis([100-5 400+5 0 1])
ylabel('Probability of fade-out','fontsize',14,'Interpreter','latex')
xlabel('Time (days)','fontsize',14,'Interpreter','latex')


legend('$$\alpha = 0.70$$','$$\alpha = 0.75$$','$$\alpha = 0.80$$',...
    '$$\alpha = 0.85$$','$$\alpha = 0.90$$','$$\alpha = 0.95$$',...
    '$$\alpha = 1.00$$','location','northwest','fontsize',14,'Interpreter','latex')

x = [0.75 0.75];
y = [0.525 0.925];
annotation('textarrow',x,y,'String','Decreasing $$\alpha$$','fontsize',12,'Interpreter','latex')

hold off
