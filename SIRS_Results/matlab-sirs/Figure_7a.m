close all;
clear all;

JSF_colour = 1/255*[131, 96, 150];
Tau_colour = 1/255*[237, 123, 123];
Gil_colour = 1/255*[240, 184, 110];

T_v = [0:0.25:2 3:1:5];
T_Vector = round(10.^T_v);

HowManySimsToDo = 10000;

dtObserve = 10^(-1); tFinal = 4*365;
TimeObserve = 0:dtObserve:tFinal;

JSF_T_Ext = zeros(length(T_v),1);
JSF_T_Fad = zeros(length(T_v),1);
JSF_T_End = zeros(length(T_v),1);

Gil_Ext = 0; Gil_Fad = 0; Gil_End = 0;


%%

master_path = 'FinalRuns/N0_5_10000/';

N0ii = 5;
for ii=1:HowManySimsToDo
    parfor Ti= 1:length(T_v)
        
        file_path = [master_path,'JumpSwitchFlow/T_',num2str(T_v(Ti)),'/JSF_N0Index_', num2str(N0ii) , '_iiSeed_', num2str(ii) ,'.csv'];        
        DATA = readmatrix(file_path);
        
        if(max(DATA(:,3)) < 50 )
            JSF_T_Ext(Ti) = JSF_T_Ext(Ti) + 1;
        elseif ((DATA(end,3)) == 0 )
            JSF_T_Fad(Ti) = JSF_T_Fad(Ti) + 1;
        else
            JSF_T_End(Ti) = JSF_T_End(Ti) + 1;
        end
    end
    
    file_path = [master_path,'Gillespie/GIL_N0Index_', num2str(N0ii) , '_iiSeed_', num2str(ii) ,'.csv'];        
    DATA = readmatrix(file_path);
    if(max(DATA(:,3)) < 50 )
        Gil_Ext = Gil_Ext + 1;
    elseif ((DATA(end,3)) == 0 )
        Gil_Fad = Gil_Fad + 1;
    else
        Gil_End = Gil_End + 1;
    end
    
    ii
end

%%

f4 = figure;
hold on;
% plot([-1],[-1],'kx','linewidth',1)
% plot([-1],[-1],'k^','linewidth',1)
% plot([-1],[-1],'ko','linewidth',1)

confidenceLevel = 0.95;
zValue = norminv(1 - (1 - confidenceLevel) / 2);


Extinction_prob = [JSF_T_Ext' Gil_Ext]./HowManySimsToDo;
FadeOut_prob = [JSF_T_Fad' Gil_Fad]./HowManySimsToDo;
Endemic_prob = [JSF_T_End' Gil_End]./HowManySimsToDo;

Ext_Error_U = Extinction_prob + zValue * sqrt(Extinction_prob.*(1-Extinction_prob)) ./ sqrt(HowManySimsToDo);
Ext_Error_L = Extinction_prob - zValue * sqrt(Extinction_prob.*(1-Extinction_prob)) ./ sqrt(HowManySimsToDo);

Fad_Error_U = FadeOut_prob + zValue * sqrt(FadeOut_prob.*(1-FadeOut_prob)) ./ sqrt(HowManySimsToDo);
Fad_Error_L = FadeOut_prob - zValue * sqrt(FadeOut_prob.*(1-FadeOut_prob)) ./ sqrt(HowManySimsToDo);

y = [Extinction_prob' FadeOut_prob' Endemic_prob'];
b = bar(y,0.7,'stacked','LineWidth',1.0);
xticks([0.5 1:13 ])
xticklabels({'\Omega^{ }= ','10^0','10^{0.25}','10^{0.5}','10^{0.75}','10^{1}','10^{1.25}','10^{1.5}','10^{1.75}','10^{2}','10^{3}','10^{4}','10^{5}','Gillespie^{ }'})
box on;
axis([0.5 13.5 0 1])

JSF_colour = 1/255*[131, 96, 150];
Tau_colour = 1/255*[237, 123, 123];
Gil_colour = 1/255*[240, 184, 110];
b(1).FaceColor = JSF_colour;
b(2).FaceColor = Tau_colour;
b(3).FaceColor = Gil_colour;
ylabel('Probability','FontSize',13, 'Interpreter','Latex')
title('Sensitivity of scenario','FontSize',13, 'Interpreter','Latex')
legend('Extinct probability','Fade-out probability','Endemic probability','location','northwest','FontSize',13, 'Interpreter','Latex')

f4.Position = [162   433   475   219];

%%

confidenceLevel = 0.10;
zValue = norminv(1 - (1 - confidenceLevel) / 2);
