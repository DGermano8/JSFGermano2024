close all;
clear all;

JSF_colour = 1/255*[131, 96, 150];
Tau_colour = 1/255*[237, 123, 123];
Gil_colour = 1/255*[240, 184, 110];

T_v = [0:0.25:2 3:1:5];
T_Vector = round(10.^T_v);

HowManySimsToDo = 50;

dtObserve = 10^(-1); tFinal = 4*365;
TimeObserve = 0:dtObserve:tFinal;

IData_JFS_T_Ext = zeros(length(T_v),HowManySimsToDo,length(TimeObserve));
IData_JFS_T_Fad = zeros(length(T_v),HowManySimsToDo,length(TimeObserve));
IData_JFS_T_End = zeros(length(T_v),HowManySimsToDo,length(TimeObserve));
JSF_T_Ext = zeros(length(T_v),1);
JSF_T_Fad = zeros(length(T_v),1);
JSF_T_End = zeros(length(T_v),1);

IData_Gil_Ext = zeros(HowManySimsToDo,length(TimeObserve));
IData_Gil_Fad = zeros(HowManySimsToDo,length(TimeObserve));
IData_Gil_End = zeros(HowManySimsToDo,length(TimeObserve));
Gil_Ext = 0; Gil_Fad = 0; Gil_End = 0;


%%

master_path = 'FinalRuns/N0_5_10000/';

N0ii = 5;
for ii=1:HowManySimsToDo
    for Ti= 1:length(T_v)
        
        file_path = [master_path,'JumpSwitchFlow/T_',num2str(T_v(Ti)),'/JSF_N0Index_', num2str(N0ii) , '_iiSeed_', num2str(ii) ,'.csv'];        
        DATA = readmatrix(file_path);
        
        if(max(DATA(:,3)) < 50 )
            JSF_T_Ext(Ti) = JSF_T_Ext(Ti) + 1;
            IData_JFS_T_Ext(Ti,JSF_T_Ext(Ti),:) = DATA(:,3);
        elseif ((DATA(end,3)) == 0 )
            JSF_T_Fad(Ti) = JSF_T_Fad(Ti) + 1;
            IData_JFS_T_Fad(Ti,JSF_T_Fad(Ti),:) = DATA(:,3);
        else
            JSF_T_End(Ti) = JSF_T_End(Ti) + 1;
            IData_JFS_T_End(Ti,JSF_T_End(Ti),:) = DATA(:,3);
        end
    end
    
    file_path = [master_path,'Gillespie/GIL_N0Index_', num2str(N0ii) , '_iiSeed_', num2str(ii) ,'.csv'];        
    DATA = readmatrix(file_path);
    if(max(DATA(:,3)) < 50 )
        Gil_Ext = Gil_Ext + 1;
        IData_Gil_Ext(Gil_Ext,:) = DATA(:,3);
    elseif ((DATA(end,3)) == 0 )
        Gil_Fad = Gil_Fad + 1;
        IData_Gil_Fad(Gil_Fad,:) = DATA(:,3);
    else
        Gil_End = Gil_End + 1;
        IData_Gil_End(Gil_End,:) = DATA(:,3);
    end
    
    ii
end

%%
[Gil_peak,pos] = max([IData_Gil_Fad(1:Gil_Fad,:)']);
Gil_Times = TimeObserve(pos);

JSF_Peaks = struct();
JSF_Times = struct();

for Ti= 1:length(T_v)
    if(JSF_T_Fad(Ti) > 0)
        holder_vals = IData_JFS_T_Fad(Ti,1:JSF_T_Fad(Ti),:);
        helding = reshape(holder_vals,[size(holder_vals,2),size(holder_vals,3)]);
        [JSF_Tmp_peak,pos] = max(helding');
        JSF_Tmp_Times = TimeObserve(pos);
        
        JSF_Peaks.(['t',num2str(Ti)]) = JSF_Tmp_peak;
        JSF_Times.(['t',num2str(Ti)]) = JSF_Tmp_Times;
    else
        % Nothing to do - there are no entries!
        JSF_Peaks.(['t',num2str(Ti)]) = [];
        JSF_Times.(['t',num2str(Ti)]) = [];
    end
end

%%

CumInf_Gil_End = zeros(Gil_End,length(TimeObserve));
for ii=1:Gil_End
    [y, x] = IntegrateByTraps(IData_Gil_End(ii,:),TimeObserve);
    CumInf_Gil_End(ii,:) = y;
end
FinalCum_Gil = CumInf_Gil_End(1:Gil_End,end);

FinalCum_JSF = struct();

for Ti= 1:length(T_v)
    if(JSF_T_End(Ti) > 0)
        cum_vals = [];
        holder_vals = IData_JFS_T_End(Ti,1:JSF_T_End(Ti),:);
        helding = reshape(holder_vals,[size(holder_vals,2),size(holder_vals,3)]);
        for ii=1:JSF_T_End(Ti)
        	[y, x] = IntegrateByTraps(helding(ii,:),TimeObserve);
            cum_vals = [cum_vals y(end)];
        end
        FinalCum_JSF.(['t',num2str(Ti)]) =cum_vals;
    else
        % Nothing to do - there are no entries!
        FinalCum_JSF.(['t',num2str(Ti)]) = [];
    end
end
%%
field_names = fieldnames(JSF_Peaks);
figure;
v = violinplot(FinalCum_JSF, field_names,...
    'QuartileStyle','boxplot', 'HalfViolin','right','Showdata',false,'GroupOrder', field_names);
ylabel({'Cumulative infections'})
title('Endemic scenario')

figure;
v = violinplot(JSF_Peaks, ['t5','t6','t7','t8','t9','t10','t11','t12'],...
    'QuartileStyle','boxplot', 'HalfViolin','right','Showdata',false,'GroupOrder', {'t5','t6','t7','t8','t9','t10','t11','t12'});
ylabel({'Peak number of infective'})
title('Epidemic fade-out scenario')

figure;
v = violinplot(JSF_Times, ['t5','t6','t7','t8','t9','t10','t11','t12'],...
    'QuartileStyle','boxplot', 'HalfViolin','right','Showdata',false,'GroupOrder', {'t5','t6','t7','t8','t9','t10','t11','t12'});
ylabel({'Peak infection time (days)'})
title('Epidemic fade-out scenario')

%%

JSF_T1_colour = 1/255*[68, 50, 78];
JSF_T2_colour = 1/255*[99, 72, 113];
JSF_T3_colour = JSF_colour;
JSF_T4_colour = 1/255*[187, 141, 78];
JSF_T5_colour = Gil_colour;
Gil_colour;

f1=figure;
% subplot(2,2,2)
Peak_Data_All = struct();
for ti = 1:length(T_v)
    Peak_Data_All.(field_names{ti}) = rmoutliers(JSF_Peaks.(field_names{ti}),"mean");
%     Peak_Data_All.(field_names{ti}) = JSF_Peaks.(field_names{ti});
    if(isempty(JSF_Peaks.(field_names{ti})))
        Peak_Data_All.(field_names{ti}) = 0;
    end
end
Peak_Data_All.Gillespie = Gil_peak;
peak_names = fieldnames(Peak_Data_All);
v = violinplot(Peak_Data_All,peak_names, 'QuartileStyle','boxplot', 'HalfViolin','right','Showdata',false,'GroupOrder',peak_names,'width',0.5);
% v = violinplot(Peak_Data_All,peak_names, 'QuartileStyle','boxplot','Showdata',false,'GroupOrder',peak_names,'ShowWhiskers',false);
axis([0.5 13.55 1.485*10^4 1.64*10^4])

ylabel({'Peak number of infective'},'FontSize',13, 'Interpreter','Latex')
title('Epidemic fade-out scenario','FontSize',13, 'Interpreter','Latex')
xticks([0.5 1:13 ])
xticklabels({'\Omega^{ }= ','10^0','10^{0.25}','10^{0.5}','10^{0.75}','10^{1}','10^{1.25}','10^{1.5}','10^{1.75}','10^{2}','10^{3}','10^{4}','10^{5}','G^{ }'})
grid on

my_colours = myColour2Gradient(13,[0 0 0], JSF_T3_colour);
for ti = 1:length(T_v)
    v(1,ti).ViolinPlot.FaceColor = my_colours(ti,:);
    v(1,ti).ViolinPlot.EdgeColor = my_colours(ti,:);
    v(1,ti).ViolinPlot.FaceAlpha = 0.6;
end
v(1,13).ViolinPlot.FaceColor = Gil_colour;
v(1,13).ViolinPlot.EdgeColor = Gil_colour;
v(1,13).ViolinPlot.FaceAlpha = 0.6;

f1.Position = [162   433   475   219];
%%

f2 = figure;
% subplot(2,2,2)
Peak_Time_All = struct();
for ti = 1:length(T_v)
    
    Peak_Time_All.(field_names{ti}) = rmoutliers(JSF_Times.(field_names{ti}),"mean");
%     Peak_Time_All.(field_names{ti}) = JSF_Times.(field_names{ti});
    if(isempty(JSF_Times.(field_names{ti})))
        Peak_Time_All.(field_names{ti}) = 0;
    end
end
Peak_Time_All.Gillespie = Gil_Times;
peak_names = fieldnames(Peak_Time_All);
v = violinplot(Peak_Time_All,peak_names, 'QuartileStyle','boxplot', 'HalfViolin','right','Showdata',false,'GroupOrder',peak_names,'width',0.5);
% v = violinplot(Peak_Time_All,peak_names, 'QuartileStyle','boxplot','Showdata',false,'GroupOrder',peak_names,'ShowWhiskers',true)
axis([0.5 13.55 58 120])

ylabel({'Peak infection time (days)'},'FontSize',13, 'Interpreter','Latex')
title('Epidemic fade-out scenario','FontSize',13, 'Interpreter','Latex')
xticks([0.5 1:13 ])
xticklabels({'\Omega^{ }= ','10^0','10^{0.25}','10^{0.5}','10^{0.75}','10^{1}','10^{1.25}','10^{1.5}','10^{1.75}','10^{2}','10^{3}','10^{4}','10^{5}','G^{ }'})
grid on


my_colours = myColour2Gradient(13,[0 0 0], JSF_T3_colour);
for ti = 1:length(T_v)
    v(1,ti).ViolinPlot.FaceColor = my_colours(ti,:);
    v(1,ti).ViolinPlot.EdgeColor = my_colours(ti,:);
    v(1,ti).ViolinPlot.FaceAlpha = 0.6;
end
v(1,13).ViolinPlot.FaceColor = Gil_colour;
v(1,13).ViolinPlot.EdgeColor = Gil_colour;
v(1,13).ViolinPlot.FaceAlpha = 0.6;

f2.Position = [162   433   475   219];

%%

f3 = figure;
% subplot(2,2,2)
Cummulative_Data_All = struct();
for ti = 1:length(T_v)
    Cummulative_Data_All.(field_names{ti}) = rmoutliers(FinalCum_JSF.(field_names{ti}),"mean");
%     Cummulative_Data_All.(field_names{ti}) = FinalCum_JSF.(field_names{ti});
    if(isempty(FinalCum_JSF.(field_names{ti})))
        Cummulative_Data_All.(field_names{ti}) = 0;
%     else
%     	Cummulative_Data_All.(field_namxes{ti}) = rmoutliers(FinalCum_JSF.(field_names{ti}),"mean");
        %     Cummulative_Data_All.(field_names{ti}) = FinalCum_JSF.(field_names{ti});
    end
end
Cummulative_Data_All.Gillespie = FinalCum_Gil;
peak_names = fieldnames(Cummulative_Data_All);
v = violinplot(Cummulative_Data_All,peak_names, 'QuartileStyle','boxplot', 'HalfViolin','right','Showdata',false,'GroupOrder',peak_names,'width',0.5);
% v = violinplot(Cummulative_Data_All,peak_names, 'QuartileStyle','boxplot','Showdata',false,'GroupOrder',peak_names);
axis([0.5 13.55 1.60*10^6 1.734*10^6])

ylabel({'Cummulative infections'},'FontSize',13, 'Interpreter','Latex')
title('Endemic scenario','FontSize',13, 'Interpreter','Latex')
xticks([0.5 1:13 ])
xticklabels({'\Omega^{ }= ','10^0','10^{0.25}','10^{0.5}','10^{0.75}','10^{1}','10^{1.25}','10^{1.5}','10^{1.75}','10^{2}','10^{3}','10^{4}','10^{5}','G^{ }'})
grid on


my_colours = myColour2Gradient(13,[0 0 0], JSF_T3_colour);
for ti = 1:length(T_v)
    v(1,ti).ViolinPlot.FaceColor = my_colours(ti,:);
    v(1,ti).ViolinPlot.EdgeColor = my_colours(ti,:);
    v(1,ti).ViolinPlot.FaceAlpha = 0.6;
end
v(1,13).ViolinPlot.FaceColor = Gil_colour;
v(1,13).ViolinPlot.EdgeColor = Gil_colour;
v(1,13).ViolinPlot.FaceAlpha = 0.6;

f3.Position = [162   433   475   219];

%%
% subplot(2,2,1)
f4 = figure;
hold on;
% plot([-1],[-1],'kx','linewidth',1)
% plot([-1],[-1],'k^','linewidth',1)
% plot([-1],[-1],'ko','linewidth',1)

Extinction_prob = [JSF_T_Ext' Gil_Ext]./HowManySimsToDo;
FadeOut_prob = [JSF_T_Fad' Gil_Fad]./HowManySimsToDo;
Endemic_prob = [JSF_T_End' Gil_End]./HowManySimsToDo;

y = [Extinction_prob' FadeOut_prob' Endemic_prob'];
b = bar(y,0.55,'stacked','LineWidth',1.5);
xticks([0.5 1:13 ])
xticklabels({'\Omega^{ }= ','10^0','10^{0.25}','10^{0.5}','10^{0.75}','10^{1}','10^{1.25}','10^{1.5}','10^{1.75}','10^{2}','10^{3}','10^{4}','10^{5}','G^{ }'})

axis([0.5 13.5 0 1])

JSF_colour = 1/255*[131, 96, 150];
Tau_colour = 1/255*[237, 123, 123];
Gil_colour = 1/255*[240, 184, 110];
b(1).FaceColor = JSF_colour;
b(2).FaceColor = Tau_colour;
b(3).FaceColor = Gil_colour;
ylabel('Probability','FontSize',13, 'Interpreter','Latex')
% title('Sensitivity of scenario','FontSize',13, 'Interpreter','Latex')
legend('Extinct','Fade-out','Endemic','location','northwest','NumColumns',3,'FontSize',13, 'Interpreter','Latex')
grid on
box on

f4.Position = [162   433   475   219];

%%

function [Integral_X, Times] = IntegrateByTraps(X,t)
    Integral_X = zeros(length(t)-1,1);
    Times = zeros(length(t)-1,1);
    for ii=2:length(t)
        Integral_X(ii) = Integral_X(ii-1) + (t(ii)-t(ii-1))*0.5*(X(ii)+X(ii-1));
        Times(ii) = t(ii);
    end

end

function PlotTracesWithDists(Times,Peaks,colours,TimeObserve,NUMBER,IData,number2plot,HowManySimsToDo)
    x = TimeObserve(1:number2plot:end);
    z = zeros(size(x));
    
    subplot(3,3,[4,5,7,8]);
    hold on;
    for ii=1:NUMBER
        hold on;

        y = IData(ii,(1:number2plot:end));
        S = surface([x;x],[y;y],[z;z],...
                'facecol','no',...
                'edgecol','interp',...
                'linew',2,...
                'edgealpha',0.01,...
                'edgecolor',colours);
    end
    plot(TimeObserve, mean(IData(1:NUMBER,:)), 'color', colours, 'LineWidth', 2);
    axis([0 365/2 0 16500 ])
    hold off;
    xlabel('Time (days)','FontSize',14)
    ylabel('Number of Infections','FontSize',14)
    title(['Fade-Out Probability = ', (num2str(NUMBER/HowManySimsToDo))],'FontSize',14)

    subplot(3,3,[1,2]);
    v = violinplot(Times,ones(size(Times)), 'QuartileStyle','boxplot', 'HalfViolin','left','Showdata',false);
    v(1,1).ViolinPlot.FaceColor = colours;
    v(1,1).ViolinPlot.EdgeColor = colours;
    v(1,1).ViolinPlot.FaceAlpha = 0.6;
    axis([0.5 1 0 365/2])
    view(90,90)
    axis off

    subplot(3,3,[6,9]);
    v = violinplot(Peaks,ones(size(Peaks)), 'QuartileStyle','boxplot', 'HalfViolin','right','Showdata',false);
    v(1,1).ViolinPlot.FaceColor = colours;
    v(1,1).ViolinPlot.EdgeColor = colours;
    v(1,1).ViolinPlot.FaceAlpha = 0.6;
    axis([1 2 0 16500 ])
    view(0,90)

    axis off
end

function PlotCUMWithDists(Peaks,colours,TimeObserve,NUMBER,IData,number2plot,HowManySimsToDo)
    [y, x] = IntegrateByTraps(IData(1,:),TimeObserve);
    z = zeros(size(x(1:number2plot:end)));
    
    subplot(1,3,[1,2]);
    hold on;
    for ii=1:NUMBER
        y = IData(ii,:)';
        S = surface([x(1:number2plot:end)';x(1:number2plot:end)'],[y(1:number2plot:end)';y(1:number2plot:end)'],[z';z'],...
        'facecol','no','edgecol','interp','linew',2,'edgealpha',0.01,'edgecolor',colours);
    end
    title(['Endemic Probability = ', (num2str(NUMBER/HowManySimsToDo))],'FontSize',14)
   
    plot(TimeObserve, mean(IData(1:NUMBER,:)), 'color', colours, 'LineWidth', 2);
    axis([0 4*365 0 1.8*10^6 ])
    hold off;
    xlabel('Time (days)','FontSize',14)
    ylabel('Cummulative Infections','FontSize',14)

    subplot(1,3,3);
    v = violinplot(Peaks,ones(size(Peaks)), 'QuartileStyle','boxplot', 'HalfViolin','right','Showdata',false);
    v(1,1).ViolinPlot.FaceColor = colours;
    v(1,1).ViolinPlot.EdgeColor = colours;
    v(1,1).ViolinPlot.FaceAlpha = 0.6;
    axis([1 2 0 1.8*10^6])
    view(0,90)

    axis off
end

