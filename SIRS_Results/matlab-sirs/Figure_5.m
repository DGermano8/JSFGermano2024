close all;
clear all;
master_path = 'FinalRuns/N0_5_10000/';

HowManySimsToDo = 20;

JSF_colour = 1/255*[131, 96, 150];
Tau_colour = 1/255*[237, 123, 123];
Gil_colour = 1/255*[240, 184, 110];

dtObserve = 10^(-1); tFinal = 4*365;
TimeObserve = 0:dtObserve:tFinal;

IData_JFS_Ext = zeros(HowManySimsToDo,length(TimeObserve));
IData_JFS_Fad = zeros(HowManySimsToDo,length(TimeObserve));
IData_JFS_End = zeros(HowManySimsToDo,length(TimeObserve));
SData_JFS_Ext = zeros(HowManySimsToDo,length(TimeObserve));
SData_JFS_Fad = zeros(HowManySimsToDo,length(TimeObserve));
SData_JFS_End = zeros(HowManySimsToDo,length(TimeObserve));
JSF_Ext = 0; JSF_Fad = 0; JSF_End = 0;

IData_Gil_Ext = zeros(HowManySimsToDo,length(TimeObserve));
IData_Gil_Fad = zeros(HowManySimsToDo,length(TimeObserve));
IData_Gil_End = zeros(HowManySimsToDo,length(TimeObserve));
SData_Gil_Ext = zeros(HowManySimsToDo,length(TimeObserve));
SData_Gil_Fad = zeros(HowManySimsToDo,length(TimeObserve));
SData_Gil_End = zeros(HowManySimsToDo,length(TimeObserve));
Gil_Ext = 0; Gil_Fad = 0; Gil_End = 0;


N0ii = 5;
for ii=1:HowManySimsToDo
    
    DATA = readmatrix([master_path,'JumpSwitchFlow/T_3/JSF_N0Index_', num2str(N0ii) , '_iiSeed_', num2str(ii) ,'.csv']);
    if(max(DATA(:,3)) < 50 )
        JSF_Ext = JSF_Ext + 1;
        IData_JFS_Ext(JSF_Ext,:) = DATA(:,3);
        SData_JFS_Ext(JSF_Ext,:) = DATA(:,2);
        
    elseif ((DATA(end,3)) == 0 )
        JSF_Fad = JSF_Fad + 1;
        IData_JFS_Fad(JSF_Fad,:) = DATA(:,3);
        SData_JFS_Fad(JSF_Fad,:) = DATA(:,2);
    else
        JSF_End = JSF_End + 1;
        IData_JFS_End(JSF_End,:) = DATA(:,3);
        SData_JFS_End(JSF_End,:) = DATA(:,2);
    end
    

    DATA = readmatrix([master_path,'Gillespie/GIL_N0Index_', num2str(N0ii) , '_iiSeed_', num2str(ii) ,'.csv']);
    if(max(DATA(:,3)) < 50 )
        Gil_Ext = Gil_Ext + 1;
        IData_Gil_Ext(Gil_Ext,:) = DATA(:,3);
        SData_Gil_Ext(Gil_Ext,:) = DATA(:,2);
    elseif ((DATA(end,3)) == 0 )
        Gil_Fad = Gil_Fad + 1;
        IData_Gil_Fad(Gil_Fad,:) = DATA(:,3);
        SData_Gil_Fad(Gil_Fad,:) = DATA(:,2);
    else
        Gil_End = Gil_End + 1;
        IData_Gil_End(Gil_End,:) = DATA(:,3);
        SData_Gil_End(Gil_End,:) = DATA(:,2);
    
    end
    
    
    
    ii
end

%%

[JSF_peak,pos] = max([IData_JFS_Fad(1:JSF_Fad,:)']);
JSF_Times = TimeObserve(pos);

[Gil_peak,pos] = max([IData_Gil_Fad(1:Gil_Fad,:)']);
Gil_Times = TimeObserve(pos);

%% PLOT THE TRACES WITH THE TWO VIOLIN PLOTS

%%
% PlotTracesWithDists(Times,Peaks,colours,TimeObserve,NUMBER,IData,number2plot,HowManySimsToDo)
num2plot = 10;

x = TimeObserve(1:num2plot:end);
z = zeros(size(x));

f1=figure;
hold on;
for ii=1:JSF_Fad
    hold on;

    y = IData_JFS_Fad(ii,(1:num2plot:end))/(10^3);
    S = surface([x;x],[y;y],[z;z],...
            'facecol','no',...
            'edgecol','interp',...
            'linew',2,...
            'edgealpha',0.01,...
            'edgecolor',JSF_colour);
end
plot(TimeObserve, mean(IData_JFS_Fad(1:JSF_Fad,:))/(10^3), 'color', JSF_colour, 'LineWidth', 2);
axis([0 365/2 0 16.5 ])
hold off;
xlabel('Time (days)','FontSize',14, 'Interpreter','Latex')
ylabel('Number of Infections $$(\times 10^{3})$$','FontSize',14, 'Interpreter','Latex')
%     title(['Fade-Out Probability = ', (num2str(NUMBER/HowManySimsToDo))],'FontSize',14)
f1.Position = [326   201   311   273];

f2=figure;
v = violinplot(JSF_Times,ones(size(JSF_Times)), 'QuartileStyle','boxplot', 'HalfViolin','left','Showdata',false);
v(1,1).ViolinPlot.FaceColor = JSF_colour;
v(1,1).ViolinPlot.EdgeColor = JSF_colour;
v(1,1).ViolinPlot.FaceAlpha = 0.6;
axis([0.65 1 0 365/2])
view(90,90)
% axis off
ylabel('Peak infetion time (days)','FontSize',14, 'Interpreter','Latex')
box off;
xlabel('-','FontSize',14, 'Interpreter','Latex','color','white')
ax = gca;
ax.XAxis.Color = 'white';
f2.Position = [327   548   311    86];

f3=figure;
v = violinplot(JSF_peak/10^3,ones(size(JSF_peak)), 'QuartileStyle','boxplot', 'HalfViolin','right','Showdata',false);
v(1,1).ViolinPlot.FaceColor = JSF_colour;
v(1,1).ViolinPlot.EdgeColor = JSF_colour;
v(1,1).ViolinPlot.FaceAlpha = 0.6;
axis([1 1.4 0 16.5 ])
view(0,90)
ylabel('Peak number of infective $$(\times 10^{3})$$','FontSize',14, 'Interpreter','Latex')
xlabel('-','FontSize',14, 'Interpreter','Latex','color','white')
ax = gca;
ax.XAxis.Color = 'white';
box off;
f3.Position = [638   201    80   274];

%%


x = TimeObserve(1:num2plot:end);
z = zeros(size(x));

f1=figure;
hold on;
for ii=1:Gil_Fad
    hold on;

    y = IData_Gil_Fad(ii,(1:num2plot:end))/(10^3);
    S = surface([x;x],[y;y],[z;z],...
            'facecol','no',...
            'edgecol','interp',...
            'linew',2,...
            'edgealpha',0.01,...
            'edgecolor',Gil_colour);
%         plot_darkmode
end
plot(TimeObserve, mean(IData_Gil_Fad(1:Gil_Fad,:))/(10^3), 'color', Gil_colour, 'LineWidth', 2);
axis([0 365/2 0 16.5 ])
hold off;
xlabel('Time (days)','FontSize',14, 'Interpreter','Latex')
ylabel('Number of Infections $$(\times 10^{3})$$','FontSize',14, 'Interpreter','Latex')
%     title(['Fade-Out Probability = ', (num2str(NUMBER/HowManySimsToDo))],'FontSize',14)
f1.Position = [326   201   311   273];

f2=figure;
v = violinplot(Gil_Times,ones(size(Gil_Times)), 'QuartileStyle','boxplot', 'HalfViolin','left','Showdata',false);
v(1,1).ViolinPlot.FaceColor = Gil_colour;
v(1,1).ViolinPlot.EdgeColor = Gil_colour;
v(1,1).ViolinPlot.FaceAlpha = 0.6;
axis([0.65 1 0 365/2])
view(90,90)
% axis off
ylabel('Peak infetion time (days)','FontSize',14, 'Interpreter','Latex')
xlabel('-','FontSize',14, 'Interpreter','Latex','color','white')
ax = gca;
ax.XAxis.Color = 'white';
box off;
f2.Position = [327   548   311    86];

f3=figure;
v = violinplot(Gil_peak/10^3,ones(size(Gil_peak)), 'QuartileStyle','boxplot', 'HalfViolin','right','Showdata',false);
v(1,1).ViolinPlot.FaceColor = Gil_colour;
v(1,1).ViolinPlot.EdgeColor = Gil_colour;
v(1,1).ViolinPlot.FaceAlpha = 0.6;
axis([1 1.4 0 16.5 ])
view(0,90)
ylabel('Peak number of infective $$(\times 10^{3})$$','FontSize',14, 'Interpreter','Latex')
xlabel('-','FontSize',14, 'Interpreter','Latex','color','white')
ax = gca;
ax.XAxis.Color = 'white';
box off;
f3.Position = [638   201    80   274];


%%
CumInf_JFS_End = zeros(JSF_End,length(TimeObserve));
CumInf_Gil_End = zeros(Gil_End,length(TimeObserve));

for ii=1:JSF_End
    [y, x] = IntegrateByTraps(IData_JFS_End(ii,:),TimeObserve);
    CumInf_JFS_End(ii,:) = y;
end
for ii=1:Gil_End
    [y, x] = IntegrateByTraps(IData_Gil_End(ii,:),TimeObserve);
    CumInf_Gil_End(ii,:) = y;
end

%% PLOT THE CUMULATIVE WITH SINGLE VIOLIN

%%
% PlotCUMWithDists(Peaks,colours,TimeObserve,NUMBER,IData,number2plot,HowManySimsToDo)
FinalCum_JSF = CumInf_JFS_End(1:JSF_End,end);

[y, x] = IntegrateByTraps(CumInf_JFS_End(1,:),TimeObserve);
z = zeros(size(x(1:num2plot:end)));

f1=figure;
hold on;
for ii=1:JSF_End
    y = CumInf_JFS_End(ii,:)';
    S = surface([x(1:num2plot:end)';x(1:num2plot:end)'],[y(1:num2plot:end)';y(1:num2plot:end)']/(10^5),[z';z'],...
    'facecol','no','edgecol','interp','linew',2,'edgealpha',0.01,'edgecolor',JSF_colour);
end
%     title(['Endemic Probability = ', (num2str(NUMBER/HowManySimsToDo))],'FontSize',14)

plot(TimeObserve, mean(CumInf_JFS_End(1:JSF_End,:))/(10^5), 'color', JSF_colour, 'LineWidth', 2);
axis([0 4*365 0 18 ])
hold off;
xlabel('Time (days)','FontSize',14, 'Interpreter','Latex')
ylabel('Cumulative Infections $$(\times 10^{5})$$','FontSize',14, 'Interpreter','Latex')
f1.Position = [326   201   311   273];

f3=figure;
v = violinplot(FinalCum_JSF,ones(size(FinalCum_JSF)), 'QuartileStyle','boxplot', 'HalfViolin','right','Showdata',false);
v(1,1).ViolinPlot.FaceColor = JSF_colour;
v(1,1).ViolinPlot.EdgeColor = JSF_colour;
v(1,1).ViolinPlot.FaceAlpha = 0.6;
axis([1 1.5 0 1.8*10^6])
view(0,90)
ylabel('Cumulative Infections $$(\times 10^{5})$$','FontSize',14, 'Interpreter','Latex')
ax = gca;
ax.XAxis.Color = 'white';
box off;
xlabel('-','FontSize',14, 'Interpreter','Latex','color','white')
f3.Position = [638   201    80   274];


%%
FinalCum_Gil = CumInf_Gil_End(1:Gil_End,end);

[y, x] = IntegrateByTraps(CumInf_Gil_End(1,:),TimeObserve);
z = zeros(size(x(1:num2plot:end)));

f1=figure;
hold on;
for ii=1:Gil_End
    y = CumInf_Gil_End(ii,:)';
    S = surface([x(1:num2plot:end)';x(1:num2plot:end)'],[y(1:num2plot:end)';y(1:num2plot:end)']/(10^5),[z';z'],...
    'facecol','no','edgecol','interp','linew',2,'edgealpha',0.01,'edgecolor',Gil_colour);
end
%     title(['Endemic Probability = ', (num2str(NUMBER/HowManySimsToDo))],'FontSize',14)

plot(TimeObserve, mean(CumInf_Gil_End(1:Gil_End,:))/(10^5), 'color', Gil_colour, 'LineWidth', 2);
axis([0 4*365 0 18 ])
hold off;
xlabel('Time (days)','FontSize',14, 'Interpreter','Latex')
ylabel('Cumulative Infections $$(\times 10^{5})$$','FontSize',14, 'Interpreter','Latex')
f1.Position = [326   201   311   273];

f3=figure;
v = violinplot(FinalCum_Gil,ones(size(FinalCum_Gil)), 'QuartileStyle','boxplot', 'HalfViolin','right','Showdata',false);
v(1,1).ViolinPlot.FaceColor = Gil_colour;
v(1,1).ViolinPlot.EdgeColor = Gil_colour;
v(1,1).ViolinPlot.FaceAlpha = 0.6;
axis([1 1.5 0 1.8*10^6])
view(0,90)
ylabel('Cumulative Infections $$(\times 10^{5})$$','FontSize',14, 'Interpreter','Latex')
xlabel('-','FontSize',14, 'Interpreter','Latex','color','white')
ax = gca;
ax.XAxis.Color = 'white';
box off;
f3.Position = [638   201    80   274];


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

        y = IData(ii,(1:number2plot:end))/(10^3);
        S = surface([x;x],[y;y],[z;z],...
                'facecol','no',...
                'edgecol','interp',...
                'linew',2,...
                'edgealpha',0.01,...
                'edgecolor',colours);
    end
    plot(TimeObserve, mean(IData(1:NUMBER,:))/(10^3), 'color', colours, 'LineWidth', 2);
    axis([0 365/2 0 16.5 ])
    hold off;
    xlabel('Time (days)','FontSize',14, 'Interpreter','Latex')
    ylabel('Number of Infections $$(\times 10^{3})$$','FontSize',14, 'Interpreter','Latex')
%     title(['Fade-Out Probability = ', (num2str(NUMBER/HowManySimsToDo))],'FontSize',14)

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
        S = surface([x(1:number2plot:end)';x(1:number2plot:end)'],[y(1:number2plot:end)';y(1:number2plot:end)']/(10^5),[z';z'],...
        'facecol','no','edgecol','interp','linew',2,'edgealpha',0.01,'edgecolor',colours);
    end
%     title(['Endemic Probability = ', (num2str(NUMBER/HowManySimsToDo))],'FontSize',14)
   
    plot(TimeObserve, mean(IData(1:NUMBER,:))/(10^5), 'color', colours, 'LineWidth', 2);
    axis([0 4*365 0 18 ])
    hold off;
    xlabel('Time (days)','FontSize',14, 'Interpreter','Latex')
    ylabel('Cummulative Infections $$(\times 10^{5})$$','FontSize',14, 'Interpreter','Latex')

    subplot(1,3,3);
    v = violinplot(Peaks,ones(size(Peaks)), 'QuartileStyle','boxplot', 'HalfViolin','right','Showdata',false);
    v(1,1).ViolinPlot.FaceColor = colours;
    v(1,1).ViolinPlot.EdgeColor = colours;
    v(1,1).ViolinPlot.FaceAlpha = 0.6;
    axis([1 2 0 1.8*10^6])
    view(0,90)

    axis off
end

