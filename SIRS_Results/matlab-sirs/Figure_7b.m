master_path = 'FinalRuns/CPU_Time_plots/';

JSF_T1_Time = readmatrix([master_path,'JSF_T1_META.csv']);
% JSF_T2_Time = readmatrix(['T2/','JSF_META.csv']);
JSF_T3_Time = readmatrix([master_path,'JSF_T3_META.csv']);
% JSF_T4_Time = readmatrix(['T4/','JSF_META.csv']);
JSF_T5_Time = readmatrix([master_path,'JSF_T5_META.csv']);

% JSF_T3_Time = readmatrix([master_path,'JSF_META.csv']);
Tau_Time = readmatrix([master_path,'TAU_META.csv']);
Gil_Time = readmatrix([master_path,'GIL_META.csv']);

JSF_colour = 1/255*[131, 96, 150];
Tau_colour = 1/255*[237, 123, 123];
Gil_colour = 1/255*[240, 184, 110];

%%
Gil_Time(9001:end,:) = [];
Tau_Time(9001:end,:) = [];
JSF_T1_Time(9001:end,:) = [];
JSF_T3_Time(9001:end,:) = [];
JSF_T5_Time(9001:end,:) = [];
%%

N0_Vector = [2:0.5:6.0];
HowManySimsToDo = 1000;


JSF_T1_Average = zeros(length(N0_Vector),1);
JSF_T1_Error = zeros(length(N0_Vector),1);
% JSF_T2_Average = zeros(length(N0_Vector),1);
% JSF_T2_Error = zeros(length(N0_Vector),1);
JSF_T3_Average = zeros(length(N0_Vector),1);
JSF_T3_Error = zeros(length(N0_Vector),1);
% JSF_T4_Average = zeros(length(N0_Vector),1);
% JSF_T4_Error = zeros(length(N0_Vector),1);
JSF_T5_Average = zeros(length(N0_Vector),1);
JSF_T5_Error = zeros(length(N0_Vector),1);

% JSF_Average = zeros(length(N0_Vector),1);
% JSF_Error = zeros(length(N0_Vector),1);
Tau_Average = zeros(length(N0_Vector),1);
Tau_Error = zeros(length(N0_Vector),1);
Gil_Average = zeros(length(N0_Vector),1);
Gil_Error = zeros(length(N0_Vector),1);

confidenceLevel = 0.95;
zValue = norminv(1 - (1 - confidenceLevel) / 2);

ii=0;
for N0ii=N0_Vector
    ii=ii+1;
    strt=(ii-1)*HowManySimsToDo+1;
    fin = (ii)*HowManySimsToDo;
    
    JSF_T1_Average(ii) = mean(JSF_T1_Time(strt:fin,3));
    JSF_T1_Error(ii) = zValue * std(JSF_T1_Time(strt:fin,3)) / sqrt(HowManySimsToDo);
%     JSF_T2_Average(ii) = mean(JSF_T2_Time(strt:fin,3));
%     JSF_T2_Error(ii) = zValue * std(JSF_T2_Time(strt:fin,3)) / sqrt(HowManySimsToDo);
    JSF_T3_Average(ii) = mean(JSF_T3_Time(strt:fin,3));
    JSF_T3_Error(ii) = zValue * std(JSF_T3_Time(strt:fin,3)) / sqrt(HowManySimsToDo);
%     JSF_T4_Average(ii) = mean(JSF_T4_Time(strt:fin,3));
%     JSF_T4_Error(ii) = zValue * std(JSF_T4_Time(strt:fin,3)) / sqrt(HowManySimsToDo);
    JSF_T5_Average(ii) = mean(JSF_T5_Time(strt:fin,3));
    JSF_T5_Error(ii) = zValue * std(JSF_T5_Time(strt:fin,3)) / sqrt(HowManySimsToDo);
    
    %     JSF_Average(ii) = mean(JSF_T3_Time(strt:fin,3));
%     JSF_Error(ii) = zValue * std(JSF_T3_Time(strt:fin,3)) / sqrt(HowManySimsToDo);
    
    Tau_Average(ii) = mean(Tau_Time(strt:fin,3));
    Tau_Error(ii) = zValue * std(Tau_Time(strt:fin,3)) / sqrt(HowManySimsToDo);

    Gil_Average(ii) = mean(Gil_Time(strt:fin,3));
    Gil_Error(ii) = zValue * std(Gil_Time(strt:fin,3)) / sqrt(HowManySimsToDo);
    
end

%%

f1=figure;

hold on;

plot(10.^(N0_Vector), JSF_T1_Average, 'o:','color', JSF_colour,'linewidth',2);
% plot(10.^(N0_Vector), JSF_T2_Average, '^-','color', JSF_colour,'linewidth',2);
plot(10.^(N0_Vector), JSF_T3_Average, 'o--','color', JSF_colour,'linewidth',2);
% plot(10.^(N0_Vector), JSF_T4_Average, '^-','color', JSF_colour,'linewidth',2);
plot(10.^(N0_Vector), JSF_T5_Average, 'o-','color', JSF_colour,'linewidth',2);

plot(10.^(N0_Vector), Tau_Average, 'o-','color', Tau_colour,'linewidth',2);

plot(10.^(N0_Vector), Gil_Average, 'o-','color', Gil_colour,'linewidth',2);



errorbar(10.^(N0_Vector), JSF_T1_Average, JSF_T1_Error,'color', JSF_colour, 'LineStyle', 'none', 'CapSize', 10,'linewidth',2);
errorbar(10.^(N0_Vector), JSF_T3_Average, JSF_T3_Error,'color', JSF_colour, 'LineStyle', 'none', 'CapSize', 10,'linewidth',2);
errorbar(10.^(N0_Vector), JSF_T5_Average, JSF_T5_Error,'color', JSF_colour, 'LineStyle', 'none', 'CapSize', 10,'linewidth',2);
errorbar(10.^(N0_Vector), Tau_Average, Tau_Error,'color', Tau_colour, 'LineStyle', 'none', 'CapSize', 10,'linewidth',2);
errorbar(10.^(N0_Vector), Gil_Average, Gil_Error,'color', Gil_colour, 'LineStyle', 'none', 'CapSize', 10,'linewidth',2);
xlabel('$$N_0$$ - log scale','fontsize',14, 'Interpreter','Latex')
ylabel('CPU run time (seconds) - log scale','fontsize',14, 'Interpreter','Latex')
legend('Jump-switch-flow ($$\Omega=10^{1}$$)','Jump-switch-flow ($$\Omega=10^{3}$$)','Jump-switch-flow ($$\Omega=10^{5}$$)','Tau-Leaping','Gilllespie',...
    'fontsize',12,'location','northwest', 'Interpreter','Latex')
hold off;
ax = gca;
ax.YAxis.Scale ="log";
ax.XAxis.Scale ="log";
axis([10^1.9 10^6.1 8*10^-3 10^2 ])
% exportgraphics(gca,'RunTime.png')clear all
f1.Position = [162   308   443   344];



%%
f2 = figure;

my_colours = myColour2Gradient(12,[0 0 0], JSF_colour);

HowManySimsToDo = 10000;

T_v = [0:0.25:2 3:1:5];
T_Vector = round(10.^T_v);

master_path = 'FinalRuns/N0_5_10000/';

N0ii = 5;

confidenceLevel = 0.95;
zValue = norminv(1 - (1 - confidenceLevel) / 2);


JSF_Time = zeros(length(T_v),10000,3);
JSF_Aver = zeros(length(T_v),1);
JSF_Err = zeros(length(T_v),1);

for ii=1:length(T_v)
    JSF_Time(ii,:,:) = readmatrix([master_path,'JSF_META_T_',num2str(ii),'.csv']);
    tmp = readmatrix([master_path,'JSF_META_T_',num2str(ii),'.csv']);
    
    JSF_Aver(ii) = mean(tmp(:,3));
    
    JSF_Err(ii) = zValue * std(tmp(:,3)) / sqrt(HowManySimsToDo);
end
Gil_Time = readmatrix([master_path,'GIL_META.csv']);
Gil_Aver = mean(Gil_Time(:,3));
Gil_Err = zValue * std(Gil_Time(:,3)) / sqrt(HowManySimsToDo);

hold on;
for ii=1:length(T_v)-1
    semilogy([ii ii+1],[JSF_Aver(ii) JSF_Aver(ii+1)],'-o','color',  my_colours(ii+1,:),'linewidth',2)
    errorbar([ii ii+1],[JSF_Aver(ii) JSF_Aver(ii+1)], [JSF_Err(ii) JSF_Err(ii+1)],'color', my_colours(ii+1,:), 'LineStyle', 'none', 'CapSize', 10,'linewidth',2);
end

semilogy([ 12 13 ],[JSF_Aver(end) Gil_Aver],'--','color', 0.5*[1 1 1],'linewidth',1.5)

semilogy([ 13 ],[Gil_Aver],'-o','color', Gil_colour,'linewidth',2)
errorbar([ 13 ], [Gil_Aver], [Gil_Err],'color', Gil_colour, 'LineStyle', 'none', 'CapSize', 10,'linewidth',2);

ssi = computureSteadyStateI(10^5);
plot([9+log(ssi)-6; 9+log(ssi)-6], [10^-2; 10^2], '--', 'color' , 0.5*[1 1 1], 'linewidth',2)
hold off;

xticks([ 1:13 ])
xticklabels({'10^0','10^{0.25}','10^{0.5}','10^{0.75}','10^{1}','10^{1.25}','10^{1.5}','10^{1.75}','10^{2}','10^{3}','10^{4}','10^{5}','G^{ }'})

ylabel('CPU run time (seconds) - log scale','fontsize',14, 'Interpreter','Latex')
xlabel('$$\Omega$$ - log scale','fontsize',14, 'Interpreter','Latex')

ax = gca;
ax.YAxis.Scale ="log";

axis([0.5 13.5 10^(-2) 10^2])

f2.Position = [162   308   443   344];



%%


%%

function SSI = computureSteadyStateI(N0)
    % These define the rates of the system
    mBeta = 2.0/7; % Infect "___" people a week
    mGamma = 1.0/7; % infecion for "___" weeks
    mDeath = 1/(85.0*365); %lifespan
    mBirth = mDeath;
    mWane = 1/(1.0*365);


    % These are the initial conditions
%     N0 = 10^5;
    I0 = 2;
    R0 = 0;
    S0 = N0-I0-R0;

    % How long to simulate for
    tFinal = 10*365;

    % These are solver options
    dt = 0.01;
    SwitchingThreshold = [10^4, 10^3, 10^3];

    % kinetic rate parameters
    X0 = [S0;I0;R0];

    % reactant stoichiometries
    nuReactant = [1,1,0;
                  0,1,0;
                  1,0,0;
                  0,1,0;
                  0,0,1;
                  1,0,0;
                  0,1,0;
                  0,0,1;
                  0,0,1];

    % product stoichiometries
    nuProduct = [0,2,0;
                 0,0,1;
                 2,0,0;
                 1,1,0;
                 1,0,1;
                 0,0,0;
                 0,0,0;
                 0,0,0;
                 1,0,0];
    % stoichiometric matrix
    nu = nuProduct - nuReactant;

    rates = @(X,t) [mBeta*(X(1)*X(2))/(X(1)+X(2)+X(3));
                    mGamma*X(2);
                    mBirth*X(1);
                    mBirth*X(2);
                    mBirth*X(3);
                    mDeath*X(1);
                    mDeath*X(2);
                    mDeath*X(3);
                    mWane*X(3)];

    solTimes = 0:dt:tFinal;
    x = zeros(3,length(solTimes));
    x(:,1) = X0';
    for ii=2:length(solTimes)
        Props = rates(x(:,ii-1),solTimes(ii));
        dXdt = (Props'*nu)';
        x(:,ii) = x(:,ii-1) + dt*dXdt;
    end
    SSI = x(2,end);
end