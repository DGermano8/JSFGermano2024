clear all;
% close all;
addpath('Solver');

%
%
% mBirth*N -----   mBeta*I*S/N     -----      mGamma*I     -----
%     ---> | S |       --->        | I |       --->        | R |
%          -----                   -----                   -----
%          | mDeath*S              | mDeath*I              | mDeath*R
%          V                       V                       V
%

% These define the rates of the system
mBeta = 2.0/7; % Infect "___" people a week
mGamma = 1.0/7; % infecion for "___" weeks
mDeath = 1/(85.0*365); %lifespan
mBirth = mDeath;
mWane = 1/(1.0*365);


% These are the initial conditions
N0 = 10^5;
I0 = 2;
R0 = 0;
S0 = N0-I0-R0;

% How long to simulate for
tFinal = 4*365;

% These are solver options
dt = 0.01;
SwitchingThreshold = [1,1,1];

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

% propensity function
k = [mBeta; mGamma; mBirth; mBirth; mBirth; mDeath; mDeath; mDeath; mWane];

rates = @(X,t) [mBeta*(X(1)*X(2))/(X(1)+X(2)+X(3));
                mGamma*X(2);
                mBirth*X(1);
                mBirth*X(2);
                mBirth*X(3);
                mDeath*X(1);
                mDeath*X(2);
                mDeath*X(3);
                mWane*X(3)];
             
% identify which reactions are discrete and which are continuous
DoDisc = [0; 0; 0];
% allow S and I to switch, but force R to be continuous
% EnforceDo = [1; 1; 1];
EnforceDo = [0; 0; 0];
% allow I to switch, but force S and R to be continuous
% EnforceDo = [1; 0; 1];

%%
stoich = struct();
stoich.nu = nu;
stoich.nuReactant = nuReactant;
stoich.DoDisc = DoDisc;
solTimes = 0:dt:tFinal;
myOpts = struct();
myOpts.EnforceDo = EnforceDo;
myOpts.dt = dt;
myOpts.SwitchingThreshold = SwitchingThreshold;

% These are the initial conditions

I0 = 2;
R0 = 0;


% N0_Vector = [2:0.5:7];
N0_Vector = [7];
HowManySimsToDo = 1000;

TimerData = [];
dtObserve = 10^(-1);
TimeObserve = 0:dtObserve:tFinal;
XObserve = zeros(3,length(TimeObserve));

%%
masterpath = 'FinalRuns/CPU_Time_plots';
% mex_WriteMatrix([masterpath,'/JSF_T1_META.csv'],[] ,'%10.10f',',','w+');
% clear mex_WriteMatrix;
% 
% mex_WriteMatrix([masterpath,'/JSF_T3_META.csv'],[] ,'%10.10f',',','w+');
% clear mex_WriteMatrix;
% 
% mex_WriteMatrix([masterpath,'/JSF_T5_META.csv'],[] ,'%10.10f',',','w+');
% clear mex_WriteMatrix;
% 
% mex_WriteMatrix([masterpath,'/GIL_META.csv'],[] ,'%10.10f',',','w+');
% clear mex_WriteMatrix;
% 
% mex_WriteMatrix([masterpath,'/TAU_META.csv'],[] ,'%10.10f',',','w+');
% clear mex_WriteMatrix;

for N0ii = N0_Vector
%     for ii=325:HowManySimsToDo
    for ii=1:HowManySimsToDo
        N0 = round(10^(N0ii));
        S0 = N0-I0-R0;
        X0 = [S0;I0;R0];
        
        rng(ii)
        tic;
        [XTau,TauArrTau] = TauLeapingMethod(X0, rates, stoich, solTimes, myOpts);
        TauTimer = toc;
        
        [UniqueT, UniqueIndecies] = unique(TauArrTau);
        XObserve = zeros(3,length(TimeObserve));
        for jj=1:3
            UniqueX = XTau(jj,UniqueIndecies);
            XObserve(jj,:) = interp1(UniqueT,UniqueX, TimeObserve,'previous');
        end
        
        mex_WriteMatrix([masterpath,'/Tau/Tau_N0Index_', num2str(N0ii) , '_iiSeed_', num2str(ii) ,'.csv'],[TimeObserve' XObserve'],'%10.4f',',','w+');

        clear mex_WriteMatrix;
        
        mex_WriteMatrix([masterpath,'/TAU_META.csv'],[N0ii ii TauTimer],'%10.10f',',','a+');
        clear mex_WriteMatrix;
%         [N0ii ii TauTimer]
        


        myOpts.SwitchingThreshold = (10^1)*SwitchingThreshold;
        rng(ii)
        tic;
        [XJSF,TauArrJSG] = JumpSwitchFlowSimulator(X0, rates, stoich, solTimes, myOpts);
        JumpSwitchFlowTimer_T1 = toc;
        
        [UniqueT, UniqueIndecies] = unique(TauArrJSG);
        XObserve = zeros(3,length(TimeObserve));
        for jj=1:3
            UniqueX = XJSF(jj,UniqueIndecies);
            XObserve(jj,:) = interp1(UniqueT,UniqueX, TimeObserve,'previous');
        end
        mex_WriteMatrix([masterpath,'/JumpSwitchFlow/T1/JSF_N0Index_', num2str(N0ii) , '_iiSeed_', num2str(ii) ,'.csv'],[TimeObserve' XObserve'],'%10.4f',',','w+');
        clear mex_WriteMatrix;
        
        mex_WriteMatrix([masterpath,'/JSF_T1_META.csv'],[N0ii ii JumpSwitchFlowTimer_T1],'%10.10f',',','a+');
        clear mex_WriteMatrix;
%         [N0ii ii JumpSwitchFlowTimer]


        myOpts.SwitchingThreshold = (10^3)*SwitchingThreshold;
        rng(ii)
        tic;
        [XJSF,TauArrJSG] = JumpSwitchFlowSimulator(X0, rates, stoich, solTimes, myOpts);
        JumpSwitchFlowTimer_T3 = toc;
        
        [UniqueT, UniqueIndecies] = unique(TauArrJSG);
        XObserve = zeros(3,length(TimeObserve));
        for jj=1:3
            UniqueX = XJSF(jj,UniqueIndecies);
            XObserve(jj,:) = interp1(UniqueT,UniqueX, TimeObserve,'previous');
        end
        mex_WriteMatrix([masterpath,'/JumpSwitchFlow/T3/JSF_N0Index_', num2str(N0ii) , '_iiSeed_', num2str(ii) ,'.csv'],[TimeObserve' XObserve'],'%10.4f',',','w+');
        clear mex_WriteMatrix;
        
        mex_WriteMatrix([masterpath,'/JSF_T3_META.csv'],[N0ii ii JumpSwitchFlowTimer_T3],'%10.10f',',','a+');
        clear mex_WriteMatrix;
        
        
        
        myOpts.SwitchingThreshold = (10^5)*SwitchingThreshold;
        rng(ii)
        tic;
        [XJSF,TauArrJSG] = JumpSwitchFlowSimulator(X0, rates, stoich, solTimes, myOpts);
        JumpSwitchFlowTimer_T5 = toc;
        
        [UniqueT, UniqueIndecies] = unique(TauArrJSG);
        XObserve = zeros(3,length(TimeObserve));
        for jj=1:3
            UniqueX = XJSF(jj,UniqueIndecies);
            XObserve(jj,:) = interp1(UniqueT,UniqueX, TimeObserve,'previous');
        end
        mex_WriteMatrix([masterpath,'/JumpSwitchFlow/T5/JSF_N0Index_', num2str(N0ii) , '_iiSeed_', num2str(ii) ,'.csv'],[TimeObserve' XObserve'],'%10.4f',',','w+');
        clear mex_WriteMatrix;
        
        mex_WriteMatrix([masterpath,'/JSF_T5_META.csv'],[N0ii ii JumpSwitchFlowTimer_T5],'%10.10f',',','a+');
        clear mex_WriteMatrix;
        
        
        
%         [N0ii ii JumpSwitchFlowTimer]

        

%         rng(ii)
%         tic;
%         [XG,TauArrG] = GillespieDirectMethod(X0, rates, stoich, solTimes, myOpts);
%         GillespieTimer = toc;
%         
% 
%         [UniqueT, UniqueIndecies] = unique(TauArrG);
%         XObserve = zeros(3,length(TimeObserve));
%         for jj=1:3
%             UniqueX = XG(jj,UniqueIndecies);
%             XObserve(jj,:) = interp1(UniqueT,UniqueX, TimeObserve,'previous');
%         end
%         mex_WriteMatrix([masterpath,'/Gillespie/GIL_N0Index_', num2str(N0ii) , '_iiSeed_', num2str(ii) ,'.csv'],[TimeObserve' XObserve'],'%10.4f',',','w+');
%         clear mex_WriteMatrix;
% 
%         mex_WriteMatrix([masterpath,'/GIL_META.csv'],[N0ii ii GillespieTimer],'%10.10f',',','a+');
%         clear mex_WriteMatrix;
%         [N0ii ii GillespieTimer]
        
%         [N0ii ii TauTimer JumpSwitchFlowTimer_T1 JumpSwitchFlowTimer_T3 JumpSwitchFlowTimer_T5 GillespieTimer]
        [N0ii ii TauTimer JumpSwitchFlowTimer_T1 JumpSwitchFlowTimer_T3 JumpSwitchFlowTimer_T5]

    end

end





