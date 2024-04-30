clear all;
close all;

% These are the initial conditions
Y = 10;
X = 5;

% How long to simulate for
tFinal = 100;

% These are solver options
dt = 0.001;
SwitchingThreshold = 10*[1,1];

% kinetic rate parameters
X0 = [Y; X];


% reactant stoichiometries
nuReactant = [1,0;
              0,1;
              1,1];

% product stoichiometries
nuProduct = [0,0;
             0,2;
             2,0];
% stoichiometric matrix
nu = nuProduct - nuReactant;

% propensity function
mA = 2; mB = 1.5; mC = 0.2;
k = [mA; mB; mC];

rates = @(X,t) [mA*X(1);
                mB*X(2);
                mC*X(1)*X(2)];
             
% identify which reactions are discrete and which are continuous
DoDisc = [0; 0];

EnforceDo = [1; 1];
EnforceDo = [0; 0];


%% 6 17
% rng(1);

stoich = struct();
stoich.nu = nu;
stoich.nuReactant = nuReactant;
stoich.DoDisc = DoDisc;
solTimes = 0:dt:tFinal;
myOpts = struct();
myOpts.EnforceDo = EnforceDo;
myOpts.dt = dt;
myOpts.SwitchingThreshold = SwitchingThreshold;

[XJSF,TauArrJSG] = JumpSwitchFlowSimulator(X0, rates, stoich, solTimes, myOpts);

figure;
hold on;
plot(TauArrJSG,XJSF(1,:),'.','LineWidth',2)
plot(TauArrJSG,XJSF(2,:),'.','LineWidth',2)
plot([0 tFinal],SwitchingThreshold, '--','color',[0.5 0.5 0.5],'LineWidth',2)
legend('Preditor','Prey')
hold off;
axis([0 10 0 40])




