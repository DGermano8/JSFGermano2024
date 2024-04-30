
%JUMPSWITCHFLOWSIMULATOR  Sample from CD-switching process.
%   [X,TAUARR] = JumpSwitchFlowSimulator(X0, RATES, STOICH, TIMES, OPTIONS)
%   simulates from the continuous-discrete-switching process with flow
%   RATES and stoichiometry STOICH, starting from initial condition
%   X0, returning the process states at the times in TIMES. The output
%   X is a matrix with as many columns as there are species, and as
%   many rows as there are time points in the output. The output
%   TAUARR is a vector of time points at which the output is sampled.
%
%   OPTIONS is a structure with the following fields:
%   - dt: the time step used for the Euler-Maruyama discretisation of
%     the continuous dynamics. Default: 0.01.
%   - EnforceDo: a boolean indicating whether to enforce the discrete
%     dynamics to be active at all times. Default: false.
%   - SwitchingThreshold: a threshold for switching between discrete
%     and continuous dynamics. Default: 0.1.
%
% TODO Currently the times are only not actually used beyond using
% TIMES(END) as the final time point.
%
% TODO Currently the documentation is confused about the name of the
% function provided and the default values for options are not
% respected.
%
% Author: Domenic P.J. Germano (2023).
function [X,TauArr] = JumpSwitchFlowSimulator(x0, rates, stoich, times, options)

%%%%%%%%%%%%%%%%% Initilise %%%%%%%%%%%%%%%%%

X0 = x0;
nu = stoich.nu;
try nuReactant = stoich.nuReactant; catch nuReactant=zeros(size(nu));
end

[nRates,nCompartments] = size(nu);

tFinal = times(end);

% Set default ODE time step
try dt = options.dt; catch dt = 10^(-3);
end

% Set default switching thresholds
try SwitchingThreshold = options.SwitchingThreshold; catch SwitchingThreshold = 1000*ones(1,nCompartments);
end

% Set default compartments to discrete
try DoDisc = stoich.DoDisc; catch DoDisc = ones(nCompartments,1);
end

% Set default dynamics to all switching
try EnforceDo = options.EnforceDo; catch EnforceDo = zeros(nCompartments,1);
end
DoCont = ~DoDisc;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% identify which compartment is in which reaction:
compartInNu = (nu~=0 | nuReactant~=0);
% compartInNu = (nu~=0);

discCompartment = zeros(size(compartInNu,1),1);
contCompartment = ones(size(compartInNu,1),1);
for ii=1:length(X0)
    for jj = 1:size(compartInNu,1)
        if(DoDisc(ii) && compartInNu(jj,ii))
            discCompartment(jj) = 1;
        end
    end
end
contCompartment = contCompartment - discCompartment;

% initialise discrete sum compartments
integralOfFiringTimes = zeros(nRates,1);
RandTimes = rand(nRates,1);
tauArray = zeros(nRates,1);

TimeMesh = 0:dt:tFinal;
overFlowAllocation = round(2.5*length(TimeMesh));

% initialise solution arrays
X = zeros(nCompartments,overFlowAllocation);
X(:,1) = X0;
TauArr = zeros(1,overFlowAllocation);
iters = 1;

% Track Absolute time
AbsT = 0;
ContT = 0;

Xprev = X0;
Xcurr = zeros(nCompartments,1);

NewDiscCompartmemt = zeros(nCompartments,1);
while ContT < tFinal
    Dtau = dt;

    Xprev = X(:,iters);
    Props = rates(Xprev,ContT);

    [Dtau, correctInteger, DoDisc, DoCont, discCompartment, contCompartment,...
        NewDoDisc, NewDoCont, NewdiscCompartment, NewcontCompartment,...
        NewDiscCompartmemt] =...
        UpdateCompartmentRegime(dt, Xprev, Dtau, Props, nu, SwitchingThreshold,...
        DoDisc, DoCont, EnforceDo, discCompartment, contCompartment, compartInNu,...
        nCompartments, NewDiscCompartmemt);
    
    % Perform the Forward Euler Step
    dXdt = ComputedXdt(Xprev, Props, nu, contCompartment, nCompartments);
    Xcurr = X(:,iters) + Dtau*(dXdt.*DoCont);

    % switch a continuous compartment to integer if needed
    OriginalDoCont = DoCont;
    if(correctInteger)
        contCompartment = NewcontCompartment;
        discCompartment = NewdiscCompartment;
        DoCont = NewDoCont;
        DoDisc = NewDoDisc;
    end

    % Dont bother doing anything discrete if its all continuous
    stayWhile = (true)*(sum(DoCont)~=length(DoCont));
    
    TimePassed = 0;
    % Perform the Stochastic Loop
    AbsT = ContT;
    
    DtauContStep = Dtau;

    while stayWhile
        
        Props = rates(Xcurr,AbsT);
%         Props = rates(Xprev,AbsT);
        integralStep = ComputeIntegralOfFiringTimes(Dtau, Props, rates, Xprev, Xcurr, AbsT);
        integralOfFiringTimes = integralOfFiringTimes + integralStep;
        
        % if a compartment has /just/ become discrete, make sure it has
        % zero sumTimes and rest the randTimes
        if(sum(NewDiscCompartmemt) == 1)
            for ii=1:length(Xprev)
                if(NewDiscCompartmemt(ii)==1 &&EnforceDo(ii)==0)
                    for jj = 1:size(compartInNu,1)
                        if(compartInNu(jj,ii)==1)
                            discCompartment(jj) = 1;
                            integralOfFiringTimes(jj) = 0.0;
                            RandTimes(jj) = rand;
                        end
                    end
                end
            end
        end

        % Identify which reactions have fired
        firedReactions = (0 > RandTimes - (1 - exp(-integralOfFiringTimes))) .* discCompartment;
        
        if( sum(firedReactions) > 0)
            tauArray = ComputeFiringTimes(firedReactions,integralOfFiringTimes,RandTimes,Props,dt,nRates,integralStep);
            
            % identify which reaction occurs first
            if(sum(tauArray) > 0)
                
                % Update the discrete compartments
                [Xcurr, Xprev, integralOfFiringTimes, integralStep, RandTimes, TimePassed, AbsT, Dtau1] = ImplementFiredReaction(tauArray ,integralOfFiringTimes,RandTimes,Props,rates,integralStep,TimePassed, AbsT, X, iters, nu, dXdt, OriginalDoCont);
                
                iters = iters + 1;
                X(:,iters) = Xcurr;
                TauArr(iters) = AbsT;
                Dtau = Dtau - Dtau1;
                                
%                 stayWhile = false;
                
            else
%                 TimePassed = DtauContStep;
                stayWhile = false;
            end
        else
%             TimePassed = DtauContStep;
            stayWhile = false;
        end

        if(TimePassed >= DtauContStep)
%             TimePassed = DtauContStep;
            stayWhile = false;
        end


    end

    iters = iters + 1;
    ContT = ContT + (DtauContStep);
    X(:,iters) = X(:,iters-1) + (DtauContStep-TimePassed)*(dXdt.*DoCont);
    
    
    
%     iters = iters + 1;
%     ContT = ContT + Dtau*(1-(TimePassed>0)) + TimePassed;
%     X(:,iters) = Xcurr;
    TauArr(iters) = ContT;
    
    if(sum(NewDiscCompartmemt)==1)
        [~,pos] = max(NewDiscCompartmemt);
        X(pos,iters) = round(X(pos,iters));
        
        for jj = 1:size(compartInNu,1)
            if(NewDiscCompartmemt(pos) && compartInNu(jj,pos))
                discCompartment(jj) = 1;
                integralOfFiringTimes(jj) = 0.0;
                RandTimes(jj) = rand;
            end
        end
        NewDiscCompartmemt(pos) = 0;
    end
    
%     if(iters >= length(TauArr))
%          X = [X zeros(nCompartments,length(TauArr))];
%          TauArr = [TauArr zeros(1,length(TauArr))];
%     end
     
end

if(iters < overFlowAllocation)
    X(:,(iters+1:end)) = [];
    TauArr((iters+1:end)) = [];
end

end

%%
function [DoDisc, DoCont, discCompartmentTmp, contCompartmentTmp] = IsDiscrete(Xprev,nu,Props,dt,SwitchingThreshold,DoDisc,DoCont,discCompartment,contCompartment, EnforceDo, compartInNu)

    OriginalDoDisc = DoDisc;
    OriginalDoCont = DoCont;
  
    DoDisc = (Xprev  <= SwitchingThreshold').*(~EnforceDo) + (EnforceDo).*OriginalDoDisc;
    DoCont = (DoDisc==0).*(~EnforceDo) + (EnforceDo).*OriginalDoCont;

    if( sum((OriginalDoDisc == DoDisc))==length(DoDisc) || (sum(EnforceDo)==length(EnforceDo)) )
        discCompartmentTmp = discCompartment;
        contCompartmentTmp = contCompartment;
    else
        discCompartmentTmp = zeros(size(compartInNu,1),1);
        contCompartmentTmp = ones(size(compartInNu,1),1);

        for ii=1:length(Xprev)
            for jj = 1:size(compartInNu,1)
                if(~EnforceDo(ii))
                    if(DoDisc(ii) && compartInNu(jj,ii))
                        discCompartmentTmp(jj) = 1;
                    end
                else
                    if(OriginalDoDisc(ii) && compartInNu(jj,ii))
                        discCompartmentTmp(jj) = 1;
                    end
                end
            end
        end
        contCompartmentTmp = contCompartmentTmp - discCompartmentTmp;
    end

end

function [Dtau, correctInteger, DoDisc, DoCont, discCompartment, contCompartment, NewDoDisc, NewDoCont, NewdiscCompartment, NewcontCompartment, NewDiscCompartmemt] = UpdateCompartmentRegime(dt, Xprev, Dtau, Props, nu, SwitchingThreshold, DoDisc, DoCont, EnforceDo, discCompartment, contCompartment, compartInNu, nCompartments, NewDiscCompartmemt)
%     NewDiscCompartmemt = zeros(nCompartments,1);
    correctInteger = 0;
    
    % identify which compartment is to be modelled with Discrete and continuous dynamics
    [NewDoDisc, NewDoCont, NewdiscCompartment, NewcontCompartment] = IsDiscrete(Xprev,nu,Props,dt,SwitchingThreshold,DoDisc,DoCont,discCompartment,contCompartment, EnforceDo, compartInNu);

    % move Euler mesh to ensure the new distcrete compartment is integer
    if(nnz(NewDoDisc) > nnz(DoDisc))

        % this ^ identifies a state has switched to discrete
        % identify which compartment is the switch
        [~,pos] = max(NewDoDisc-DoDisc);

        % Perform the Forward Euler Step
        dXdt = sum(Props.*(contCompartment.*nu),1)';

        Dtau = min(dt,abs((round(Xprev(pos)) - Xprev(pos))/dXdt(pos)));

        if(Dtau < dt)
            % used to implement the remainder of the step, and then
            % switch to discrete later
            NewDiscCompartmemt = ((NewDoDisc - DoDisc) ==1);
            correctInteger = 1;
        end
    else
        contCompartment = NewcontCompartment;
        discCompartment = NewdiscCompartment;
        DoCont = NewDoCont;
        DoDisc = NewDoDisc;
    end
end

function dXdt = ComputedXdt(Xprev, Props, nu, contCompartment, nCompartments)
%     dXdt = sum(Props.*(contCompartment.*nu),1)';
    dXdt = (Props'*(contCompartment.*nu))';
end

function integralStep = ComputeIntegralOfFiringTimes(Dtau, Props, rates, Xprev, Xcurr, AbsT)
    % Integrate the cumulative wait times using trapezoid method
    integralStep = Dtau*0.5*(Props + rates(Xcurr, AbsT+Dtau));
end

function tauArray = ComputeFiringTimes(firedReactions,integralOfFiringTimes,RandTimes,Props,dt,nRates,integralStep)
    tauArray = zeros(nRates,1);
    for kk=1:length(firedReactions)
        tauArray(kk) = Inf;
        if(firedReactions(kk))
            % calculate time tau until event using linearisation of integral:
            ExpInt = exp(-(integralOfFiringTimes(kk)-integralStep(kk)));
            % This method explicitly assumes that the rates are
            % varying slowly in time (more specifically, time independant)
            Integral = -1*log((1-RandTimes(kk))/ExpInt);
            
            tau_val_1 = Integral/((Props(kk)));

            % tau_val_2 = -1;
            tau_val =  tau_val_1;
            % were doing a linear approximation to solve this, so it may be
	        % off, in which case we just fix it to a small order here
            % if(tau_val < 0)
            %     if(abs(tau_val_1) < dt^(2))
            %         tau_val_2 = abs(tau_val_1);
            %     end
            %     tau_val_1 = 0;
            %     tau_val = max(tau_val_1,tau_val_2);
            % 
            % end
            tauArray(kk) = tau_val;
        end
    end

end

function [Xcurr, Xprev, integralOfFiringTimes, integralStep, RandTimes, TimePassed, AbsT, Dtau1] = ImplementFiredReaction(tauArray ,integralOfFiringTimes,RandTimes,Props,rates,integralStep,TimePassed, AbsT, X, iters, nu, dXdt, OriginalDoCont);

%     tauArray(tauArray==0.0) = Inf;
    % Identify first event occurance time and type of event
    [Dtau1,pos] = min(tauArray);

    TimePassed = TimePassed + Dtau1;
    AbsT = AbsT + Dtau1;

    % implement first reaction

    Xcurr = X(:,iters) + nu(pos,:)' + Dtau1*(dXdt.*OriginalDoCont);
    Xprev = X(:,iters);

    Props = rates(Xcurr,AbsT);

    % Bring compartments up to date
    integralOfFiringTimes = integralOfFiringTimes - integralStep;
    integralStep = ComputeIntegralOfFiringTimes(Dtau1, Props, rates, Xprev, Xcurr, AbsT);

    integralOfFiringTimes = integralOfFiringTimes+integralStep;

    % reset timers and sums
    RandTimes(pos) = rand;
    integralOfFiringTimes(pos) = 0.0;
                
end

