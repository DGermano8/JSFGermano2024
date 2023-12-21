function [X,TauArr] = GillespieDirectMethod(x0, rates, stoich, times, options)

%%%%%%%%%%%%%%%%% Initilise %%%%%%%%%%%%%%%%%
X0 = x0;
nu = stoich.nu;

tFinal = times(end);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[nRates,nCompartments] = size(nu);


% initialise solution arrays
X = zeros(nCompartments,2*length(times));
X(:,1) = X0;
TauArr = zeros(1,2*length(times));
iters = 1;


while TauArr(iters) <= tFinal
    % compute propensities
    Props = rates(X(:,iters),TauArr(iters));
    % sample exponential waiting time
    dt = exprnd(1/sum(Props));
    % check if the simulation is finished
%     if TauArr(iters) + dt <= tFinal
        iters = iters + 1;
        % sample the next reaction event
        j = randsample(nRates,1,true,Props);
        % update copy numbers  
        X(:,iters) = X(:,iters-1) + nu(j,:)';
        TauArr(iters) = TauArr(iters-1) + dt;
%     else
%         return;
%     end
     if(iters >= length(TauArr))
         X = [X zeros(nCompartments,length(TauArr))];
         TauArr = [TauArr zeros(1,length(TauArr))];
     end
end

X(:,(iters+1:end)) = [];
TauArr((iters+1:end)) = [];
