close all;

Y0 = 8;
X0 = 5;
tspan = [0:0.01:10];

[t,y] = ode15s(@simple_ode,tspan,[Y0; X0]);

figure;
hold on;
plot(t,y(:,1),'LineWidth',2)
plot(t,y(:,2),'LineWidth',2)
legend('Preditor','Prey')
hold off;


function dydt = simple_ode(t,y)
    a = 2.0;
    b = 1.5;
    c = 0.2;

    Y = y(1); % preditor
    X = y(2); % prey

    dYdt = c*X*Y - a*Y;
    dXdt = b*X - c*Y*X;


    dydt = [dYdt; dXdt];
end