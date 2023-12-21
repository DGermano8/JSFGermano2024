
num2plot = 1;
linwidth=2.5;
markertype_jsf = '-';
markertype_gil = ':';
x = TimeObserve(1:num2plot:end);
z = zeros(size(x));
fsz = 12;
tsz = 14;

figure;
hold on;
y = IData_JFS_Ext(4,(1:num2plot:end));
plot(x,y,markertype_jsf,'color',JSF_colour, 'LineWidth', linwidth);
y = IData_Gil_Ext(5,(1:num2plot:end));
plot(x,y,markertype_gil,'color',Gil_colour, 'LineWidth', linwidth);
hold off;
axis([0 10 0 5 ])
% title('Extiction','fontsize',tsz)
ylabel('I','fontsize',fsz, 'Interpreter','Latex')
xlabel('Time (days)','fontsize',fsz, 'Interpreter','Latex')
cf = gcf;
cf.Position=[440 378 200 180];


figure;
hold on;
y = IData_JFS_Fad(8,(1:num2plot:end));
plot(x,y/10^4,markertype_jsf,'color',JSF_colour, 'LineWidth', linwidth);
y = IData_Gil_Fad(2,(1:num2plot:end));
plot(x,y/10^4,markertype_gil,'color',Gil_colour, 'LineWidth', linwidth);
axis([0 365/2 0 1.6])
plot(x,10^3*ones(length(x),1)/(10^4),'k--')
hold off;
% title('Fade-Out','fontsize',tsz)
ylabel('I $$(\times 10^{4})$$','fontsize',fsz, 'Interpreter','Latex')
xlabel('Time (days)','fontsize',fsz, 'Interpreter','Latex')
cf = gcf;
cf.Position=[440 378 200 180];

figure;
hold on;
y = IData_JFS_End(1,(1:num2plot:end));
plot(x,y/10^4,markertype_jsf,'color',JSF_colour, 'LineWidth', linwidth);
y = IData_Gil_End(1,(1:num2plot:end));
plot(x,y/10^4,markertype_gil,'color',Gil_colour, 'LineWidth', linwidth);
axis([0 3*365 0 1.6 ])
plot(x,10^3*ones(length(x),1)/(10^4),'k--')
hold off;
% title('Endemic','fontsize',tsz)
ylabel('I $$(\times 10^{4})$$','fontsize',fsz, 'Interpreter','Latex')
xlabel('Time (days)','fontsize',fsz, 'Interpreter','Latex')
cf = gcf;
cf.Position=[440 378 200 180];

figure;
hold on;
y = IData_JFS_Ext(4,(1:num2plot:end));
x = SData_JFS_Ext(4,(1:num2plot:end));
plot(x,y,markertype_jsf,'color',JSF_colour, 'LineWidth', linwidth);
y = IData_Gil_Ext(5,(1:num2plot:end));
x = SData_Gil_Ext(5,(1:num2plot:end));
plot(x,y,markertype_gil,'color',Gil_colour, 'LineWidth', linwidth);
axis([9.998*10^4 10^5 1 5])
plot(0:1:10^5,10^3*ones(length(0:1:10^5),1),'k--')
hold off;
ax = gca;
ylabel('I','fontsize',fsz, 'Interpreter','Latex')
xlabel('S','fontsize',fsz, 'Interpreter','Latex')
cf = gcf;
cf.Position=[440 378 200 180];

figure;
hold on;
y = IData_JFS_Fad(8,(1:num2plot:end));
x = SData_JFS_Fad(8,(1:num2plot:end));
plot(x,y,markertype_jsf,'color',JSF_colour, 'LineWidth', linwidth);
y = IData_Gil_Fad(2,(1:num2plot:end));
x = SData_Gil_Fad(2,(1:num2plot:end));
plot(x,y,markertype_gil,'color',Gil_colour, 'LineWidth', linwidth);
axis([2.4*10^4 10^5 1 2*10^4 ])
yticks([10^0,10^2,10^4])
yticklabels({'10^0','10^2','10^4'})
plot(0:1:10^5,10^3*ones(length(0:1:10^5),1),'k--')
hold off;
ax = gca;
ax.YAxis.Scale ="log";
ax.XAxis.Scale ="log";
ylabel('I','fontsize',fsz, 'Interpreter','Latex')
xlabel('S','fontsize',fsz, 'Interpreter','Latex')
cf = gcf;
cf.Position=[440 378 200 180];

figure;
hold on;
y = IData_JFS_End(1,(1:num2plot:end));
x = SData_JFS_End(1,(1:num2plot:end));
plot(x,y,markertype_jsf,'color',JSF_colour, 'LineWidth', linwidth);
y = IData_Gil_End(1,(1:num2plot:end));
x = SData_Gil_End(1,(1:num2plot:end));
plot(x,y,markertype_gil,'color',Gil_colour, 'LineWidth', linwidth);
axis([2.4*10^4 10^5 1 2*10^4 ])
yticks([10^0,10^2,10^4])
yticklabels({'10^0','10^2','10^4'})
plot(0:1:10^5,10^3*ones(length(0:1:10^5),1),'k--')
hold off;
ax = gca;
ax.YAxis.Scale ="log";
ax.XAxis.Scale ="log";
xlabel('S','fontsize',fsz, 'Interpreter','Latex')
ylabel('I','fontsize',fsz, 'Interpreter','Latex')
cf = gcf;
cf.Position=[440 378 200 180];
