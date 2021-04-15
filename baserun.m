close all
clear all
%% Run the model and get results in the res-struct:
res = NPD(param);

plotNPD(res, param)

%% Ds
figure
y = -[res10.zAtIndex res30.zAtIndex res50.zAtIndex];
x = [res10.maxvalue res30.maxvalue res50.maxvalue];
plot(x,y,'b-','HandleVisibility','off')
hold on
plot(x(1),y(1),'b*',x(2),y(2),'r*',x(3),y(3),'g*')
legend('D = 10 [m^2/d]', 'D = 30 [m^2/d]', 'D = 50 [m^2/d]', 'Location', 'southwest')
ylabel('Depth of maximum P conc. [m]')
xlabel('Maximum P conc. [cells/m^3]')

%% D
HN = 0.0425;
HI = 30*86400;
N = res.N(3833,:);
P = res.P(3833,:);
t = res.t(3833);

I0 = 300*86400; %param.I0;
Iamplitude = 300*86400; %param.Iamplitude;
kw = 0.045; %param.kw;
kp = 6*10^-10; %param.kp;
dz = 300/100; %param.dz;
z = dz/2:dz:(300-dz/2); %param.z;

P = res.P(3833,:);
I = calclight_1(P, t, I0, Iamplitude, kw, kp, z, dz);

for i = 1:length(N)
    N_lim(i) = N(i)/(HN+N(i));
    I_lim(i) = I(i)/(HI+I(i));
end
figure(5)
hold on
plot(N_lim, -res.z, 'g--')
plot(I_lim, -res.z, 'r--')
%line(N_lim, -res.z, 'Color','b')
%line(I_lim, -res.z, 'Color','r')
ax1 = gca;
ax1_pos = ax1.Position;
ax2 = axes('Position', ax1_pos,'XAxisLocation', 'top', 'YAxisLocation', 'right');
line(P, -res.z, 'Color', 'k')
hold off
%plot(N_lim, -res.z, 'g--')
%hold on
%plot(I_lim, -res.z, 'r--')
%xline((0.01*24)/(0.04*24), 'k--')
%hold off