format compact
clear
clc
%close all
clf reset

fliplength = 17;

a = 1:fliplength;

k = 0.5;
numerator = [1,k];
denominator = [3,1.5,k];
sys = tf(numerator,denominator);
[stepsys,time]=step(sys, fliplength);
stepsys = stepsys.*180;
time = time+2;

b = interp1(time,stepsys,a,"linear");
b(end) = 180;
b(1) = 0;
b(2) = mean([b(1),b(3)]);

b = 180 - b;
stepsys = 180-stepsys;

hold on
grid on
plot(time,stepsys,"k",LineWidth=2)
scatter(a,b,"r","filled")
xlim([0,fliplength])