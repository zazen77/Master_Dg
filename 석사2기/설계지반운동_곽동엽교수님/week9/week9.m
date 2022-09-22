clc
clear

x1 = 0 : 0.001 : 0.5 * 18.5;
y1 = ((1.6 / 18.5) * x1 +0.2 ) * 43.65572606;

x2 = 0.5 * 18.5 : 0.001 : 2 * 18.5;
y2 = ((-0.6 / 18.5) * x2 + 1.3) * 43.65572606;

x3 = 2 * 18.5 : 0.001 : 4 * 18.5;
y3 = ((-0.05 / 18.5) * x3 +0.2)*43.65572606 ;

p = plot(x1,y1,x2,y2,x3,y3);
set(gca,'YDir','reverse');
xlabel('d = distance from the wall') 
ylabel('\deltav = vertical settlement at the distance d') 
legend({'\deltav = ((1.6 / 18.5) * d +0.2 ) * 43.65572606',...
    '\deltav = ((-0.6 / 18.5) * d + 1.3) * 43.65572606',...
    '\deltav = (-0.05 / 18.5) * d +0.2) * 43.65572606'},...
    'Location','southeast');