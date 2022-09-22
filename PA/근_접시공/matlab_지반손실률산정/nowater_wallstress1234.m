clc
clear

% 실제 데이터
x = [-5
-4
-3
-2
-1
0
1
3
4
5
6
7
8
10
11
12
13
14
15
16
17
18];
y = [0.006792495
0
0.006792495
0.02716998
0.01358499
-0.006792495
-0.203774851
-0.400757206
-0.427927186
-0.448304671
-0.468682157
-0.495852137
-0.509437127
-0.516229622
-0.516229622
-0.516229622
-0.523022117
-0.523022117
-0.529814612
-0.529814612
-0.536607107
-0.529814612]; 

%보간
xq1 = -5:0.1:15;

%선형보간
y0 = interp1(x,y,xq1);
p = pchip(x,y,xq1);
s = spline(x,y,xq1);
m = makima(x,y,xq1);

plot(x,y,'o',xq1,p,'-',xq1,s,'-.',xq1,m,'--')
grid on
xlabel('Elapsed Time from Arrival Date (day)') 
ylabel('Volume Loss (%)') 
legend('Sample Points','pchip',...
    'spline','makima','Location','North East')


% 
% figure(1) % 원래 데이터
% plot(x1,y1,'o')
% grid on
% xlabel('Elapsed Time from Arrival Date (day)') 
% ylabel('Volume Loss (%)') 
% legend({'Volume Loss'},...
%     'Location','northeast');
% % ylim([0 8])
%  
% 
% figure(2) % 선형 보간된 데이터
% plot(x0,y0)
% grid on
% xlabel('Elapsed Time from Arrival Date (day)') 
% ylabel('Volume Loss (%)') 
% legend({'Volume Loss'},...
%     'Location','northeast');
% 
% 
% % p = plot(x1,y1);
% % set(gca,'YDir','reverse');
% % xlabel('Elapsed Time from Arrival Date (day)') 
% % ylabel('Volume Loss (%)') 
% 
% legend({'Volume Loss'},...
%     'Location','northeast');