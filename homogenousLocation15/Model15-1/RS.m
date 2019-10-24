%%%%%%% Calculate the rotational period of skew bridge models. 
%%%% clear all command 
close all; clear all; clc ;

%%%% Define the size of the figure 
Lox=100;
Loy=100;
width=750;
height=460;

% for i=23:23
filename1=sprintf('center_top_acc.txt')
filename2=sprintf('x_acc_SW4_node_455.txt')

Acc_top=load(filename1);

acc1=386.4*Acc_top(1,:)/9.81;
acc1=acc1';
acc2=386.4*importdata(filename2)/9.81;
[T1,Spa1,Spv1,Sd1]=SPEC(0.002820874,0.01,acc1,0.05,386.4,4);
[T2,Spa2,Spv2,Sd2]=SPEC(0.002820874,0.01,acc2,0.05,386.4,4);

figure(1)
set(gcf,'Position',[Lox, Loy, width, height]);
plot (T1,Spa1,'k-',T2,Spa2,'r-','LineWidth',2.0)
set(gca, 'gridlinestyle','--','GridAlpha',1,'FontSize',14,'FontWeight', 'Normal','FontName','Arial');
xlabel('Period (s)')
ylabel('Acceleration response spectrum (g)')
title('x-acc response spectrum for center top surface node (5% damping)','FontName','Arial','FontSize',14,'FontWeight','Normal')
legend1=legend('ESSI','SW4');
set(legend1,'FontSize',14,'FontName','Arial','FontWeight','Normal')
grid on
savefig(gcf,'1-x-acc response spectrum for center top surface node');
saveas(gcf,'1-x-acc response spectrum for center top surface node','meta');
saveas(gcf,'1-x-acc response spectrum for center top surface node','png');

% for i=23:23

filename4=sprintf('z_acc_SW4_node_455.txt')

acc3=386.4*Acc_top(2,:)/9.81;
acc3=acc3';
n1=length(acc3);
acc4=-386.4*importdata(filename4)/9.81;
[T3,Spa3,Spv3,Sd3]=SPEC(0.002820874,0.01,acc3(2:n1-2),0.05,386.4,4);
[T4,Spa4,Spv4,Sd4]=SPEC(0.002820874,0.01,acc4,0.05,386.4,4);

figure(2)
set(gcf,'Position',[Lox, Loy, width, height]);
plot (T3,Spa3,'k-',T4,Spa4,'r-','LineWidth',2.0)
set(gca, 'gridlinestyle','--','GridAlpha',1,'FontSize',14,'FontWeight', 'Normal','FontName','Arial');
xlabel('Period (s)')
ylabel('Acceleration response spectrum (g)')
title('y-acc response spectrum for center top surface node (5% damping)','FontName','Arial','FontSize',14,'FontWeight','Normal')
legend1=legend('ESSI','SW4');
set(legend1,'FontSize',14,'FontName','Arial','FontWeight','Normal')
grid on
savefig(gcf,'2-y-acc response spectrum for center top surface node');
saveas(gcf,'2-y-acc response spectrum for center top surface node','meta');
saveas(gcf,'2-y-acc response spectrum for center top surface node','png');



vel_SW4_x=importdata('120_120_0_xvel.csv');
vel_SW4_y=-importdata('120_120_0_zvel.csv');

Disp_top=importdata('Center_top_disp.txt');

disp_ESSI_x=Disp_top(1,:);
disp_ESSI_y=Disp_top(2,:);
acc_ESSI_x=Acc_top(1,:);
acc_ESSI_y=Acc_top(2,:);

dt=0.002820874;

disp_SW4_x(1)=0;
acc_SW4_x(1)=0;
n=length(vel_SW4_x);
for i=1:n-1
    disp_SW4_x(i+1)=disp_SW4_x(i)+vel_SW4_x(i)*dt;
    acc_SW4_x(i+1)=(vel_SW4_x(i+1)-vel_SW4_x(i))/dt;
end

disp_SW4_y(1)=0;
acc_SW4_y(1)=0;
for i=1:n-1
    disp_SW4_y(i+1)=disp_SW4_y(i)+vel_SW4_y(i)*dt;
    acc_SW4_y(i+1)=(vel_SW4_y(i+1)-vel_SW4_y(i))/dt;
end

for i=1:n
    t(i)=(i-1)*dt;
end


figure(3)
set(gcf,'Position',[Lox, Loy, width, height]);
plot (t,disp_ESSI_x(2:n+1)-disp_ESSI_x(1),'k-',t,disp_SW4_x,'r-','LineWidth',2.0)
set(gca, 'gridlinestyle','--','GridAlpha',1,'FontSize',14,'FontWeight', 'Normal','FontName','Arial');
xlabel('Time (s)')
ylabel('Displacement (m)')
title('x-displacement response history of center top surface node','FontName','Arial','FontSize',14,'FontWeight','Normal')
legend1=legend('ESSI','SW4');
set(legend1,'FontSize',14,'FontName','Arial','FontWeight','Normal')
grid on
savefig(gcf,'3-x-displacement response history of center top surface node');
saveas(gcf,'3-x-displacement response history of center top surface node','meta');
saveas(gcf,'3-x-displacement response history of center top surface node','png');

figure(4)
set(gcf,'Position',[Lox, Loy, width, height]);
plot (t,disp_ESSI_y(2:n+1)-disp_ESSI_y(1),'k-',t,disp_SW4_y,'r-','LineWidth',2.0)
set(gca, 'gridlinestyle','--','GridAlpha',1,'FontSize',14,'FontWeight', 'Normal','FontName','Arial');
xlabel('Time (s)')
ylabel('Displacement (m)')
title('y-displacement response history of center top surface node','FontName','Arial','FontSize',14,'FontWeight','Normal')
legend1=legend('ESSI','SW4');
set(legend1,'FontSize',14,'FontName','Arial','FontWeight','Normal')
grid on
savefig(gcf,'4-y-displacement response history of center top surface node');
saveas(gcf,'4-y-displacement response history of center top surface node','meta');
saveas(gcf,'4-y-displacement response history of center top surface node','png');


figure(5)
set(gcf,'Position',[Lox, Loy, width, height]);
plot (t,acc_ESSI_x(2:n+1),'k-',t,acc_SW4_x,'r-','LineWidth',2.0)
set(gca, 'gridlinestyle','--','GridAlpha',1,'FontSize',14,'FontWeight', 'Normal','FontName','Arial');
xlabel('Time (s)')
xlim([0 20]);
ylabel('Acceleration (m/s^2)')
title('x-acceleration response history of center top surface node','FontName','Arial','FontSize',14,'FontWeight','Normal')
legend1=legend('ESSI','SW4');
set(legend1,'FontSize',14,'FontName','Arial','FontWeight','Normal')
grid on
savefig(gcf,'5-x-acceleration response history of center top surface node');
saveas(gcf,'5-x-acceleration response history of center top surface node','meta');
saveas(gcf,'5-x-acceleration response history of center top surface node','png');

figure(6)
set(gcf,'Position',[Lox, Loy, width, height]);
plot (t,acc_ESSI_y(2:n+1),'k-',t,acc_SW4_y,'r-','LineWidth',2.0)
set(gca, 'gridlinestyle','--','GridAlpha',1,'FontSize',14,'FontWeight', 'Normal','FontName','Arial');
xlabel('Time (s)')
xlim([0,20]);
ylabel('Acceleration (m/s^2)')
title('y-acceleration response history of center top surface node','FontName','Arial','FontSize',14,'FontWeight','Normal')
legend1=legend('ESSI','SW4');
set(legend1,'FontSize',14,'FontName','Arial','FontWeight','Normal')
grid on
savefig(gcf,'6-y-acceleration response history of center top surface node');
saveas(gcf,'6-y-acceleration response history of center top surface node','meta');
saveas(gcf,'6-y-acceleration response history of center top surface node','png');


vel_SW4_left_x=importdata('left_xvel.csv');
vel_SW4_right_x=importdata('right_xvel.csv');
vel_SW4_bottom_x=importdata('bot_xvel.csv');

vel_SW4_left_y=-importdata('left_zvel.csv');
vel_SW4_right_y=-importdata('right_zvel.csv');
vel_SW4_bottom_y=-importdata('bot_zvel.csv');

[disp_SW4_left_x,acc_SW4_left_x]=Dispandacc(dt, vel_SW4_left_x); 
[disp_SW4_right_x,acc_SW4_right_x]=Dispandacc(dt, vel_SW4_right_x); 
[disp_SW4_bottom_x,acc_SW4_bottom_x]=Dispandacc(dt, vel_SW4_bottom_x); 

[disp_SW4_left_y,acc_SW4_left_y]=Dispandacc(dt, vel_SW4_left_y); 
[disp_SW4_right_y,acc_SW4_right_y]=Dispandacc(dt, vel_SW4_right_y); 
[disp_SW4_bottom_y,acc_SW4_bottom_y]=Dispandacc(dt, vel_SW4_bottom_y); 

disp_ESSI_left=importdata('Left_node_disp.txt');
disp_ESSI_right=importdata('Right_node_disp.txt');
disp_ESSI_bottom=importdata('Bottom_node_disp.txt');

acc_ESSI_left=importdata('Left_node_acc.txt');
acc_ESSI_right=importdata('Right_node_acc.txt');
acc_ESSI_bottom=importdata('Bottom_node_acc.txt');

disp_ESSI_left_x=disp_ESSI_left(1,:);
disp_ESSI_left_y=disp_ESSI_left(2,:);

disp_ESSI_right_x=disp_ESSI_right(1,:);
disp_ESSI_right_y=disp_ESSI_right(2,:);

disp_ESSI_bottom_x=disp_ESSI_bottom(1,:);
disp_ESSI_bottom_y=disp_ESSI_bottom(2,:);

acc_ESSI_left_x=acc_ESSI_left(1,:);
acc_ESSI_left_y=acc_ESSI_left(2,:);

acc_ESSI_right_x=acc_ESSI_right(1,:);
acc_ESSI_right_y=acc_ESSI_right(2,:);

acc_ESSI_bottom_x=acc_ESSI_bottom(1,:);
acc_ESSI_bottom_y=acc_ESSI_bottom(2,:);

figure(7)
set(gcf,'Position',[Lox, Loy, width, height]);
plot (t,disp_ESSI_left_x(2:n+1)-disp_ESSI_left_x(1),'k-',t,disp_SW4_left_x,'r-','LineWidth',2.0)
set(gca, 'gridlinestyle','--','GridAlpha',1,'FontSize',14,'FontWeight', 'Normal','FontName','Arial');
xlabel('Time (s)')
ylabel('Displacement (m)')
title('Comparison of x-displacement response history of left node','FontName','Arial','FontSize',14,'FontWeight','Normal')
legend1=legend('ESSI','SW4');
set(legend1,'FontSize',14,'FontName','Arial','FontWeight','Normal')
grid on
savefig(gcf,'7-x-displacement response history of left node');
saveas(gcf,'7-x-displacement response history of left node','meta');
saveas(gcf,'7-x-displacement response history of left node','png');


figure(8)
set(gcf,'Position',[Lox, Loy, width, height]);
plot (t,disp_ESSI_right_x(2:n+1)-disp_ESSI_right_x(1),'k-',t,disp_SW4_right_x,'r-','LineWidth',2.0)
set(gca, 'gridlinestyle','--','GridAlpha',1,'FontSize',14,'FontWeight', 'Normal','FontName','Arial');
xlabel('Time (s)')
ylabel('Displacement (m)')
title('Comparison of x-displacement response history of right node','FontName','Arial','FontSize',14,'FontWeight','Normal')
legend1=legend('ESSI','SW4');
set(legend1,'FontSize',14,'FontName','Arial','FontWeight','Normal')
grid on
savefig(gcf,'8-x-displacement response history of right node');
saveas(gcf,'8-x-displacement response history of right node','meta');
saveas(gcf,'8-x-displacement response history of right node','png');

figure(9)
set(gcf,'Position',[Lox, Loy, width, height]);
plot (t,disp_ESSI_bottom_x(2:n+1)-disp_ESSI_bottom_x(1),'k-',t,disp_SW4_bottom_x,'r-','LineWidth',2.0)
set(gca, 'gridlinestyle','--','GridAlpha',1,'FontSize',14,'FontWeight', 'Normal','FontName','Arial');
xlabel('Time (s)')
ylabel('Displacement (m)')
title('Comparison of x-displacement response history of bottom node','FontName','Arial','FontSize',14,'FontWeight','Normal')
legend1=legend('ESSI','SW4');
set(legend1,'FontSize',14,'FontName','Arial','FontWeight','Normal')
grid on
savefig(gcf,'9-x-displacement response history of bottom node');
saveas(gcf,'9-x-displacement response history of bottom node','meta');
saveas(gcf,'9-x-displacement response history of bottom node','png');


figure(10)
set(gcf,'Position',[Lox, Loy, width, height]);
plot (t,disp_ESSI_left_y(2:n+1)-disp_ESSI_left_y(1),'k-',t,disp_SW4_left_y,'r-','LineWidth',2.0)
set(gca, 'gridlinestyle','--','GridAlpha',1,'FontSize',14,'FontWeight', 'Normal','FontName','Arial');
xlabel('Time (s)')
ylabel('Displacement (m)')
title('Comparison of y-displacement response history of left node','FontName','Arial','FontSize',14,'FontWeight','Normal')
legend1=legend('ESSI','SW4');
set(legend1,'FontSize',14,'FontName','Arial','FontWeight','Normal')
grid on
savefig(gcf,'10-y-displacement response history of left node');
saveas(gcf,'10-y-displacement response history of left node','meta');
saveas(gcf,'10-y-displacement response history of left node','png');


figure(11)
set(gcf,'Position',[Lox, Loy, width, height]);
plot (t,disp_ESSI_right_y(2:n+1)-disp_ESSI_right_y(1),'k-',t,disp_SW4_right_y,'r-','LineWidth',2.0)
set(gca, 'gridlinestyle','--','GridAlpha',1,'FontSize',14,'FontWeight', 'Normal','FontName','Arial');
xlabel('Time (s)')
ylabel('Displacement (m)')
title('Comparison of y-displacement response history of right node','FontName','Arial','FontSize',14,'FontWeight','Normal')
legend1=legend('ESSI','SW4');
set(legend1,'FontSize',14,'FontName','Arial','FontWeight','Normal')
grid on
savefig(gcf,'11-y-displacement response history of right node');
saveas(gcf,'11-y-displacement response history of right node','meta');
saveas(gcf,'11-y-displacement response history of right node','png');

figure(12)
set(gcf,'Position',[Lox, Loy, width, height]);
plot (t,disp_ESSI_bottom_y(2:n+1)-disp_ESSI_bottom_y(1),'k-',t,disp_SW4_bottom_y,'r-','LineWidth',2.0)
set(gca, 'gridlinestyle','--','GridAlpha',1,'FontSize',14,'FontWeight', 'Normal','FontName','Arial');
xlabel('Time (s)')
ylabel('Displacement (m)')
title('Comparison of y-displacement response history of bottom node','FontName','Arial','FontSize',14,'FontWeight','Normal')
legend1=legend('ESSI','SW4');
set(legend1,'FontSize',14,'FontName','Arial','FontWeight','Normal')
grid on
savefig(gcf,'12-y-displacement response history of bottom node');
saveas(gcf,'12-y-displacement response history of bottom node','meta');
saveas(gcf,'12-y-displacement response history of bottom node','png');

figure(13)
set(gcf,'Position',[Lox, Loy, width, height]);
plot (t(1:n-2),acc_ESSI_left_x(2:n-1),'k-',t,acc_SW4_left_x,'r-','LineWidth',2.0)
set(gca, 'gridlinestyle','--','GridAlpha',1,'FontSize',14,'FontWeight', 'Normal','FontName','Arial');
xlabel('Time (s)')
xlim([0 20]);
ylabel('Acceleration (m/s^2)')
title('Comparison of x-acceleration response history of left node','FontName','Arial','FontSize',14,'FontWeight','Normal')
legend1=legend('ESSI','SW4');
set(legend1,'FontSize',14,'FontName','Arial','FontWeight','Normal')
grid on
savefig(gcf,'13-x-acceleration response history of left node');
saveas(gcf,'13-x-acceleration response history of left node','meta');
saveas(gcf,'13-x-acceleration response history of left node','png');

figure(14)
set(gcf,'Position',[Lox, Loy, width, height]);
plot (t(1:n-2),acc_ESSI_right_x(2:n-1),'k-',t,acc_SW4_right_x,'r-','LineWidth',2.0)
set(gca, 'gridlinestyle','--','GridAlpha',1,'FontSize',14,'FontWeight', 'Normal','FontName','Arial');
xlabel('Time (s)')
xlim([0 20]);
ylabel('Acceleration (m/s^2)')
title('Comparison of x-acceleration response history of right node','FontName','Arial','FontSize',14,'FontWeight','Normal')
legend1=legend('ESSI','SW4');
set(legend1,'FontSize',14,'FontName','Arial','FontWeight','Normal')
grid on
savefig(gcf,'14-x-acceleration response history of right node');
saveas(gcf,'14-x-acceleration response history of right node','meta');
saveas(gcf,'14-x-acceleration response history of right node','png');

figure(15)
set(gcf,'Position',[Lox, Loy, width, height]);
plot (t(1:n-2),acc_ESSI_bottom_x(2:n-1),'k-',t,acc_SW4_bottom_x,'r-','LineWidth',2.0)
set(gca, 'gridlinestyle','--','GridAlpha',1,'FontSize',14,'FontWeight', 'Normal','FontName','Arial');
xlabel('Time (s)')
xlim([0 20]);
ylabel('Acceleration (m/s^2)')
title('Comparison of x-acceleration response history of bottom node','FontName','Arial','FontSize',14,'FontWeight','Normal')
legend1=legend('ESSI','SW4');
set(legend1,'FontSize',14,'FontName','Arial','FontWeight','Normal')
grid on
savefig(gcf,'15-x-acceleration response history of bottom node');
saveas(gcf,'15-x-acceleration response history of bottom node','meta');
saveas(gcf,'15-x-acceleration response history of bottom node','png');

figure(16)
set(gcf,'Position',[Lox, Loy, width, height]);
plot (t(1:n-2),acc_ESSI_left_y(2:n-1),'k-',t,acc_SW4_left_y,'r-','LineWidth',2.0)
set(gca, 'gridlinestyle','--','GridAlpha',1,'FontSize',14,'FontWeight', 'Normal','FontName','Arial');
xlabel('Time (s)')
xlim([0,20]);
ylabel('Acceleration (m/s^2)')
title('Comparison of y-acceleration response history of left node','FontName','Arial','FontSize',14,'FontWeight','Normal')
legend1=legend('ESSI','SW4');
set(legend1,'FontSize',14,'FontName','Arial','FontWeight','Normal')
grid on
savefig(gcf,'16-y-acceleration response history of left node');
saveas(gcf,'16-y-acceleration response history of left node','meta');
saveas(gcf,'16-y-acceleration response history of left node','png');

figure(17)
set(gcf,'Position',[Lox, Loy, width, height]);
plot (t(1:n-2),acc_ESSI_right_y(2:n-1),'k-',t,acc_SW4_right_y,'r-','LineWidth',2.0)
set(gca, 'gridlinestyle','--','GridAlpha',1,'FontSize',14,'FontWeight', 'Normal','FontName','Arial');
xlabel('Time (s)')
xlim([0,20]);
ylabel('Acceleration (m/s^2)')
title('Comparison of y-acceleration response history of right node','FontName','Arial','FontSize',14,'FontWeight','Normal')
legend1=legend('ESSI','SW4');
set(legend1,'FontSize',14,'FontName','Arial','FontWeight','Normal')
grid on
savefig(gcf,'17-y-acceleration response history of right node');
saveas(gcf,'17-y-acceleration response history of right node','meta');
saveas(gcf,'17-y-acceleration response history of right node','png');

figure(18)
set(gcf,'Position',[Lox, Loy, width, height]);
plot (t(1:n),acc_ESSI_bottom_y(2:n+1),'k-',t,acc_SW4_bottom_y,'r-','LineWidth',2.0)
set(gca, 'gridlinestyle','--','GridAlpha',1,'FontSize',14,'FontWeight', 'Normal','FontName','Arial');
xlabel('Time (s)')
xlim([0,20]);
ylabel('Acceleration (m/s^2)')
title('Comparison of y-acceleration response history of bottom node','FontName','Arial','FontSize',14,'FontWeight','Normal')
legend1=legend('ESSI','SW4');
set(legend1,'FontSize',14,'FontName','Arial','FontWeight','Normal')
grid on
savefig(gcf,'18-y-acceleration response history of bottom node');
saveas(gcf,'18-y-acceleration response history of bottom node','meta');
saveas(gcf,'18-y-acceleration response history of bottom node','png');


[T5,Spa5,Spv5,Sd5]=SPEC(0.002820874,0.01,acc_ESSI_left_x(2:n-2),0.05,386.4,4);
[T6,Spa6,Spv6,Sd6]=SPEC(0.002820874,0.01,acc_SW4_left_x,0.05,386.4,4);

[T7,Spa7,Spv7,Sd7]=SPEC(0.002820874,0.01,acc_ESSI_left_y(2:n-2),0.05,386.4,4);
[T8,Spa8,Spv8,Sd8]=SPEC(0.002820874,0.01,acc_SW4_left_y,0.05,386.4,4);

[T9,Spa9,Spv9,Sd9]=SPEC(0.002820874,0.01,acc_ESSI_right_x(2:n-2),0.05,386.4,4);
[T10,Spa10,Spv10,Sd10]=SPEC(0.002820874,0.01,acc_SW4_right_x,0.05,386.4,4);

[T11,Spa11,Spv11,Sd11]=SPEC(0.002820874,0.01,acc_ESSI_right_y(2:n-2),0.05,386.4,4);
[T12,Spa12,Spv12,Sd12]=SPEC(0.002820874,0.01,acc_SW4_right_y,0.05,386.4,4);

[T13,Spa13,Spv13,Sd13]=SPEC(0.002820874,0.01,acc_ESSI_bottom_x(2:n-2),0.05,386.4,4);
[T14,Spa14,Spv14,Sd14]=SPEC(0.002820874,0.01,acc_SW4_bottom_x,0.05,386.4,4);

[T15,Spa15,Spv15,Sd15]=SPEC(0.002820874,0.01,acc_ESSI_bottom_y(2:n-2),0.05,386.4,4);
[T16,Spa16,Spv16,Sd16]=SPEC(0.002820874,0.01,acc_SW4_bottom_y,0.05,386.4,4);

figure(19)
set(gcf,'Position',[Lox, Loy, width, height]);
plot (T5,Spa5,'k-',T6,Spa6,'r-','LineWidth',2.0)
set(gca, 'gridlinestyle','--','GridAlpha',1,'FontSize',14,'FontWeight', 'Normal','FontName','Arial');
xlabel('Period (s)')
ylabel('Acceleration response spectrum (g)')
title('x-acc response spectrum for left node (5% damping)','FontName','Arial','FontSize',14,'FontWeight','Normal')
legend1=legend('ESSI','SW4');
set(legend1,'FontSize',14,'FontName','Arial','FontWeight','Normal')
grid on
savefig(gcf,'19-x-acc response spectrum for left node');
saveas(gcf,'19-x-acc response spectrum for left node','meta');
saveas(gcf,'19-x-acc response spectrum for left node','png');


figure(20)
set(gcf,'Position',[Lox, Loy, width, height]);
plot (T7,Spa7,'k-',T8,Spa8,'r-','LineWidth',2.0)
set(gca, 'gridlinestyle','--','GridAlpha',1,'FontSize',14,'FontWeight', 'Normal','FontName','Arial');
xlabel('Period (s)')
ylabel('Acceleration response spectrum (g)')
title('y-acc response spectrum for left node (5% damping)','FontName','Arial','FontSize',14,'FontWeight','Normal')
legend1=legend('ESSI','SW4');
set(legend1,'FontSize',14,'FontName','Arial','FontWeight','Normal')
grid on
savefig(gcf,'20-y-acc response spectrum for left node');
saveas(gcf,'20-y-acc response spectrum for left node','meta');
saveas(gcf,'20-y-acc response spectrum for left node','png');

figure(21)
set(gcf,'Position',[Lox, Loy, width, height]);
plot (T9,Spa9,'k-',T10,Spa10,'r-','LineWidth',2.0)
set(gca, 'gridlinestyle','--','GridAlpha',1,'FontSize',14,'FontWeight', 'Normal','FontName','Arial');
xlabel('Period (s)')
ylabel('Acceleration response spectrum (g)')
title('x-acc response spectrum for right node (5% damping)','FontName','Arial','FontSize',14,'FontWeight','Normal')
legend1=legend('ESSI','SW4');
set(legend1,'FontSize',14,'FontName','Arial','FontWeight','Normal')
grid on
savefig(gcf,'21-x-acc response spectrum for right node');
saveas(gcf,'21-x-acc response spectrum for right node','meta');
saveas(gcf,'21-x-acc response spectrum for right node','png');

figure(22)
set(gcf,'Position',[Lox, Loy, width, height]);
plot (T11,Spa11,'k-',T12,Spa12,'r-','LineWidth',2.0)
set(gca, 'gridlinestyle','--','GridAlpha',1,'FontSize',14,'FontWeight', 'Normal','FontName','Arial');
xlabel('Period (s)')
ylabel('Acceleration response spectrum (g)')
title('y-acc response spectrum for right node (5% damping)','FontName','Arial','FontSize',14,'FontWeight','Normal')
legend1=legend('ESSI','SW4');
set(legend1,'FontSize',14,'FontName','Arial','FontWeight','Normal')
grid on
savefig(gcf,'22-y-acc response spectrum for right node');
saveas(gcf,'22-y-acc response spectrum for right node','meta');
saveas(gcf,'22-y-acc response spectrum for right node','png');

figure(23)
set(gcf,'Position',[Lox, Loy, width, height]);
plot (T13,Spa13,'k-',T14,Spa14,'r-','LineWidth',2.0)
set(gca, 'gridlinestyle','--','GridAlpha',1,'FontSize',14,'FontWeight', 'Normal','FontName','Arial');
xlabel('Period (s)')
ylabel('Acceleration response spectrum (g)')
title('x-acc response spectrum for bottom node (5% damping)','FontName','Arial','FontSize',14,'FontWeight','Normal')
legend1=legend('ESSI','SW4');
set(legend1,'FontSize',14,'FontName','Arial','FontWeight','Normal')
grid on
savefig(gcf,'23-x-acc response spectrum for bottom node');
saveas(gcf,'23-x-acc response spectrum for bottom node','meta');
saveas(gcf,'23-x-acc response spectrum for bottom node','png');

figure(24)
set(gcf,'Position',[Lox, Loy, width, height]);
plot (T15,Spa15,'k-',T16,Spa16,'r-','LineWidth',2.0)
set(gca, 'gridlinestyle','--','GridAlpha',1,'FontSize',14,'FontWeight', 'Normal','FontName','Arial');
xlabel('Period (s)')
ylabel('Acceleration response spectrum (g)')
title('y-acc response spectrum for bottom node (5% damping)','FontName','Arial','FontSize',14,'FontWeight','Normal')
legend1=legend('ESSI','SW4');
set(legend1,'FontSize',14,'FontName','Arial','FontWeight','Normal')
grid on
savefig(gcf,'24-y-acc response spectrum for bottom node');
saveas(gcf,'24-y-acc response spectrum for bottom node','meta');
saveas(gcf,'24-y-acc response spectrum for bottom node','png');

disp_ESSI_node_190=importdata('Node_190_disp.txt');
acc_ESSI_node_190=importdata('Node_190_acc.txt');

% disp_ESSI_node_12=importdata('disp_node_12.txt');
% acc_ESSI_node_12=importdata('acc_node_12.txt');
% disp_ESSI_node_12_x=disp_ESSI_node_12(1,:);
% disp_ESSI_node_12_y=disp_ESSI_node_12(2,:);
% acc_ESSI_node_12_x=acc_ESSI_node_12(1,:);
% acc_ESSI_node_12_y=acc_ESSI_node_12(2,:);

disp_ESSI_node_190_x=disp_ESSI_node_190(1,:);
disp_ESSI_node_190_y=disp_ESSI_node_190(2,:);
acc_ESSI_node_190_x=acc_ESSI_node_190(1,:);
acc_ESSI_node_190_y=acc_ESSI_node_190(2,:);

figure(25)
set(gcf,'Position',[Lox, Loy, width, height]);
plot (t,disp_ESSI_node_190_x(2:n+1)-disp_ESSI_node_190_x(1),'k-',t(1:n-2),disp_ESSI_left_x(2:n-1)-disp_ESSI_left_x(1),'-r','LineWidth',2.0)
set(gca, 'gridlinestyle','--','GridAlpha',1,'FontSize',14,'FontWeight', 'Normal','FontName','Arial');
xlabel('Time (s)')
ylabel('Displacement (m)')
title('x-displacement response history of node 190 and left node','FontName','Arial','FontSize',14,'FontWeight','Normal')
legend1=legend('Node 190','Left node');
set(legend1,'FontSize',14,'FontName','Arial','FontWeight','Normal')
grid on
savefig(gcf,'25-x-displacement response history of node 190 and left node');
saveas(gcf,'25-x-displacement response history of node 190 and left node','meta');
saveas(gcf,'25-x-displacement response history of node 190 and left node','png');

figure(26)
set(gcf,'Position',[Lox, Loy, width, height]);
plot (t,disp_ESSI_node_190_y(2:n+1)-disp_ESSI_node_190_y(1),'k-',t(1:n-2),disp_ESSI_left_y(2:n-1)-disp_ESSI_left_y(1),'-r','LineWidth',2.0)
set(gca, 'gridlinestyle','--','GridAlpha',1,'FontSize',14,'FontWeight', 'Normal','FontName','Arial');
xlabel('Time (s)')
ylabel('Displacement (m)')
title('y-displacement response history of node 190 and left node','FontName','Arial','FontSize',14,'FontWeight','Normal')
legend1=legend('Node 190','Left node');
set(legend1,'FontSize',14,'FontName','Arial','FontWeight','Normal')
grid on
savefig(gcf,'26-y-displacement response history of node 190 and left node');
saveas(gcf,'26-y-displacement response history of node 190 and left node','meta');
saveas(gcf,'26-y-displacement response history of node 190 and left node','png');

figure(27)
set(gcf,'Position',[Lox, Loy, width, height]);
plot (t,acc_ESSI_node_190_x(2:n+1),'k-',t(1:n-2),acc_ESSI_left_x(2:n-1),'-r','LineWidth',2.0)
set(gca, 'gridlinestyle','--','GridAlpha',1,'FontSize',14,'FontWeight', 'Normal','FontName','Arial');
xlabel('Time (s)')
xlim([0 20]);
ylabel('Acceleration (m/s^2)')
title('x-acceleration response history of node 190 and left node','FontName','Arial','FontSize',14,'FontWeight','Normal')
legend1=legend('Node 190','Left node');
set(legend1,'FontSize',14,'FontName','Arial','FontWeight','Normal')
grid on
savefig(gcf,'27-x-acceleration response history of node 190 and left node');
saveas(gcf,'27-x-acceleration response history of node 190 and left node','meta');
saveas(gcf,'27-x-acceleration response history of node 190 and left node','png');

figure(28)
set(gcf,'Position',[Lox, Loy, width, height]);
plot (t,acc_ESSI_node_190_y(2:n+1),'k-',t(1:n-2),acc_ESSI_left_y(2:n-1),'-r','LineWidth',2.0)
set(gca, 'gridlinestyle','--','GridAlpha',1,'FontSize',14,'FontWeight', 'Normal','FontName','Arial');
xlabel('Time (s)')
xlim([0 20]);
ylabel('Acceleration (m/s^2)')
title('y-acceleration response history of node 190 and left node','FontName','Arial','FontSize',14,'FontWeight','Normal')
legend1=legend('Node 190','Left node');
set(legend1,'FontSize',14,'FontName','Arial','FontWeight','Normal')
grid on
savefig(gcf,'28-y-acceleration response history of node 190 and left node');
saveas(gcf,'28-y-acceleration response history of node 190 and left node','meta');
saveas(gcf,'28-y-acceleration response history of node 190 and left node','png');

figure(29)
set(gcf,'Position',[Lox, Loy, width, height]);
plot (t,disp_ESSI_left_x(2:n+1)-disp_ESSI_left_x(1),'k-',t,disp_ESSI_right_x(2:n+1)-disp_ESSI_right_x(1),'r-',t,disp_ESSI_x(2:n+1),'m-',t,disp_ESSI_bottom_x(2:n+1),'b-','LineWidth',2.0)
set(gca, 'gridlinestyle','--','GridAlpha',1,'FontSize',14,'FontWeight', 'Normal','FontName','Arial');
xlabel('Time (s)')
ylabel('Displacement (m)')
title('Comparison of x-displacement response history of ESSI model','FontName','Arial','FontSize',14,'FontWeight','Normal')
legend1=legend('Left','Right','Top','Bottom');
set(legend1,'FontSize',14,'FontName','Arial','FontWeight','Normal')
grid on
savefig(gcf,'29-x-displacement response history of ESSI model');
saveas(gcf,'29-x-displacement response history of ESSI model','meta');
saveas(gcf,'29-x-displacement response history of ESSI model','png');


figure(291)
set(gcf,'Position',[Lox, Loy, width, height]);
plot (t,disp_ESSI_left_y(2:n+1)-disp_ESSI_left_y(1),'k-',t,disp_ESSI_right_y(2:n+1)-disp_ESSI_right_y(1),'r-',t,disp_ESSI_y(2:n+1)-disp_ESSI_y(1),'m-',t,disp_ESSI_bottom_y(2:n+1)-disp_ESSI_bottom_y(1),'b-','LineWidth',2.0)
set(gca, 'gridlinestyle','--','GridAlpha',1,'FontSize',14,'FontWeight', 'Normal','FontName','Arial');
xlabel('Time (s)')
ylabel('Displacement (m)')
title('Comparison of y-displacement response history of ESSI model','FontName','Arial','FontSize',14,'FontWeight','Normal')
legend1=legend('Left','Right','Top','Bottom');
set(legend1,'FontSize',14,'FontName','Arial','FontWeight','Normal')
grid on
savefig(gcf,'291-y-displacement response history of ESSI model');
saveas(gcf,'291-y-displacement response history of ESSI model','meta');
saveas(gcf,'291-y-displacement response history of ESSI model','png');


figure(30)
set(gcf,'Position',[Lox, Loy, width, height]);
plot (t,disp_SW4_left_x,'k-',t,disp_SW4_right_x,'r-',t,disp_SW4_x,'m-',t,disp_SW4_bottom_x,'b-','LineWidth',2.0)
set(gca, 'gridlinestyle','--','GridAlpha',1,'FontSize',14,'FontWeight', 'Normal','FontName','Arial');
xlabel('Time (s)')
ylabel('Displacement (m)')
title('Comparison of x-displacement response history of SW4 model','FontName','Arial','FontSize',14,'FontWeight','Normal')
legend1=legend('Left','Right','Top','Bottom');
set(legend1,'FontSize',14,'FontName','Arial','FontWeight','Normal')
grid on
savefig(gcf,'30-x-displacement response history of SW4 model');
saveas(gcf,'30-x-displacement response history of SW4 model','meta');
saveas(gcf,'30-x-displacement response history of SW4 model','png');

figure(31)
set(gcf,'Position',[Lox, Loy, width, height]);
plot (t,acc_ESSI_left_x(2:n+1)-acc_ESSI_left_x(1),'k-',t,acc_ESSI_right_x(2:n+1)-acc_ESSI_right_x(1),'r-',t,acc_ESSI_x(2:n+1),'m-',t,acc_ESSI_bottom_x(2:n+1),'b-','LineWidth',2.0)
set(gca, 'gridlinestyle','--','GridAlpha',1,'FontSize',14,'FontWeight', 'Normal','FontName','Arial');
xlabel('Time (s)')
xlim([0 20]);
ylabel('acc (m/^2)')
title('Comparison of x-acc response history of ESSI model','FontName','Arial','FontSize',14,'FontWeight','Normal')
legend1=legend('Left','Right','Top','Bottom');
set(legend1,'FontSize',14,'FontName','Arial','FontWeight','Normal')
grid on
savefig(gcf,'31-x-acc response history of ESSI model');
saveas(gcf,'31-x-acc response history of ESSI model','meta');
saveas(gcf,'31-x-acc response history of ESSI model','png');

figure(311)
set(gcf,'Position',[Lox, Loy, width, height]);
plot (t,acc_ESSI_left_y(2:n+1)-acc_ESSI_left_y(1),'k-',t,acc_ESSI_right_y(2:n+1)-acc_ESSI_right_y(1),'r-',t,acc_ESSI_y(2:n+1),'m-',t,acc_ESSI_bottom_y(2:n+1),'b-','LineWidth',2.0)
set(gca, 'gridlinestyle','--','GridAlpha',1,'FontSize',14,'FontWeight', 'Normal','FontName','Arial');
xlabel('Time (s)')
xlim([0 20]);
ylabel('acc (m/^2)')
title('Comparison of y-acc response history of ESSI model','FontName','Arial','FontSize',14,'FontWeight','Normal')
legend1=legend('Left','Right','Top','Bottom');
set(legend1,'FontSize',14,'FontName','Arial','FontWeight','Normal')
grid on
savefig(gcf,'311-y-acc response history of ESSI model');
saveas(gcf,'311-y-acc response history of ESSI model','meta');
saveas(gcf,'311-y-acc response history of ESSI model','png');


figure(32)
set(gcf,'Position',[Lox, Loy, width, height]);
plot (t,acc_SW4_left_x,'k-',t,acc_SW4_right_x,'r-',t,acc_SW4_x,'m-',t,acc_SW4_bottom_x,'b-','LineWidth',2.0)
set(gca, 'gridlinestyle','--','GridAlpha',1,'FontSize',14,'FontWeight', 'Normal','FontName','Arial');
xlabel('Time (s)')
ylabel('Acceleration (m/s^2)')
title('Comparison of x-acc response history of SW4 model','FontName','Arial','FontSize',14,'FontWeight','Normal')
legend1=legend('Left','Right','Top','Bottom');
set(legend1,'FontSize',14,'FontName','Arial','FontWeight','Normal')
grid on
savefig(gcf,'32-x-acc response history of SW4 model');
saveas(gcf,'32-x-acc response history of SW4 model','meta');
saveas(gcf,'32-x-acc response history of SW4 model','png');