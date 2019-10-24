i%%%%%%% Calculate the rotational period of skew bridge models. 
%%%% clear all command 
close all; clear all; clc ;
 
%%%% Define the size of the figure 
Lox=50;
Loy=50;
width=1100;
height=650;
 
dt=0.002821;



for i=1:11
    filename1=sprintf('x_vel_sw4_X%d.csv',i)
    vel_SW4_L_x(i,:)=importdata(filename1);
    filename2=sprintf('x_vel_sw4_X%d.csv',i+20)
    vel_SW4_R_x(i,:)=importdata(filename2);
    filename3=sprintf('x_vel_sw4_Z%d.csv',i)
    vel_SW4_M_x(i,:)=importdata(filename3);
        
    filename1=sprintf('z_vel_sw4_X%d.csv',i)
    vel_SW4_L_z(i,:)=importdata(filename1);
    filename2=sprintf('z_vel_sw4_X%d.csv',i+20)
    vel_SW4_R_z(i,:)=importdata(filename2);
    filename3=sprintf('z_vel_sw4_Z%d.csv',i)
    vel_SW4_M_z(i,:)=importdata(filename3);             
 
    end
for i=1:17
    filename4=sprintf('z_vel_sw4_T%d.csv',i)
    vel_SW4_T_z(i,:)=importdata(filename4);
end 


% for i=1:11
% [disp_SW4_L_x(i,:),acc_SW4_L_x(i,:)]=Dispandacc(dt, vel_SW4_L_x(i,:)); 
% [disp_SW4_R_x(i,:),acc_SW4_R_x(i,:)]=Dispandacc(dt, vel_SW4_R_x(i,:)); 
% [disp_SW4_M_x(i,:),acc_SW4_M_x(i,:)]=Dispandacc(dt, vel_SW4_M_x(i,:)); 
% 
% [disp_SW4_L_z(i,:),acc_SW4_L_z(i,:)]=Dispandacc(dt, vel_SW4_L_z(i,:)); 
% [disp_SW4_R_z(i,:),acc_SW4_R_z(i,:)]=Dispandacc(dt, vel_SW4_R_z(i,:)); 
% [disp_SW4_M_z(i,:),acc_SW4_M_z(i,:)]=Dispandacc(dt, vel_SW4_M_z(i,:)); 
% end 
% 
% for i=1:17
% [disp_SW4_T_z(i,:),acc_SW4_T_z(i,:)]=Dispandacc(dt, vel_SW4_T_z(i,:));  
% end


for i=1:11
    
    filename1=sprintf('x_disp_sw4_X%d.csv',i)
    disp_SW4_L_x(i,:)=importdata(filename1);
    filename2=sprintf('x_disp_sw4_X%d.csv',i+20)
    disp_SW4_R_x(i,:)=importdata(filename2);
    filename3=sprintf('x_disp_sw4_Z%d.csv',i)
    disp_SW4_M_x(i,:)=importdata(filename3);
        
    filename1=sprintf('z_disp_sw4_X%d.csv',i)
    disp_SW4_L_z(i,:)=importdata(filename1);
    filename2=sprintf('z_disp_sw4_X%d.csv',i+20)
    disp_SW4_R_z(i,:)=importdata(filename2);
    filename3=sprintf('z_disp_sw4_Z%d.csv',i)
    disp_SW4_M_z(i,:)=importdata(filename3);
    
    
    filename1=sprintf('x_acc_sw4_X%d.csv',i)
    acc_SW4_L_x(i,:)=importdata(filename1);
    filename2=sprintf('x_acc_sw4_X%d.csv',i+20)
    acc_SW4_R_x(i,:)=importdata(filename2);
    filename3=sprintf('x_acc_sw4_Z%d.csv',i)
    acc_SW4_M_x(i,:)=importdata(filename3);
        
    filename1=sprintf('z_acc_sw4_X%d.csv',i)
    acc_SW4_L_z(i,:)=importdata(filename1);
    filename2=sprintf('z_acc_sw4_X%d.csv',i+20)
    acc_SW4_R_z(i,:)=importdata(filename2);
    filename3=sprintf('z_acc_sw4_Z%d.csv',i)
    acc_SW4_M_z(i,:)=importdata(filename3);
    
end 

for i=1:17
    filename1=sprintf('z_disp_sw4_T%d.csv',i)
    disp_SW4_T_z(i,:)=importdata(filename1);   
end



for i=1:11
filename1=sprintf('X%d_disp_2.5.txt',i)
disp_ESSI_L(((3*i-2):3*i),:)=importdata(filename1);
filename2=sprintf('X%d_disp_2.5.txt',i+20)
disp_ESSI_R(((3*i-2):3*i),:)=importdata(filename2);
filename3=sprintf('Z%d_disp_2.5.txt',i)
disp_ESSI_M(((3*i-2):3*i),:)=importdata(filename3);

filename1=sprintf('X%d_acc_2.5.txt',i)
acc_ESSI_L(((3*i-2):3*i),:)=importdata(filename1);
filename2=sprintf('X%d_acc_2.5.txt',i+20)
acc_ESSI_R(((3*i-2):3*i),:)=importdata(filename2);
filename3=sprintf('Z%d_acc_2.5.txt',i)
acc_ESSI_M(((3*i-2):3*i),:)=importdata(filename3);

end 
for i=1:17
filename4=sprintf('T%d_disp_2.5.txt',i)
disp_ESSI_T(((3*i-2):3*i),:)=importdata(filename4);
end

for i=1:11  
disp_ESSI_L_x(i,:)=disp_ESSI_L(3*i-2,:);
disp_ESSI_R_x(i,:)=disp_ESSI_R(3*i-2,:);
disp_ESSI_M_x(i,:)=disp_ESSI_M(3*i-2,:);

disp_ESSI_L_z(i,:)=disp_ESSI_L(3*i,:);
disp_ESSI_R_z(i,:)=disp_ESSI_R(3*i,:);
disp_ESSI_M_z(i,:)=disp_ESSI_M(3*i,:);

acc_ESSI_L_x(i,:)=acc_ESSI_L(3*i-2,:);
acc_ESSI_R_x(i,:)=acc_ESSI_R(3*i-2,:);
acc_ESSI_M_x(i,:)=acc_ESSI_M(3*i-2,:);

acc_ESSI_L_z(i,:)=acc_ESSI_L(3*i,:);
acc_ESSI_R_z(i,:)=acc_ESSI_R(3*i,:);
acc_ESSI_M_z(i,:)=acc_ESSI_M(3*i,:);
end 

for i=1:17
disp_ESSI_T_z(i,:)=disp_ESSI_T(3*i,:);
end

for i=1:8
 theta_ESSI_T(i,:)=(disp_ESSI_T_z(i,:)-disp_ESSI_T_z(18-i,:))/(20*(9-i));
 theta_sw4_T(i,:)=(disp_SW4_T_z(i,:)-disp_SW4_T_z(18-i,:))/(20*(9-i));
end 


N=length(acc_ESSI_L_x(1,:));
n=N-5;
 
for i=1:N
    t(i)=(i-1)*dt;
end

N1=length(acc_SW4_L_x(1,:));
n1=N1-5;
 
for i=1:N1
    t1(i)=(i-1)*dt;
end

t0=0.0;
t10=0.0;
t11=20.0;
a=2.0;
b=1.0;
c=0.2;


figure(11)
set(gcf,'Position',[Lox, Loy, width, height]);
hold on
p=plot (t,1.0*disp_ESSI_L_x(11,:),'b-',t1-t0,1.0*disp_SW4_L_x(11,:),'r-')
p(1).LineWidth=a;
p(2).LineWidth=b;
set(gca, 'gridlinestyle','--','GridAlpha',1,'FontSize',14,'FontWeight', 'Normal','FontName','Arial');
xlabel('Time (s)')
ylabel('Horizontal displacement (m)')
xlim([t10 t11]);
title('Center top surface','FontName','Arial','FontSize',14,'FontWeight','Normal')
legend1=legend('ESSI','SW4');
set(legend1,'FontSize',14,'FontName','Arial','FontWeight','Normal','Location','Southeast')
grid on
savefig(gcf,'11-comprison of disp history in x(horizontal) direction of center top surface');
saveas(gcf,'11-comprison of disp history in x(horizontal) direction of center top surface','meta');
saveas(gcf,'11-comprison of disp history in x(horizontal) direction of center top surface','png');

figure(12)
set(gcf,'Position',[Lox, Loy, width, height]);
hold on
p=plot (t,1.0*disp_ESSI_L_z(11,:),'b-',t1-t0,1.0*disp_SW4_L_z(11,:),'r-')
p(1).LineWidth=a;
p(2).LineWidth=b;
set(gca, 'gridlinestyle','--','GridAlpha',1,'FontSize',14,'FontWeight', 'Normal','FontName','Arial');
xlabel('Time (s)')
ylabel('Vertical displacement (m)')
xlim([t10 t11]);
title('Center top surface','FontName','Arial','FontSize',14,'FontWeight','Normal')
legend1=legend('ESSI','SW4');
set(legend1,'FontSize',14,'FontName','Arial','FontWeight','Normal','Location','Southeast')
grid on
savefig(gcf,'12-comprison of disp history in z(vertical) direction of center top surface');
saveas(gcf,'12-comprison of disp history in z(vertical) direction of center top surface','meta');
saveas(gcf,'12-comprison of disp history in z(vertical) direction of center top surface','png');

figure(13)
set(gcf,'Position',[Lox, Loy, width, height]);
hold on
p=plot (t,1.0*acc_ESSI_L_x(11,:),'b-',t1-t0,1.0*acc_SW4_L_x(11,:),'r-')
p(1).LineWidth=a;
p(2).LineWidth=b;
set(gca, 'gridlinestyle','--','GridAlpha',1,'FontSize',14,'FontWeight', 'Normal','FontName','Arial');
xlabel('Time (s)')
ylabel('Horizontal acceleration (m/s^2)')
xlim([t10 t11]);
title('Center top surface','FontName','Arial','FontSize',14,'FontWeight','Normal')
legend1=legend('ESSI','SW4');
set(legend1,'FontSize',14,'FontName','Arial','FontWeight','Normal','Location','Southeast')
grid on
savefig(gcf,'13-comprison of acc history in x(horizontal) direction of center top surface');
saveas(gcf,'13-comprison of acc response history in x(horizontal) direction of center top surface','meta');
saveas(gcf,'13-comprison of acc response history in x(horizontal) direction of center top surface','png');

figure(14)
set(gcf,'Position',[Lox, Loy, width, height]);
hold on
p=plot (t,1.0*acc_ESSI_L_z(11,:),'b-',t1-t0,1.0*acc_SW4_L_z(11,:),'r-')
p(1).LineWidth=a;
p(2).LineWidth=b;
set(gca, 'gridlinestyle','--','GridAlpha',1,'FontSize',14,'FontWeight', 'Normal','FontName','Arial');
xlabel('Time (s)')
ylabel('Vertical acceleration (m/s^2)')
xlim([t10 t11]);
title('Center top surface','FontName','Arial','FontSize',14,'FontWeight','Normal')
legend1=legend('ESSI','SW4');
set(legend1,'FontSize',14,'FontName','Arial','FontWeight','Normal','Location','Southeast')
grid on
savefig(gcf,'14-comprison of acc history in z(vertical) direction of center top surface');
saveas(gcf,'14-comprison of acc response history in z(vertical) direction of center top surface','meta');
saveas(gcf,'14-comprison of acc response history in z(vertical) direction of center top surface','png');


figure(15)
set(gcf,'Position',[Lox, Loy, width, height]);
hold on
p=plot (t,1.0*theta_ESSI_T(6,:),'b-',t1-t0,1.0*theta_sw4_T(6,:),'r-')
p(1).LineWidth=a;
p(2).LineWidth=b;
set(gca, 'gridlinestyle','--','GridAlpha',1,'FontSize',14,'FontWeight', 'Normal','FontName','Arial');
xlabel('Time (s)')
ylabel('Roation (rad)')
xlim([t10 t11]);
title('Surface rotation','FontName','Arial','FontSize',14,'FontWeight','Normal')
legend1=legend('ESSI','SW4');
set(legend1,'FontSize',14,'FontName','Arial','FontWeight','Normal','Location','Southeast')
grid on
savefig(gcf,'15-comprison of surface rotation response history');
saveas(gcf,'15-comprison of surface rotation response history','meta');
saveas(gcf,'15-comprison of surface rotation response history','png');





% figure(1)
% set(gcf,'Position',[Lox, Loy, width, height]);
% hold on
% for i=1:11
% p=plot (t,1.0*disp_ESSI_L_x(i,:)+c*(i+1),'b-',t1-t0,1.0*disp_SW4_L_x(i,:)+c*(i+1),'r-')
% p(1).LineWidth=a;
% p(2).LineWidth=b;
% end
% set(gca, 'gridlinestyle','--','GridAlpha',1,'FontSize',14,'FontWeight', 'Normal','FontName','Arial');
% xlabel('Time (s)')
% ylabel('Horizontal displacement (m)')
% xlim([t10 t11]);
% title('Left column','FontName','Arial','FontSize',14,'FontWeight','Normal')
% legend1=legend('ESSI','SW4');
% set(legend1,'FontSize',14,'FontName','Arial','FontWeight','Normal','Location','Southeast')
% grid on
% savefig(gcf,'1-comprison of displacement response history in x direction of left column');
% saveas(gcf,'1-comprison of displacement response history in x direction of left column','meta');
% saveas(gcf,'1-comprison of displacement response history in x direction of left column','png');
% 
% figure(2)
% set(gcf,'Position',[Lox, Loy, width, height]);
% hold on
% for i=1:11
% p=plot (t,1.0*disp_ESSI_L_z(i,:)+c*(i+1),'b-',t1-t0,1.0*disp_SW4_L_z(i,:)+c*(i+1),'r-')
% p(1).LineWidth=a;
% p(2).LineWidth=b;
% end
% set(gca, 'gridlinestyle','--','GridAlpha',1,'FontSize',14,'FontWeight', 'Normal','FontName','Arial');
% xlabel('Time (s)')
% ylabel('Vertical displacement (m)')
% xlim([t10 t11]);
% title('Left column','FontName','Arial','FontSize',14,'FontWeight','Normal')
% legend1=legend('ESSI','SW4');
% set(legend1,'FontSize',14,'FontName','Arial','FontWeight','Normal','Location','Southeast')
% grid on
% savefig(gcf,'2-comprison of displacement response history in z direction of left column');
% saveas(gcf,'2-comprison of displacement response history in z direction of left column','meta');
% saveas(gcf,'2-comprison of displacement response history in z direction of left column','png');
% 
% figure(3)
% set(gcf,'Position',[Lox, Loy, width, height]);
% hold on
% for i=1:11
% p=plot (t,1.0*disp_ESSI_M_x(i,:)+c*(i+1),'b-',t1-t0,1.0*disp_SW4_M_x(i,:)+c*(i+1),'r-')
% p(1).LineWidth=a;
% p(2).LineWidth=b;
% end
% set(gca, 'gridlinestyle','--','GridAlpha',1,'FontSize',14,'FontWeight', 'Normal','FontName','Arial');
% xlabel('Time (s)')
% ylabel('Horizontal displacement (m)')
% xlim([t10 t11]);
% title('Center column','FontName','Arial','FontSize',14,'FontWeight','Normal')
% legend1=legend('ESSI','SW4');
% set(legend1,'FontSize',14,'FontName','Arial','FontWeight','Normal','Location','Southeast')
% grid on
% savefig(gcf,'3-comprison of displacement response history in x direction of center column');
% saveas(gcf,'3-comprison of displacement response history in x direction of center column','meta');
% saveas(gcf,'3-comprison of displacement response history in x direction of center column','png');
% 
% figure(4)
% set(gcf,'Position',[Lox, Loy, width, height]);
% hold on
% for i=1:11
% p=plot (t,1.0*disp_ESSI_R_z(i,:)+c*(i+1),'b-',t1-t0,1.0*disp_SW4_R_z(i,:)+c*(i+1),'r-')
% p(1).LineWidth=a;
% p(2).LineWidth=b;
% end
% set(gca, 'gridlinestyle','--','GridAlpha',1,'FontSize',14,'FontWeight', 'Normal','FontName','Arial');
% xlabel('Time (s)')
% ylabel('Vertical displacement (m)')
% xlim([t10 t11]);
% title('Center column','FontName','Arial','FontSize',14,'FontWeight','Normal')
% legend1=legend('ESSI','SW4');
% set(legend1,'FontSize',14,'FontName','Arial','FontWeight','Normal','Location','Southeast')
% grid on
% savefig(gcf,'4-comprison of displacement response history in x direction of center column');
% saveas(gcf,'4-comprison of displacement response history in x direction of center column','meta');
% saveas(gcf,'4-comprison of displacement response history in x direction of center column','png');
% 
% figure(5)
% set(gcf,'Position',[Lox, Loy, width, height]);
% hold on
% for i=1:11
% p=plot (t,1.0*disp_ESSI_R_x(i,:)+c*(i+1),'b-',t1-t0,1.0*disp_SW4_R_x(i,:)+c*(i+1),'r-')
% p(1).LineWidth=a;
% p(2).LineWidth=b;
% end
% set(gca, 'gridlinestyle','--','GridAlpha',1,'FontSize',14,'FontWeight', 'Normal','FontName','Arial');
% xlabel('Time (s)')
% ylabel('Horizontal displacement (m)')
% xlim([t10 t11]);
% title('Right column','FontName','Arial','FontSize',14,'FontWeight','Normal')
% legend1=legend('ESSI','SW4');
% set(legend1,'FontSize',14,'FontName','Arial','FontWeight','Normal','Location','Southeast')
% grid on
% savefig(gcf,'1-comprison of displacement response history in x direction of right column');
% saveas(gcf,'1-comprison of displacement response history in x direction of right column','meta');
% saveas(gcf,'1-comprison of displacement response history in x direction of right column','png');
% 
% figure(6)
% set(gcf,'Position',[Lox, Loy, width, height]);
% hold on
% for i=1:11
% p=plot (t,1.0*disp_ESSI_R_z(i,:)+c*(i+1),'b-',t1-t0,1.0*disp_SW4_R_z(i,:)+c*(i+1),'r-')
% p(1).LineWidth=a;
% p(2).LineWidth=b;
% end
% set(gca, 'gridlinestyle','--','GridAlpha',1,'FontSize',14,'FontWeight', 'Normal','FontName','Arial');
% xlabel('Time (s)')
% ylabel('Vertical displacement (m)')
% xlim([t10 t11]);
% title('Right column','FontName','Arial','FontSize',14,'FontWeight','Normal')
% legend1=legend('ESSI','SW4');
% set(legend1,'FontSize',14,'FontName','Arial','FontWeight','Normal','Location','Southeast')
% grid on
% savefig(gcf,'6-comprison of displacement response history in z direction of right column');
% saveas(gcf,'6-comprison of displacement response history in z direction of right column','meta');
% saveas(gcf,'6-comprison of displacement response history in z direction of right column','png');
% 
% figure(7)
% set(gcf,'Position',[Lox, Loy, width, height]);
% hold on
% for i=1:11
% p=plot (t,1.0*acc_ESSI_L_x(i,:)+10*(i+1),'b-',t1-t0,1.0*acc_SW4_L_x(i,:)+10*(i+1),'r-')
% p(1).LineWidth=a;
% p(2).LineWidth=b;
% end
% set(gca, 'gridlinestyle','--','GridAlpha',1,'FontSize',14,'FontWeight', 'Normal','FontName','Arial');
% xlabel('Time (s)')
% ylabel('Horizontal acceleration (m/s^2)')
% xlim([t10 t11]);
% title('Left column','FontName','Arial','FontSize',14,'FontWeight','Normal')
% legend1=legend('ESSI','SW4');
% set(legend1,'FontSize',14,'FontName','Arial','FontWeight','Normal','Location','Northeast')
% grid on
% savefig(gcf,'7-comprison of acceleration response history in x direction of left column');
% saveas(gcf,'7-comprison of acceleration response history in x direction of left column','meta');
% saveas(gcf,'7-comprison of acceleration response history in x direction of left column','png');
%  
% figure(8)
% set(gcf,'Position',[Lox, Loy, width, height]);
% hold on
% for i=1:11
% p=plot (t,1.0*acc_ESSI_L_z(i,:)+10*(i+1),'b-',t1-t0,1.0*acc_SW4_L_z(i,:)+10*(i+1),'r-')
% p(1).LineWidth=a;
% p(2).LineWidth=b;
% end
% set(gca, 'gridlinestyle','--','GridAlpha',1,'FontSize',14,'FontWeight', 'Normal','FontName','Arial');
% xlabel('Time (s)')
% ylabel('Vertical acceleration (m/s^2)')
% xlim([t10 t11]);
% title('Left column','FontName','Arial','FontSize',14,'FontWeight','Normal')
% legend1=legend('ESSI','SW4');
% set(legend1,'FontSize',14,'FontName','Arial','FontWeight','Normal','Location','Northeast')
% grid on
% savefig(gcf,'8-comprison of acceleration response history in z direction of left column');
% saveas(gcf,'8-comprison of acceleration response history in z direction of left column','meta');
% saveas(gcf,'8-comprison of acceleration response history in z direction of left column','png');
%  
% figure(9)
% set(gcf,'Position',[Lox, Loy, width, height]);
% hold on
% for i=1:11
% p=plot (t,1.0*acc_ESSI_M_x(i,:)+10*(i+1),'b-',t1-t0,1.0*acc_SW4_M_x(i,:)+10*(i+1),'r-')
% p(1).LineWidth=a;
% p(2).LineWidth=b;
% end
% set(gca, 'gridlinestyle','--','GridAlpha',1,'FontSize',14,'FontWeight', 'Normal','FontName','Arial');
% xlabel('Time (s)')
% ylabel('Horizontal acceleration (m/s^2)')
% xlim([t10 t11]);
% title('Center column','FontName','Arial','FontSize',14,'FontWeight','Normal')
% legend1=legend('ESSI','SW4');
% set(legend1,'FontSize',14,'FontName','Arial','FontWeight','Normal','Location','Northeast')
% grid on
% savefig(gcf,'9-comprison of acceleration response history in x direction of center column');
% saveas(gcf,'9-comprison of acceleration response history in x direction of center column','meta');
% saveas(gcf,'9-comprison of acceleration response history in x direction of center column','png');
%  
% figure(10)
% set(gcf,'Position',[Lox, Loy, width, height]);
% hold on
% for i=1:11
% p=plot (t,1.0*acc_ESSI_R_z(i,:)+10*(i+1),'b-',t1-t0,1.0*acc_SW4_R_z(i,:)+10*(i+1),'r-')
% p(1).LineWidth=a;
% p(2).LineWidth=b;
% end
% set(gca, 'gridlinestyle','--','GridAlpha',1,'FontSize',14,'FontWeight', 'Normal','FontName','Arial');
% xlabel('Time (s)')
% ylabel('Vertical acceleration (m/s^2)')
% xlim([t10 t11]);
% title('Center column','FontName','Arial','FontSize',14,'FontWeight','Normal')
% legend1=legend('ESSI','SW4');
% set(legend1,'FontSize',14,'FontName','Arial','FontWeight','Normal','Location','Northeast')
% grid on
% savefig(gcf,'10-comprison of acceleration response history in x direction of center column');
% saveas(gcf,'10-comprison of acceleration response history in x direction of center column','meta');
% saveas(gcf,'10-comprison of acceleration response history in x direction of center column','png');
%  
% figure(11)
% set(gcf,'Position',[Lox, Loy, width, height]);
% hold on
% for i=1:11
% p=plot (t,1.0*acc_ESSI_R_x(i,:)+10*(i+1),'b-',t1-t0,1.0*acc_SW4_R_x(i,:)+10*(i+1),'r-')
% p(1).LineWidth=a;
% p(2).LineWidth=b;
% end
% set(gca, 'gridlinestyle','--','GridAlpha',1,'FontSize',14,'FontWeight', 'Normal','FontName','Arial');
% xlabel('Time (s)')
% ylabel('Horizontal acceleration (m/s^2)')
% xlim([t10 t11]);
% title('Right column','FontName','Arial','FontSize',14,'FontWeight','Normal')
% legend1=legend('ESSI','SW4');
% set(legend1,'FontSize',14,'FontName','Arial','FontWeight','Normal','Location','Northeast')
% grid on
% savefig(gcf,'11-comprison of acceleration response history in x direction of right column');
% saveas(gcf,'11-comprison of acceleration response history in x direction of right column','meta');
% saveas(gcf,'11-comprison of acceleration response history in x direction of right column','png');
%  
% figure(12)
% set(gcf,'Position',[Lox, Loy, width, height]);
% hold on
% for i=1:11
% p=plot (t,1.0*acc_ESSI_R_z(i,:)+10*(i+1),'b-',t1-t0,1.0*acc_SW4_R_z(i,:)+10*(i+1),'r-')
% p(1).LineWidth=a;
% p(2).LineWidth=b;
% end
% set(gca, 'gridlinestyle','--','GridAlpha',1,'FontSize',14,'FontWeight', 'Normal','FontName','Arial');
% xlabel('Time (s)')
% ylabel('Vertical acceleration (m/s^2)')
% xlim([t10 t11]);
% title('Right column','FontName','Arial','FontSize',14,'FontWeight','Normal')
% legend1=legend('ESSI','SW4');
% set(legend1,'FontSize',14,'FontName','Arial','FontWeight','Normal','Location','Northeast')
% grid on
% savefig(gcf,'12-comprison of acceleration response history in z direction of right column');
% saveas(gcf,'12-comprison of acceleration response history in z direction of right column','meta');
% saveas(gcf,'12-comprison of acceleration response history in z direction of right column','png');
% 
% 
% % for i=1:8
% % deltaT1=0.002821;
% % order1=6;
% % cutout1=2.0;
% % [e1,f1]=butter(order1,cutout1/(1/(2*deltaT1)),'low');
% % theta_ESSI_T(i,:)=filtfilt(e1,f1,theta_ESSI_T(i,:));
% % theta_sw4_T(i,:)=filtfilt(e1,f1,theta_sw4_T(i,:));
% % end
% 
% 
% figure(13)
% set(gcf,'Position',[Lox, Loy, width, height]);
% hold on
% for i=1:8
% p=plot (t,1.0*theta_ESSI_T(i,:)+0.0001*(i+1),'b-',t1-t0,1.0*theta_sw4_T(i,:)+0.0001*(i+1),'r-')
% p(1).LineWidth=a;
% p(2).LineWidth=b;
% end
% 
% set(gca, 'gridlinestyle','--','GridAlpha',1,'FontSize',14,'FontWeight', 'Normal','FontName','Arial');
% ylabel('Time (s)')
% xlabel('Rotation (rad/s^2)')
% xlim([t10 t11]);
% title('Right column','FontName','Arial','FontSize',14,'FontWeight','Normal')
% legend1=legend('ESSI','SW4');
% set(legend1,'FontSize',14,'FontName','Arial','FontWeight','Normal','Location','Northeast')
% grid on
% savefig(gcf,'12-comprison of rotation response history of right column');
% saveas(gcf,'12-comprison of rotation response history of right column','meta');
% saveas(gcf,'12-comprison of rotation response history of right column','png');
% 
% 
% % figure(20)
% % set(gcf,'Position',[Lox, Loy, width, height]);
% % subplot(2,2,1)
% % p=plot (t(1:n),Rot_ESSI_y_1(1:n),'b-',t1(1:N1-1)-t0,Rot_SW4_y_1(1:N1-1),'r-')
% % p(1).LineWidth=a;
% % p(2).LineWidth=b;
% % set(gca, 'gridlinestyle','--','GridAlpha',1,'FontSize',12,'FontWeight', 'Normal','FontName','Arial');
% % xlabel('Time (s)')
% % xlim([t10 t11]);
% % ylabel('Rotation about y axis (rad)')
% % title('Rotation 1 about y axis','FontName','Arial','FontSize',12,'FontWeight','Normal')
% % legend1=legend('ESSI','SW4');
% % set(legend1,'FontSize',12,'FontName','Arial','FontWeight','Normal','Location','Northwest','orientation','Horizontal')
% % grid on
% %  
% % subplot(2,2,2)
% % p=plot (t(1:n),Rot_ESSI_y_2(1:n),'b-',t1(1:N1-1)-t0,Rot_SW4_y_2(1:N1-1),'r-')
% % p(1).LineWidth=a;
% % p(2).LineWidth=b;
% % set(gca, 'gridlinestyle','--','GridAlpha',1,'FontSize',12,'FontWeight', 'Normal','FontName','Arial');
% % xlabel('Time (s)')
% % xlim([t10 t11]);
% % ylabel('Rotation about y axis (rad)')
% % title('Rotation 2 about y axis','FontName','Arial','FontSize',12,'FontWeight','Normal')
% % grid on
% %  
% % subplot(2,2,3)
% % p=plot (t(1:n),Rot_ESSI_y_3(1:n),'b-',t1(1:N1-1)-t0,Rot_SW4_y_3(1:N1-1),'r-')
% % p(1).LineWidth=a;
% % p(2).LineWidth=b;
% % set(gca, 'gridlinestyle','--','GridAlpha',1,'FontSize',12,'FontWeight', 'Normal','FontName','Arial');
% % xlabel('Time (s)')
% % xlim([t10 t11]);
% % ylabel('Rotation about y axis (rad)')
% % title('Rotation 3 about y axis','FontName','Arial','FontSize',12,'FontWeight','Normal')
% % grid on
% %  
% % subplot(2,2,4)
% % p=plot (t(1:n),Rot_ESSI_y_4(1:n),'b-',t1(1:N1-1)-t0,Rot_SW4_y_4(1:N1-1),'r-')
% % p(1).LineWidth=a;
% % p(2).LineWidth=b;
% % set(gca, 'gridlinestyle','--','GridAlpha',1,'FontSize',12,'FontWeight', 'Normal','FontName','Arial');
% % xlabel('Time (s)')
% % xlim([t10 t11]);
% % ylabel('Rotation about y axis (rad)')
% % title('Rotation 4 about y axis','FontName','Arial','FontSize',12,'FontWeight','Normal')
% % grid on
% % savefig(gcf,'20-comprison of rotation response history about y axis');
% % saveas(gcf,'20-comprison of rotation response history about y axis','meta');
% % saveas(gcf,'20-comprison of rotation response history about y axis','png');
% %  
% %  
% % figure(21)
% % set(gcf,'Position',[Lox, Loy, width, height]);
% % subplot(2,2,1)
% % p=plot (t(1:n),Rot_ESSI_x_1(1:n),'b-',t1(1:N1-1),Rot_SW4_x_1(1:N1-1),'r-')
% % p(1).LineWidth=a;
% % p(2).LineWidth=b;
% % set(gca, 'gridlinestyle','--','GridAlpha',1,'FontSize',12,'FontWeight', 'Normal','FontName','Arial');
% % xlabel('Time (s)')
% % xlim([t10 t11]);
% % ylabel('Rotation about x axis (rad)')
% % title('Rotation 1 about x axis','FontName','Arial','FontSize',12,'FontWeight','Normal')
% % legend1=legend('ESSI','SW4');
% % set(legend1,'FontSize',12,'FontName','Arial','FontWeight','Normal','Location','Northwest','orientation','Horizontal')
% % grid on
% %  
% % subplot(2,2,2)
% % p=plot (t(1:n),Rot_ESSI_x_2(1:n),'b-',t1(1:N1-1),Rot_SW4_x_2(1:N1-1),'r-')
% % p(1).LineWidth=a;
% % p(2).LineWidth=b;
% % set(gca, 'gridlinestyle','--','GridAlpha',1,'FontSize',12,'FontWeight', 'Normal','FontName','Arial');
% % xlabel('Time (s)')
% % xlim([t10 t11]);
% % ylabel('Rotation about x axis (rad)')
% % title('Rotation 2 about x axis','FontName','Arial','FontSize',12,'FontWeight','Normal')
% % grid on
% %  
% % subplot(2,2,3)
% % p=plot (t(1:n),Rot_ESSI_x_3(1:n),'b-',t1(1:N1-1),Rot_SW4_x_3(1:N1-1),'r-')
% % p(1).LineWidth=a;
% % p(2).LineWidth=b;
% % set(gca, 'gridlinestyle','--','GridAlpha',1,'FontSize',12,'FontWeight', 'Normal','FontName','Arial');
% % xlabel('Time (s)')
% % xlim([t10 t11]);
% % ylabel('Rotation about x axis (rad)')
% % title('Rotation 3 about x axis','FontName','Arial','FontSize',12,'FontWeight','Normal')
% % grid on
% %  
% % subplot(2,2,4)
% % p=plot (t(1:n),Rot_ESSI_x_4(1:n),'b-',t1(1:N1-1),Rot_SW4_x_4(1:N1-1),'r-')
% % p(1).LineWidth=a;
% % p(2).LineWidth=b;
% % set(gca, 'gridlinestyle','--','GridAlpha',1,'FontSize',12,'FontWeight', 'Normal','FontName','Arial');
% % xlabel('Time (s)')
% % xlim([t10 t11]);
% % ylabel('Rotation about x axis (rad)')
% % title('Rotation 4 about x axis','FontName','Arial','FontSize',12,'FontWeight','Normal')
% % grid on
% % savefig(gcf,'21-comprison of rotation response history about x axis');
% % saveas(gcf,'21-comprison of rotation response history about x axis','meta');
% % saveas(gcf,'21-comprison of rotation response history about x axis','png');
% % 
% % 
% % t10=0.0;
% % t11=20.0;
% % figure(1)
% % set(gcf,'Position',[Lox, Loy, width, height]);
% % subplot(3,2,1)
% % p=plot (t,disp_ESSI_x1_x,'b-',t1(1:N1-1)-t0,disp_SW4_x1_x(1:N1-1),'r-')
% % p(1).LineWidth=a;
% % p(2).LineWidth=b;
% % set(gca, 'gridlinestyle','--','GridAlpha',1,'FontSize',14,'FontWeight', 'Normal','FontName','Arial');
% % % xlabel('Time (s)')
% % % ylabel('Displacement (m)')
% % xlim([t10 t11]);
% % title('x1','FontName','Arial','FontSize',14,'FontWeight','Normal')
% % legend1=legend('ESSI','SW4');
% % set(legend1,'FontSize',14,'FontName','Arial','FontWeight','Normal','Location','Northeast')
% % grid on
% %  
% % subplot(3,2,2)
% % p=plot (t,disp_ESSI_x2_x,'b-',t1(1:N1-1)-t0,disp_SW4_x2_x(1:N1-1),'r-')
% % p(1).LineWidth=a;
% % p(2).LineWidth=b;
% % set(gca, 'gridlinestyle','--','GridAlpha',1,'FontSize',14,'FontWeight', 'Normal','FontName','Arial');
% % % xlabel('Time (s)')
% % % ylabel('Displacement (m)')
% % xlim([t10 t11]);
% % title('x2','FontName','Arial','FontSize',14,'FontWeight','Normal')
% % grid on
% %  
% %  
% % subplot(3,2,3)
% % p=plot (t,disp_ESSI_y1_x,'b-',t1(1:N1-1)-t0,disp_SW4_y1_x(1:N1-1),'r-')
% % p(1).LineWidth=a;
% % p(2).LineWidth=b;
% % set(gca, 'gridlinestyle','--','GridAlpha',1,'FontSize',14,'FontWeight', 'Normal','FontName','Arial');
% % % xlabel('Time (s)')
% % ylabel('Displacement (m)')
% % xlim([t10 t11]);
% % title('y1','FontName','Arial','FontSize',14,'FontWeight','Normal')
% % grid on
% %  
% % subplot(3,2,4)
% % p=plot (t,disp_ESSI_y2_x,'b-',t1(1:N1-1)-t0,disp_SW4_y2_x(1:N1-1),'r-')
% % p(1).LineWidth=a;
% % p(2).LineWidth=b;
% % set(gca, 'gridlinestyle','--','GridAlpha',1,'FontSize',14,'FontWeight', 'Normal','FontName','Arial');
% % % xlabel('Time (s)')
% % ylabel('Displacement (m)')
% % xlim([t10 t11]);
% % title('y2','FontName','Arial','FontSize',14,'FontWeight','Normal')
% % grid on
% %  
% % subplot(3,2,5)
% % p=plot (t,disp_ESSI_z1_x,'b-',t1(1:N1-1)-t0,disp_SW4_z1_x(1:N1-1),'r-')
% % p(1).LineWidth=a;
% % p(2).LineWidth=b;
% % set(gca, 'gridlinestyle','--','GridAlpha',1,'FontSize',14,'FontWeight', 'Normal','FontName','Arial');
% % xlabel('Time (s)')
% % % ylabel('Displacement (m)')
% % xlim([t10 t11]);
% % title('z1','FontName','Arial','FontSize',14,'FontWeight','Normal')
% % grid on
% %  
% % subplot(3,2,6)
% % p=plot (t,disp_ESSI_z2_x,'b-',t1(1:N1-1)-t0,disp_SW4_z2_x(1:N1-1),'r-')
% % p(1).LineWidth=a;
% % p(2).LineWidth=b;
% % set(gca, 'gridlinestyle','--','GridAlpha',1,'FontSize',14,'FontWeight', 'Normal','FontName','Arial');
% % xlabel('Time (s)')
% % % ylabel('Displacement (m)')
% % xlim([t10 t11]);
% % title('z2','FontName','Arial','FontSize',14,'FontWeight','Normal')
% % grid on
% % savefig(gcf,'1-comprison of displacement response history in x direction');
% % saveas(gcf,'1-comprison of displacement response history in x direction','meta');
% % saveas(gcf,'1-comprison of displacement response history in x direction','png');
% %  
% %  
% % figure(2)
% % set(gcf,'Position',[Lox, Loy, width, height]);
% % subplot(3,2,1)
% % p=plot (t,disp_ESSI_x1_y,'b-',t1(1:N1-1)-t0,disp_SW4_x1_y(1:N1-1),'r-')
% % p(1).LineWidth=a;
% % p(2).LineWidth=b;
% % set(gca, 'gridlinestyle','--','GridAlpha',1,'FontSize',14,'FontWeight', 'Normal','FontName','Arial');
% % % xlabel('Time (s)')
% % % ylabel('Displacement (m)')
% % xlim([t10 t11]);
% % title('x1','FontName','Arial','FontSize',14,'FontWeight','Normal')
% % legend1=legend('ESSI','SW4');
% % set(legend1,'FontSize',14,'FontName','Arial','FontWeight','Normal','Location','Northeast')
% % grid on
% %  
% % subplot(3,2,2)
% % p=plot (t,disp_ESSI_x2_y,'b-',t1(1:N1-1)-t0,disp_SW4_x2_y(1:N1-1),'r-')
% % p(1).LineWidth=a;
% % p(2).LineWidth=b;
% % set(gca, 'gridlinestyle','--','GridAlpha',1,'FontSize',14,'FontWeight', 'Normal','FontName','Arial');
% % % xlabel('Time (s)')
% % % ylabel('Displacement (m)')
% % xlim([t10 t11]);
% % title('x2','FontName','Arial','FontSize',14,'FontWeight','Normal')
% % grid on
% %  
% %  
% % subplot(3,2,3)
% % p=plot (t,disp_ESSI_y1_y,'b-',t1(1:N1-1)-t0,disp_SW4_y1_y(1:N1-1),'r-')
% % p(1).LineWidth=a;
% % p(2).LineWidth=b;
% % set(gca, 'gridlinestyle','--','GridAlpha',1,'FontSize',14,'FontWeight', 'Normal','FontName','Arial');
% % % xlabel('Time (s)')
% % ylabel('Displacement (m)')
% % xlim([t10 t11]);
% % title('y1','FontName','Arial','FontSize',14,'FontWeight','Normal')
% % grid on
% %  
% % subplot(3,2,4)
% % p=plot (t,disp_ESSI_y2_y,'b-',t1(1:N1-1)-t0,disp_SW4_y2_y(1:N1-1),'r-')
% % p(1).LineWidth=a;
% % p(2).LineWidth=b;
% % set(gca, 'gridlinestyle','--','GridAlpha',1,'FontSize',14,'FontWeight', 'Normal','FontName','Arial');
% % % xlabel('Time (s)')
% % ylabel('Displacement (m)')
% % xlim([t10 t11]);
% % title('y2','FontName','Arial','FontSize',14,'FontWeight','Normal')
% % grid on
% %  
% % subplot(3,2,5)
% % p=plot (t,disp_ESSI_z1_y,'b-',t1(1:N1-1)-t0,disp_SW4_z1_y(1:N1-1),'r-')
% % p(1).LineWidth=a;
% % p(2).LineWidth=b;
% % set(gca, 'gridlinestyle','--','GridAlpha',1,'FontSize',14,'FontWeight', 'Normal','FontName','Arial');
% % xlabel('Time (s)')
% % % ylabel('Displacement (m)')
% % xlim([t10 t11]);
% % title('z1','FontName','Arial','FontSize',14,'FontWeight','Normal')
% % grid on
% %  
% % subplot(3,2,6)
% % p=plot (t,disp_ESSI_z2_y,'b-',t1(1:N1-1)-t0,disp_SW4_z2_y(1:N1-1),'r-')
% % p(1).LineWidth=a;
% % p(2).LineWidth=b;
% % set(gca, 'gridlinestyle','--','GridAlpha',1,'FontSize',14,'FontWeight', 'Normal','FontName','Arial');
% % xlabel('Time (s)')
% % % ylabel('Displacement (m)')
% % title('z2','FontName','Arial','FontSize',14,'FontWeight','Normal')
% % xlim([t10 t11]);
% % grid on
% % savefig(gcf,'2-comprison of displacement response history in y direction');
% % saveas(gcf,'2-comprison of displacement response history in y direction','meta');
% % saveas(gcf,'2-comprison of displacement response history in y direction','png');
% %  
% %  
% % figure(3)
% % set(gcf,'Position',[Lox, Loy, width, height]);
% % subplot(3,2,1)
% % p=plot (t,disp_ESSI_x1_z,'b-',t1(1:N1-1)-t0,disp_SW4_x1_z(1:N1-1),'r-')
% % p(1).LineWidth=a;
% % p(2).LineWidth=b;
% % set(gca, 'gridlinestyle','--','GridAlpha',1,'FontSize',14,'FontWeight', 'Normal','FontName','Arial');
% % % xlabel('Time (s)')
% % % ylabel('Displacement (m)')
% % xlim([t10 t11]);
% % title('x1','FontName','Arial','FontSize',14,'FontWeight','Normal')
% % legend1=legend('ESSI','SW4');
% % set(legend1,'FontSize',14,'FontName','Arial','FontWeight','Normal','Location','Southeast')
% % grid on
% %  
% % subplot(3,2,2)
% % p=plot (t,disp_ESSI_x2_z,'b-',t1(1:N1-1)-t0,disp_SW4_x2_z(1:N1-1),'r-')
% % p(1).LineWidth=a;
% % p(2).LineWidth=b;
% % set(gca, 'gridlinestyle','--','GridAlpha',1,'FontSize',14,'FontWeight', 'Normal','FontName','Arial');
% % % xlabel('Time (s)')
% % % ylabel('Displacement (m)')
% % xlim([t10 t11]);
% % title('x2','FontName','Arial','FontSize',14,'FontWeight','Normal')
% % grid on
% %  
% %  
% % subplot(3,2,3)
% % p=plot (t,disp_ESSI_y1_z,'b-',t1(1:N1-1)-t0,disp_SW4_y1_z(1:N1-1),'r-')
% % p(1).LineWidth=a;
% % p(2).LineWidth=b;
% % set(gca, 'gridlinestyle','--','GridAlpha',1,'FontSize',14,'FontWeight', 'Normal','FontName','Arial');
% % % xlabel('Time (s)')
% % ylabel('Displacement (m)')
% % xlim([t10 t11]);
% % title('y1','FontName','Arial','FontSize',14,'FontWeight','Normal')
% % grid on
% %  
% % subplot(3,2,4)
% % p=plot (t,disp_ESSI_y2_z,'b-',t1(1:N1-1)-t0,disp_SW4_y2_z(1:N1-1),'r-')
% % p(1).LineWidth=a;
% % p(2).LineWidth=b;
% % set(gca, 'gridlinestyle','--','GridAlpha',1,'FontSize',14,'FontWeight', 'Normal','FontName','Arial');
% % % xlabel('Time (s)')
% % ylabel('Displacement (m)')
% % xlim([t10 t11]);
% % title('y2','FontName','Arial','FontSize',14,'FontWeight','Normal')
% % grid on
% %  
% % subplot(3,2,5)
% % p=plot (t,disp_ESSI_z1_z,'b-',t1(1:N1-1)-t0,disp_SW4_z1_z(1:N1-1),'r-')
% % p(1).LineWidth=a;
% % p(2).LineWidth=b;
% % set(gca, 'gridlinestyle','--','GridAlpha',1,'FontSize',14,'FontWeight', 'Normal','FontName','Arial');
% % xlabel('Time (s)')
% % % ylabel('Displacement (m)')
% % xlim([t10 t11]);
% % title('z1','FontName','Arial','FontSize',14,'FontWeight','Normal')
% % grid on
% %  
% % subplot(3,2,6)
% % p=plot (t,disp_ESSI_z2_z,'b-',t1(1:N1-1)-t0,disp_SW4_z2_z(1:N1-1),'r-')
% % p(1).LineWidth=a;
% % p(2).LineWidth=b;
% % set(gca, 'gridlinestyle','--','GridAlpha',1,'FontSize',14,'FontWeight', 'Normal','FontName','Arial');
% % xlabel('Time (s)')
% % % ylabel('Displacement (m)')
% % xlim([t10 t11]);
% % title('z2','FontName','Arial','FontSize',14,'FontWeight','Normal')
% % grid on
% % savefig(gcf,'3-comprison of displacement response history in z direction');
% % saveas(gcf,'3-comprison of displacement response history in z direction','meta');
% % saveas(gcf,'3-comprison of displacement response history in z direction','png');
% % 
% % 
% % figure(4)
% % set(gcf,'Position',[Lox, Loy, width, height]);
% % subplot(3,2,1)
% % p=plot (t(1:n),acc_ESSI_x1_x(1:n),'b-',t1(1:N1-1)-t0,acc_SW4_x1_x(1:N1-1),'r-')
% % p(1).LineWidth=a;
% % p(2).LineWidth=b;
% % set(gca, 'gridlinestyle','--','GridAlpha',1,'FontSize',14,'FontWeight', 'Normal','FontName','Arial');
% % % xlabel('Time (s)')
% % % ylabel('Acceleration (m/s^2)')
% % xlim([t10 t11]);
% % title('x1','FontName','Arial','FontSize',14,'FontWeight','Normal')
% % legend1=legend('ESSI','SW4');
% % set(legend1,'FontSize',14,'FontName','Arial','FontWeight','Normal','Location','Northeast')
% % grid on
% %  
% % subplot(3,2,2)
% % p=plot (t(1:n),acc_ESSI_x2_x(1:n),'b-',t1(1:N1-1)-t0,acc_SW4_x2_x(1:N1-1),'r-')
% % p(1).LineWidth=a;
% % p(2).LineWidth=b;
% % set(gca, 'gridlinestyle','--','GridAlpha',1,'FontSize',14,'FontWeight', 'Normal','FontName','Arial');
% % % xlabel('Time (s)')
% % % ylabel('Acceleration (m/s^2)')
% % xlim([t10 t11]);
% % title('x2','FontName','Arial','FontSize',14,'FontWeight','Normal')
% % grid on
% %  
% %  
% % subplot(3,2,3)
% % p=plot (t(1:n),acc_ESSI_y1_x(1:n),'b-',t1(1:N1-1)-t0,acc_SW4_y1_x(1:N1-1),'r-')
% % p(1).LineWidth=a;
% % p(2).LineWidth=b;
% % set(gca, 'gridlinestyle','--','GridAlpha',1,'FontSize',14,'FontWeight', 'Normal','FontName','Arial');
% % % xlabel('Time (s)')
% % ylabel('Acceleration (m/s^2)')
% % xlim([t10 t11]);
% % title('y1','FontName','Arial','FontSize',14,'FontWeight','Normal')
% % grid on
% %  
% % subplot(3,2,4)
% % p=plot (t(1:n),acc_ESSI_y2_x(1:n),'b-',t1(1:N1-1)-t0,acc_SW4_y2_x(1:N1-1),'r-')
% % p(1).LineWidth=a;
% % p(2).LineWidth=b;
% % set(gca, 'gridlinestyle','--','GridAlpha',1,'FontSize',14,'FontWeight', 'Normal','FontName','Arial');
% % % xlabel('Time (s)')
% % ylabel('Acceleration (m/s^2)')
% % xlim([t10 t11]);
% % title('y2','FontName','Arial','FontSize',14,'FontWeight','Normal')
% % grid on
% %  
% % subplot(3,2,5)
% % p=plot (t(1:n),acc_ESSI_z1_x(1:n),'b-',t1(1:N1-1)-t0,acc_SW4_z1_x(1:N1-1),'r-')
% % p(1).LineWidth=a;
% % p(2).LineWidth=b;
% % set(gca, 'gridlinestyle','--','GridAlpha',1,'FontSize',14,'FontWeight', 'Normal','FontName','Arial');
% % xlabel('Time (s)')
% % % ylabel('Acceleration (m/s^2)')
% % xlim([t10 t11]);
% % title('z1','FontName','Arial','FontSize',14,'FontWeight','Normal')
% % grid on
% %  
% % subplot(3,2,6)
% % p=plot (t(1:n),acc_ESSI_z2_x(1:n),'b-',t1(1:N1-1)-t0,acc_SW4_z2_x(1:N1-1),'r-')
% % p(1).LineWidth=a;
% % p(2).LineWidth=b;
% % set(gca, 'gridlinestyle','--','GridAlpha',1,'FontSize',14,'FontWeight', 'Normal','FontName','Arial');
% % xlabel('Time (s)')
% % % ylabel('Acceleration (m/s^2)')
% % xlim([t10 t11]);
% % title('z2','FontName','Arial','FontSize',14,'FontWeight','Normal')
% % grid on
% % savefig(gcf,'4-comprison of acceleration response history in x direction');
% % saveas(gcf,'4-comprison of acceleration response history in x direction','meta');
% % saveas(gcf,'4-comprison of acceleration response history in x direction','png');
% % 
% % for i=2:n-1
% %     acc_interpolate_ESSI_x1_x(i)=(disp_ESSI_x1_x(i+1)+disp_ESSI_x1_x(i-1)-2*disp_ESSI_x1_x(i))/dt^2;
% % end
% % acc_interpolate_ESSI_x1_x(1)=0.0;
% % acc_interpolate_ESSI_x1_x(n)=0.0;
% % 
% % for i=2:N1-2
% %     acc_interpolate_SW4_x1_x(i)=(disp_SW4_x1_x(i+1)+disp_SW4_x1_x(i-1)-2*disp_SW4_x1_x(i))/dt^2;
% % end
% % acc_interpolate_SW4_x1_x(1)=0.0;
% % acc_interpolate_SW4_x1_x(N1-1)=0.0;
% % 
% % 
% % figure(41)
% % set(gcf,'Position',[Lox, Loy, width, height]);
% % p=plot (t(1:n),acc_ESSI_x1_x(1:n),'b-',t1(1:N1-1)-t0,acc_SW4_x1_x(1:N1-1),'r-',t(1:n),acc_interpolate_ESSI_x1_x(1:n),'k-',t1(1:N1-1)-t0,acc_interpolate_SW4_x1_x(1:N1-1),'g-')
% % p(1).LineWidth=a;
% % p(2).LineWidth=b;
% % p(3).LineWidth=b;
% % p(4).LineWidth=b;
% % set(gca, 'gridlinestyle','--','GridAlpha',1,'FontSize',14,'FontWeight', 'Normal','FontName','Arial');
% % % xlabel('Time (s)')
% % % ylabel('Acceleration (m/s^2)')
% % xlim([t10 t11]);
% % title('x1','FontName','Arial','FontSize',14,'FontWeight','Normal')
% % legend1=legend('ESSI','SW4','ESSI-interpolate','SW4-interpolate');
% % set(legend1,'FontSize',14,'FontName','Arial','FontWeight','Normal','Location','Northeast')
% % grid on
% % savefig(gcf,'41-comprison of acceleration response history in x direction');
% % saveas(gcf,'41-comprison of acceleration response history in x direction','meta');
% % saveas(gcf,'41-comprison of acceleration response history in x direction','png');
% % 
% % 
% % figure(5)
% % set(gcf,'Position',[Lox, Loy, width, height]);
% % subplot(3,2,1)
% % p=plot (t(1:n),acc_ESSI_x1_y(1:n),'b-',t1(1:N1-1)-t0,acc_SW4_x1_y(1:N1-1),'r-')
% % p(1).LineWidth=a;
% % p(2).LineWidth=b;
% % set(gca, 'gridlinestyle','--','GridAlpha',1,'FontSize',14,'FontWeight', 'Normal','FontName','Arial');
% % % xlabel('Time (s)')
% % % ylabel('Acceleration (m/s^2)')
% % xlim([t10 t11]);
% % title('x1','FontName','Arial','FontSize',14,'FontWeight','Normal')
% % legend1=legend('ESSI','SW4');
% % set(legend1,'FontSize',14,'FontName','Arial','FontWeight','Normal','Location','Northeast')
% % grid on
% %  
% % subplot(3,2,2)
% % p=plot (t(1:n),acc_ESSI_x2_y(1:n),'b-',t1(1:N1-1)-t0,acc_SW4_x2_y(1:N1-1),'r-')
% % p(1).LineWidth=a;
% % p(2).LineWidth=b;
% % set(gca, 'gridlinestyle','--','GridAlpha',1,'FontSize',14,'FontWeight', 'Normal','FontName','Arial');
% % % xlabel('Time (s)')
% % % ylabel('Acceleration (m/s^2)')
% % xlim([t10 t11]);
% % title('x2','FontName','Arial','FontSize',14,'FontWeight','Normal')
% % grid on
% %  
% %  
% % subplot(3,2,3)
% % p=plot (t(1:n),acc_ESSI_y1_y(1:n),'b-',t1(1:N1-1)-t0,acc_SW4_y1_y(1:N1-1),'r-')
% % p(1).LineWidth=a;
% % p(2).LineWidth=b;
% % set(gca, 'gridlinestyle','--','GridAlpha',1,'FontSize',14,'FontWeight', 'Normal','FontName','Arial');
% % % xlabel('Time (s)')
% % ylabel('Acceleration (m/s^2)')
% % xlim([t10 t11]);
% % title('y1','FontName','Arial','FontSize',14,'FontWeight','Normal')
% % grid on
% %  
% % subplot(3,2,4)
% % p=plot (t(1:n),acc_ESSI_y2_y(1:n),'b-',t1(1:N1-1)-t0,acc_SW4_y2_y(1:N1-1),'r-')
% % p(1).LineWidth=a;
% % p(2).LineWidth=b;
% % set(gca, 'gridlinestyle','--','GridAlpha',1,'FontSize',14,'FontWeight', 'Normal','FontName','Arial');
% % % xlabel('Time (s)')
% % ylabel('Acceleration (m/s^2)')
% % xlim([t10 t11]);
% % title('y2','FontName','Arial','FontSize',14,'FontWeight','Normal')
% % grid on
% %  
% % subplot(3,2,5)
% % p=plot (t(1:n),acc_ESSI_z1_y(1:n),'b-',t1(1:N1-1)-t0,acc_SW4_z1_y(1:N1-1),'r-')
% % p(1).LineWidth=a;
% % p(2).LineWidth=b;
% % set(gca, 'gridlinestyle','--','GridAlpha',1,'FontSize',14,'FontWeight', 'Normal','FontName','Arial');
% % xlabel('Time (s)')
% % % ylabel('Acceleration (m/s^2)')
% % xlim([t10 t11]);
% % title('z1','FontName','Arial','FontSize',14,'FontWeight','Normal')
% % grid on
% %  
% % subplot(3,2,6)
% % p=plot (t(1:n),acc_ESSI_z2_y(1:n),'b-',t1(1:N1-1)-t0,acc_SW4_z2_y(1:N1-1),'r-')
% % p(1).LineWidth=a;
% % p(2).LineWidth=b;
% % set(gca, 'gridlinestyle','--','GridAlpha',1,'FontSize',14,'FontWeight', 'Normal','FontName','Arial');
% % xlabel('Time (s)')
% % % ylabel('Acceleration (m/s^2)')
% % xlim([t10 t11]);
% % title('z2','FontName','Arial','FontSize',14,'FontWeight','Normal')
% % grid on
% % savefig(gcf,'5-comprison of acceleration response history in y direction');
% % saveas(gcf,'5-comprison of acceleration response history in y direction','meta');
% % saveas(gcf,'5-comprison of acceleration response history in y direction','png');
% %  
% %  
% % figure(6)
% % set(gcf,'Position',[Lox, Loy, width, height]);
% % subplot(3,2,1)
% % p=plot (t(1:n),acc_ESSI_x1_z(1:n),'b-',t1(1:N1-1)-t0,acc_SW4_x1_z(1:N1-1),'r-')
% % p(1).LineWidth=a;
% % p(2).LineWidth=b;
% % set(gca, 'gridlinestyle','--','GridAlpha',1,'FontSize',14,'FontWeight', 'Normal','FontName','Arial');
% % % xlabel('Time (s)')
% % % ylabel('Acceleration (m/s )')
% % xlim([t10 t11]);
% % title('x1','FontName','Arial','FontSize',14,'FontWeight','Normal')
% % legend1=legend('ESSI','SW4');
% % set(legend1,'FontSize',14,'FontName','Arial','FontWeight','Normal','Location','Northeast','Orientation','Horizontal')
% % grid on
% %  
% % subplot(3,2,2)
% % p=plot (t(1:n),acc_ESSI_x2_z(1:n),'b-',t1(1:N1-1)-t0,acc_SW4_x2_z(1:N1-1),'r-')
% % p(1).LineWidth=a;
% % p(2).LineWidth=b;
% % set(gca, 'gridlinestyle','--','GridAlpha',1,'FontSize',14,'FontWeight', 'Normal','FontName','Arial');
% % % xlabel('Time (s)')
% % % ylabel('Acceleration (m/s^2)')
% % xlim([t10 t11]);
% % title('x2','FontName','Arial','FontSize',14,'FontWeight','Normal')
% % grid on
% %  
% %  
% % subplot(3,2,3)
% % p=plot (t(1:n),acc_ESSI_y1_z(1:n),'b-',t1(1:N1-1)-t0,acc_SW4_y1_z(1:N1-1),'r-')
% % p(1).LineWidth=a;
% % p(2).LineWidth=b;
% % set(gca, 'gridlinestyle','--','GridAlpha',1,'FontSize',14,'FontWeight', 'Normal','FontName','Arial');
% % % xlabel('Time (s)')
% % ylabel('Acceleration (m/s^2)')
% % xlim([t10 t11]);
% % title('y1','FontName','Arial','FontSize',14,'FontWeight','Normal')
% % grid on
% %  
% % subplot(3,2,4)
% % p=plot (t(1:n),acc_ESSI_y2_z(1:n),'b-',t1(1:N1-1)-t0,acc_SW4_y2_z(1:N1-1),'r-')
% % p(1).LineWidth=a;
% % p(2).LineWidth=b;
% % set(gca, 'gridlinestyle','--','GridAlpha',1,'FontSize',14,'FontWeight', 'Normal','FontName','Arial');
% % % xlabel('Time (s)')
% % ylabel('Acceleration (m/s^2)')
% % xlim([t10 t11]);
% % title('y2','FontName','Arial','FontSize',14,'FontWeight','Normal')
% % grid on
% %  
% % subplot(3,2,5)
% % p=plot (t(1:n),acc_ESSI_z1_z(1:n),'b-',t1(1:N1-1)-t0,acc_SW4_z1_z(1:N1-1),'r-')
% % p(1).LineWidth=a;
% % p(2).LineWidth=b;
% % set(gca, 'gridlinestyle','--','GridAlpha',1,'FontSize',14,'FontWeight', 'Normal','FontName','Arial');
% % xlabel('Time (s)')
% % % ylabel('Acceleration (m/s^2)')
% % xlim([t10 t11]);
% % title('z1','FontName','Arial','FontSize',14,'FontWeight','Normal')
% % grid on
% %  
% % subplot(3,2,6)
% % p=plot (t(1:n),acc_ESSI_z2_z(1:n),'b-',t1(1:N1-1)-t0,acc_SW4_z2_z(1:N1-1),'r-')
% % p(1).LineWidth=a;
% % p(2).LineWidth=b;
% % set(gca, 'gridlinestyle','--','GridAlpha',1,'FontSize',14,'FontWeight', 'Normal','FontName','Arial');
% % xlabel('Time (s)')
% % % ylabel('Acceleration (m/s^2)')
% % xlim([t10 t11]);
% % title('z2','FontName','Arial','FontSize',14,'FontWeight','Normal')
% % grid on
% % savefig(gcf,'6-comprison of acceleration response history in z direction');
% % saveas(gcf,'6-comprison of acceleration response history in z direction','meta');
% % saveas(gcf,'6-comprison of acceleration response history in z direction','png');
% %  
% % [T_ESSI_x1_x,Spa_ESSI_x1_x,Spv_ESSI_x1_x,Sd_ESSI_x1_x]=SPEC(1*dt,0.01,acc_ESSI_x1_x(1:n),0.05,386.4,4);
% % [T_ESSI_x1_y,Spa_ESSI_x1_y,Spv_ESSI_x1_y,Sd_ESSI_x1_y]=SPEC(1*dt,0.01,acc_ESSI_x1_y(1:n),0.05,386.4,4);
% % [T_ESSI_x1_z,Spa_ESSI_x1_z,Spv_ESSI_x1_z,Sd_ESSI_x1_z]=SPEC(1*dt,0.01,acc_ESSI_x1_z(1:n),0.05,386.4,4);
% %  
% % [T_ESSI_x2_x,Spa_ESSI_x2_x,Spv_ESSI_x2_x,Sd_ESSI_x2_x]=SPEC(1*dt,0.01,acc_ESSI_x2_x(1:n),0.05,386.4,4);
% % [T_ESSI_x2_y,Spa_ESSI_x2_y,Spv_ESSI_x2_y,Sd_ESSI_x2_y]=SPEC(1*dt,0.01,acc_ESSI_x2_y(1:n),0.05,386.4,4);
% % [T_ESSI_x2_z,Spa_ESSI_x2_z,Spv_ESSI_x2_z,Sd_ESSI_x2_z]=SPEC(1*dt,0.01,acc_ESSI_x2_z(1:n),0.05,386.4,4);
% %  
% % [T_ESSI_y1_x,Spa_ESSI_y1_x,Spv_ESSI_y1_x,Sd_ESSI_y1_x]=SPEC(1*dt,0.01,acc_ESSI_y1_x(1:n),0.05,386.4,4);
% % [T_ESSI_y1_y,Spa_ESSI_y1_y,Spv_ESSI_y1_y,Sd_ESSI_y1_y]=SPEC(1*dt,0.01,acc_ESSI_y1_y(1:n),0.05,386.4,4);
% % [T_ESSI_y1_z,Spa_ESSI_y1_z,Spv_ESSI_y1_z,Sd_ESSI_y1_z]=SPEC(1*dt,0.01,acc_ESSI_y1_z(1:n),0.05,386.4,4);
% %  
% % [T_ESSI_y2_x,Spa_ESSI_y2_x,Spv_ESSI_y2_x,Sd_ESSI_y2_x]=SPEC(1*dt,0.01,acc_ESSI_y2_x(1:n),0.05,386.4,4);
% % [T_ESSI_y2_y,Spa_ESSI_y2_y,Spv_ESSI_y2_y,Sd_ESSI_y2_y]=SPEC(1*dt,0.01,acc_ESSI_y2_y(1:n),0.05,386.4,4);
% % [T_ESSI_y2_z,Spa_ESSI_y2_z,Spv_ESSI_y2_z,Sd_ESSI_y2_z]=SPEC(1*dt,0.01,acc_ESSI_y2_z(1:n),0.05,386.4,4);
% %  
% % [T_ESSI_z1_x,Spa_ESSI_z1_x,Spv_ESSI_z1_x,Sd_ESSI_z1_x]=SPEC(1*dt,0.01,acc_ESSI_z1_x(1:n),0.05,386.4,4);
% % [T_ESSI_z1_y,Spa_ESSI_z1_y,Spv_ESSI_z1_y,Sd_ESSI_z1_y]=SPEC(1*dt,0.01,acc_ESSI_z1_y(1:n),0.05,386.4,4);
% % [T_ESSI_z1_z,Spa_ESSI_z1_z,Spv_ESSI_z1_z,Sd_ESSI_z1_z]=SPEC(1*dt,0.01,acc_ESSI_z1_z(1:n),0.05,386.4,4);
% %  
% % [T_ESSI_z2_x,Spa_ESSI_z2_x,Spv_ESSI_z2_x,Sd_ESSI_z2_x]=SPEC(1*dt,0.01,acc_ESSI_z2_x(1:n),0.05,386.4,4);
% % [T_ESSI_z2_y,Spa_ESSI_z2_y,Spv_ESSI_z2_y,Sd_ESSI_z2_y]=SPEC(1*dt,0.01,acc_ESSI_z2_y(1:n),0.05,386.4,4);
% % [T_ESSI_z2_z,Spa_ESSI_z2_z,Spv_ESSI_z2_z,Sd_ESSI_z2_z]=SPEC(1*dt,0.01,acc_ESSI_z2_z(1:n),0.05,386.4,4);
% %  
% % [T_SW4_x1_x,Spa_SW4_x1_x,Spv_SW4_x1_x,Sd_SW4_x1_x]=SPEC(dt,0.01,acc_SW4_x1_x(1:n),0.05,386.4,4);
% % [T_SW4_x1_y,Spa_SW4_x1_y,Spv_SW4_x1_y,Sd_SW4_x1_y]=SPEC(dt,0.01,acc_SW4_x1_y(1:n),0.05,386.4,4);
% % [T_SW4_x1_z,Spa_SW4_x1_z,Spv_SW4_x1_z,Sd_SW4_x1_z]=SPEC(dt,0.01,acc_SW4_x1_z(1:n),0.05,386.4,4);
% %  
% % [T_SW4_x2_x,Spa_SW4_x2_x,Spv_SW4_x2_x,Sd_SW4_x2_x]=SPEC(dt,0.01,acc_SW4_x2_x(1:n),0.05,386.4,4);
% % [T_SW4_x2_y,Spa_SW4_x2_y,Spv_SW4_x2_y,Sd_SW4_x2_y]=SPEC(dt,0.01,acc_SW4_x2_y(1:n),0.05,386.4,4);
% % [T_SW4_x2_z,Spa_SW4_x2_z,Spv_SW4_x2_z,Sd_SW4_x2_z]=SPEC(dt,0.01,acc_SW4_x2_z(1:n),0.05,386.4,4);
% %  
% % [T_SW4_y1_x,Spa_SW4_y1_x,Spv_SW4_y1_x,Sd_SW4_y1_x]=SPEC(dt,0.01,acc_SW4_y1_x(1:n),0.05,386.4,4);
% % [T_SW4_y1_y,Spa_SW4_y1_y,Spv_SW4_y1_y,Sd_SW4_y1_y]=SPEC(dt,0.01,acc_SW4_y1_y(1:n),0.05,386.4,4);
% % [T_SW4_y1_z,Spa_SW4_y1_z,Spv_SW4_y1_z,Sd_SW4_y1_z]=SPEC(dt,0.01,acc_SW4_y1_z(1:n),0.05,386.4,4);
% %  
% % [T_SW4_y2_x,Spa_SW4_y2_x,Spv_SW4_y2_x,Sd_SW4_y2_x]=SPEC(dt,0.01,acc_SW4_y2_x(1:n),0.05,386.4,4);
% % [T_SW4_y2_y,Spa_SW4_y2_y,Spv_SW4_y2_y,Sd_SW4_y2_y]=SPEC(dt,0.01,acc_SW4_y2_y(1:n),0.05,386.4,4);
% % [T_SW4_y2_z,Spa_SW4_y2_z,Spv_SW4_y2_z,Sd_SW4_y2_z]=SPEC(dt,0.01,acc_SW4_y2_z(1:n),0.05,386.4,4);
% %  
% % [T_SW4_z1_x,Spa_SW4_z1_x,Spv_SW4_z1_x,Sd_SW4_z1_x]=SPEC(dt,0.01,acc_SW4_z1_x(1:n),0.05,386.4,4);
% % [T_SW4_z1_y,Spa_SW4_z1_y,Spv_SW4_z1_y,Sd_SW4_z1_y]=SPEC(dt,0.01,acc_SW4_z1_y(1:n),0.05,386.4,4);
% % [T_SW4_z1_z,Spa_SW4_z1_z,Spv_SW4_z1_z,Sd_SW4_z1_z]=SPEC(dt,0.01,acc_SW4_z1_z(1:n),0.05,386.4,4);
% %  
% % [T_SW4_z2_x,Spa_SW4_z2_x,Spv_SW4_z2_x,Sd_SW4_z2_x]=SPEC(dt,0.01,acc_SW4_z2_x(1:n),0.05,386.4,4);
% % [T_SW4_z2_y,Spa_SW4_z2_y,Spv_SW4_z2_y,Sd_SW4_z2_y]=SPEC(dt,0.01,acc_SW4_z2_y(1:n),0.05,386.4,4);
% % [T_SW4_z2_z,Spa_SW4_z2_z,Spv_SW4_z2_z,Sd_SW4_z2_z]=SPEC(dt,0.01,acc_SW4_z2_z(1:n),0.05,386.4,4);
% %  
% %  
% % % [T_ESSI_left_x1,Spa_ESSI_left_x1,Spv_ESSI_left_x1,Sd_ESSI_left_x1]=SPEC(dt,0.01,acc_ESSI_left_x1(1:n),0.05,386.4,4);
% % % [T_ESSI_left_y1,Spa_ESSI_left_y1,Spv_ESSI_left_y1,Sd_ESSI_left_y1]=SPEC(dt,0.01,acc_ESSI_left_y1(1:n),0.05,386.4,4);
% % % [T_ESSI_left_z1,Spa_ESSI_left_z1,Spv_ESSI_left_z1,Sd_ESSI_left_z1]=SPEC(dt,0.01,acc_ESSI_left_z1(1:n),0.05,386.4,4);
% % %  
% % % [T_ESSI_right_x1,Spa_ESSI_right_x1,Spv_ESSI_right_x1,Sd_ESSI_right_x1]=SPEC(dt,0.01,acc_ESSI_right_x1(1:n),0.05,386.4,4);
% % % [T_ESSI_right_y1,Spa_ESSI_right_y1,Spv_ESSI_right_y1,Sd_ESSI_right_y1]=SPEC(dt,0.01,acc_ESSI_right_y1(1:n),0.05,386.4,4);
% % % [T_ESSI_right_z1,Spa_ESSI_right_z1,Spv_ESSI_right_z1,Sd_ESSI_right_z1]=SPEC(dt,0.01,acc_ESSI_right_z1(1:n),0.05,386.4,4);
% % % 
% % % [T_ESSI_left_x2,Spa_ESSI_left_x2,Spv_ESSI_left_x2,Sd_ESSI_left_x2]=SPEC(dt,0.01,acc_ESSI_left_x2(1:n),0.05,386.4,4);
% % % [T_ESSI_left_y2,Spa_ESSI_left_y2,Spv_ESSI_left_y2,Sd_ESSI_left_y2]=SPEC(dt,0.01,acc_ESSI_left_y2(1:n),0.05,386.4,4);
% % % [T_ESSI_left_z2,Spa_ESSI_left_z2,Spv_ESSI_left_z2,Sd_ESSI_left_z2]=SPEC(dt,0.01,acc_ESSI_left_z2(1:n),0.05,386.4,4);
% % %  
% % % [T_ESSI_right_x2,Spa_ESSI_right_x2,Spv_ESSI_right_x2,Sd_ESSI_right_x2]=SPEC(dt,0.01,acc_ESSI_right_x2(1:n),0.05,386.4,4);
% % % [T_ESSI_right_y2,Spa_ESSI_right_y2,Spv_ESSI_right_y2,Sd_ESSI_right_y2]=SPEC(dt,0.01,acc_ESSI_right_y2(1:n),0.05,386.4,4);
% % % [T_ESSI_right_z2,Spa_ESSI_right_z2,Spv_ESSI_right_z2,Sd_ESSI_right_z2]=SPEC(dt,0.01,acc_ESSI_right_z2(1:n),0.05,386.4,4);
% % % 
% % % [T_ESSI_left_x3,Spa_ESSI_left_x3,Spv_ESSI_left_x3,Sd_ESSI_left_x3]=SPEC(dt,0.01,acc_ESSI_left_x3(1:n),0.05,386.4,4);
% % % [T_ESSI_left_y3,Spa_ESSI_left_y3,Spv_ESSI_left_y3,Sd_ESSI_left_y3]=SPEC(dt,0.01,acc_ESSI_left_y3(1:n),0.05,386.4,4);
% % % [T_ESSI_left_z3,Spa_ESSI_left_z3,Spv_ESSI_left_z3,Sd_ESSI_left_z3]=SPEC(dt,0.01,acc_ESSI_left_z3(1:n),0.05,386.4,4);
% % %  
% % % [T_ESSI_right_x3,Spa_ESSI_right_x3,Spv_ESSI_right_x3,Sd_ESSI_right_x3]=SPEC(dt,0.01,acc_ESSI_right_x3(1:n),0.05,386.4,4);
% % % [T_ESSI_right_y3,Spa_ESSI_right_y3,Spv_ESSI_right_y3,Sd_ESSI_right_y3]=SPEC(dt,0.01,acc_ESSI_right_y3(1:n),0.05,386.4,4);
% % % [T_ESSI_right_z3,Spa_ESSI_right_z3,Spv_ESSI_right_z3,Sd_ESSI_right_z3]=SPEC(dt,0.01,acc_ESSI_right_z3(1:n),0.05,386.4,4);
% % % 
% % % [T_ESSI_middle_x,Spa_ESSI_middle_x,Spv_ESSI_middle_x,Sd_ESSI_middle_x]=SPEC(dt,0.01,acc_ESSI_middle_x(1:n),0.05,386.4,4);
% % % [T_ESSI_middle_y,Spa_ESSI_middle_y,Spv_ESSI_middle_y,Sd_ESSI_middle_y]=SPEC(dt,0.01,acc_ESSI_middle_y(1:n),0.05,386.4,4);
% % % [T_ESSI_middle_z,Spa_ESSI_middle_z,Spv_ESSI_middle_z,Sd_ESSI_middle_z]=SPEC(dt,0.01,acc_ESSI_middle_z(1:n),0.05,386.4,4);
% % % 
% % % [T_SW4_left_x1,Spa_SW4_left_x1,Spv_SW4_left_x1,Sd_SW4_left_x1]=SPEC(dt,0.01,acc_SW4_left_x1(1:n),0.05,386.4,4);
% % % [T_SW4_left_y1,Spa_SW4_left_y1,Spv_SW4_left_y1,Sd_SW4_left_y1]=SPEC(dt,0.01,acc_SW4_left_y1(1:n),0.05,386.4,4);
% % % [T_SW4_left_z1,Spa_SW4_left_z1,Spv_SW4_left_z1,Sd_SW4_left_z1]=SPEC(dt,0.01,acc_SW4_left_z1(1:n),0.05,386.4,4);
% % %  
% % % [T_SW4_right_x1,Spa_SW4_right_x1,Spv_SW4_right_x1,Sd_SW4_right_x1]=SPEC(dt,0.01,acc_SW4_right_x1(1:n),0.05,386.4,4);
% % % [T_SW4_right_y1,Spa_SW4_right_y1,Spv_SW4_right_y1,Sd_SW4_right_y1]=SPEC(dt,0.01,acc_SW4_right_y1(1:n),0.05,386.4,4);
% % % [T_SW4_right_z1,Spa_SW4_right_z1,Spv_SW4_right_z1,Sd_SW4_right_z1]=SPEC(dt,0.01,acc_SW4_right_z1(1:n),0.05,386.4,4);
% % %  
% % % [T_SW4_left_x2,Spa_SW4_left_x2,Spv_SW4_left_x2,Sd_SW4_left_x2]=SPEC(dt,0.01,acc_SW4_left_x2(1:n),0.05,386.4,4);
% % % [T_SW4_left_y2,Spa_SW4_left_y2,Spv_SW4_left_y2,Sd_SW4_left_y2]=SPEC(dt,0.01,acc_SW4_left_y2(1:n),0.05,386.4,4);
% % % [T_SW4_left_z2,Spa_SW4_left_z2,Spv_SW4_left_z2,Sd_SW4_left_z2]=SPEC(dt,0.01,acc_SW4_left_z2(1:n),0.05,386.4,4);
% % %  
% % % [T_SW4_right_x2,Spa_SW4_right_x2,Spv_SW4_right_x2,Sd_SW4_right_x2]=SPEC(dt,0.01,acc_SW4_right_x2(1:n),0.05,386.4,4);
% % % [T_SW4_right_y2,Spa_SW4_right_y2,Spv_SW4_right_y2,Sd_SW4_right_y2]=SPEC(dt,0.01,acc_SW4_right_y2(1:n),0.05,386.4,4);
% % % [T_SW4_right_z2,Spa_SW4_right_z2,Spv_SW4_right_z2,Sd_SW4_right_z2]=SPEC(dt,0.01,acc_SW4_right_z2(1:n),0.05,386.4,4);
% % %  
% % % [T_SW4_left_x3,Spa_SW4_left_x3,Spv_SW4_left_x3,Sd_SW4_left_x3]=SPEC(dt,0.01,acc_SW4_left_x3(1:n),0.05,386.4,4);
% % % [T_SW4_left_y3,Spa_SW4_left_y3,Spv_SW4_left_y3,Sd_SW4_left_y3]=SPEC(dt,0.01,acc_SW4_left_y3(1:n),0.05,386.4,4);
% % % [T_SW4_left_z3,Spa_SW4_left_z3,Spv_SW4_left_z3,Sd_SW4_left_z3]=SPEC(dt,0.01,acc_SW4_left_z3(1:n),0.05,386.4,4);
% % %  
% % % [T_SW4_right_x3,Spa_SW4_right_x3,Spv_SW4_right_x3,Sd_SW4_right_x3]=SPEC(dt,0.01,acc_SW4_right_x3(1:n),0.05,386.4,4);
% % % [T_SW4_right_y3,Spa_SW4_right_y3,Spv_SW4_right_y3,Sd_SW4_right_y3]=SPEC(dt,0.01,acc_SW4_right_y3(1:n),0.05,386.4,4);
% % % [T_SW4_right_z3,Spa_SW4_right_z3,Spv_SW4_right_z3,Sd_SW4_right_z3]=SPEC(dt,0.01,acc_SW4_right_z3(1:n),0.05,386.4,4);
% % %  
% % % [T_SW4_middle_x,Spa_SW4_middle_x,Spv_SW4_middle_x,Sd_SW4_middle_x]=SPEC(dt,0.01,acc_SW4_middle_x(1:n),0.05,386.4,4);
% % % [T_SW4_middle_y,Spa_SW4_middle_y,Spv_SW4_middle_y,Sd_SW4_middle_y]=SPEC(dt,0.01,acc_SW4_middle_y(1:n),0.05,386.4,4);
% % % [T_SW4_middle_z,Spa_SW4_middle_z,Spv_SW4_middle_z,Sd_SW4_middle_z]=SPEC(dt,0.01,acc_SW4_middle_z(1:n),0.05,386.4,4);
% % %  
% %  
% % figure(7)
% % set(gcf,'Position',[Lox, Loy, width, height]);
% % subplot(3,2,1)
% % p=plot (T_ESSI_x1_x,Spa_ESSI_x1_x,'b-',T_SW4_x1_x,Spa_SW4_x1_x,'r-')
% % p(1).LineWidth=a;
% % p(2).LineWidth=b;
% % set(gca, 'gridlinestyle','--','GridAlpha',1,'FontSize',14,'FontWeight', 'Normal','FontName','Arial');
% % % xlabel('Time (s)')
% % % ylabel('Accleration response spectrum (g)')
% % title('x1','FontName','Arial','FontSize',14,'FontWeight','Normal')
% % legend1=legend('ESSI','SW4');
% % set(legend1,'FontSize',14,'FontName','Arial','FontWeight','Normal','Location','Northeast')
% % grid on
% %  
% % subplot(3,2,2)
% % p=plot (T_ESSI_x2_x,Spa_ESSI_x2_x,'b-',T_SW4_x2_x,Spa_SW4_x2_x,'r-')
% % p(1).LineWidth=a;
% % p(2).LineWidth=b;
% % set(gca, 'gridlinestyle','--','GridAlpha',1,'FontSize',14,'FontWeight', 'Normal','FontName','Arial');
% % % xlabel('Time (s)')
% % % ylabel('Accleration response spectrum (g)')
% % title('x2','FontName','Arial','FontSize',14,'FontWeight','Normal')
% % grid on
% %  
% %  
% % subplot(3,2,3)
% % p=plot (T_ESSI_y1_x,Spa_ESSI_y1_x,'b-',T_SW4_y1_x,Spa_SW4_y1_x,'r-')
% % p(1).LineWidth=a;
% % p(2).LineWidth=b;
% % set(gca, 'gridlinestyle','--','GridAlpha',1,'FontSize',14,'FontWeight', 'Normal','FontName','Arial');
% % % xlabel('Time (s)')
% % ylabel('Accleration response spectrum (g)')
% % title('y1','FontName','Arial','FontSize',14,'FontWeight','Normal')
% % grid on
% %  
% % subplot(3,2,4)
% % p=plot (T_ESSI_y2_x,Spa_ESSI_y2_x,'b-',T_SW4_y2_x,Spa_SW4_y2_x,'r-')
% % p(1).LineWidth=a;
% % p(2).LineWidth=b;
% % set(gca, 'gridlinestyle','--','GridAlpha',1,'FontSize',14,'FontWeight', 'Normal','FontName','Arial');
% % % Accleration response spectrum (g)')
% % ylabel('Accleration response spectrum (g)')
% % title('y2','FontName','Arial','FontSize',14,'FontWeight','Normal')
% % grid on
% %  
% % subplot(3,2,5)
% % p=plot (T_ESSI_z1_x,Spa_ESSI_z1_x,'b-',T_SW4_z1_x,Spa_SW4_z1_x,'r-')
% % p(1).LineWidth=a;
% % p(2).LineWidth=b;
% % set(gca, 'gridlinestyle','--','GridAlpha',1,'FontSize',14,'FontWeight', 'Normal','FontName','Arial');
% % xlabel('Period (s)')
% % % ylabel('Accleration response spectrum (g)')
% % title('z1','FontName','Arial','FontSize',14,'FontWeight','Normal')
% % grid on
% %  
% % subplot(3,2,6)
% % p=plot (T_ESSI_z2_x,Spa_ESSI_z2_x,'b-',T_SW4_z2_x,Spa_SW4_z2_x,'r-')
% % p(1).LineWidth=a;
% % p(2).LineWidth=b;
% % set(gca, 'gridlinestyle','--','GridAlpha',1,'FontSize',14,'FontWeight', 'Normal','FontName','Arial');
% % xlabel('Period (s)')
% % % Accleration response spectrum (g)')
% % title('z2','FontName','Arial','FontSize',14,'FontWeight','Normal')
% % grid on
% % savefig(gcf,'7-comprison of acceleration response spectrum in x direction');
% % saveas(gcf,'7-comprison of acceleration response spectrum in x direction','meta');
% % saveas(gcf,'7-comprison of acceleration response spectrum in x direction','png');
% %  
% %  
% % figure(8)
% % set(gcf,'Position',[Lox, Loy, width, height]);
% % subplot(3,2,1)
% % p=plot (T_ESSI_x1_y,Spa_ESSI_x1_y,'b-',T_SW4_x1_y,Spa_SW4_x1_y,'r-')
% % p(1).LineWidth=a;
% % p(2).LineWidth=b;
% % set(gca, 'gridlinestyle','--','GridAlpha',1,'FontSize',14,'FontWeight', 'Normal','FontName','Arial');
% % % xlabel('Time (s)')
% % % ylabel('Accleration response spectrum (g)')
% % title('x1','FontName','Arial','FontSize',14,'FontWeight','Normal')
% % legend1=legend('ESSI','SW4');
% % set(legend1,'FontSize',14,'FontName','Arial','FontWeight','Normal','Location','Northeast')
% % grid on
% %  
% % subplot(3,2,2)
% % p=plot (T_ESSI_x2_y,Spa_ESSI_x2_y,'b-',T_SW4_x2_y,Spa_SW4_x2_y,'r-')
% % p(1).LineWidth=a;
% % p(2).LineWidth=b;
% % set(gca, 'gridlinestyle','--','GridAlpha',1,'FontSize',14,'FontWeight', 'Normal','FontName','Arial');
% % % xlabel('Time (s)')
% % % ylabel('Accleration response spectrum (g)')
% % title('x2','FontName','Arial','FontSize',14,'FontWeight','Normal')
% % grid on
% %  
% %  
% % subplot(3,2,3)
% % p=plot (T_ESSI_y1_y,Spa_ESSI_y1_y,'b-',T_SW4_y1_y,Spa_SW4_y1_y,'r-')
% % p(1).LineWidth=a;
% % p(2).LineWidth=b;
% % set(gca, 'gridlinestyle','--','GridAlpha',1,'FontSize',14,'FontWeight', 'Normal','FontName','Arial');
% % % xlabel('Time (s)')
% % ylabel('Accleration response spectrum (g)')
% % title('y1','FontName','Arial','FontSize',14,'FontWeight','Normal')
% % grid on
% %  
% % subplot(3,2,4)
% % p=plot (T_ESSI_y2_y,Spa_ESSI_y2_y,'b-',T_SW4_y2_y,Spa_SW4_y2_y,'r-')
% % p(1).LineWidth=a;
% % p(2).LineWidth=b;
% % set(gca, 'gridlinestyle','--','GridAlpha',1,'FontSize',14,'FontWeight', 'Normal','FontName','Arial');
% % % Accleration response spectrum (g)')
% % ylabel('Accleration response spectrum (g)')
% % title('y2','FontName','Arial','FontSize',14,'FontWeight','Normal')
% % grid on
% %  
% % subplot(3,2,5)
% % p=plot (T_ESSI_z1_y,Spa_ESSI_z1_y,'b-',T_SW4_z1_y,Spa_SW4_z1_y,'r-')
% % p(1).LineWidth=a;
% % p(2).LineWidth=b;
% % set(gca, 'gridlinestyle','--','GridAlpha',1,'FontSize',14,'FontWeight', 'Normal','FontName','Arial');
% % xlabel('Period (s)')
% % % ylabel('Accleration response spectrum (g)')
% % title('z1','FontName','Arial','FontSize',14,'FontWeight','Normal')
% % grid on
% %  
% % subplot(3,2,6)
% % p=plot (T_ESSI_z2_y,Spa_ESSI_z2_y,'b-',T_SW4_z2_y,Spa_SW4_z2_y,'r-')
% % p(1).LineWidth=a;
% % p(2).LineWidth=b;
% % set(gca, 'gridlinestyle','--','GridAlpha',1,'FontSize',14,'FontWeight', 'Normal','FontName','Arial');
% % xlabel('Period (s)')
% % % Accleration response spectrum (g)')
% % title('z2','FontName','Arial','FontSize',14,'FontWeight','Normal')
% % grid on
% % savefig(gcf,'8-comprison of acceleration response spectrum in y direction');
% % saveas(gcf,'8-comprison of acceleration response spectrum in y direction','meta');
% % saveas(gcf,'8-comprison of acceleration response spectrum in y direction','png');
% %  
% %  
% % figure(9)
% % set(gcf,'Position',[Lox, Loy, width, height]);
% % subplot(3,2,1)
% % p=plot (T_ESSI_x1_z,Spa_ESSI_x1_z,'b-',T_SW4_x1_z,Spa_SW4_x1_z,'r-')
% % p(1).LineWidth=a;
% % p(2).LineWidth=b;
% % set(gca, 'gridlinestyle','--','GridAlpha',1,'FontSize',14,'FontWeight', 'Normal','FontName','Arial');
% % % xlabel('Time (s)')
% % xlim([0 2]);
% % % ylabel('Accleration response spectrum (g)')
% % title('x1','FontName','Arial','FontSize',14,'FontWeight','Normal')
% % legend1=legend('ESSI','SW4');
% % set(legend1,'FontSize',14,'FontName','Arial','FontWeight','Normal','Location','Northeast')
% % grid on
% %  
% % subplot(3,2,2)
% % p=plot (T_ESSI_x2_z,Spa_ESSI_x2_z,'b-',T_SW4_x2_z,Spa_SW4_x2_z,'r-')
% % p(1).LineWidth=a;
% % p(2).LineWidth=b;
% % set(gca, 'gridlinestyle','--','GridAlpha',1,'FontSize',14,'FontWeight', 'Normal','FontName','Arial');
% % % xlabel('Time (s)')
% % xlim([0 2]);
% % % ylabel('Accleration response spectrum (g)')
% % title('x2','FontName','Arial','FontSize',14,'FontWeight','Normal')
% % grid on
% %  
% %  
% % subplot(3,2,3)
% % p=plot (T_ESSI_y1_z,Spa_ESSI_y1_z,'b-',T_SW4_y1_z,Spa_SW4_y1_z,'r-')
% % p(1).LineWidth=a;
% % p(2).LineWidth=b;
% % set(gca, 'gridlinestyle','--','GridAlpha',1,'FontSize',14,'FontWeight', 'Normal','FontName','Arial');
% % % xlabel('Time (s)')
% % xlim([0 2]);
% % ylabel('Accleration response spectrum (g)')
% % title('y1','FontName','Arial','FontSize',14,'FontWeight','Normal')
% % grid on
% %  
% % subplot(3,2,4)
% % p=plot (T_ESSI_y2_z,Spa_ESSI_y2_z,'b-',T_SW4_y2_z,Spa_SW4_y2_z,'r-')
% % p(1).LineWidth=a;
% % p(2).LineWidth=b;
% % set(gca, 'gridlinestyle','--','GridAlpha',1,'FontSize',14,'FontWeight', 'Normal','FontName','Arial');
% % % Accleration response spectrum (g)')
% % xlim([0 2]);
% % ylabel('Accleration response spectrum (g)')
% % title('y2','FontName','Arial','FontSize',14,'FontWeight','Normal')
% % grid on
% %  
% % subplot(3,2,5)
% % p=plot (T_ESSI_z1_z,Spa_ESSI_z1_z,'b-',T_SW4_z1_z,Spa_SW4_z1_z,'r-')
% % p(1).LineWidth=a;
% % p(2).LineWidth=b;
% % set(gca, 'gridlinestyle','--','GridAlpha',1,'FontSize',14,'FontWeight', 'Normal','FontName','Arial');
% % xlabel('Period (s)')
% % xlim([0 2]);
% % % ylabel('Accleration response spectrum (g)')
% % title('z1','FontName','Arial','FontSize',14,'FontWeight','Normal')
% % grid on
% %  
% % subplot(3,2,6)
% % p=plot (T_ESSI_z2_z,Spa_ESSI_z2_z,'b-',T_SW4_z2_z,Spa_SW4_z2_z,'r-')
% % p(1).LineWidth=a;
% % p(2).LineWidth=b;
% % set(gca, 'gridlinestyle','--','GridAlpha',1,'FontSize',14,'FontWeight', 'Normal','FontName','Arial');
% % xlabel('Period(s)')
% % xlim([0 2]);
% % % Accleration response spectrum (g)')
% % title('z2','FontName','Arial','FontSize',14,'FontWeight','Normal')
% % grid on
% % savefig(gcf,'9-comprison of acceleration response spectrum in z direction');
% % saveas(gcf,'9-comprison of acceleration response spectrum in z direction','meta');
% % saveas(gcf,'9-comprison of acceleration response spectrum in z direction','png');
% %  
% %  
% %  


