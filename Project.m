%% Flight conditions
H=0; g=9.81;
rho=1.226;
Uo=70.1;
M=0.2;
theta0=0.2041;
% Geometric data
S=49.24;
b=11.8;
c_=4.88;
A=2.83;
e=0.89;
% Inertial data
m=15060;
Ixx=32146;
Iyy=159375;
Izz=181349;
Ixz=2170;
% Steady state conditions:
Clo=1;
Cdo=0.2;
Ctxo=0.2;
Cmo=0;
Cmto=0;
% Longitudinal Aerodynamic derivatives
Cmu=0;
Cmalpha=-0.098;
Cmalphapoint=-0.95;
Cmq=-2;
Cmtu=0;
Cmtalpha=0;
Clu=0;
Clalpha=2.8;
Clalphapoint=0;
Clq=0;
Cdu=0;
Cdalpha=0.555;
Cdalphapoint=0;
Cdq=0;
Ctu=0;
Cldeltae=0.24;
Cddeltae=-0.14;
Cmdeltae=-0.322;
% Lateral aerodynamic derivatives
Clbeta=-0.156;
Clp=-0.272;
Clr=0.205;
Cldeltaa=0.057;
Cldeltar=0.0009;
Cnbeta=0.199;
Cnp=0.013;
Cnr=-0.32;
Cndeltaa=-0.0041;
Cndeltar=-0.072;
Cybeta=-0.655;
Cyp=0;
Cyr=0;
Cydeltaa=-0.0355;
Cydeltar=0.124;
Cddeltat=0;Cldeltat=0;Cmdeltat=0;
%% Longitudinal Dynamic Stability of Airplane
% Matrix A and B
q_=1/2*rho*Uo^2;
Xu=-q_*S/(m*Uo)*(2*Cdo+Cdu);
Zu=-q_*S/(m*Uo)*(2*Clo+Clu);
Mu=q_*S*c_/(Iyy*Uo)*Cmu;
Xw=q_*S/(m*Uo)*(Clo*(1-2/(pi*e*A)*Clalpha));
Zw=-q_*S/(m*Uo)*(Cdo+Clalpha);
Mw=q_*S*c_/(Iyy*Uo)*Cmalpha;
Xwpoint=0;
Zwpoint=q_*S*c_/(2*m*Uo^2)*(Cdo+Clalphapoint); %???????? 3.1064*10-5
Mwpoint=q_*S*c_^2/(2*Iyy*Uo^2)*Cmalphapoint;
Xq=0;
Zq=q_*S*c_/(2*m*Uo)*Clq;
Mq=q_*S*c_^2/(2*Iyy*Uo)*Cmq;
%Matrix A
Along=[Xu Xw 0 -g*cos(theta0);Zu Zw Uo -g*sin(theta0);Mu+Zu*Mwpoint Mw+Zw*Mwpoint Mq+Uo*Mwpoint 0;0 0 1 0];
%Elevator derivatives
Xdeltae=q_*S*Cddeltae/(m*Uo);
Zdeltae=q_*S*Cldeltae/(m*Uo);
Mdeltae=q_*S*c_*Cmdeltae/(Iyy*Uo);
Xdeltat=q_*S*Cddeltat/(m*Uo);
Zdeltat=q_*S*Cldeltat/(m*Uo);
Mdeltat=q_*S*c_*Cmdeltat/(Iyy*Uo);
% Matrix B
Blong=[Xdeltae Xdeltat;Zdeltae Zdeltat;Mdeltae+Zdeltae*Mwpoint Mdeltat+Zdeltat*Mwpoint;0 0];
%% Characteristic equation
eqlong=poly(Along);
%% Eigenvalues
[eigvectorslong, eigvalueslong]=eig(Along);
eigvalueslong=[eigvalueslong(1,1);eigvalueslong(2,2);eigvalueslong(3,3);eigvalueslong(4,4)];
%% Different modes of longitudinal stability
% Short period mode
syms lambda;
Sp(lambda)=(lambda-eigvalueslong(1))*(lambda-eigvalueslong(2));
Wsp=round((Sp(0))^(1/2),5);
Epssp(lambda)=diff(Sp(lambda));
Epssp=round(Epssp(0)/(2*Wsp),5);
% Phugoid mode
Ph(lambda)=(lambda-eigvalueslong(3))*(lambda-eigvalueslong(4));
Wph=round((Ph(0))^(1/2),5);
Epsph(lambda)=diff(Ph(lambda));
Epsph=round(Epsph(0)/(2*Wph),5);
%% Curves of longitudinal motions
tlong=0:1:600;
deltau=(eigvectorslong(1,1)*exp(eigvalueslong(1)*tlong))+(eigvectorslong(1,2)*exp(eigvalueslong(2)*tlong))+(eigvectorslong(1,3)*exp(eigvalueslong(3)*tlong))+(eigvectorslong(1,4)*exp(eigvalueslong(4)*tlong));
deltaw=(eigvectorslong(2,1)*exp(eigvalueslong(1)*tlong))+(eigvectorslong(2,2)*exp(eigvalueslong(2)*tlong))+(eigvectorslong(2,3)*exp(eigvalueslong(3)*tlong))+(eigvectorslong(2,4)*exp(eigvalueslong(4)*tlong));
deltaq=(eigvectorslong(3,1)*exp(eigvalueslong(1)*tlong))+(eigvectorslong(3,2)*exp(eigvalueslong(2)*tlong))+(eigvectorslong(3,3)*exp(eigvalueslong(3)*tlong))+(eigvectorslong(3,4)*exp(eigvalueslong(4)*tlong));
deltatheta=(eigvectorslong(4,1)*exp(eigvalueslong(1)*tlong))+(eigvectorslong(4,2)*exp(eigvalueslong(2)*tlong))+(eigvectorslong(4,3)*exp(eigvalueslong(3)*tlong))+(eigvectorslong(4,4)*exp(eigvalueslong(4)*tlong));
figure(4)
subplot(2,2,1)
plot(tlong, deltau);
subtitle("Phugoid mode delta u")
subplot(2,2,2)
plot(tlong, deltatheta);
subtitle("Phugoid mode delta teta");
subplot(2,2,3)
plot(tlong, deltaw);
subtitle("short period mode delta w");
subplot(2,2,4)
plot(tlong, deltaq);
subtitle("short period mode delta q");

%% Transfer functions of each variable
s=tf('s'); 
I=[1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1];
tflong=minreal((s*I-Along)\Blong(:,1),1e-4); 
u=zpk(tflong(1,1));
w=zpk(tflong(2,1));
q=zpk(tflong(3,1));
theta=zpk(tflong(4,1));

%% Lateral Dynamic Stability of Airplane
% Matrix A and B
%Lateral Stability derivatives
Yv=q_*S/(m*Uo)*Cybeta;
Yp=q_*S*b/(2*m*Uo)*Cyp;
Yr=q_*S*b/(2*m*Uo)*Cyr;
Lv=q_*S*b/(Ixx*Uo)*Clbeta;
Lp=q_*S*b^2/(2*Ixx*Uo)*Clp;
Lr=q_*S*b^2/(2*Ixx*Uo)*Clr;
Nv= q_*S*b/(Izz*Uo)*Cnbeta;
Np= q_*S*b^2/(2*Izz*Uo)*Cnp;
Nr= q_*S*b^2/(2*Izz*Uo)*Cnr;
% Matrix A
Alat=[Yv Yp -(Uo-Yr) g*cos(theta0);Lv Lp Lr 0;Nv Np Nr 0;0 1 0 0];
% Control Matrix Elements
Ydeltar=q_*S/(m*Uo)*Cydeltar;
Ydeltaa=q_*S/(m*Uo)*Cydeltaa;
Ldeltar=q_*S*b/(Ixx*Uo)*Cldeltar;
Ldeltaa=q_*S*b/(Ixx*Uo)*Cldeltaa;
Ndeltar=q_*S*b/(Izz*Uo)*Cndeltar;
Ndeltaa=q_*S*b/(Izz*Uo)*Cndeltaa;
% Matrix B
Blat=[Ydeltar Ydeltaa;Ldeltar Ldeltaa;Ndeltar Ndeltaa;0 0];
%% Characteristic equation
eqlat=poly(Alat);
%% Eigenvalues
[eigvectorslat,eigvalueslat]=eig(Alat);
eigvalueslat=[eigvalueslat(1,1);eigvalueslat(2,2);eigvalueslat(3,3);eigvalueslat(4,4)];
%% Modes of lateral stability
% Rolling mode
Lroll=eigvalueslat(3);
% Spiral mode
Lspiral=eigvalueslat(4);
% Dutch roll mode
Wdr=(eigvalueslat(1)*eigvalueslat(2))^(1/2);
Epsdr=(eigvalueslat(1)+eigvalueslat(2))/(-2*Wdr);

%% Curves of longitudinal motions
tlat=0:1:300;
deltavspiral=(eigvectorslat(1,2)*exp(eigvalueslat(4)*tlat));
deltapspiral=(eigvectorslat(2,2)*exp(eigvalueslat(4)*tlat));
deltarspiral=(eigvectorslat(3,2)*exp(eigvalueslat(4)*tlat));
deltaphispiral=(eigvectorslat(4,2)*exp(eigvalueslat(4)*tlat));
figure(1)
hold on
plot(tlat, deltavspiral);
plot(tlat, deltapspiral);
plot(tlat, deltarspiral);
plot(tlat, deltaphispiral);
title('Spiral mode')
xlabel('Time (s)')
legend('Roll angle','Yaw rate','Side slip','Roll rate');

tlat=0:0.1:6;
deltavroll=(eigvectorslat(1,1)*exp(eigvalueslat(3)*tlat));
deltaproll=(eigvectorslat(2,1)*exp(eigvalueslat(3)*tlat));
deltarroll=(eigvectorslat(3,1)*exp(eigvalueslat(3)*tlat));
deltaphiroll=(eigvectorslat(4,1)*exp(eigvalueslat(3)*tlat));
figure(2)
hold on
plot(tlat, deltavroll);
plot(tlat, deltaproll);
plot(tlat, deltarroll);
plot(tlat, deltaphiroll);
title('Rolling mode')
xlabel('Time (s)')
legend('Side velocity','Roll rate','Yaw rate','Slide angle');

tlat=0:0.1:100;
deltavdutch=(eigvectorslat(1,3)*exp(eigvalueslat(1)*tlat))+(eigvectorslat(1,4)*exp(eigvalueslat(2)*tlat));
deltapdutch=(eigvectorslat(2,3)*exp(eigvalueslat(1)*tlat))+(eigvectorslat(2,4)*exp(eigvalueslat(2)*tlat));
deltardutch=(eigvectorslat(3,3)*exp(eigvalueslat(1)*tlat))+(eigvectorslat(3,4)*exp(eigvalueslat(2)*tlat));
deltaphidutch=(eigvectorslat(4,3)*exp(eigvalueslat(1)*tlat))+(eigvectorslat(4,4)*exp(eigvalueslat(2)*tlat));
figure(3)
hold on
plot(tlat, deltavdutch);
plot(tlat, deltapdutch);
plot(tlat, deltardutch);
plot(tlat, deltaphidutch);
title('Dutch roll mode')
xlabel('Time (s)')
legend('Side velocity','Roll rate','Yaw rate','Slide angle');

%% Transfer functions of each variable
tflat=minreal((s*I-Alat)\Blat,1e-4); 
v=zpk(tflat(1,1));
p=zpk(tflat(2,1));
r=zpk(tflat(3,1));
phi=zpk(tflat(4,1));




