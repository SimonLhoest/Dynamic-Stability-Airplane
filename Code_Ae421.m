%{
------------------------------------------------------
Flight Conditions 
------------------------------------------------------
%}
 
H = 0;
Rho = 1.226;
U0 = 70.1;
M = 0.2;
GC = 0;
Theta0 = 0.2041;
 
%{
------------------------------------------------------
Geometric data
------------------------------------------------------
%}
 
S = 49.24;
b = 11.8;
Cbar = 4.88;
A = 2.83;
e = 0.89;
 
%{
------------------------------------------------------
Inertial data 
------------------------------------------------------
%}
 
m = 15060;
Ixx = 32146;
Iyy = 159375;
Izz = 181349;
Ixz = 2170;
 
%{
------------------------------------------------------
Steady state conditions 
------------------------------------------------------
%}
 
CL0 = 1;
CD0 = 0.2;
CTx0 = 0.2;
Cm0 = 0;
Cmt0 = 0;
 
%{
------------------------------------------------------
Longitudinal Aerodynamic deravatives  
------------------------------------------------------
%}
 
Cmu = 0;
Cm_alpha = -0.098;
Cm_alpha_dot = -0.95;
Cmq = -2;
Cm_t_u = 0;
Cm_t_alpha = 0;
C_L_u = 0;
C_L_alpha = 2.8;
C_L_alpha_dot = 0;
CLq = 0;
C_D_u = 0;
C_D_alpha = 0.555;
C_D_alpha_dot = 0;
C_Dq = 0;
C_T_u = 0;
C_L_delta_e = 0.24;
C_D_delta_e = -0.14;
C_m_delta_e = -0.322;
 
%{
------------------------------------------------------
Lateral Aerodynamic derivatives 
------------------------------------------------------
%}
 
C_I_beta = -0.156;
C_Ip = -0.272;
C_Ir = 0.205;
C_I_delta_a = 0.057;
C_I_delta_r = 0.0009;
C_n_beta = 0.199;
C_np = 0.013;
C_nr = -0.32;
C_n_delta_a = -0.0041;
C_n_delta_r = -0.072;
C_y_beta = -0.655;
C_yp = 0;
C_yr = 0;
C_y_delta_a = -0.0355;
C_y_delta_r = 0.124;
C_T_delta_T = 0;
C_L_delta_T = 0;
C_m_delta_T = 0;
 
qbar = (1/2)*Rho*U0^2;
 
g = 9.81;
 
%{
------------------------------------------------------
Matrix A 
------------------------------------------------------
%}

 
Xu = (-qbar * S)*(2*CD0 + C_D_u)/(m*U0);
Xw = ((qbar * S)/(m * U0))*(CL0*(1-(2*C_L_alpha)/(pi*e*A)));
 
Zu = -((qbar*S)/(m*U0))*(2*CL0 + C_L_u);
Zw = (-qbar * S)*(CD0 + C_L_alpha)/(m*U0);
 
Mq = (qbar*S*(Cbar^2) * Cmq) / (2 * Iyy * U0);
Mu = ((qbar*S*Cbar)/(Iyy*U0))* Cmu;
Mw = ((qbar*S*Cbar)/(Iyy*U0))* Cm_alpha;
%{
------------------------------------------------------
Matrix B 
------------------------------------------------------
%}
 
 
X_delta_e = (qbar * S * C_D_delta_e) / (m * U0);
X_delta_T = (qbar * S * C_T_delta_T) / (m * U0);
 
Z_delta_e = (qbar * S * C_L_delta_e) / (m * U0);
Z_delta_T = (qbar * S * C_L_delta_T) / (m * U0);
 
M_delta_e = (qbar * S * Cbar * C_m_delta_e) / (Iyy * U0);
M_delta_T = (qbar * S * Cbar * C_m_delta_T) / (Iyy * U0);
 
Mw_point = (qbar * S * Cbar^2 * Cm_alpha_dot) / (2 * Iyy * U0^2);
 
 
 
A = [Xu Xw 0 -g*cos(Theta0); Zu Zw U0 -g*sin(Theta0); Mu+Zu*Mw_point Mw+Zw*Mw_point Mq + U0*Mw_point 0; 0 0 1 0];
B = [X_delta_e X_delta_T;Z_delta_e Z_delta_T;M_delta_e+Z_delta_e*Mw_point M_delta_T+Z_delta_T*Mw_point;0 0];



%{
------------------------------------------------------
Characteristic Equation
------------------------------------------------------
%}

eq_characteristic = poly(A);



%{
------------------------------------------------------
Eigen_values
------------------------------------------------------
%}
[Eigen_vector, Eigen_values] = eig(A);

lambda1 = Eigen_values(1,1);
lambda2 = Eigen_values(2,2);
lambda3 = Eigen_values(3,3);
lambda4 = Eigen_values(4,4);


%{
------------------------------------------------------
Short period Mode
------------------------------------------------------
%}
wnsp = sqrt(lambda1*lambda2);
epsilon_sp = (-lambda1-lambda2)/(wnsp * 2);


%{
------------------------------------------------------
Phugoid Mode
------------------------------------------------------
%}
wnp = sqrt(lambda3*lambda4);
epsilon_p = (-lambda3-lambda4)/(wnp * 2);


%{
------------------------------------------------------
Curves of longitudinal motion
------------------------------------------------------
%}
t_long = 0:1:700;
delta_u = (Eigen_vector(1,1)*exp(lambda1*t_long))+(Eigen_vector(1,2)*exp(lambda2*t_long))+(Eigen_vector(1,3)*exp(lambda3*t_long))+(Eigen_vector(1,4)*exp(lambda4*t_long));
delta_w = (Eigen_vector(2,1)*exp(lambda1*t_long))+(Eigen_vector(2,2)*exp(lambda2*t_long))+(Eigen_vector(2,3)*exp(lambda3*t_long))+(Eigen_vector(2,4)*exp(lambda4*t_long));
delta_q = (Eigen_vector(3,1)*exp(lambda1*t_long))+(Eigen_vector(3,2)*exp(lambda2*t_long))+(Eigen_vector(3,3)*exp(lambda3*t_long))+(Eigen_vector(3,4)*exp(lambda4*t_long));
delta_teta = (Eigen_vector(4,1)*exp(lambda1*t_long))+(Eigen_vector(4,2)*exp(lambda2*t_long))+(Eigen_vector(4,3)*exp(lambda3*t_long))+(Eigen_vector(4,4)*exp(lambda4*t_long));
%%

subplot(2,2,1)
plot(t_long, delta_u);
subtitle("Phugoid mode delta u")
subplot(2,2,2)
plot(t_long, delta_w);
subtitle("short period mode delta w");
subplot(2,2,3)
plot(t_long, delta_q);
subtitle("short period mode delta q");
subplot(2,2,4)
plot(t_long, delta_teta);
subtitle("Phugoid mode delta teta");
%%


%{
------------------------------------------------------
Transfer Function of Each variable
------------------------------------------------------
%}
I_long = [1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1];
s1=tf('s'); 

tf_long = minreal(inv(s1*I_long-A)*B,1e-4); 
vda1_long = zpk(tf_long(1,1));
vda2_long = zpk(tf_long(2,1));
vda3_long = zpk(tf_long(3,1));
vda4_long = zpk(tf_long(4,1));


%{
------------------------------------------------------
II) Lateral Dynamic Stability Airplane 
------------------------------------------------------
%}


%{
------------------------------------------------------
Matrix A lateral
------------------------------------------------------
%}
 
Yv = (qbar*S*C_y_beta)/(m*U0);
Lv = (qbar*S*b*C_I_beta)/(Ixx*U0);
Nv = (qbar*S*b*C_n_beta)/(Izz*U0);
 
Yp = (qbar*S*b*C_yp)/(2*m*U0);
Lp = (qbar*S*b^2*C_Ip)/(2*Ixx*U0);
Np = (qbar*S*b^2*C_np)/(2*Izz*U0);
 
Yr = (qbar*S*b*C_yr)/(2*m*U0);
Lr = (qbar*S*b^2*C_Ir)/(2*Ixx*U0);
Nr = (qbar*S*b^2*C_nr)/(2*Izz*U0);

A_Lat = [Yv Yp -(U0-Yr) g*cos(Theta0); Lv Lp Lr 0; Nv Np Nr 0; 0 1 0 0];

%{
------------------------------------------------------
Matrix B lateral 
------------------------------------------------------
%}
 
Y_delta_r = ((qbar*S)/(m*U0))*C_y_delta_r;
L_delta_r = ((qbar*S*b)/(Ixx*U0))*C_I_delta_r;
N_delta_r = ((qbar*S*b)/(Izz*U0))*C_n_delta_r;
 
Y_delta_a = ((qbar*S)/(m*U0))*C_y_delta_a;
L_delta_a = ((qbar*S*b)/(Ixx*U0))*C_I_delta_a;
N_delta_a = ((qbar*S*b)/(Izz*U0))*C_n_delta_a;
 

B_Lat = [Y_delta_r Y_delta_a; L_delta_r L_delta_a; N_delta_r N_delta_a; 0 0];


%{
------------------------------------------------------
Characteristic Equation
------------------------------------------------------
%}

eq_characteristic_Lat = poly(A_Lat);


%{
------------------------------------------------------
Eigen_values
------------------------------------------------------
%}
[Eigen_vector_Lat, Eigen_values_Lat] = eig(A_Lat);

lambda1_Lat = Eigen_values_Lat(1,1);
lambda2_Lat = Eigen_values_Lat(2,2);
lambda3_Lat = Eigen_values_Lat(3,3);
lambda4_Lat = Eigen_values_Lat(4,4);


%{
------------------------------------------------------
Spiral Mode
------------------------------------------------------
%}
t_Lat = 0:1:100;
delta_v_spiral = (Eigen_vector_Lat(1,2)*exp(lambda2_Lat*t_Lat));
delta_p_spiral = (Eigen_vector_Lat(2,2)*exp(lambda2_Lat*t_Lat));
delta_r_spiral = (Eigen_vector_Lat(3,2)*exp(lambda2_Lat*t_Lat));
delta_phi_spiral = (Eigen_vector_Lat(4,2)*exp(lambda2_Lat*t_Lat));

plot(t_Lat, delta_v_spiral);

plot(t_Lat, delta_p_spiral);

plot(t_Lat, delta_r_spiral);

plot(t_Lat, delta_phi_spiral);



%{
------------------------------------------------------
Rolling Mode
------------------------------------------------------
%}
delta_v_rolling = (Eigen_vector_Lat(1,1)*exp(lambda1_Lat*t_Lat));
delta_p_rolling = (Eigen_vector_Lat(2,1)*exp(lambda1_Lat*t_Lat));
delta_r_rolling = (Eigen_vector_Lat(3,1)*exp(lambda1_Lat*t_Lat));
delta_phi_rolling = (Eigen_vector_Lat(4,1)*exp(lambda1_Lat*t_Lat));


plot(t_Lat, delta_v_rolling);

plot(t_Lat, delta_p_rolling);

plot(t_Lat, delta_r_rolling);

plot(t_Lat, delta_phi_rolling);




%{
------------------------------------------------------
Dutch Mode
------------------------------------------------------
%}
delta_v_Dutch = (Eigen_vector_Lat(1,3)*exp(lambda3_Lat*t_Lat))+(Eigen_vector_Lat(1,4)*exp(lambda4_Lat*t_Lat));
delta_p_Dutch = (Eigen_vector_Lat(2,3)*exp(lambda3_Lat*t_Lat))+(Eigen_vector_Lat(2,4)*exp(lambda4_Lat*t_Lat));
delta_r_Dutch = (Eigen_vector_Lat(3,3)*exp(lambda3_Lat*t_Lat))+(Eigen_vector_Lat(3,4)*exp(lambda4_Lat*t_Lat));
delta_phi_Dutch = (Eigen_vector_Lat(4,3)*exp(lambda3_Lat*t_Lat))+(Eigen_vector_Lat(4,4)*exp(lambda4_Lat*t_Lat));

subplot(1,3,3)
plot(t_Lat, delta_v_Dutch);
subplot(1,3,3)
plot(t_Lat, delta_p_Dutch);
subplot(1,3,3)
plot(t_Lat, delta_r_Dutch);
subplot(1,3,3)
plot(t_Lat, delta_phi_Dutch);


%{
------------------------------------------------------
Transfer Function of Each variable
------------------------------------------------------
%}

%TF
I_lat = [1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1];
s=tf('s'); 

tf_lat = minreal(inv(s*I_lat-A_Lat)*B_Lat,1e-4); 
vda1_lat = zpk(tf_lat(1,1));
vda2_lat = zpk(tf_lat(2,1));
vda3_lat = zpk(tf_lat(3,1));
vda4_lat = zpk(tf_lat(4,1));

