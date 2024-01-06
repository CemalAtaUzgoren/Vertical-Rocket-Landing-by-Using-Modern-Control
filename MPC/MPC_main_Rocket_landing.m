clc; 
clear
close all;

%% Constant parameters of the rocket and state space modeling
l=70;%meter
g=9.81;%gravity
m=505.806;%mass of the rocket as kg
j=(m*(l^2))/12;%inertia kgm^2
u1_avg=100000/2;%side motors thrust force as kg
u2_avg=154674.99817/2;%average thrust force as kg

% state space matrix
Ac=[0 1 0 0 0 0;
   0 0 0 0 -u2_avg/m 0;
   0 0 0 1 0 0;
   0 0 0 0 u1_avg/m 0;
   0 0 0 0 0 1;
   0 0 0 0 0 0];
% control matrix
Bc=[0 0;
   1/m 0;
   0 0;
   0 ((u2_avg-g)*m)/(u2_avg-(g*m));
   0 0;
   l/(2*j) 0];
% out-put matrix
Cc=[1 0 0 0 0 0 0 0 0;
    0 0 1 0 0 0 0 0 0;
    0 0 0 0 1 0 0 0 0];

%% constructing augmented system matrises for mpc control
%desired system response (x_tilda=a*x) 
a=-2; % for x position
b=-2; % for y position
c=-2; % for tehta

% Augmented A matrix
size_ex_Ac=size(Ac,2); % size of the not augmented Ac matrix
Ac(size(Ac,1)+1,size(Ac,2)+1)=a;
Ac(size(Ac,1)+1,size(Ac,2)+1)=b;
Ac(size(Ac,1)+1,size(Ac,2)+1)=c;


% Augmented B matrix
Bc(size(Bc,1)+1,:)=0;
Bc(size(Bc,1)+1,:)=0;
Bc(size(Bc,1)+1,:)=0;

sys=ss(Ac,Bc,Cc,zeros(3,2));
eig(sys)
%% error and Q matrix for cost function
% Error matrix
gama=1;
e = [-1 0 -1 0 -1 0 1 1*gama 1];
% Q matrix
Q=zeros(size(Ac,1),size(Ac,2));
for i=1:1:size(Ac,2)
    if e(:,i)~=0 && i<=size_ex_Ac 
        Q(:,i)=e*(-1);
    elseif i>size_ex_Ac 
        Q(:,i)=e;
    end
end
%Q(:,8)=Q(:,8)*gama;


%% solving the matrix riccati equation reverse in time domain to find optimum R matrix

Ro=10^(-14);

init_for_riccati=zeros(45,1);
tspan = [0 -(10^3)];
options=odeset('RelTol',1e-9,'AbsTol',1e-9);
[t,x] = ode15s(@(t,x) matrix_riccati(Ac,Bc,Q,Ro,x,t), tspan, init_for_riccati,options); % solving riccati

%constructing legend matrix
r= sym('r', [9 9]);
A_2 = tril(r,0) + tril(r,-1).';
count=1;
for l=1:1:size(Ac,1) %sutun
    for k=l:1:size(Ac,2)%satir
        legends(count,1)=A_2(k,l);
        count=count+1;
    end
end
leg_str=cellstr(string(legends));
%ploting r values to visualize converged values
i=1;
j=1;
for g=1:1:45
    if mod(g,9)==0
       figure(i)
       i=i+1;
       j=1;
    end
    subplot(3,3,j);
    plot(t,x(:,g),'LineWidth',2,'Color','g');
    legend (leg_str(g,1))
    grid on
    j=j+1;
end
%constructing R matrix with riccati solutions
count=1;
R=zeros(9,9);
for l=1:1:size(Ac,1) %sutun
    for k=l:1:size(Ac,2)%satir
        R(k,l)=x(size(x,1),count);
        count=count+1;
    end
end
R = tril(R,0) + tril(R,-1).';

%% calculating optimal feed back gain matrix K

K_T= (1/Ro) * Bc' * R;

%% system response
z0=[200;20000;2000;400;30;2;200;2000;30];
tspan = [0 200];
[t,z] = ode15s(@(t,z) response(K_T,Ac,Bc,z,t),tspan,z0);
figure
hold on
grid on
for u=1:2:size(Ac,2)-4
plot(t,z(:,u))
end




