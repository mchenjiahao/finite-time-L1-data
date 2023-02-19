clear
%Known value
A=[0 0.001;0.006 0]; Ar=[-1.2 0.1;0.2 -1.8];
E=[0.02;0.01]; Er=[0.03;0.03];
C=[0.01 0.03]; Cr=[0.02 0.03];
B=[0.03;0.02];
%Given  value
T=10; h=2; l=[1;1];
gamma=1;
p=0.1;
%for p=0.1;0.1;1;
    ka=sdpvar(1,1); %赋予初始值
    kb=sdpvar(1,1);
    K1=[ka kb];
    kc=sdpvar(1,1); %赋予初始值
    kd=sdpvar(1,1);
    K2=[kc kd];
    
    v1=sdpvar(2,1);
    v2=sdpvar(2,1);
    v=[v1;v2];
    F=[-0.001-0.03*kb<=0,-0.006-0.02*ka<=0];
    F=[B*K2>=0,v1>=0,v2>=0];
    F=[0<=K2<=0.1];
    F=[F,K2'*B'*v1+Ar'*v2-p.*v2<=0];
    F=[F,A'*v1+K1'*B'*v1+C'- p.*v1<=0];
    F=[F,E'*v1-gamma<=0];
    F=[F,Er'*v2-gamma<=0];
    
    k11=value(ka);
    k12=value(kb);
    k21=value(kc);
    k22=value(kd);

f=@(t,x)[-1.2*x(3)+0.1*x(4)+0.003;0.2*x(3)-1.8*x(4)+0.003;...
     0.03*k11*x(1)+0.001*x(2)+0.03*k12*x(2)+0.03*k21*x(3)+0.03*k22*x(3)+0.02*exp(-t);...
   (0.006+0.02*k11)*x(1)+0.02*k12*x(2)+0.02*k21*x(3)+0.02*k22*x(4)+0.01*exp(-t)];
[t,z]=ode45(f,[0,10],[0;0;2;1]);

x1=z(:,1);
x2=z(:,2);
xr1=z(:,3);
xr2=z(:,4);
y=0.01*x1+0.03*x2;
yr=0.02*xr1+0.03*xr2;
e=y-yr;

objective=norm(e,1);
solvesdp(F,objective);
norm(e,1)
   