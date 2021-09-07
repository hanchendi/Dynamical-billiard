clear
clc

l_R=2;
A=pi;
R=sqrt(A/(2*l_R+pi));
l=R*l_R;
N=1000;%Boundary discretisize
M_total=1000;%total poins choose
dx=l/N;
dtheta=pi/N;
N_total=10;%How many bounce in one time
%B=zeros(M_total,N_total);
Length_B=zeros(M_total,N_total);
Theta_B=zeros(M_total,N_total);

xx1=-l/2+dx:dx:l/2;
xx2=l/2-dx:-dx:-l/2;
theta=dtheta-pi/2:dtheta:pi/2;
x_ex(1:N)=xx1;
y_ex(1:N)=-ones(1,N)*R;
theta_ex(1:N)=zeros(1,N);
kappa(1:N)=99999;

x_ex(N+1:2*N)=l/2+R*cos(theta);
y_ex(N+1:2*N)=R*sin(theta);
theta_ex(N+1:2*N)=dtheta:dtheta:pi;
kappa(N+1:2*N)=R;

x_ex(2*N+1:3*N)=xx2;
y_ex(2*N+1:3*N)=ones(1,N)*R;
theta_ex(2*N+1:3*N)=ones(1,N)*pi;
kappa(2*N+1:3*N)=99999;

x_ex(3*N+1:4*N)=-l/2+R*cos(theta+pi);
y_ex(3*N+1:4*N)=R*sin(theta+pi);
theta_ex(3*N+1:4*N)=dtheta+pi:dtheta:2*pi;
kappa(3*N+1:4*N)=R;

Lya_final=zeros(1,M_total);

for o=1:1
    B(1)=randperm(4*N,1);
    angle_begin=rand*pi+theta_ex(B(1));
    
    x_test=x_ex-x_ex(B(1));
    y_test=y_ex-y_ex(B(1));
    sin_test=y_test./sqrt(x_test.^2+y_test.^2);
    cos_test=x_test./sqrt(x_test.^2+y_test.^2);
    
    B_test=abs(cos_test-cos(angle_begin))+abs(sin_test-sin(angle_begin));
    B_test_min=find(B_test==min(B_test));
    B(2)=B_test_min(1);
    Length_B(o,1)=sqrt((x_ex(B(2))-x_ex(B(1)))^2+(y_ex(B(2))-y_ex(B(1)))^2);

    for t=2:N_total
    
        t1=(y_ex(B(t))-y_ex(B(t-1)))/(x_ex(B(t))-x_ex(B(t-1)));
        t2=x_ex(B(t))-x_ex(B(t-1));
        if t1>=0 && t2>=0
            angle0=atan(t1);
        elseif t1<0 && t2<0
            angle0=atan(t1)+pi;
        elseif t1>=0 && t2<0
            angle0=atan(t1)+pi;
        elseif t1<0 && t2>=0
            angle0=atan(t1)+2*pi;
        end
    
        x_test=x_ex-x_ex(B(t));
        y_test=y_ex-y_ex(B(t));
        sin_test=y_test./sqrt(x_test.^2+y_test.^2);
        cos_test=x_test./sqrt(x_test.^2+y_test.^2);
    
        angle_reflect=2*pi-angle0+2*theta_ex(B(t));
        B_test=abs(cos_test-cos(angle_reflect))+abs(sin_test-sin(angle_reflect));
        B_test_min=find(B_test==min(B_test));
        B(t+1)=B_test_min(1);
        Length_B(o,t)=Length_B(o,t-1)+sqrt((x_ex(B(t))-x_ex(B(t-1)))^2+(y_ex(B(t))-y_ex(B(t-1)))^2);
        Theta_B(o,t-1)=abs(pi-abs(mod(angle_reflect,2*pi)-angle0))/2;
    end
    disp(o)
end

plot(x_ex,y_ex,'k','linewidth',1.5);hold on;plot(x_ex(B(1)),y_ex(B(1)),'r*');hold on
plot([x_ex(B(1)),x_ex(B(2))],[y_ex(B(1)),y_ex(B(2))],'r');hold on
for i=2:N_total-1
    plot([x_ex(B(i)),x_ex(B(i+1))],[y_ex(B(i)),y_ex(B(i+1))],'b');hold on
end
axis equal