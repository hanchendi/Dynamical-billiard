clear
clc

r1=1;
r2=0.5;
deviate=0.3;
A=pi;
N=10000;%Boundary discretisize
M_total=200;%total poins choose
dtheta=2*pi/N;
N_total=100;%How many bounce in one time
N_line=200;%test inside point

theta1=asin((r2+deviate)/r1);
theta2=asin((r2-deviate)/r1);

Length_B=zeros(M_total,N_total);
Theta_B=zeros(M_total,N_total);
Point_B=zeros(M_total,N_total);

theta=dtheta:dtheta:2*pi;
x_ex(1:N)=r1*cos(theta);
y_ex(1:N)=r1*sin(theta);
theta_ex(1:N)=mod(theta+pi/2,2*pi);

x_ex(N+1:2*N)=r2*cos(theta)+deviate;
y_ex(N+1:2*N)=r2*sin(theta);
theta_ex(N+1:2*N)=mod(theta+pi/2,2*pi);

x_ex_test=r2*cos(linspace(0,2*pi,500))+deviate;
y_ex_test=r2*sin(linspace(0,2*pi,500));
theta_choose=asin(linspace(10^(-5),1-10^(-5),M_total))*2;
for o=89:89
    B(1:N_total)=zeros(1,N_total);
    %B(1)=randperm(N,1);
    %angle_begin=rand*pi+theta_ex(B(1));
    B(1)=5000;
    angle_begin=theta_choose(o)+theta_ex(B(1));
    angle_test=angle_begin-theta_ex(B(1));
    
    %if abs(angle_test-pi/2)>theta1
    %    theta_0=abs(angle_test-pi/2);
    %    Length_B(o,1)=1./(2*r1*cos(theta_0));
    %    Theta_B(o,1:N_total)=ones(1,N_total)*theta_0;
    %    for i=2:N_total
    %        Length_B(o,i)=Length_B(o,i-1)+Length_B(o,1);
    %    end
    %    
    %else
        x_test=x_ex-x_ex(B(1));
        y_test=y_ex-y_ex(B(1));
        sin_test=y_test./sqrt(x_test.^2+y_test.^2);
        cos_test=x_test./sqrt(x_test.^2+y_test.^2);
    
        B_test=abs(cos_test-cos(angle_begin))+abs(sin_test-sin(angle_begin));
            test=0;
            while test==0
                a=find(B_test==min(B_test));
                for i=1:length(a)
                    xx_test=linspace(x_ex(B(1)),x_ex(a(i)),N_line);
                    yy_test=linspace(y_ex(B(1)),y_ex(a(i)),N_line);
                    for k=1:N_line
                        P=inpolygon(xx_test(k),yy_test(k),x_ex_test,y_ex_test);
                        if P==1
                            test=0;
                            B_test(a(i))=NaN;
                            break
                        else
                            test=1;
                        end
                    end
                end
            end
            a=find(B_test==min(B_test));   
            B(2)=a(1);
        Length_B(o,1)=sqrt((x_ex(B(2))-x_ex(B(1)))^2+(y_ex(B(2))-y_ex(B(1)))^2);

         for t=2:N_total   
            t1=(y_ex(B(t))-y_ex(B(t-1)))/(x_ex(B(t))-x_ex(B(t-1)));
            t2=x_ex(B(t))-x_ex(B(t-1));
            Length_B(o,1)=sqrt((x_ex(B(2))-x_ex(B(1)))^2+(y_ex(B(2))-y_ex(B(1)))^2);
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
            test=0;
            while test==0
                a=find(B_test==min(B_test));
                for i=1:length(a)
                    xx_test=linspace(x_ex(B(t)),x_ex(a(i)),N_line);
                    yy_test=linspace(y_ex(B(t)),y_ex(a(i)),N_line);
                    for k=1:N_line
                        P=inpolygon(xx_test(k),yy_test(k),x_ex_test,y_ex_test);
                        if P==1
                            test=0;
                            B_test(a(i))=NaN;
                            break
                        else
                            test=1;
                        end
                    end
                end
            end
            a=find(B_test==min(B_test));   
            B(t+1)=a(1);
            Length_B(o,t)=Length_B(o,t-1)+sqrt((x_ex(B(t))-x_ex(B(t-1)))^2+(y_ex(B(t))-y_ex(B(t-1)))^2);
            Theta_B(o,t-1)=abs(pi-abs(mod(angle_reflect,2*pi)-angle0))/2;
          end
    %end
    Point_B(o,:)=B(1:N_total);
    disp(o)
end

for i=1:N_total
    plot([x_ex(B(i)) x_ex(B(i+1))],[y_ex(B(i)) y_ex(B(i+1))],'b');hold on
end
plot(x_ex(1:N),y_ex(1:N),'r','linewidth',1.5);hold on
plot(x_ex(N+1:2*N),y_ex(N+1:2*N),'r','linewidth',1.5);hold on
axis equal;axis off
%save([pwd,'/Theta_B_',num2str(deviate),'.mat'],'Theta_B')
%save([pwd,'/Length_B_',num2str(deviate),'.mat'],'Length_B')
%save([pwd,'/Point_B_',num2str(deviate),'.mat'],'Point_B')