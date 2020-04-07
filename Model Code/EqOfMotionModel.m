function [DEuler]=EqOfMotionModel()
      clear all
      m=1;  %Disk mass (g)
      g=9.8;      %Gravitation constant (m/s/s)
      a=1.0;  %Disk radius (m)
      h=2;    %Disk height (mm)
      C=[0.03, 0.073,0.185]';
      [I,J,rG]=CreateConstantMatrices(m,g,h,a);
      IntegrationTimes=0:0.0001:0.05;
      
      Euler0=deg2rad([65;0;0]);
      Euler1=deg2rad([64.999;0;0]);
      DEuler0=[0;0;80];
      DEuler1=[0;0;80];
      DE2o=@(theta)[1,0,0;0,1,sin(theta);0,0,cos(theta)];
      omega0=DE2o(Euler0(1))*DEuler0;
      omega1=DE2o(Euler1(1))*DEuler1;
      position0=[0;0;a*sin(Euler0(1))+(h/2)*cos(Euler0(1))];
      position1=[cos(Euler1(3))*(cos(Euler1(1))-cos(Euler0(1)));...
          sin(Euler1(3))*(cos(Euler1(1))-cos(Euler0(1)));...
          a*sin(Euler1(1))+(h/2)*cos(Euler1(1))];
      
      %ODE45 Method
      opts=odeset('RelTol',1e-4,'AbsTol',1e-5,'MaxStep',1e-1,...
          'InitialStep',1e-10,'Events',@(t,x)FallingEventFcn(t,x,a,h));
      [t,x]=ode45(@(t,x)func(t,x,I,J,rG,g,m,h,C),IntegrationTimes,...
          [position0;Euler0;omega0],opts);
      [t2,x2]=ode45(@(t,x)func(t,x,I,J,rG,g,m,h,C),IntegrationTimes,...
          [position1;Euler1;omega1],opts);
      
      %Position of Centre Of Mass
      figure(1)
      plot3(x(:,1),x(:,2),x(:,3))
      xlabel('x');
      ylabel('y');
      zlabel('z');
      %saveas(figure(1),'CentreOfMassMotion.png')
      
      %Omega vector
      figure(2)
      plot3(x(:,7),x(:,8),x(:,9))
      xlabel('\omega_1');
      ylabel('\omega_2');
      zlabel('\omega_3');
      %saveas(figure(2),'AngularMomentumVector.png')
      
      %Theta vs Time
      figure(3)
      plot(t,x(:,4)./pi)
      xlabel('Time, t');
      ylabel('\theta/\pi');
      %saveas(figure(3),'Theta.png')
      
      %Lyapunov exponent plotter
      LyapunovDiff=sqrt((x(:,1)-x2(:,1)).^2+(x(:,2)-x2(:,2)).^2+(x(:,3)-...
          x2(:,3)).^2);
      figure(4)
      plot(IntegrationTimes,log(LyapunovDiff))
      P = polyfit(IntegrationTimes(74:187),log(LyapunovDiff(74:187)'),1);
      LyapunovExponent = P(1)
      LogInitialSeparation = log(LyapunovDiff(1))
      YIntercept = P(2)
      xlabel('Time, t');
      ylabel('ln(d_t)');
      saveas(figure(4),'Lyapunov(65v0.001,0,0,0,0,80).png')
end

function [I,J,rG]=CreateConstantMatrices(m,g,h,a)
  I=diag(m.*[(3*a^2+h^2)/12,(a^2)/2,(3*a^2+h^2)/12]);
  J=I+m.*[a^2+(h^2)/4,0,0;0,a^2,-a*h/2;0,-a*h/2,(h^2)/4];
  rG=[0;h/2;a];
end

function [E2C]=Euler2Cartesian(Euler)
  theta=Euler(1);
  psi=Euler(2);
  phi=Euler(3);
  R1=@(angle)[1,0,0;0,cos(angle),sin(angle);0,-sin(angle),cos(angle)];
  R2=@(angle)[cos(angle),0,-sin(angle);0,1,0;sin(angle),0,cos(angle)];
  R3=@(angle)[cos(angle),sin(angle),0;-sin(angle),cos(angle),0;0,0,1];
  
  C2E=R1(theta)*R3(phi);
  E2C=inv(C2E);
end

function [Dx]=func(t,x,I,J,rG,g,m,h,C)
  Euler=x(4:6);
  omega=x(7:9);
  theta=Euler(1);
  o2DE=[1,0,0;0,1,-tan(theta);0,0,1/cos(theta)];
  Omega=omega-[0;omega(2)-omega(3)*tan(theta);0];
  LG=I*omega;
  gvec=-g*[0;sin(theta);cos(theta)];
  vc=cross(omega,rG);
  taudrag=-C.*omega.*abs(omega);
  Domega=J\(-cross(Omega,LG)-m*cross(rG,cross(Omega,vc))+m*...
      cross(rG,gvec))+taudrag;
  DEuler=o2DE*omega;
  Dposition=Euler2Cartesian(Euler)*vc;
  Dx=[Dposition;DEuler;Domega];
end

function [value,isterminal,direction] = FallingEventFcn(t,x,a,h)
  isterminal=[1,1];
  direction=[0,0];
  value=[x(4); x(4)-pi/2];
end
