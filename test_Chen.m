
a = 35.;
b = 3.;
c = 28.;

Chen = @(xv, t) [a*(xv(2) - xv(1)); ...
                  (c - a)*xv(1) - xv(1)*xv(3) + c*xv(2); ... 
                  xv(1)*xv(2) - b*xv(3)];

% integrate the system
x0 = [-8; 7; 27];
T = 80.;
DT = 1e-4;
IOSTEP = 20;
[X,t] = AB3(Chen,x0,T,DT,IOSTEP);
X = X';



Xnoisy = X +  1.*randn(size(X));


% compute the real derivative
N = length(X);
Xdot = zeros(N,3);
for i=1:N
    Xdot(i,:) = Chen(X(i,:), 0.)';
end

iters = 20;
alpha = 10.;
% compute TVR derivative
Xdot_TVRa(:,1) = TVRegDiff( Xnoisy(:,1), iters, alpha, [], 'small', 1e12, DT*IOSTEP, 0, 0);
Xdot_TVRa(:,2) = TVRegDiff( Xnoisy(:,2), iters, alpha, [], 'small', 1e12, DT*IOSTEP, 0, 0);
Xdot_TVRa(:,3) = TVRegDiff( Xnoisy(:,3), iters, alpha, [], 'small', 1e12, DT*IOSTEP, 0, 0);
Xdot_TVRb(:,1) = Xdot_TVRa(2:end,1);
Xdot_TVRb(:,2) = Xdot_TVRa(2:end,2);
Xdot_TVRb(:,3) = Xdot_TVRa(2:end,3);



fontsize=18.0;
figure(1)
clf
hold
plot(t, Xdot_TVRb(:,1), 'r-');
plot(t, Xdot(:,1), 'b-');
set(gca,'Fontsize',fontsize-6.0);
xlabel('$t$', 'Interpreter','Latex','Fontsize',fontsize);
ylabel('$x_1''(t) $', 'Interpreter','Latex','Fontsize',fontsize);
figure(2)
clf
hold
plot(t, Xdot_TVRb(:,2), 'r-');
plot(t, Xdot(:,2), 'b-');
set(gca,'Fontsize',fontsize-6.0);
xlabel('$t$', 'Interpreter','Latex','Fontsize',fontsize);
ylabel('$x_2''(t) $', 'Interpreter','Latex','Fontsize',fontsize);
figure(3)
clf
hold
plot(t, Xdot_TVRb(:,3), 'r-');
plot(t, Xdot(:,3), 'b-');
set(gca,'Fontsize',fontsize-6.0);
xlabel('$t$', 'Interpreter','Latex','Fontsize',fontsize);
ylabel('$x_3''(t) $', 'Interpreter','Latex','Fontsize',fontsize);



Xnew1(:,1) = cumsum(Xdot_TVRb(:,1))*IOSTEP*DT;
Xnew1(:,2) = cumsum(Xdot_TVRb(:,2))*IOSTEP*DT;
Xnew1(:,3) = cumsum(Xdot_TVRb(:,3))*IOSTEP*DT;
Xnew1(:,1) = Xnew1(:,1) - (mean(Xnew1(:,1)) - mean(Xnoisy(:,1)));
Xnew1(:,2) = Xnew1(:,2) - (mean(Xnew1(:,2)) - mean(Xnoisy(:,2)));
Xnew1(:,3) = Xnew1(:,3) - (mean(Xnew1(:,3)) - mean(Xnoisy(:,3)));


order = 3;
Theta = polyTheta(Xnew1, order);
lambda = 0.45;
Xi = SparseRegression(Theta,Xdot_TVRb,lambda,3)

newfunc = @(x,t) (polyTheta(x',order)*Xi)';
[Xnew,tnew] = AB3(newfunc,x0,T,DT,IOSTEP);
Xnew = Xnew';



figure(4)
clf
plot3(Xnoisy(:,1), Xnoisy(:,2), Xnoisy(:,3));
view(57,23);
axis([-50 50 -40 40 0 50]);
set(gca,'Fontsize',fontsize-6.0);
xlabel('$x_1$', 'Interpreter','Latex','Fontsize',fontsize);
ylabel('$x_2$', 'Interpreter','Latex','Fontsize',fontsize);
zlabel('$x_3$', 'Interpreter','Latex','Fontsize',fontsize);



figure(5)
clf
plot3(Xnew(:,1), Xnew(:,2), Xnew(:,3));
view(57,23);
axis([-50 50 -40 40 0 50]);
set(gca,'Fontsize',fontsize-6.0);
xlabel('$x_1$', 'Interpreter','Latex','Fontsize',fontsize);
ylabel('$x_2$', 'Interpreter','Latex','Fontsize',fontsize);
zlabel('$x_3$', 'Interpreter','Latex','Fontsize',fontsize);


