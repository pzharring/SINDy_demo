fontsize=18.0;
noise = 0.003;

iter = 5;
alpha = 2;

X = [];
Xdot = [];
T = 60.;
DT = 1e-3;
IOSTEP = 20;

% generate mu < 1/2 data
for mu_val= 0.38:0.025:0.48
    clear Xmu Xdot_mu t;
    Func = @(xv, t) [ (xv(1)*xv(1))*(1 - xv(1)) - xv(1)*xv(2); ...
                      xv(1)*xv(2) - mu_val*xv(2)];
    x0 = [0.5; 0.575]; % outside limit cycle
    [Xmu,t] = AB3(Func,x0,T,DT,IOSTEP);
    Xmu = Xmu';
    tlength = size(Xmu, 1);
    Xmu = [Xmu mu_val*ones(size(Xmu, 1), 1)];
    Xmu = Xmu + noise*randn(size(Xmu));
    Xdot_mu(:,1) = TVRegDiff( Xmu(:,1), iter, alpha, [], 'small', 1e+2, DT*IOSTEP, 0, 0 );
    Xdot_mu(:,2) = TVRegDiff( Xmu(:,2), iter, alpha, [], 'small', 1e+2, DT*IOSTEP, 0, 0 );
    Xdot_mu(:,3) = zeros(size(Xmu, 1) + 1, 1);
    Xdot_mu = Xdot_mu(2:end,:);
    X = [X; Xmu];
    Xdot = [Xdot; Xdot_mu];
    
    clear Xmu Xdot_mu t;
    x0 = [0.5; 0.375]; % inside limit cycle
    [Xmu,t] = AB3(Func,x0,T,DT,IOSTEP);
    Xmu = Xmu';
    tlength = size(Xmu, 1);
    Xmu = [Xmu mu_val*ones(size(Xmu, 1), 1)];
    Xmu = Xmu + noise*randn(size(Xmu));
    Xdot_mu(:,1) = TVRegDiff( Xmu(:,1), iter, alpha, [], 'small', 1e+2, DT*IOSTEP, 0, 0 );
    Xdot_mu(:,2) = TVRegDiff( Xmu(:,2), iter, alpha, [], 'small', 1e+2, DT*IOSTEP, 0, 0 );
    Xdot_mu(:,3) = zeros(size(Xmu, 1) + 1, 1);
    Xdot_mu = Xdot_mu(2:end,:);
    X = [X; Xmu];
    Xdot = [Xdot; Xdot_mu];
end

% generate mu > 1/2 data
for mu_val= 0.505:0.025:0.555
    clear Xmu Xdot_mu t;
    Func = @(xv, t) [ (xv(1)*xv(1))*(1 - xv(1)) - xv(1)*xv(2); ...
                      xv(1)*xv(2) - mu_val*xv(2)];
    x0 = [0.5; 0.375];
    [Xmu,t] = AB3(Func,x0,T,DT,IOSTEP);
    Xmu = Xmu';
    tlength = size(Xmu, 1);
    Xmu = [Xmu mu_val*ones(size(Xmu, 1), 1)];
    Xmu = Xmu + noise*randn(size(Xmu));
    Xdot_mu(:,1) = TVRegDiff( Xmu(:,1), iter, alpha, [], 'small', 1e+2, DT*IOSTEP, 0, 0 );
    Xdot_mu(:,2) = TVRegDiff( Xmu(:,2), iter, alpha, [], 'small', 1e+2, DT*IOSTEP, 0, 0 );
    Xdot_mu(:,3) = zeros(size(Xmu, 1) + 1, 1);
    Xdot_mu = Xdot_mu(1:end-1,:);
    X = [X; Xmu];
    Xdot = [Xdot; Xdot_mu];
end


% plot
figure('Position', [10 10 300 600])
for k=0:12
    if k <=9  && mod(k,2)==1
        str = 'r-';
    else
        str = 'b-';
    end
    plot3(X(k*tlength + 1:(k+1)*tlength,1), ...
          X(k*tlength + 1:(k+1)*tlength,2), ...
          X(k*tlength + 1:(k+1)*tlength,3), str)
    hold on
end
xlim([0 1])
ylim([0 1])
zlim([0.35 0.6])
view(-144,10)
grid on
set(gca,'Fontsize',fontsize-6.0);
xlabel('$x$', 'Interpreter','Latex','Fontsize',fontsize);
ylabel('$y$', 'Interpreter','Latex','Fontsize',fontsize);
zlabel('$\mu$', 'Interpreter','Latex','Fontsize',fontsize);


order = 3;
Theta = polyTheta(X, order);
lambda = 0.45;
Xi = SparseRegression(Theta,Xdot,lambda,3)

newfunc = @(x,t) (polyTheta(x',order)*Xi)';

X_SINDy = [];

% generate mu < 1/2 data for learned system
for mu_val= 0.38:0.025:0.48
    clear Xmu  t;
    x0 = [0.5; 0.575; mu_val]; % outside limit cycle
    [Xmu,t] = AB3(newfunc,x0,T,DT,IOSTEP);
    Xmu = Xmu';
    X_SINDy = [X_SINDy; Xmu];
    
    clear Xmu t;
    x0 = [0.5; 0.375; mu_val]; % inside limit cycle
    DT = 1e-3;
    [Xmu,t] = AB3(newfunc,x0,T,DT,IOSTEP);
    Xmu = Xmu';
    tlength = size(Xmu, 1);
    X_SINDy = [X_SINDy; Xmu];
end

% generate mu > 1/2 data for learned system
for mu_val= 0.505:0.025:0.555
    clear Xmu t;
    x0 = [0.5; 0.375; mu_val];
    [Xmu,t] = AB3(newfunc,x0,T,DT,IOSTEP);
    Xmu = Xmu';
    tlength = size(Xmu, 1);
    X_SINDy = [X_SINDy; Xmu];
end





% plot
figure('Position', [10 10 300 600])
for k=0:12
    if k <= 9 && mod(k,2)==1
        str = 'r-';
    else
        str = 'b-';
    end
    plot3(X_SINDy(k*tlength + 1:(k+1)*tlength,1), ...
          X_SINDy(k*tlength + 1:(k+1)*tlength,2), ...
          X_SINDy(k*tlength + 1:(k+1)*tlength,3), str)
    hold on
end
xlim([0 1])
ylim([0 1])
zlim([0.35 0.6])
view(-144,10)
grid on
set(gca,'Fontsize',fontsize-6.0);
xlabel('$x$', 'Interpreter','Latex','Fontsize',fontsize);
ylabel('$y$', 'Interpreter','Latex','Fontsize',fontsize);
zlabel('$\mu$', 'Interpreter','Latex','Fontsize',fontsize);


