a = -1.;
b = 1.;
num = 100;

x = linspace(a, b, num)';
dx = (b - a)/(num-1);
f = @(y) abs(y);
df_true = @(y) sign(y);

func = f(x) + 0.05*normrnd(0,1,num,1);


df = TVRegDiff( func , 7000, 0.2, [0; diff(func) ; 0], 'small', 1e-6, dx, 0, 0);
df = df(2:end);

x_FD = (x + 0.5*dx);
x_FD = x_FD(2:end);
df_FD = diff(func)/dx;




fontsize=18.0;
figure(1)
clf
plot(x, func, 'DisplayName', 'TVR');
hold
set(gca,'Fontsize',fontsize-6.0);
xlabel('$x$', 'Interpreter','Latex','Fontsize',fontsize);
ylabel('$f(x)$', 'Interpreter','Latex','Fontsize',fontsize);





figure(2)
clf
plot(x, df, 'DisplayName', 'TVR');
hold
plot(x_FD, df_FD, 'DisplayName', 'FD');
plot(x, df_true(x), 'DisplayName', 'df/dx')
legend('location', 'northeast')
set(gca,'Fontsize',fontsize-6.0);
xlabel('$x$', 'Interpreter','Latex','Fontsize',fontsize);
ylabel('$f''(x)$', 'Interpreter','Latex','Fontsize',fontsize);
axis([a b -max(df_FD) max(df_FD)])


figure(3)
clf
semilogy(x, abs(df_true(x)-df), 'DisplayName', 'TVR')
hold
semilogy(x_FD, abs(df_true(x_FD)-df_FD), 'DisplayName', 'FD')
legend()
set(gca,'Fontsize',fontsize-6.0);
xlabel('$x$', 'Interpreter','Latex','Fontsize',fontsize);
ylabel('$|f''(x) - \tilde f''(x)|$', 'Interpreter','Latex','Fontsize',fontsize);
axis([a b min(abs(df_true(x)-df)) max(abs(df_true(x_FD)-df_FD))])
