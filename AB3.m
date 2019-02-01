function [y,t] = AB3(fun,y0,T,DT,IOSTEP)
% Solves a system of ODEs using the third-order Adams-Bashforth Scheme, 
% and the Heun scheme for the first two steps.
%
% INPUTS: fun    -- function handle representing 1st derivative of y
%         y0     -- column vector (length n) of intial values y(0)
%         T      -- period of integration
%         DT     -- timestep
%         IOSTEP -- input/output step
%
% OUTPUTS: y -- n x S matrix collecting S snapshots of the solution
%          t -- vector collecting the time values of each output timestep

NSTEPS = ceil(T/DT);
y = y0;
t = 0.;
ts = 0.;

% Heun scheme for startup
y1 = y0+0.5*DT*(fun(y0+DT*fun(y0,ts),DT)+fun(y0,ts));
ts = DT;
y2 = y1+0.5*DT*(fun(y1+DT*fun(y1,ts),ts+DT)+fun(y1,ts));

% Adams-Bashforth Scheme, 3rd order
for i=2:NSTEPS
    ts = (i-1)*DT;
    y3 = y2 + (DT/12.0)*(23*fun(y2,ts-DT)-16*fun(y1,ts-2*DT)+5*fun(y0,ts-3*DT));
    if mod(i,IOSTEP)==0
        y =[ y y3];
        t = [ t ts+DT];
    end
    y0 = y1;
    y1 = y2;
    y2 = y3;
end

end