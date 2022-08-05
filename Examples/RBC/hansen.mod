@#define indivisible_labor=1

var c 
    w 
    r 
    y 
    h 
    k 
    i
    lam
    z; 
    
varexo eps_a;

parameters beta 
    delta 
    theta 
    gamma 
    A 
    h_0 
    sigma_eps 
    B ;

//Calibration, p. 319
beta = 0.99;
delta = 0.025;
theta = 0.36;
gamma = 0.95;
A = 2;
sigma_eps=0.00712;
h_0=0.53;

model;
//1. Euler Equation
1/c = beta*((1/c(+1))*(r(+1) +(1-delta)));
//2. Labor FOC
@#if indivisible_labor
    (1-theta)*(y/h) = B*c;
@#else
    (1-theta)*(y/h) = A/(1-h)*c;
@#endif
//3. Resource constraint
c = y +(1-delta)*k(-1) - k;
//4. LOM capital
k= (1-delta)*k(-1) + i;
//5. Production function
y = lam*k(-1)^(theta)*h^(1-theta);
//6. Real wage
r = theta*(y/k(-1));
//7. Real interest rate
w = (1-theta)*(y/h);
//8. LOM TFP
log(lam)=gamma*log(lam(-1))+eps_a;
//9. Productivity
z= y/h;
end;

steady_state_model;
B=-A*(log(1-h_0))/h_0; 
lam = 1;
@#if indivisible_labor
    h = (1-theta)*(1/beta -(1-delta))/(B*(1/beta -(1-delta)-theta*delta));
@#else
    h = (1+(A/(1-theta))*(1 - (beta*delta*theta)/(1-beta*(1-delta))))^(-1);
@#endif
k = h*((1/beta -(1-delta))/(theta*lam))^(1/(theta-1));
i = delta*k;
y = lam*k^(theta)*h^(1-theta);
c = y-delta*k;
r =  1/beta - (1-delta);
w = (1-theta)*(y/h);
z = y/h;
end;

steady;

shocks;
var eps_a; 
stderr sigma_eps;
end;

check;
steady;