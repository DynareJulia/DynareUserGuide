var pi 
    y_gap 
    y_nat     
    y 
    r_nat 
    r_real 
    i 
    n 
    m_real 
    m_growth_ann 
    nu
    a  
    r_real_ann 
    i_ann 
    r_nat_ann 
    pi_ann 
    ;     

varexo eps_a 
       eps_nu 
       ;

parameters alppha 
    betta 
    rho_a 
    rho_nu   
    siggma 
    phi 
    phi_pi 
    phi_y 
    eta 
    epsilon 
    theta 
    ;

siggma = 1;
phi=1;
phi_pi = 1.5;
phi_y  = .5/4;
theta=2/3;
rho_nu =0.5;
rho_a  = 0.9;
betta = 0.99;
eta  =4;
alppha=1/3;
epsilon=6;



model(linear); 
//Composite parameters
#Omega=(1-alppha)/(1-alppha+alppha*epsilon);  //defined on page 47
#psi_n_ya=(1+phi)/(siggma*(1-alppha)+phi+alppha); //defined on page 48
#lambda=(1-theta)*(1-betta*theta)/theta*Omega; //defined on page 47
#kappa=lambda*(siggma+(phi+alppha)/(1-alppha));  //defined on page 49

//1. New Keynesian Phillips Curve eq. (21)
pi=betta*pi(+1)+kappa*y_gap;
//2. Dynamic IS Curve eq. (22)
y_gap=-1/siggma*(i-pi(+1)-r_nat)+y_gap(+1);
//3. Interest Rate Rule eq. (25)
i=phi_pi*pi+phi_y*y_gap+nu;
//4. Definition natural rate of interest eq. (23)
r_nat=siggma*psi_n_ya*(a(+1)-a);
//5. Definition real interest rate
r_real=i-pi(+1);
//6. Definition natural output, eq. (19)
y_nat=psi_n_ya*a;
//7. Definition output gap
y_gap=y-y_nat;
//8. Monetary policy shock
nu=rho_nu*nu(-1)+eps_nu;
//9. TFP shock
a=rho_a*a(-1)+eps_a;
//10. Production function (eq. 13)
y=a+(1-alppha)*n;
//11. Money growth (derived from eq. (4))
m_growth_ann=4*(y-y(-1)-eta*(i-i(-1))+pi);
//12. Real money demand (eq. 4)
m_real=y-eta*i;


//13. Annualized nominal interest rate
i_ann=4*i;
//14. Annualized real interest rate
r_real_ann=4*r_real;
//15. Annualized natural interest rate
r_nat_ann=4*r_nat;
//16. Annualized inflation
pi_ann=4*pi;
end;



shocks;
var eps_nu = 0.25^2; 
end;


stoch_simul(order=1);

