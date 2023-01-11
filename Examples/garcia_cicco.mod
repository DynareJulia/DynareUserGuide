
@#define RBC =0
//set to 1 for RBC model and to 0 for Financial Frictions Model

var c       $c$ 
    k       $k$ 
    a       $a$ 
    h       $h$ 
    d       $d$ 
    y       $y$ 
    invest  $i$  
    tb      $tb$ 
    mu_c    ${MU_C}$ 
    tb_y    ${\frac{TB}{Y}}$ 
    g_y     ${\Delta Y}$
    g_c     ${\Delta C}$
    g_invest ${\Delta I}$
    g       ${g}$
    r       ${r}$
    mu      ${\mu}$
    nu      ${\nu}$
    @#if RBC == 0
    s       ${s}$
    @# endif
; 

predetermined_variables k d;

%Define parameters
parameters beta     ${\beta}$ 
        gamma       ${\gamma}$ 
        delta       ${\delta}$
        alpha       ${\alpha}$
        psi         ${\psi}$
        omega       ${\omega}$
        theta       ${\theta}$
        phi         ${\phi}$
        dbar        ${\bar d}$
        gbar        ${\bar g}$
        rho_a       ${\rho_a}$
        rho_g       ${\rho_g}$
        rho_nu      ${\rho_\nu}$
        rho_mu      ${\rho_\mu}$
        rho_s       ${\rho_s}$
    @#if RBC == 0
        s_share     ${sshare}$
        S           ${S}$
    @# endif
;

varexo eps_a ${\varepsilon_a}$
        eps_g ${\varepsilon_g}$ 
        eps_nu ${\varepsilon_\nu}$
        eps_mu ${\varepsilon_\mu}$
    @#if RBC == 0
        eps_s ${\varepsilon_s}$
    @# endif
;

        
@#if RBC == 1
    gbar  = 1.0050; %Gross growth rate of output
    rho_g = 0.8280; %Serial correlation of innovation in permanent technology shock
    rho_a = 0.7650; %Serial correlation of transitory technology shock
    phi   = 3.3000; %Adjustment cost parameter
@# else
    gbar  = 1.009890776104921;
    rho_g = 0.323027844166870;
    rho_a = 0.864571930755821;
    phi   = 4.810804146604144;
@# endif
            
%parameters only used for Financial frictions model, irrelevant for RBC where volatilities are 0        
rho_nu =    0.850328786147732;
rho_s  = 0.205034667802314;
rho_mu = 0.906802888826967;

%From Table 2, except for psi, which is estimated for Financial Frictions model

gamma = 2; %intertemporal elasticity of substitution
delta = 1.03^4-1;%0.03; %Depreciation rate
alpha = 0.32; %Capital elasticity of the production function
omega = 1.6; %exponent of labor in utility function
theta = 1.4*omega;
beta = 0.98^4;%0.98;%discount factor
dbar = 0.007;

@#if RBC == 1
    psi = 0.001;
@# else
    psi = 2.867166241970346; %parameter governing the debt elasticity of the interest rate. 
    s_share = 0.10; %Share of public spending in GDP
@# endif
        
        
model;
#RSTAR = 1/beta * gbar^gamma; %World interest rate
%1. Interest Rate
r = RSTAR + psi*(exp(d-dbar) - 1)+exp(mu-1)-1;

%2. Marginal utility of consumption
mu_c = nu * (c - theta/omega*h^omega)^(-gamma);

%3. Resource constraint (see the remark on the fixed typo in the preamble)
@#if RBC == 1
    y= log(tb) + c + invest + phi/2 * (k(+1)/k*g -gbar)^2*k;
@# else
    y= log(tb) + c + s + invest + phi/2 * (k(+1)/k*g -gbar)^2*k;
@# endif

%4. Trade balance
log(tb)= d - d(+1)*g/r;

%5. Definition output
y= a*k^alpha*(g*h)^(1-alpha);

%6. Definition investment
invest= k(+1)*g - (1-delta) *k;

%7. Euler equation
mu_c= beta/g^gamma*r*mu_c(+1);

%8. First order condition labor
theta*h^(omega-1)=(1-alpha)*a*g^(1-alpha)*(k/h)^alpha;

%9. First order condition investment
mu_c*(1+phi*(k(+1)/k*g-gbar))= beta/g^gamma*mu_c(+1)*(1-delta+alpha*a(+1)*(g(+1)*h(+1)/k(+1))^(1-alpha) +phi*k(+2)/k(+1)*g(+1)*(k(+2)/k(+1)*g(+1)-gbar) - phi/2*(k(+2)/k(+1)*g(+1)-gbar)^2);

%10. Definition trade-balance to output ratio
log(tb_y) = log(tb)/y; 

%11. Output growth
g_y= y/y(-1)*g(-1);

%12. Consumption growth
g_c = c/c(-1)*g(-1);

%13. Investment growth
g_invest = invest/invest(-1)*g(-1);

%14. LOM temporary TFP
log(a)=rho_a * log(a(-1))+eps_a; 

%15. LOM permanent TFP Growth
log(g/gbar)=rho_g*log(g(-1)/gbar)+eps_g; 

%16. Preference shock
log(nu) =rho_nu * log(nu(-1))+eps_nu;

%17. exogenous stochastic country premium shock
log(mu)= rho_mu * log(mu(-1))+eps_mu;

@#if RBC == 1
    
@# else
    %18. Exogenous spending shock
    log(s/S)= rho_s * log(s(-1)/S) + eps_s;
@# endif

end;

steady_state_model;
    r   = 1/beta*gbar^gamma; %World interest rate
    d   = dbar; %foreign debt
    k_over_gh =((gbar^gamma/beta-1+delta)/alpha)^(1/(alpha-1)); %k/(g*h)
    h   = ((1-alpha)*gbar*k_over_gh^alpha/theta)^(1/(omega-1)); %hours
    k   = k_over_gh*gbar*h; %capital
    invest = (gbar-1+delta)*k; %investment
    y   = k^alpha*(h*gbar)^(1-alpha); %output
    @#if RBC == 1
        s = 0;    
    @# else
        s = y*s_share;
        S = s;
    @# endif
    c   = (gbar/r-1)*d +y-s-invest; %Consumption
    tb  = y - c - s - invest; %Trade balance
    tb_y = tb /y;
    mu_c = (c - theta/omega*h^omega)^(-gamma); %marginal utility of wealth
    a   = 1; %productivity shock 
    g   = gbar; %Growth rate of nonstationary productivity shock
    g_c = g;
    g_invest = g;
    g_y = g;
    nu  = 1;
    mu  = 1;
    tb  = exp(tb);
    tb_y= exp(tb_y);
end;

shocks;
@#if RBC == 1
    var eps_a; stderr 0.0270;
    var eps_g; stderr 0.0300;
    var eps_nu; stderr  0;
    var eps_mu; stderr  0;
@# else
    var eps_a; stderr 0.033055089525252;
    var eps_g; stderr 0.010561526060797;
    var eps_nu; stderr  0.539099453618175;
    var eps_s; stderr   0.018834174505537;
    var eps_mu; stderr  0.057195449717680;
@# endif
end;


stoch_simul(loglinear,order=1,irf=0) g_y g_c g_invest tb_y;
