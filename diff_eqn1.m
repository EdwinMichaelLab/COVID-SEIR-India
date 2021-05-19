%% diff_eqn1.m:

% FUNCTION NAME:
%   diff_eqn1
%
% DESCRIPTION:
%   This function defines our system of ODEs that we are integrating
%   via ode45().
%
% INPUTS:
%   t: Array of timesteps.
%   pop0: Array of initial class values
%   P: Array of parameters.
%
% OUTPUT:
%   None.

function dPop = diff_eqn1(t,pop0,P,mr,vacdata)

%% Define variables
beta = P(1);
alpha = P(2);
sigma = P(3);
rho = P(4);
gammaA = P(5);
gammaM = P(6);
gammaH = P(7);
gammaC = P(8);
delta1 = P(9);
delta2 = P(10);
delta3 = P(11);
m  = P(12);
lockdown_ratio = P(13);
epsilon = P(14);
x1 = P(15);
x2 = P(16);
x3 = P(17);
d = P(18);
q = 0; %start at zero, switch on after date. P(19);
quarantine_start = P(20);
z = P(21);
gammaQ = P(22);
deltaQ = P(23);
p = P(24);
beta2 = P(25);
vac_start = P(26);
movement_range = P(27);

S = pop0(1);
E = pop0(2);
IA = pop0(3);
IP = pop0(4);
IM = pop0(5);
IH = pop0(6);
IC = pop0(7);
D = pop0(8);
R1 = pop0(9);
R2 = pop0(10);
V = pop0(11);
S2 = pop0(12);
R3 = pop0(13);
B = pop0(14);

% 2ND VARIANT
E2 = pop0(15);
IA2 = pop0(16);
IP2 = pop0(17);
IM2 = pop0(18);
IH2 = pop0(19);
IC2 = pop0(20);
R4 = pop0(21);
D2 = pop0(22);


% Define Fixed Quarantine Parameters
f = 0.10;
vac_rate = 0; %(5% / day?)
boost_rate = 0.025;%0.025; now directly moving using data!
vac_efficacy = 0.75;
boost_efficacy = 0.90;
vac_efficacy2 = 0.70;
boost_efficacy2 = 0.85;

% vac_efficacy = 1.0;
% boost_efficacy = 1.0;
% vac_efficacy2 = 1.0;
% boost_efficacy2 = 1.0;

%% Make intervention modifications
% if after quarantine start, switch
% if t > quarantine_start
%     q = P(19) * c2switch(t - quarantine_start, 4, false);
%     %q = 0;
% end
q=0;

%lambda = alpha/mr(floor(t));

% firstmr = round(mr(1), 4);
% if firstmr == 0.5178 % ALACHUA
%     rels = [409 440 470 501 532];
% elseif firstmr == 0.4336 % HILLSBOROUGH
%     rels = [410 441 471 502 533];
% elseif firstmr == 0.5523 % Miami
%     rels = [413 444 474 505 536];
% end
IExtra = 0;
% % %411 May1st, 442 jun1, 472 July1, 503 aug1
% if t >= rels(5)
%     d = 1;
%     lockdown_ratio = mr(floor(t));
%     lockdown_fraction = lockdown_ratio/(1+lockdown_ratio);
%     lockdown_fraction = lockdown_fraction * 0.1;
%     lockdown_ratio = lockdown_fraction/(1-lockdown_fraction);
%     lambda = alpha/lockdown_ratio;
%     IExtra = 0/3725688;
% % % else
% if t < 427+30
%     d = d*0.50;
% end
% % 
%     IExtra = 0/3725688;
% end

% if t < (vac_start + 28)
%     boost_rate = 0;
% else
%     boost_rate = 0.0015/2;
% end

if S2 ~= 0
    alpha = alpha*2;
    lambda = lambda*2;
end
beta2 = 0;
E2 = 0;
alpha=0;
lambda=0;


%% Define differential equations for each compartment
dPop = zeros(length(pop0),1);
ci = 1-(R2*0.1);

dPop(1) = -d*beta*S*(IA+IP+IM+IH+IC+IExtra) - ci*d*beta2*S*(IA2+IP2+IM2+IH2+IC2+IExtra)  - alpha*S +lambda*R1 -  S*vac_rate; % S
dPop(2) = d*beta*(S+S2)*(IA+IP+IM+IH+IC+IExtra) - sigma*rho*E - sigma*(1-rho)*E + beta*d*V*(1-vac_efficacy)*(IA+IP+IM+IH+IC)  + beta*d*B*(1-boost_efficacy)*(IA+IP+IM+IH+IC); % E
dPop(3) = sigma*rho*E - gammaA*IA; % IA
dPop(4) = sigma*(1-rho)*E - delta1*IP; % IP
dPop(5) = delta1*IP - x1*delta2*IM - (1-x1)*gammaM*IM; % IM
dPop(6) = x1*delta2*IM - x2*delta3*IH - (1-x2)*gammaH*IH; % IH
dPop(7) = x2*delta3*IH - (1-x3)*gammaC*IC - x3*m*IC; % IC
dPop(8) = x3*m*IC; % D
dPop(9) = alpha*S - lambda*R1; % R1
dPop(10) = gammaA*IA + (1-x1)*gammaM*IM + (1-x2)*gammaH*IH + (1-x3)*gammaC*IC; % R2
dPop(11) = S*vac_rate - beta*d*V*(1-vac_efficacy)*(IA+IP+IM+IH+IC) - ci*beta2*d*V*(1-vac_efficacy2)*(IA2+IP2+IM2+IH2+IC2) - V*boost_rate;
dPop(12) = -d*beta*(S2)*(IA+IP+IM+IH+IC+IExtra) -ci*d*beta2*(S2)*(IA2+IP2+IM2+IH2+IC2+IExtra) - alpha*S2 +lambda*R3;
dPop(13) = alpha*S2 -lambda*R3;
dPop(14) = V*boost_rate - beta*d*B*(1-boost_efficacy)*(IA+IP+IM+IH+IC+IExtra) - ci*beta2*d*B*(1-boost_efficacy2)*(IA2+IP2+IM2+IH2+IC2); % B


% 2ND VARIANT:
dPop(15) = ci*d*beta2*(S+S2)*(IA2+IP2+IM2+IH2+IC2+IExtra) - sigma*rho*E2 - sigma*(1-rho)*E2 + ci*beta2*d*V*(1-vac_efficacy2)*(IA2+IP2+IM2+IH2+IC2)  + ci*beta2*d*B*(1-boost_efficacy2)*(IA2+IP2+IM2+IH2+IC2); % E
dPop(16) = sigma*rho*E2 - gammaA*IA2; % IA
dPop(17) = sigma*(1-rho)*E2 - delta1*IP2; % IP
dPop(18) = delta1*IP2 - x1*delta2*IM2 - (1-x1)*gammaM*IM2; % IM
dPop(19) = x1*delta2*IM2 - x2*delta3*IH2 - (1-x2)*gammaH*IH2; % IH
dPop(20) = x2*delta3*IH2 - (1-x3)*gammaC*IC2 - x3*m*IC2; % IC
dPop(21) = gammaA*IA2 + (1-x1)*gammaM*IM2 + (1-x2)*gammaH*IH2 + (1-x3)*gammaC*IC2; % R2
dPop(22) = x3*m*IC2; % D4

end
