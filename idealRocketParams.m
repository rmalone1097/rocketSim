% Function that finds parameters of rocket consisting of conical nose (D), tube (L + D), and triangular fins (D)
%  assuming ideal gas properties, linear burn rate, and no drag

%% Inputs Key

% M_L = payload mass (kg)
% h_max = maximum height (apogee) (ft)
% a_max = maximum acceleration (unitless, in terms of g)
% SM = static margin ((X_CP - X_CG)/D)
% rho_s = density of structure (kg/m^3)
% rho_p = density of propelant (kg/m^3)
% sig_s = working shell stress (Pa)
% N = number of fins

%% Outputs Key

% R = mass ratio
% W_eq = exit velocity (m/s)
% t_b = burn time (s)
% p0_pa = pressure ratio
% delta_D = thickness(delta)/D
% D = diameter (m)
% L = lemgth of tube up to fins (m)
% L_D = L/D
% X_CG = center of gravity (m)
% X_CP = center of pressure (m)
% M_p = mass of propelant (kg)
% M_s = mass of structure (kg)
% M_o = total mass of vehicle (kg)
% lam = payload ratio
% eps = structural ratio
% L_p = length of propelant tank

%% Funciton
function val_vector = idealRocketParams(M_L,h_max,a_max,SM,rho_s,rho_p,sig_s,N)

%Known constants
g = 9.81;
gamma = 1.4;
a = 346;
p_a = 101.325e3;    
h_max = h_max*0.3048;

%R,W_eq,t_b assuming linear burn rate
R = 1 + a_max;
W_eq = sqrt((h_max*g)/((log(R)^2)/2 - (1/(1 + a_max))*(R*log(R) - (R - 1))));
t_b = ((R - 1)*W_eq)/(g*R);

%Mach number and stagnation pressure from isentropic relations
M = W_eq/a;
p0_pa = (1 + ((gamma - 1)/2)*M^2)^(gamma/(gamma - 1));
p_0 = p0_pa*p_a;

%Initialization of D and Lambda where D has solutions
D = (M_L/16.67):0.0001:(M_L/6.67);
lam = 0;

% Loop to find D that maximizes payload ratio
for i = 1:length(D)
    delta = (D(1,i)*p_0)/(2*sig_s);
    M_nose = delta*rho_s*pi*D(1,i)*(D(1,i) + sqrt(D(1,i)^2 + (D(1,i)/2)^2));
    M_fins = ((1/2)*D(i)^2)*3*delta*rho_s;
    M_fb = 0.7854*(D(i))*(D(i)^2 - (D(i)-2*delta)^2)*rho_s;
    % Initial L = L_p
    L = ((R - 1)*(M_L + M_nose + M_fins + M_fb) - (pi*D(i)^3)/4*rho_p)/((pi*D(i)^2/4)*rho_p - (R - 1)*pi*D(i)*rho_s*delta);
    error = 1;
    lam_log = lam;
    % Loop to find L that meets SM and is less than L_p
    while (error > 0.0001)
        M_tube = 0.7854*(L + D(i))*(D(i)^2 - (D(i)-2*delta)^2)*rho_s;
        M_s = M_nose + M_tube + M_fins;
        M_p = (R - 1)*(M_s + M_L);
        M_o = M_s + M_p + M_L;
        L_p = M_p/((pi*D(i)^2*rho_p)/4);
            if (L_p > (L + D(i)))
                L = L + 0.0001;
                continue
            end
        % X_CP calculation
        barX_n = (2/3)*D(i);
        CNa_n = 2;
        l = sqrt((D(i)/2)^2 + D(i)^2);
        m = D(i);
        a = D(i);
        b = 0;
        CNa_f = (4*N)/(1 + sqrt(1 + (2*l/(a + b))^2));
        K_fb = 1 + (D(i)/2)/(1.5*D(i));
        CNa_fb = K_fb*CNa_f;
        X_f = (D(i) + L);
        barX_f = X_f + (m*(a + 2*b)/(3*(a + b))) + (1/6)*(a + b - (a*b/(a + b)));
        CNa_tot = CNa_n + CNa_fb;
        X_CP = (CNa_n*barX_n + CNa_fb*barX_f)/CNa_tot;
        %X_CG calculation
        X_CG = (2/3*D(i)*M_nose + 2/3*D(i)*M_L + (L + D(i))/2*M_tube + (L + 2*D(i) - 0.5*L_p)*M_p + 2/3*D(i)*M_fins) / M_o;
        error = X_CP - X_CG - D(i)*SM;
        L = L + 0.0001;
    end
    lam = M_L/(M_p + M_s);
    if (lam - lam_log < 0)
        D_max = D(i-1);
        break
    end
end

% Assembly of parameter vector
D = D_max;
L_D = L/D;
delta_D = delta/D;
eps = M_s/(M_p + M_s);
val_vector = [R W_eq t_b p0_pa delta_D D L L_D X_CG X_CP M_p M_s M_o lam eps L_p]';
end
