clc
clear all
clear figure
hold on

%%%%%%%%%  Note's: Beta is the Angle of oblique Shock 



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%  Sector 1 : Upstream Flow information  %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

R = 287;
gamma = 1.405;
C_p = (gamma*R)/(gamma-1);
M_inf = 2.5;
Temp = 27 + 273.15;
a_inf = sqrt(gamma*R*Temp);
V_inf = M_inf * a_inf;
P_inf = 101320;
V_max = sqrt(2 * C_p * Temp * (1 + (gamma-1)/2 * M_inf^2)) ;
miu = asin(1/M_inf) ;
miu2 = miu *(180/pi);
%note that First guess of Beta_first_guess must be greater than miu2
tetha_cone = 20;                                                  
tetha_cone_radian = tetha_cone *(pi/180);
Beta_first_guess = 60 ;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%  Sector 2 : Solving Taylor-Mccoll for each gueesed Beta_first_guess %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


Beta = Beta_first_guess * (pi/180) ; 
counter1 = 0;

best_err = inf;
for ANGLE_iterator=Beta:-0.001:miu

    %Caluculate Flow parameter right after oblique shock
    ebarat = 2*cot(ANGLE_iterator) * ( ( (M_inf^2 * (sin(ANGLE_iterator))^2) - 1 ) / ( (M_inf^2*(gamma + cos(2*ANGLE_iterator))) + 2 ) );   
    tetha_shock = atan(ebarat);
    M_n1 = M_inf * sin(ANGLE_iterator);
    M_n2 = sqrt( ( M_n1^2 + (2/(gamma-1)) ) / ( ((2*gamma/(gamma-1)) * M_n1^2) -1 ) );   
    M_2 = M_n2 / sin(ANGLE_iterator - tetha_shock);
    a_2_bibod = sqrt((1 + (gamma-1)/2 * M_inf^2) / (1 + (gamma-1)/2 * M_2^2) * (a_inf^2 / V_max^2));
    

    % nondimminssionalizing the Velocitie's before Solving Taylor-Mccoll eq
    V_2_bibod  = M_2 * a_2_bibod;
    u_2_bibod = V_2_bibod * cos(ANGLE_iterator - tetha_shock);
    v_2_bibod = -V_2_bibod * sin(ANGLE_iterator - tetha_shock);
    % Note that this nondimminssionalizing has been done with V_max

    % Assigning initial Condition's and Calculating domain
    tspan = [ANGLE_iterator 0.00001];
    y0 = [u_2_bibod ; v_2_bibod];
    options = odeset('Events',@coneEvent);  
    [t,y,te,ye] = ode45(@(t,y) Taylor(t,y,gamma),tspan,y0,options);
    if isempty(te)
        continue   
    end

    Sol_tetha_radian = te;
    err = abs(Sol_tetha_radian - tetha_cone_radian);

    if err < best_err
        best_err = err;
        Beta_best = ANGLE_iterator;
        t_best = t;
        y_best = y;
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%  Sector 3 : Plotting Result of Calculation's  %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Beta_shock = Beta_best*180/pi;
theta_plot = t_best*180/pi;

%Claculating flow parameter on the Surface of the Cone
u_surface = y_best(:,1);
v_surface = y_best(:,2);

V_surface = sqrt(u_surface.^2 + v_surface.^2);
a_local = sqrt((gamma-1)/2 .* (1 - V_surface.^2));
Mach = V_surface ./ a_local;

P_ratio = ((1 + (gamma-1)/2 * M_inf^2) ./ ...
           (1 + (gamma-1)/2 * Mach.^2)).^(gamma/(gamma-1));
p0inf_over_pinf = (1 + (gamma-1)/2*M_inf^2)^(gamma/(gamma-1));

P_over_p0inf = P_ratio / p0inf_over_pinf;

figure

% Plotting Mach and ratio of locall Static Presure to totall presuare of
% upstream flow , between shock and cone
subplot(1,2,1)
plot(theta_plot,Mach,'LineWidth',2)
set(gca,'XDir','reverse')
xlabel('\theta (deg)')
xlim([min(theta_plot) max(theta_plot)])
ylabel('Mach')
title('Mach distribution')
grid on

info_txt = sprintf('M_\\infty = %.2f   |   \\beta = %.2f^\\circ   |   \\theta_c = %.2f^\\circ',M_inf, Beta_shock, tetha_cone);
annotation('textbox',[0.2 0.93 0.6 0.05], 'String',info_txt, 'HorizontalAlignment','center', 'EdgeColor','none', 'FontSize',11, 'FontWeight','bold');

subplot(1,2,2)
plot(theta_plot,P_over_p0inf,'LineWidth',2)
set(gca,'XDir','reverse')
xlabel('\theta (deg)')
ylabel('p / p_{0\infty}')
title('Pressure ratio')
grid on




u_cone_over_c = u_surface(end);              % u/c on Cone Surface
a2_cone_over_c2 = a_local(end)^2;            % a^2/c^2 on Cone Surface
M_cone = Mach(end);                          % Mach on Cone Surface
p_cone_p0inf = P_over_p0inf(end);            % p/p∞ on Cone Surface

% Note that c is the V_max and it's the parameter that we used to
% non-dimminsionalized the Velociti's with


fprintf('\n=== Flow Parameter on Cone Surface ===\n');
fprintf('u/c        = %.10f\n',u_cone_over_c);
fprintf('a^2/c^2    = %.10f\n',a2_cone_over_c2);
fprintf('M_cone     = %.10f\n',M_cone);
fprintf('p_cone/p0∞  = %.10f\n', p_cone_p0inf);
fprintf('\n=== Angles that we wanted ===\n');
fprintf('tetha_cone  = %.10f\n', tetha_cone);
fprintf('Beta_shock  = %.10f\n', Beta_shock);



function [value,isterminal,direction] = coneEvent(~,y)

value = y(2) + 1e-4;    
isterminal = 1;         
direction  = +1;      

end

