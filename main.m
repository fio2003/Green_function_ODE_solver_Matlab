clear;
clc;
format shortg
global a b; %bounds
a = 0; % lower integration limit
b = 1; % upper integration limit
h = 0.00001; % runge-kutta step
% below are our boundary coefficients
alpha0 = 1; % for lower limit near first derivative
alpha1 = -1; % for lower limit near zero derivative
beta0 = 0; % for upper limit near first derivative
beta1 = 1; % for upper limit near zero derivative

str_RHS = inputdlg('Enter right-hand side');

fileID_RHS = fopen('RHS_fx_func.m','W');
fprintf(fileID_RHS,'function out = RHS_fx_func(x)\n    global a b;\n    out = %s;\nend',str_RHS{1});
fclose(fileID_RHS);

const_num = 0;

str_p0 = inputdlg('Enter function p0(x)');
check_empty = strfind(str_p0, 'x');
if isempty(check_empty{1}) == 1
    const_num = const_num + 1;
end

str_p1 = inputdlg('Enter function p1(x)');
check_empty = strfind(str_p1, 'x');
if isempty(check_empty{1}) == 1
    const_num = const_num + 1;
end

str_p2 = inputdlg('Enter function p2(x)');
check_empty = strfind(str_p2, 'x');
if isempty(check_empty{1}) == 1
    const_num = const_num + 1;
end

fileID_p0 = fopen('p0.m','W');
fprintf(fileID_p0,'function out = p0(x)\n    out = %s;\nend',str_p0{1});
fclose(fileID_p0);

fileID_p1 = fopen('p1.m','W');
fprintf(fileID_p1,'function out = p1(x)\n    out = %s;\nend',str_p1{1});
fclose(fileID_p1);

fileID_p2 = fopen('p2.m','W');
fprintf(fileID_p2,'function out = p2(x)\n    out = %s;\nend',str_p2{1});
fclose(fileID_p2);

debug = 10;
% -5.25 + (3.75 + 0.5 * x) * x | -4 | 4 | 1 | -5.25 + (3.75 + 0.5 x) x
%  2.75 + (3.75 + 0.5 * x) * x |  4 | 4 | 1 |
%  6.75 + (3.75 + 0.5 * x) * x |  8 | 4 | 1 |
tic
if const_num == 3 && debug == 1
    fprintf('You entered three constant arguments. Using fast path\n');
    coef_0 = eval(str_p0{1});
    coef_1 = eval(str_p1{1});
    coef_2 = eval(str_p2{1});

    D = coef_1 * coef_1 - 4 * coef_0 * coef_2;
    if D > 0
        m_real(1) = ( -coef_1 + sqrt(D) ) / (2 * coef_0);
        m_real(2) = ( -coef_1 - sqrt(D) ) / (2 * coef_0);
        m_im(1) = 0;
        m_im(2) = 0;
    elseif D == 0
        m_real(1) = -coef_1	/ (2 * coef_0);
        m_real(2) = m_real(1);
        m_im(1) = 0;
        m_im(2) = 0;
    else
        alpha = -coef_1	/ (2 * coef_0);
        beta = sqrt(-D)/ (2 * coef_0);
        m_real(1) = -coef_1	/ (2 * coef_0);
        m_real(2) = -coef_1  / (2 * coef_0);
        m_im(1) =    sqrt(-D)/ (2 * coef_0);
        m_im(2) =   sqrt(-D)/ (2 * coef_0);
    end
    totiter = round(abs(b - a)/ h);
    x = zeros(1,totiter);
    u = zeros(2, totiter);
    u_prime = zeros(2, totiter);

    x(1) = a;
    for i = 1:totiter-1
        x(i+1) = x(i) + h;
    end

    if m_real(1) ~= m_real(2)
        for j = 1:2
            for i=1:totiter
                u(j, i) = exp(x(i)  * m_real(j) );
                u_prime(j, i) = u(j, i) * m_real(j);
            end
        end
    elseif m_im(1) == 0
        for i=1:totiter
            u(1,i) = exp(m_real(1)*x(i));
            u_prime(1,i) = m_real(1)*exp(m_real(1)*x(i));
            u(2,i) = x(i)*exp(m_real(1)*x(i));
            u_prime(2,i) = exp(m_real(1)*x(i)) + m_real(1)*x(i)*exp(m_real(1)*x(i));
        end
    else
        for i=1:totiter
            u(1,i) = exp(alpha*x(i))*cos(beta*x(i));
            u_prime(1,i) = alpha*exp(alpha*x(i))*cos(beta*x(i)) - beta*exp(alpha*x(i))*sin(beta*x(i));
            u(2,i) = exp(alpha*x(i))*sin(beta*x(i));
            u_prime(2,i) = alpha*exp(alpha*x(i))*sin(beta*x(i)) + beta*exp(alpha*x(i))*cos(beta*x(i));
        end
%         disp('HERE!');
    end
        
    
else
 
    [x, u, u_s_prime, u_d_prime, u_t_prime, totiter] = rungekutta(a, b, h);
end

C1 = zeros(1,totiter);
C2 = zeros(1,totiter);
C3 = zeros(1,totiter);
C4 = zeros(1,totiter);

W = zeros(1,totiter);
for i=1:totiter
    
    
    W(i) = u(2,i) * u_s_prime(4,i) * u_d_prime(3,i) * u_t_prime(1,i) ...
         - u(2,i) * u_s_prime(3,i) * u_d_prime(4,i) * u_t_prime(1,i) ...
         - u(1,i) * u_s_prime(4,i) * u_d_prime(3,i) * u_t_prime(2,i) ...
         + u(1,i) * u_s_prime(3,i) * u_d_prime(4,i) * u_t_prime(2,i) ...
         - u(2,i) * u_s_prime(4,i) * u_d_prime(1,i) * u_t_prime(3,i) ...
         + u(1,i) * u_s_prime(4,i) * u_d_prime(2,i) * u_t_prime(3,i) ...
         + u(2,i) * u_s_prime(1,i) * u_d_prime(4,i) * u_t_prime(3,i) ...
         - u(1,i) * u_s_prime(2,i) * u_d_prime(4,i) * u_t_prime(3,i) ...
         + u(4,i) * ( ...
                      u_d_prime(3,i) * ( - u_s_prime(2,i) * u_t_prime(1,i) + u_s_prime(1,i) * u_t_prime(2,i) ) ...
                    + u_s_prime(3,i) * ( u_d_prime(2,i) * u_t_prime(1,i) - u_d_prime(1,i) * u_t_prime(2,i) ) ...
                    + u_t_prime(3,i) * ( u_s_prime(2,i) * u_d_prime(1,i) - u_s_prime(1,i) * u_d_prime(2,i) ) ...
                    ) ...
         + u_t_prime(4,i)*( ...
                      u_s_prime(3,i) * ( u(2,i) * u_d_prime(1,i) - u(1,i) * u_d_prime(2,i) ) ...
                    + u_d_prime(3,i) * (-u(2,i) * u_s_prime(1,i) + u(1,i) * u_s_prime(2,i) )...
                    ) ...
         + u(3,i) * ( ...
                      u_d_prime(4,i) * ( u_s_prime(2,i) * u_t_prime(1,i) - u_s_prime(1,i) * u_t_prime(2,i) ) ...
                    + u_s_prime(4,i) * (-u_d_prime(2,i) * u_t_prime(1,i) + u_d_prime(1,i) * u_t_prime(2,i) ) ...
                    + u_t_prime(4,i) * (-u_s_prime(2,i) * u_d_prime(1,i) + u_s_prime(1,i) * u_d_prime(2,i) ) ...
                    );
end



for i = 1:totiter - 1
    
    C1(i+1) = C1(i) + 0.5*h*RHS_fx_func(x(i) )/( p0( x(i) )   * W(i) )...
        *(u_s_prime(4,i) * (u(3,i)*u_d_prime(2,i) - u(2,i)*u_d_prime(3,i) )...
        + u(4,i)*(-u_s_prime(3,i)*u_d_prime(2,i) + u_s_prime(2,i)*u_d_prime(3,i) )...
        + u_d_prime(4,i)*(-u(3,i)*u_s_prime(2,i) + u(2,i)*u_s_prime(3,i) ) )...
        ...
        + 0.5*h*RHS_fx_func(x(i+1) )/( p0( x(i+1) )   * W(i+1) )...
        *(u_s_prime(4,i+1) * (u(3,i+1)*u_d_prime(2,i+1) - u(2,i+1)*u_d_prime(3,i+1) )...
        + u(4,i+1)*(-u_s_prime(3,i+1)*u_d_prime(2,i+1) + u_s_prime(2,i+1)*u_d_prime(3,i+1) )...
        + u_d_prime(4,i+1)*(-u(3,i+1)*u_s_prime(2,i+1) + u(2,i+1)*u_s_prime(3,i+1) ) );

    
    C2(i+1) = C2(i) + 0.5*h*RHS_fx_func(x(i) )/( p0( x(i) )   * W(i) )...
        *(u_s_prime(4,i) * (-u(3,i)*u_d_prime(1,i) + u(1,i)*u_d_prime(3,i) )...
        + u(4,i) * (u_s_prime(3,i)*u_d_prime(1,i) - u_s_prime(1,i)*u_d_prime(3,i) )...
        + u_d_prime(4,i) * (u(3,i)*u_s_prime(1,i) - u(1,i)*u_s_prime(3,i) ) )...
        ...
        + 0.5*h*RHS_fx_func(x(i+1) )/( p0( x(i+1) )   * W(i+1) )...
        *(u_s_prime(4,i+1) * (-u(3,i+1)*u_d_prime(1,i+1) + u(1,i+1)*u_d_prime(3,i+1) )...
        + u(4,i+1) * (u_s_prime(3,i+1)*u_d_prime(1,i+1) - u_s_prime(1,i+1)*u_d_prime(3,i+1) )...
        + u_d_prime(4,i+1) * (u(3,i+1)*u_s_prime(1,i+1) - u(1,i+1)*u_s_prime(3,i+1) ) );
        
    C3(i+1) = C3(i) + 0.5*h*RHS_fx_func(x(i) )/( p0( x(i) )   * W(i) )...
        *(u_s_prime(4,i) * (u(2,i)*u_d_prime(1,i) - u(1,i)*u_d_prime(2,i) )...
        + u(4,i) * (-u_s_prime(2,i)*u_d_prime(1,i) + u_s_prime(1,i)*u_d_prime(2,i) )...
        + u_d_prime(4,i) * (-u(2,i)*u_s_prime(1,i) + u(1,i)*u_s_prime(2,i) ) )...
        ...
        + 0.5*h*RHS_fx_func(x(i+1) )/( p0( x(i+1) )   * W(i+1) )...
        *(u_s_prime(4,i+1) * (u(2,i+1)*u_d_prime(1,i+1) - u(1,i+1)*u_d_prime(2,i+1) )...
        + u(4,i+1) * (-u_s_prime(2,i+1)*u_d_prime(1,i+1) + u_s_prime(1,i+1)*u_d_prime(2,i+1) )...
        + u_d_prime(4,i+1) * (-u(2,i+1)*u_s_prime(1,i+1) + u(1,i+1)*u_s_prime(2,i+1) ) );
    
    C4(i+1) = C4(i) + 0.5*h*RHS_fx_func(x(i) )/( p0( x(i) )   * W(i) )...
        *(u_s_prime(3,i) * (-u(2,i)*u_d_prime(1,i) + u(1,i)*u_d_prime(2,i) )...
        + u(3,i)*(u_s_prime(2,i)*u_d_prime(1,i) - u_s_prime(1,i)*u_d_prime(2,i) )...
        + u_d_prime(3,i)*(u(2,i)*u_s_prime(1,i) - u(1,i)*u_s_prime(2,i) ) )...
        ...
        + 0.5*h*RHS_fx_func(x(i+1) )          /( p0( x(i+1) ) * W(i+1) )...
        * (u_s_prime(3,i+1) * (-u(2,i+1)*u_d_prime(1,i+1) + u(1,i+1)*u_d_prime(2,i+1) )...
        + u(3,i+1)*(u_s_prime(2,i+1)*u_d_prime(1,i+1) - u_s_prime(1,i+1)*u_d_prime(2,i+1) )...
        + u_d_prime(3,i+1)*(u(2,i+1)*u_s_prime(1,i+1) - u(1,i+1)*u_s_prime(2,i+1) ) );    
end

P1 = C1(totiter)*(delta1*u(1,totiter) + delta2*u_s_prime(1,totiter) + delta3*u_d_prime(1,totiter) + delta4*u_t_prime(1,totiter))+...
    C2(totiter)*(delta1*u(2,totiter) + delta2*u_s_prime(2,totiter) + delta3*u_d_prime(2,totiter) + delta4*u_t_prime(2,totiter))+...
    C3(totiter)*(delta1*u(3,totiter) + delta2*u_s_prime(3,totiter) + delta3*u_d_prime(3,totiter) + delta4*u_t_prime(3,totiter))+...
    C4(totiter)*(delta1*u(4,totiter) + delta2*u_s_prime(4,totiter) + delta3*u_d_prime(4,totiter) + delta4*u_t_prime(4,totiter));

P2 = C1(totiter)*(gamma1*u(1,totiter) + gamma2*u_s_prime(1,totiter) + gamma3*u_d_prime(1,totiter) + gamma4*u_t_prime(1,totiter))+...
    C2(totiter)*(gamma1*u(2,totiter) + gamma2*u_s_prime(2,totiter) + gamma3*u_d_prime(2,totiter) + gamma4*u_t_prime(2,totiter))+...
    C3(totiter)*(gamma1*u(3,totiter) + gamma2*u_s_prime(3,totiter) + gamma3*u_d_prime(3,totiter) + gamma4*u_t_prime(3,totiter))+...
    C4(totiter)*(gamma1*u(4,totiter) + gamma2*u_s_prime(4,totiter) + gamma3*u_d_prime(4,totiter) + gamma4*u_t_prime(4,totiter));

fprintf('C1 = %5.10f\n',C1(1,totiter));
fprintf('C2 = %5.10f\n',C2(1,totiter));
fprintf('P = %5.10f\n',P(1,totiter));


Am(1,1) = alpha4 * u(1,1) + alpha3 * u_s_prime(1,1) + alpha2 * u_d_prime(1,1) + alpha1 * u_t_prime(1,1);
Am(1,2) = alpha4 * u(2,1) + alpha3 * u_s_prime(2,1) + alpha2 * u_d_prime(2,1) + alpha1 * u_t_prime(2,1);
Am(1,3) = alpha4 * u(3,1) + alpha3 * u_s_prime(3,1) + alpha2 * u_d_prime(3,1) + alpha1 * u_t_prime(3,1);
Am(1,4) = alpha4 * u(4,1) + alpha3 * u_s_prime(4,1) + alpha2 * u_d_prime(4,1) + alpha1 * u_t_prime(4,1);

Am(2,1) = beta4 * u(1,1) + beta3 * u_s_prime(1,1) + beta2 * u_d_prime(1,1) + beta1 * u_t_prime(1,1);
Am(2,2) = beta4 * u(2,1) + beta3 * u_s_prime(2,1) + beta2 * u_d_prime(2,1) + beta1 * u_t_prime(2,1);
Am(2,3) = beta4 * u(3,1) + beta3 * u_s_prime(3,1) + beta2 * u_d_prime(3,1) + beta1 * u_t_prime(3,1);
Am(2,4) = beta4 * u(4,1) + beta3 * u_s_prime(4,1) + beta2 * u_d_prime(4,1) + beta1 * u_t_prime(4,1);

Am(3,1) = gamma4 * u(1,totiter) + gamma3 * u_s_prime(1,totiter) + gamma2 * u_d_prime(1,totiter) + gamma1 * u_t_prime(1,totiter);
Am(3,2) = gamma4 * u(2,totiter) + gamma3 * u_s_prime(2,totiter) + gamma2 * u_d_prime(2,totiter) + gamma1 * u_t_prime(2,totiter);
Am(3,3) = gamma4 * u(3,totiter) + gamma3 * u_s_prime(3,totiter) + gamma2 * u_d_prime(3,totiter) + gamma1 * u_t_prime(3,totiter);
Am(3,4) = gamma4 * u(4,totiter) + gamma3 * u_s_prime(4,totiter) + gamma2 * u_d_prime(4,totiter) + gamma1 * u_t_prime(4,totiter);

Am(4,1) = delta4 * u(1,totiter) + delta3 * u_s_prime(1,totiter) + delta2 * u_d_prime(1,totiter) + delta1 * u_t_prime(1,totiter);
Am(4,2) = delta4 * u(2,totiter) + delta3 * u_s_prime(2,totiter) + delta2 * u_d_prime(2,totiter) + delta1 * u_t_prime(2,totiter);
Am(4,3) = delta4 * u(3,totiter) + delta3 * u_s_prime(3,totiter) + delta2 * u_d_prime(3,totiter) + delta1 * u_t_prime(3,totiter);
Am(4,4) = delta4 * u(4,totiter) + delta3 * u_s_prime(4,totiter) + delta2 * u_d_prime(4,totiter) + delta1 * u_t_prime(4,totiter);


Rm(1) = 0;%usually == 0, since our integral part cancell itself because of the integration limits(from a to a)
Rm(2) = 0;%usually == 0, since our integral part cancell itself because of the integration limits(from a to a)

Rm(3) = P1;
Rm(4) = P2;
% fprintf('// %f\t%f \\\\    // H1 \\\\   __   // %f \\\\\n',  Am(1,1), Am(1,2), Rm(1));
% fprintf('\\\\ %f\t%f  //    \\\\ H2 //   --   \\\\ %f //\n\n',Am(2,1), Am(2,2), Rm(2));

% H = zeros(1,2);
H = Am\Rm';

% fprintf('// %f\t%f \\\\    // %f  \\\\   __   // %f \\\\\n',Am(1,1), Am(1,2), H(1), Rm(1));
% fprintf('\\\\ %f\t%f  //    \\\\ %f //   --   \\\\ %f //\n',Am(2,1), Am(2,2), H(2), Rm(2));


% Hprov = (-C1(totiter)*sinh(b) - C2(totiter)*cosh(b))/sinh(b);
% disp(Hprov);
solution = zeros(1,totiter);
for i= 1:totiter
    solution(i) = ( H(1) + C1(i) ) * u(1,i) + ( H(2) + C2(i) ) * u(2,i) + ( H(3) + C3(i) ) * u(3,i) + ( H(4) + C4(i) ) * u(4,i); 
end

toc

check = zeros(1,totiter);
for i = 1:totiter
    check(i) = 0.25 * (2 * x(i) * x(i) - x(i) - 1);
end

toterror = 0;

for i = 1:totiter
    toterror = toterror + abs(check(i) - solution(i));
end



% hold off
plot(x, solution, x, check );
% some_ptmp  = 0.5 - 1.5 * exp(-1) + 1.5 - 0.5*exp(1) ;
% P(totiter) - some_ptmp
% h1_tmp = -some_ptmp /(exp(1) - exp(-1));
% h2_tmp = some_ptmp /(exp(1) - exp(-1));
% plot(x, h1_tmp .* exp(x) + h2_tmp .* exp(-x) + 0.5 .* (x.*x - 3.*x + 3) + 0.5 .* (x.*x + x + 1 - exp(x)) - 1.5 * exp(-x), x, solution);
% plot(x, h1_tmp .* exp(x) + h2_tmp .* exp(-x) + 0.5 .* (x.*x - 3.*x + 3) + 0.5 .* (x.*x + x + 1), x, solution);

plot(x, solution - check );


% hold on
% plot(x, check )
% hold off





