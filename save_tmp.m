clear;
clc;
format shortg
global a b; %bounds
a = 0; % lower integration limit
b = 1; % upper integration limit
h = 0.0001; % runge-kutta step
% below are our boundary coefficients
alpha0 = 0; % for lower limit near first derivative
alpha1 = 1; % for lower limit near zero derivative
beta0 = 0; % for upper limit near first derivative
beta1 = 1; % for upper limit near zero derivative

% str_RHS = inputdlg('Enter right-hand side');
% 
% fileID_RHS = fopen('RHS_fx_func.m','W');
% fprintf(fileID_RHS,'function out = RHS_fx_func(x)\n    global a b;\n    out = %s;\nend',str_RHS{1});
% fclose(fileID_RHS);

const_num = 0;

% str_p0 = inputdlg('Enter function p0(x)');
% check_empty = strfind(str_p0, 'x')
% if isempty(check_empty{1}) == 1
%     const_num = const_num + 1;
% end
% 
% str_p1 = inputdlg('Enter function p1(x)');
% check_empty = strfind(str_p1, 'x')
% if isempty(check_empty{1}) == 1
%     const_num = const_num + 1;
% end
% 
% str_p2 = inputdlg('Enter right-hand side');
% check_empty = strfind(str_p2, 'x')
% if isempty(check_empty{1}) == 1
%     const_num = const_num + 1;
% end

% fileID_p0 = fopen('p0.m','W');
% fprintf(fileID_p0,'function out = p0(x)\n    out = %s;\nend',str_p0{1});
% fclose(fileID_p0);
% 
% fileID_p1 = fopen('p1.m','W');
% fprintf(fileID_p1,'function out = p1(x)\n    out = %s;\nend',str_p1{1});
% fclose(fileID_p1);
% 
% fileID_p2 = fopen('p2.m','W');
% fprintf(fileID_p2,'function out = p2(x)\n    out = %s;\nend',str_p2{1});
% fclose(fileID_p2);

debug = 1
% -5.25 + (3.75 + 0.5 * x) * x | -4 | 4 | 1 | -5.25 + (3.75 + 0.5 x) x
%  2.75 + (3.75 + 0.5 * x) * x |  4 | 4 | 1 |
%  6.75 + (3.75 + 0.5 * x) * x |  8 | 4 | 1 |
%
if const_num == 3 || debug == 1
    
%     coef_0 = eval(str_p0{1});
%     coef_1 = eval(str_p1{1});
%     coef_2 = eval(str_p2{1});

    coef_0 = 1;
    coef_1 = 0;
    coef_2 = -1;
    
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
        m_real(1) = -coef_1	/ (2 * coef_0);
        m_real(2) = -coef_1  / (2 * coef_0);
        m_im(1) =    sqrt(-D)/ (2 * coef_0);
        m_im(2) =   -sqrt(-D)/ (2 * coef_0);
    end
    totiter = round(abs(b - a)/ h);
    x = zeros(1,totiter);
    u = zeros(2, totiter);
    u_prime = zeros(2, totiter);

    x(1) = a;
    for i = 1:totiter-1
        x(i+1) = x(i) + h;
    end

    
    
    u(1,1) = 0;
    u_prime(1,1) = 1;
    u(2,1) = 1;
    u_prime(2,1) = 0;
    if m_im(1) == 0 && m_im(2) == 0
        for j = 1:2
            for i=1:totiter - 1
                u(j, i+1) = exp(x(i+1)  * m_real(j) );
                u_prime(j, i+1) = u(j, i+1) * m_real(j);
%                 if j == 2
%                     error = u(j,i+1) - exp(-x(i+1)) 
%                     error2 = u_prime(j, i+1) + exp(-x(i+1))
%                 end
            end
        end
    else
        for j = 1:2
            for i=1:totiter - 1
                if j == 1
                    part = cos( x(i) * m_im(j) );
                else
                    part = sin( x(i) * m_im(j) );
                end
                u(j, i+1) = exp( x(i) * m_real(j) ) * part / 2;
                if j == 1
                    part = -sin( x(i) * m_im(j) ) * m_im(j);
                else
                    part = cos( x(i) * m_im(j) )  * m_im(j);
                end
                u_prime(j, i+1) = u(j, i+1) * m_real(j) * part;
            end
        end
    end
        
    
else
 
    [x, u, u_prime, totiter] = rungekutta(a, b, h);
end
C1 = zeros(1,totiter);
C2 = zeros(1,totiter);
P = zeros(1,totiter);
W = zeros(1,totiter);
for i=1:totiter
    W(i) = u(1,i) * u_prime(2,i) - u(2,i) * u_prime(1,i);
end

partA = beta0 * u_prime(1,totiter) + beta1 * u(1,totiter) ;
partB = beta0 * u_prime(2,totiter) + beta1 * u(2,totiter);

for i = 1:totiter - 1
    C1(i + 1) = C1(i)   - 0.5 * h * u(2,i)   * RHS_fx_func(x(i) )   / ( p0( x(i) )   * W(i) )...
                        - 0.5 * h * u(2,i+1) * RHS_fx_func(x(i+1) ) / ( p0( x(i+1) ) * W(i+1) );
    C2(i + 1) = C2(i)   + 0.5 * h * u(1,i)   * RHS_fx_func(x(i) )   / ( p0( x(i) )   * W(i) )...
                        + 0.5 * h * u(1,i+1) * RHS_fx_func(x(i+1) ) / ( p0( x(i+1) ) * W(i+1) );
%     error1 = C1(i) - 0.5 .* exp(-x(i)) .* (x.*x + x + 1)
%     error2 = C2(i) - 0.5 .* exp( x(i)) .* (x.*x - 3.*x + 3) 
%    error3 =  C1(i) - exp(-x(i)) .* (-0.5 + 0.5 * exp(x(i)) - 0.5 * x(i) - 0.5 * x(i)*x(i))
error4 =  C2(i) - 1.5 + exp(x(i)) .* (1.5 - 1.5 * x(i) + 0.5 * x(i)*x(i))
                    
%     C1(i + 1) = C1(i)   + 0.5 * h * (alpha0 * u_prime(2,i)   + alpha1 * u(2,i))   * RHS_fx_func(x(i) )   / ( p0( x(i) )   * W(i) )...
%                         + 0.5 * h * (alpha0 * u_prime(2,i+1) + alpha1 * u(2,i+1)) * RHS_fx_func(x(i+1) ) / ( p0( x(i+1) ) * W(i+1) );
%     C2(i + 1) = C2(i)   - 0.5 * h * (beta0  * u_prime(1,i)   + beta1  * u(1,i))   * RHS_fx_func(x(i) )   / ( p0( x(i) )   * W(i) )...
%                         - 0.5 * h * (beta0  * u_prime(1,i+1) + beta1  * u(1,i+1)) * RHS_fx_func(x(i+1) ) / ( p0( x(i+1) ) * W(i+1) );
    
     P(i + 1) = P(i) + 0.5 * h * ( u(2, i)   * partA - u(1, i)   * partB) *  RHS_fx_func(x(i)  ) / ( p0( x(i)   ) * W(i)   ) ...
                     + 0.5 * h * ( u(2, i+1) * partA - u(1, i+1) * partB) *  RHS_fx_func(x(i+1)) / ( p0( x(i+1) ) * W(i+1) );
end

% sum = 0;
% 
% for i=1:totiter
%     term = h * ( u(2, i)   * partA - u(1, i)   * partB) *  RHS_fx_func(x(i)  ) / ( p0( x(i)   ) * W(i)   );
%     if (i==1 || i == totiter)
%         term = 0.5*term;
%     end
%     sum = sum + term;
%     
% end
% 
% P(totiter) = sum;

% CERR = C1(42) + (sin(pi*x(42)) * sinh(x(42)) - pi*cos(pi*x(42)) * cosh(x(42)) + pi)/(1+pi*pi);
% 
% fprintf('\nERROR = %f\n',CERR);
% 
% CERR = C2(totiter) - (sin(pi*b) * cosh(b) - pi*cos(pi*b) * sinh(b))/(1+pi*pi);
% 
% fprintf('\nERROR = %f\n',CERR);
% 
% Am(1,1) = u(1,1);
% Am(1,2) = u(2,1);
% Am(2,1) = u(1, totiter);
% Am(2,2) = u(2, totiter);

Am(1,1) = alpha0 * u_prime(1,1)       + alpha1 * u(1,1);
Am(1,2) = alpha0 * u_prime(2,1)       + alpha1 * u(2,1);
Am(2,1) = beta0  * u_prime(1,totiter) + beta1  * u(1,totiter);
Am(2,2) = beta0  * u_prime(2,totiter) + beta1  * u(2,totiter);


Rm(1) = 0;%usually == 0
% Rm(2) = - C1(totiter) * u(1,totiter) - C2(totiter) * u(2,totiter);
% Rm(2) = - C1(totiter) * u(1,totiter) * partA - C2(totiter) * u(2,totiter) * partB;
Rm(2) = P(totiter);
fprintf('// %f\t%f \\\\    // H1 \\\\   __   // %f \\\\\n',  Am(1,1), Am(1,2), Rm(1));
fprintf('\\\\ %f\t%f  //    \\\\ H2 //   --   \\\\ %f //\n\n',Am(2,1), Am(2,2), Rm(2));

% H = zeros(1,2);
H = Am\Rm';

fprintf('// %f\t%f \\\\    // %f  \\\\   __   // %f \\\\\n',Am(1,1), Am(1,2), H(1), Rm(1));
fprintf('\\\\ %f\t%f  //    \\\\ %f //   --   \\\\ %f //\n',Am(2,1), Am(2,2), H(2), Rm(2));


% Hprov = (-C1(totiter)*sinh(b) - C2(totiter)*cosh(b))/sinh(b);
% disp(Hprov);
solution = zeros(1,totiter);
for i= 1:totiter
    solution(i) = ( H(1) + C1(i) ) * u(1,i) + ( H(2) + C2(i) ) * u(2,i); 
end

check = zeros(1,totiter);
for i = 1:totiter
    check(i) = 0.25 * (2 * x(i) * x(i) - x(i) - 1);
end


% hold off
% plot(x, solution, x, check );
h1_tmp = (2*exp(-1) - 2) /(exp(1) - exp(-1));
h2_tmp = (-2*exp(1) + 2) /(exp(1) - exp(-1));
plot(x, h1_tmp .* exp(x) + h2_tmp .* exp(-x) + 0.5 .* (x.*x - 3.*x + 3) + 0.5 .* (x.*x + x + 1), x, solution);
% hold on
% plot(x, check )
% hold off





