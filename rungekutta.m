function [x, u, v, w, z, totiter]=rungekutta(a, b, h)

totiter = round(abs(b - a)/ h);
% h = abs(b - a)/(totiter);

x = zeros(1,totiter);
u = zeros(4, totiter);
v = zeros(4, totiter);%in our case v is u_prime
w = zeros(4, totiter);%in our case v is v_prime (u prime prime)
z = zeros(4, totiter);%in our case v is w_prime (u prime prime prime )

x(1) = a;
for i = 1:totiter-1
    x(i+1) = x(i) + h;
end

u(1,1) = 1;
v(1,1) = 0;
w(1,1) = 0;
z(1,1) = 0;

u(2,1) = 0;
v(2,1) = 1;
w(2,1) = 0;
z(2,1) = 0;

u(3,1) = 0;
v(3,1) = 0;
w(3,1) = 1;
z(3,1) = 0;

u(4,1) = 0;
v(4,1) = 0;
w(4,1) = 0;
z(4,1) = 1;

for j = 1:4
    for i=1:totiter - 1
        var_u_1 = h * f_u_prime(v(j,i));
        var_v_1 = h * f_v_prime(w(j,i));
        var_w_1 = h * f_w_prime(z(j,i));
        var_z_1 = h * f_z_prime(x(i), u(j,i) , v(j,i), w(j,i), z(j,i));
        
        
        var_u_2 = h * f_u_prime(v(j,i) + var_u_1/2);
        var_v_2 = h * f_v_prime(w(j,i) + var_v_1/2);
        var_w_2 = h * f_w_prime(z(j,i) + var_w_1/2);
        var_z_2 = h * f_z_prime(x(i) + h/2, u(j,i) + var_u_1/2, v(j,i) + var_v_1/2, w(j,i) + var_w_1/2, z(j,i) + var_z_1/2 );
        
        var_u_3 = h * f_u_prime(v(j,i) + var_u_2/2);
        var_v_3 = h * f_v_prime(w(j,i) + var_v_2/2);
        var_w_3 = h * f_w_prime(z(j,i) + var_w_2/2);
        var_z_3 = h * f_z_prime(x(i) + h/2, u(j,i) + var_u_2/2, v(j,i) + var_v_2/2, w(j,i) + var_w_2/2, z(j,i) + var_z_2/2 );

        var_u_4 = h * f_u_prime(v(j,i) + var_u_3);
        var_v_4 = h * f_v_prime(w(j,i) + var_v_3);
        var_w_4 = h * f_w_prime(z(j,i) + var_w_3);
        var_z_4 = h * f_z_prime(x(i) + h, u(j,i) + var_u_3, v(j,i) + var_v_3, w(j,i) + var_w_3, z(j,i) + var_z_3 );

          
        u(j, i+1) = u(j, i) + (var_u_1 + 2*var_u_2 + 2*var_u_3 + var_u_4)/6;
        v(j, i+1) = v(j, i) + (var_v_1 + 2*var_v_2 + 2*var_v_3 + var_v_4)/6;
        w(j, i+1) = w(j, i) + (var_w_1 + 2*var_w_2 + 2*var_w_3 + var_w_4)/6;
        z(j, i+1) = z(j, i) + (var_z_1 + 2*var_z_2 + 2*var_z_3 + var_z_4)/6;
            
        
%         if j == 1
%             fprintf('Error = %8.10f\n', u(j,i) - (exp(x(i) ) - exp(-x(i) ) )/2);
%         else
%             fprintf('Error = %8.10f\n', u(j,i) - (exp(x(i) ) + exp(-x(i) ) )/2);
%         end
    end
end


end


function val_out = f_u_prime(v)
    val_out = v;
end

function val_out = f_v_prime(w)
    val_out = w;
end

function val_out = f_w_prime(z)
    val_out = z;
end


function val_out = f_z_prime(x, u, v, w, z)
    val_out = (-p1(x)*z - p2(x)*w - p3(x)*v - p4(x)*u )/p0(x);
end
