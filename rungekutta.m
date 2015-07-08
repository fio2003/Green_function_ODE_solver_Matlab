function [x, u, v, totiter]=rungekutta(a, b, h)

totiter = round(abs(b - a)/ h);
% h = abs(b - a)/(totiter);

x = zeros(1,totiter);
u = zeros(2, totiter);
v = zeros(2, totiter);

x(1) = a;
for i = 1:totiter-1
    x(i+1) = x(i) + h;
end

u(1,1) = 0;
v(1,1) = 1;
u(2,1) = 1;
v(2,1) = 0;

for j = 1:2
    for i=1:totiter - 1
        k1 = h * f1(x(i), u(j,i) , v(j,i));
        l1 = h * f2(x(i), u(j,i) , v(j,i));
        k2 = h*f1(x(i) + h/2, u(j,i) + k1/2 , v(j,i) + l1/2);
        l2 = h*f2(x(i) + h/2, u(j,i) + k1/2 , v(j,i) + l1/2);
        k3 = h*f1(x(i) + h/2, u(j,i) + k2/2 , v(j,i) + l2/2);
        l3 = h*f2(x(i) + h/2, u(j,i) + k2/2 , v(j,i) + l2/2);
        k4 = h*f1(x(i) + h, u(j,i) + k3, v(j,i) + l3);
        l4 = h*f2(x(i) + h, u(j,i) + k3, v(j,i) + l3);
        u(j, i+1) = u(j, i) + (k1 + 2*k2 + 2*k3 + k4)/6;
        v(j, i+1) = v(j, i) + (l1 + 2*l2 + 2*l3 + l4)/6;
%         if j == 1
%             fprintf('Error = %8.10f\n', u(j,i) - (exp(x(i) ) - exp(-x(i) ) )/2);
%         else
%             fprintf('Error = %8.10f\n', u(j,i) - (exp(x(i) ) + exp(-x(i) ) )/2);
%         end
    end
end


end


function val_out = f1(x, u, v)
    val_out = v;
end


function val_out = f2(x, u, v)
    val_out = (-p1(x)*v - p2(x)*u )/p0(x);
end
