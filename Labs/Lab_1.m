eps = 0.01;
a = 0;
b = 1;
%[xRes, fRes, xi, fi, iterCount] = BitwiseSearch(a,b, eps);
%[xRes, fRes, xi, fi, iterCount] = GoldenSection(a,b, eps);
%[xRes, fRes, xi, fi, iterCount] = ParabolasMethod(a,b, eps);
[xRes, fRes, xi, fi, iterCount] = NewtonMethod(a,b, eps);

%Получение данных для построения графика целевой функции
xArr = zeros(1,iterCount);
fArr = zeros(1,iterCount);
a = -4;
step = 1;
while(step < 600)
    xArr(step) =  step*0.01 + a;
    fArr(step) = Func(step*0.01 + a);
    step = step + 1;
end

%Получение данных для построения точек, приближающихся к минимуму
xiArr = zeros(1,iterCount);
fiArr = zeros(1,iterCount);

step = 1;
for i=1:50
    if(xi(i)~=0)
        xiArr(step) = xi(i);
        fiArr(step) = fi(i);
        step = step + 1;
    end
end

fprintf('Количество вычислений целевой функции: %1.d.\n',iterCount);
fprintf('x* = %1.10f.\n',xRes);
fprintf('f(x*) = %1.10f.\n',fRes);

plot(xArr, fArr,xRes, fRes, 'ro', xiArr, fiArr, 'b*');
grid on;
title('График целевой функции f(x)');
xlabel('x');
ylabel('Значение целевой функции f(x)');
ylim([-2 15])

%Вычисление значения целевой функции в точке х
function X = Func(x)
firstPart = cosh((3*(x^3) + 2*(x^2) - 4*x + 5)/3);
secondPart = tanh((x^3 - 3*sqrt(2)*x -2)/(2*x + sqrt(2)));
X = firstPart + secondPart - 2.5;
%X = cosh((3*(x^3)+2*(x^2)-4*x+5)/3)+tanh((x^3-3*sqrt(2)*x-2)/(2*x+sqrt(2)))-2.5;
end

%Метод поразрядного поиска
function [xResult, fResult, xI,fI, iterationCount] = BitwiseSearch(a, b, eps)
% a - начало отрезка, b - конец отрезка, eps - эпсилон, точность поиска
% xResult - оптимальный x*, %fResult - значение целевой функции в x*
% xI-последовательность xi, приближающих точку искомого минимума
% fI-последовательность fi, приближающих точку искомого минимума
% iterationCount - число вычислений значения целевой функции

xi = zeros(1,50);
fi = zeros(1,50);

delta = (b - a)/4;
x0 = a;
f0 = Func(x0);
iter = 0;
while(true)
    iter = iter + 1;
    x1 = x0 + delta;
    f1 = Func(x1);
    xi(iter) = x1;
    fi(iter) = f1;
    
    if(f0 > f1)
        x0 = x1;
        f0 = f1;
        if(~(a < x0 && x0 < b))
            if(abs(delta) <= eps)
                xResult = x0;
                fResult = f0;
                break;
            else
                x0 = x1;
                f0 = f1;
                delta = - (delta/4);
            end
        end
    else
        if(abs(delta) <= eps)
            xResult = x0;
            fResult = f0;
            break;
        else
            x0 = x1;
            f0 = f1;
            delta = - (delta/4);
        end
    end
end
iterationCount = iter;
xI = xi;
fI = fi;
end

%Метод золотого сечения
function [xResult, fResult, xI,fI, iterationCount] = GoldenSection(a, b, eps)
% a - начало отрезка, b - конец отрезка, eps - эпсилон, точность поиска
% xResult - оптимальный x*, %fResult - значение целевой функции в x*
% xI-последовательность xi, приближающих точку искомого минимума
% fI-последовательность fi, приближающих точку искомого минимума
% iterationCount - число вычислений значения целевой функции

xi = zeros(1,50);
fi = zeros(1,50);
phi = (1 + sqrt(5)) / 2;

x1 = b - (b - a)/phi;
x2 = a + (b - a)/phi;
xi(1) = x1;
xi(2) = x2;

f1 = Func(x1);
f2 = Func(x2);
fi(1) = f1;
fi(2) = f2;

iter = 2;
while(abs(b - a) > eps)
    iter = iter + 1;
    if(f1 < f2)
        b = x2;
        x2 = x1;
        f2 = f1;
        x1 = a + (b - a)/phi;
        f1 = Func(x1);
    else
        a = x1;
        x1 = x2;
        f1 = f2;
        x2 = b - (b - a)/phi;
        f2 = Func(x2);
    end
    xi(iter) = (x1 + x2)/2;
    fi(iter) = Func((x1 + x2)/2);
end
iterationCount = iter;
xI = xi;
fI = fi;
xResult = (x1 + x2)/2;
fResult = Func(xResult);
end

%Метод парабол
function [xResult, fResult, xI,fI, iterationCount] = ParabolasMethod(a, b, eps)
% a - начало отрезка, b - конец отрезка, eps - эпсилон, точность поиска
% xResult - оптимальный x*, %fResult - значение целевой функции в x*
% xI-последовательность xi, приближающих точку искомого минимума
% fI-последовательность fi, приближающих точку искомого минимума
% iterationCount - число вычислений значения целевой функции

xi = zeros(1,50);
fi = zeros(1,50);

xCenter = (b - a)/2;
step = 0.001*xCenter;
xMin = 0;
iter = 0;
while(true)
    iter = iter + 1;
    xLeft = xCenter - step;
    xRight = xCenter + step;
    xMin = 0.5*(Func(xLeft)*(xRight + xCenter) - 2*(Func(xCenter)*(xRight + xLeft)) + Func(xRight)*(xCenter + xLeft))/(Func(xLeft) - 2*Func(xCenter) + Func(xRight));
    if(abs(xMin - xCenter) < eps)
        break;
    end
    xCenter = xMin;
    xi(iter) = xMin;
    fi(iter) = Func(xMin);
end
iterationCount = iter - 1;
xI = xi;
fI = fi;
xResult = xMin;
fResult = Func(xResult);
end

%Метод Ньютона
function [xResult, fResult, xI,fI, iterationCount] = NewtonMethod(a, b, eps)
% a - начало отрезка, b - конец отрезка, eps - эпсилон, точность поиска
% xResult - оптимальный x*, %fResult - значение целевой функции в x*
% xI-последовательность xi, приближающих точку искомого минимума
% fI-последовательность fi, приближающих точку искомого минимума
% iterationCount - число вычислений значения целевой функции

xi = zeros(1,50);
fi = zeros(1,50);

x0 = a;
x1 = 0;
dx = 0;
iter = 0;
while(dx < eps)
    iter = iter + 1;
    res1 = DX(x0);
    res2 = D2X(x0);
    
    x1 = x0 - res1/res2;
    dx = abs(x1 - x0);
    x0 = x1;
    
    difff = D2X(x1);
    if( difff < 0)
        break;
    end
    xi(iter) = x1;
    fi(iter) = Func(x1);
end
iterationCount = iter;
xI = xi;
fI = fi;
xResult = x1;
fResult = Func(xResult);
end

function value = DX(pointX)
syms x;
f = cosh((3*(x^3)+2*(x^2)-4*x+5)/3)+tanh((x^3-3*sqrt(2)*x-2)/(2*x+sqrt(2)))-2.5;
F = diff(f);
str1 = char(F);
x = pointX;
value=eval(str1);
end

function value = D2X(pointX)
syms x;
f = cosh((3*(x^3)+2*(x^2)-4*x+5)/3)+tanh((x^3-3*sqrt(2)*x-2)/(2*x+sqrt(2)))-2.5;
F = diff(f,2);
str1 = char(F);
x = pointX;
value=eval(str1);
end
