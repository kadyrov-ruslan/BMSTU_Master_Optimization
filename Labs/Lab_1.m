eps = 1e-6;
a = 0;
b = 1;
%[xRes, fRes, xi, fi, iterCount] = BitwiseSearch(a,b, eps);
%[xRes, fRes, xi, fi, iterCount] = GoldenSection(a,b, eps);
[xRes, fRes, xi, fi, iterCount] = ParabolasMethod(a,b, eps,1);
% x = 0:1e-4:1;
% y = Func(x);
% plot(x,y);
%Получение данных для построения графика целевой функции
xArr = zeros(1,iterCount);
fArr = zeros(1,iterCount);
a = 0;
step = 1;
while(step < 600)
    xArr(step) =  step*0.01 + a;
    fArr(step) = Func(step*0.01 + a);
    step = step + 1;
end

%Получение данных для построения точек, приближающихся к минимуму
xI_Arr = zeros(1,iterCount);
fI_Arr = zeros(1,iterCount);

step = 1;
for i=1:MaxIterationCount()
    if(xi(i)~=0)
        xI_Arr(step) = xi(i);
        fI_Arr(step) = fi(i);
        step = step + 1;
    end
end

fprintf('Количество вычислений целевой функции: %1.d.\n',iterCount);
fprintf('x* = %1.10f.\n',xRes);
fprintf('f(x*) = %1.10f.\n',fRes);
xRes - 0.113

plot(xArr, fArr,xRes, fRes, 'ro',xI_Arr,fI_Arr,'k-*');
grid on;
title('График целевой функции f(x)');
xlabel('x');
ylabel('Значение целевой функции f(x)');
ylim([-2 15])

%Вычисление значения целевой функции в точке х
function X = Func(x)
% firstPart = cosh((3*(x^3) + 2*(x^2) - 4*x + 5)/3);
% secondPart = tanh((x^3 - 3*sqrt(2)*x -2)/(2*x + sqrt(2)));
% X = firstPart + secondPart - 2.5;
%X = cosh((3*(x^3)+2*(x^2)-4*x+5)/3)+tanh((x^3-3*sqrt(2)*x-2)/(2*x+sqrt(2)))-2.5;
X =(cosh(x-0.113))^2;
%X = tanh((x-0.222)^4);
end

%Метод поразрядного поиска
function [xResult, fResult, xI,fI, iterationCount] = BitwiseSearch(a, b, eps)
% a - начало отрезка, b - конец отрезка, eps - эпсилон, точность поиска
% xResult - оптимальный x*, %fResult - значение целевой функции в x*
% xI-последовательность xi, приближающих точку искомого минимума
% fI-последовательность fi, приближающих точку искомого минимума
% iterationCount - число вычислений значения целевой функции

xi = zeros(1,MaxIterationCount());
fi = zeros(1,MaxIterationCount());

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

xi = zeros(1,MaxIterationCount());
fi = zeros(1,MaxIterationCount());

phi = (sqrt(5) - 1) / 2;
x2 = a + (b - a)*phi;
x1 = a + b - x2;

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
        x1 = a + b - x1;
        f1 = Func(x1);
        xi(iter) = x1;
        fi(iter) = f1;
    else
        a = x1;
        x1 = x2;
        f1 = f2;
        x2 = a + b - x2;
        f2 = Func(x2);
        xi(iter) = x2;
        fi(iter) = Func(x2);
    end
end
iterationCount = iter;
xI = xi;
fI = fi;
xResult = (x1 + x2)/2;
fResult = Func(xResult);
end

%Метод парабол
function [xResult, fResult, xI,fI, iterationCount] = ParabolasMethod(a, b, eps, iterNum)
% a - начало отрезка, b - конец отрезка, eps - эпсилон, точность поиска
% xResult - оптимальный x*, %fResult - значение целевой функции в x*
% xI-последовательность xi, приближающих точку искомого минимума
% fI-последовательность fi, приближающих точку искомого минимума
% iterationCount - число вычислений значения целевой функции
[p1, p2] = GetPointsByGoldenSection(a, b, iterNum);
xi = zeros(1,MaxIterationCount());
fi = zeros(1,MaxIterationCount());

x1 = p1;
x2 = (p1 + p2)/2;
x3 = p2;

f1 = Func(x1);
f2 = Func(x2);
f3 = Func(x3);

[x, f] = CalculateResult(x1, x2, x3, f1, f2, f3);
iter = 4;
pointsCount = 1;
while(true)
    if (x < x2)
        if (f > f2)
            x1 = x;
            f1 = f;
        else
            x3 = x2;
            f3 = f2;
            x2 = x;
            f2 = f;
        end
    else
        if (f >= f2)
            x3 = x;
            f3 = f;
        else
            x1 = x2;
            f1 = f2;
            x2 = x;
            f2 = f;
        end
    end   
    xtmp = x;
    [x, f] = CalculateResult(x1, x2, x3, f1, f2, f3);   
    iter = iter + 1;
    if (abs(x - xtmp) <= eps)
        break;
    end
end
iterationCount = iter - 4;
xI = xi;
fI = fi;
xResult = x;
fResult = f;
end

function num = MaxIterationCount
num = 250;
end

function [p1, p2] = GetPointsByGoldenSection(a, b, iterationNum)
t = (sqrt(5) - 1)/2;
l = b - a;

x1 = b - t * l;
x2 = a + t * l;
f1 = Func(x1);
f2 = Func(x2);

for i = 1:iterationNum
    if (f1 <= f2)
        b = x2;
        l = b - a;
        x2 = x1;
        f2 = f1;
        x1 = b - t * l;
        f1 = Func(x1);
    else
        a = x1;
        l = b - a;
        x1 = x2;
        f1 = f2;
        x2 = a + t * l;
        f2 = Func(x2);
    end
end

p1 = a;
p2 = b;
end

function [x, f] = CalculateResult(x1, x2, x3, f1, f2, f3)
a1 = (f2 - f1) / (x2 - x1);
a2 = ((f3 - f1) / (x3 - x1) - (f2 - f1) / (x2 - x1)) / (x3 - x2);
x = (x1 + x2 - a1 / a2) / 2;
f = Func(x);
end
