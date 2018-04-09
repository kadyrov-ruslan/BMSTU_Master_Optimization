eps = 0.000001;
a = 0;
b = 1;
FuncString = 'cosh((3*(x^3) + 2*(x^2) - 4*x + 5)/3) + tanh((x^3 - 3*sqrt(2)*x -2)/(2*x + sqrt(2))) - 2.5;'

[xRes, fRes, xi, fi, iterCount] = NewtonMethod(FuncString,a,b, eps);

%��������� ������ ��� ���������� ������� ������� �������
xArr = zeros(1,iterCount);
fArr = zeros(1,iterCount);
a = 0;
step = 1;
while(step < 600)
    xArr(step) =  step*0.01 + a;
    fArr(step) = Func(step*0.01 + a);
    step = step + 1;
end

%��������� ������ ��� ���������� �����, �������������� � ��������
xI_Arr = zeros(1,iterCount);
fI_Arr = zeros(1,iterCount);

step = 1;
if(iterCount > 0)
    for i=1:iterCount
        if(xi(i)~=0)
            xI_Arr(step) = xi(i);
            fI_Arr(step) = fi(i);
            step = step + 1;
        end
    end
end
fprintf('������� � ������� ������ ������� \n');
fprintf('���������� ���������� ������� �������: %1.d.\n',iterCount);
fprintf('x* = %1.10f.\n',xRes);
fprintf('f(x*) = %1.10f.\n',fRes);

[xRes, fRes, xi, fi, iterCount] = StandartFunc(FuncString,a,b);
fprintf('������� � ������� ����������� ������� fminbnd \n');
fprintf('���������� ���������� ������� �������: %1.d.\n',iterCount);
fprintf('x* = %1.10f.\n',xRes);
fprintf('f(x*) = %1.10f.\n',fRes);

plot(xArr, fArr,xRes, fRes, 'ro',xI_Arr,fI_Arr,'k*');
grid on;
title('������ ������� ������� f(x)');
xlabel('x');
ylabel('�������� ������� ������� f(x)');
ylim([-2 15])

%���������� �������� ������� ������� � ����� �
function X = Func(x)
% firstPart = cosh((3*(x^3) + 2*(x^2) - 4*x + 5)/3);
% secondPart = tanh((x^3 - 3*sqrt(2)*x -2)/(2*x + sqrt(2)));
% X = firstPart + secondPart - 2.5;
X = cosh((3*(x^3)+2*(x^2)-4*x+5)/3)+tanh((x^3-3*sqrt(2)*x-2)/(2*x+sqrt(2)))-2.5;
%X = (cosh(x-0.111))^2;
end

%����� �������
function [xResult, fResult, xI,fI, iterationCount] = NewtonMethod(FunctionString, a, b, eps)
% a - ������ �������, b - ����� �������, eps - �������, �������� ������
% xResult - ����������� x*, %fResult - �������� ������� ������� � x*
% xI-������������������ xi, ������������ ����� �������� ��������
% fI-������������������ fi, ������������ ����� �������� ��������
% iterationCount - ����� ���������� �������� ������� �������
inlineFunc = inline(FunctionString);
iterCount = 0;
xI = zeros(0, MaxIterationCount());
fI = zeros(0, MaxIterationCount());

df = char(diff(str2sym(FunctionString)));
ddf = char(diff(str2sym(df)));

dfun = inline(df);
ddfun = inline(ddf);

if (inlineFunc(a) * ddfun(a) > 0)
    xResult = b;
else
    xResult = a;
end

while (abs(dfun(xResult)) > eps)
    X0 = xResult; % X0 - �������� ����������� ����
    xResult =  X0 - (dfun(X0)/(ddfun(X0))); % ������ ������ ��������
    fResult = inlineFunc(xResult);
    iterCount = iterCount + 1;
    xI(iterCount) = xResult;
    fI(iterCount) = fResult;
end
iterationCount = iterCount;
end

% ����������� ������� fminbnd
function [xResult, fResult, xI, fI, iterCount, iterOverflow] = StandartFunc(FuncString, a, b)
[xResult, fResult, exitflag, output] = fminbnd(FuncString, a, b);
iterCount = output.funcCount;
xI = zeros(0, MaxIterationCount());
fI = zeros(0, MaxIterationCount());
end

function num = MaxIterationCount
num = 250;
end
