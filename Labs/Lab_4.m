eps = 1e-6;
a = 0;
b = 1;
%FuncString = 'cosh((3*(x^3) + 2*(x^2) - 4*x + 5)/3) + tanh((x^3 - 3*sqrt(2)*x -2)/(2*x + sqrt(2))) - 2.5;'
FuncString = '(cosh(x-0.111))^8;'
[xRes, fRes, xi, fi, iterCount] = newton(FuncString,a,b, eps);

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
xRes - 0.111

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
%X = cosh((3*(x^3)+2*(x^2)-4*x+5)/3)+tanh((x^3-3*sqrt(2)*x-2)/(2*x+sqrt(2)))-2.5;
X = (cosh(x-0.111))^8;
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

% ����� �������
function [xResult, fResult, xI, fI, iterCount] = newton(FuncString, a, b, eps)
iterCount = 0;
xI = zeros(0, MaxIterationCount());
fI = zeros(0, MaxIterationCount());
inlineFunc = inline(FuncString);

xResult = (a + b) / 2;

while (abs(funcDiff1(FuncString, xResult, eps)) > eps)
    iterCount = iterCount + 1;
    X0 = xResult;
    xI(iterCount) = X0 - (funcDiff1(FuncString, X0, eps)/(funcDiff2(FuncString, X0, eps)));
    fI(iterCount) = inlineFunc(xResult);
    xResult =  xI(iterCount);
    fResult = fI(iterCount);
end
end

function [value] = funcDiff1(FuncString, x, eps)
    inlineFunc = inline(FuncString);
    value = (inlineFunc(x + eps/2) - inlineFunc(x-eps/2)) / eps;
end 

function [value] = funcDiff2(FuncString, x, eps) 
    inlineFunc = inline(FuncString);
    value = (inlineFunc(x + eps) - 2 * inlineFunc(x) + inlineFunc(x - eps)) / eps / eps;
end

function [symdif] = derivativeRightDiff(symfun, start, ending, step)
symdif = ''; % ��������������� ������� ����������� (� ��������� �������������)
n = 10; % ������� �������� ���������������� �-���
xi = start:step:ending; % ������� ��������
xI = xi(1:length(xi)-1);
fun = str2func(['@(x)', symfun]);
dy = zeros(0, length(xi));
for i=1:length(xI)
    dy(i)=(fun(xi(i+1))-fun(xi(i)))/step;
end
p = polyfit(xI,dy,n); % ������������� ��������� ������� n
for i=1:length(p)
    symdif = strcat(symdif, num2str(p(i)), '*x^', num2str(length(p) - i));
    if i ~= length(p)
        symdif = strcat(symdif, '+');
    end
end
end

function num = MaxIterationCount
num = 250;
end
