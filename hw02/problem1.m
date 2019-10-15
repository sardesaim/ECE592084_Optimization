clc; clear all; close all;

syms x1 x2 alph;
f = @(x1, x2) x1.^2+x2.^2+x1.*x2;
X_init{1} = [0.8,0.25];
gx = gradient(f(x1, x2))';
a = []; a0 = NaN ; b0 = NaN;
epsilon = 0.075
phia = @(alph) subs(f(x1-alph*gx(1), x2-alph*gx(1)), [x1 x2], X_init(1));
a(1) = -0.5; a(2) = a(1)+epsilon; a(3) = a(2)+2*epsilon;
%bracketing procedure
while (isnan(a0) && isnan(b0))
    if(double(phia(a(1)))>double(phia(a(2)))&& double(phia(a(2)))< double(phia(a(3))))
        a0 = a(1); b0 = a(3);
    elseif(double(phia(a(1))) > double(phia(a(2))) && double(phia(a(2))) > double(phia(a(3))))
        a(4) = a(3)+4*epsilon;
        if(double(phia(a(4)))>double(phia(a(3))))
            a0=a(1); b0 = a(4);
        else
            a(4) = a(4)+8*epsilon;
            a0=a(1); b0 = a(4);
        end
    else
    end
end
%golden section
epsi = 0.075;
N = ceil(log(0.01/(b0-a0))/log(0.61803));
rho = 0.382;
a = a0;
b = b0;
s = a+rho*(b-a);
t = a+(1-rho)*(b-a);
f1=double(phia(s))
f2=double(phia(t))
for n = 1:N
    if double(phia(s))< double(phia(t))
       b = t
       t = s
       s = a+rho*(b-a)
       f2 = f1
       f1 = double(phia(s))
    else 
       a = s
       s = t
       s = a+(1-rho)*(b-a)
       f1 = f2
       f2 = double(phia(t))
    end
end
alpha = (s+t)/2;
X_init{2} = X_init{1} - alpha*double(subs(gx,[x1 x2], X_init(1)));

x = linspace(-1,1,50);
y = x;
[x1,x2] = meshgrid(x,y);
box on; mesh(x1,x2,f(x1, x2));
figure;
contour(x1,x2,f(x1,x2), 20);
for i = 1:length(X_init)
    px(i) = X_init{i}(1);
    py(i) = X_init{i}(2);
end
hold on;
plot(px, py, '^-'); %plot sequence of points starting from starting point1