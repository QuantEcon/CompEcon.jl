function demqua09

format long
fdjac(@f1,[0;1])
fdhess(@f2,[1;0])

function y = f1(x)
y = [exp(x(1))-x(2);
x(1)+x(2)^2;
(1-x(1))*log(x(2))];

function y = f2(x)
y = x(1)^2*exp(-x(2));