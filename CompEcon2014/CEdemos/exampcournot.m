%% Exampcournot

  x = 2.3;
  for it=1:50
    step = -(x^4-2)/(4*x^3);
    x = x + step
    if norm(step)<1.e-10, break, end
  end
  
  f = @(x) x^4-2;
  broyden(f,2.3)