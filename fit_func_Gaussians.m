function F = fit_func_Gaussians(x,P,t)


F = P - (...
    x(7) + ...
    x(1)*x(3)*sqrt(2*pi)*normpdf(t,x(2),x(3)) + ...
    x(4)*x(6)*sqrt(2*pi)*normpdf(t,x(5),x(6)) );