function q = vnewton(func,v,p)
    % Perform one Newton step
    f = feval(func,v,p);
    j = vjacob(func,v,p);
    q = v - j\f;
 
