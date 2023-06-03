
while 1
    n = input('Enter n: ');
    A = zeros(n,n);
    for j = 1:n
        for k = 1:n
            A(j,k) = n+1-max(j,k);
        end
    end
    A
    sindex = input('Enter indices where s_i is not +/- 1: ');
    P = A(sindex,sindex)
    u = input('Enter u values (+/-1): ');
    uindex = setdiff([1:n],sindex);
    s = -inv(P)*A(sindex,uindex)*u
end

