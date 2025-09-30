[m,n]=size(A);
itera = 0;
tol=1.0e-4;
iteramax = 3e+3;
eps=1e-8;
b_til=[b;u];
A_til = [A zeros(size(A)); eye(n), eye(n)];
[L_til,U_til]=lu(A_til*A_til');
A_a = A*A';
[L,U]=lu(A_a);
x0=A_til'*(U_til\(L_til\b_til));
s=x0(n+1:n+n);
x0=x0(1:n);
eps2=100;
eps1=max(-min(x0),eps2);
eps1=max(eps1,norm(b,1)/(eps2*norm(A,1)));
x=max(x0,eps1);
y=zeros(m,1);
z=x;
itrc=0;
fprintf('%s %15s\n','iteração','gap')

while (norm(b-A*x)/(norm(b)+1)>1e-8)||(norm(c-A'*y-z)/(norm(c)+1)>1e-8)||(x'*z/(abs(c'*x)+abs(b'*y)+1)>1e-8)
    rp = b - A*x; 
    rd = c - A'*y - z;
    gamma = x'*z;
    mu = gamma/n^2;
    rc = mu * ones(n,1) - x.*z;
    d = x./z;
    S = A*diag(sparse(d))*A';
    R=chol(S);
    dy = R'*R\(rp + A*(d.*(rd - rc./x)));
    dx = d.*(A'*dy - rd + rc./x);
    dz = (rc - z.*dx)./x;
    %dv = ru - dx;

    alfap = min(1,-.99995/min(-.99995,min(dx./x)));
    alfad = min(1,-.99995/min(-.99995,min(dz./z)));

    gamma = x'*z;
    gamma_1=(x+alfap.*dx)'*(z+alfad.*dz);
    if(gamma_1<1)
        mu = gamma/n^2;
    else
        mu = ((gamma_1/gamma)^2)*gamma_1/n;
    end
    rc = mu * ones(n,1) - x.*z;

    dy = R'*R\(rp + A*(d.*(rd - rc./x)));
    dx = d.*(A'*dy - rd + rc./x);
    dz = (rc - z.*dx)./x;

    alfap = min(1,-.99995/min(-.99995,min(dx./x)));
    alfad = min(1,-.99995/min(-.99995,min(dz./z)));
    x=x+alfap*dx;
    %disp(min(x))
    y=y+alfad*dy;
    z=z+alfad*dz;
    itrc = itrc+1;
    fprintf('%d% 12f\n',itrc,gamma)
end

disp(dot(c,x)); 
%p = symamd(S);
%S = S(p,p);
R = chol(sparse(S));
nnz_R = nnz(R); 
figure(2);
spy((R));
title('Matriz R Reordenada:')



