% min_{S>=0, S'*1=1, S*1=1, F'*F=I}  ||S - A||^2 + 2*lambda*trace(F'*Ln*F)
function [y1,y2,SS,U,V,term2] = coclustering_bipartite_fast1(A, G, c, sum_alpha,NITER, islocal)

if nargin < 6
    islocal = 1;
end

if nargin < 5
    NITER = 30;
end

lambda = 0.1;

[n,m] = size(G);
G = sparse(G);
a1 = sum(G,2);
D1a = spdiags(1./sqrt(a1),0,n,n);
a2 = sum(G,1);
D2a = spdiags(1./sqrt(a2'),0,m,m);
A1 = D1a*G*D2a;
SS2 = A1'*A1;
SS2 = full(SS2); 
[V, ev0, ~]=eig1(SS2,m);

V = V(:,1:c);
U=(A1*V)./(ones(n,1)*sqrt(ev0(1:c)'));
U = sqrt(2)/2*U; V = sqrt(2)/2*V;

A = full(A);

idxa = cell(n,1);
for i=1:n
    if islocal == 1
        idxa0 = find(A(i,:)>0);
    else
        idxa0 = 1:m;
    end;
    idxa{i} = idxa0;
end;

idxam = cell(m,1);
for i=1:m
    if islocal == 1
        idxa0 = find(A(:,i)>0);
    else
        idxa0 = 1:n;
    end;
    idxam{i} = idxa0;
end;

D1 = D1a; D2 = D2a; 
for iter = 1:NITER
    U1 = D1*U;
    V1 = D2*V;
    dist = L2_distance_1(U1',V1'); 
    S = zeros(n,m);
    for i=1:n
        idxa0 = idxa{i};
        ai = A(i,idxa0);
        di = dist(i,idxa0);
        ad = (ai-0.5*lambda*di)/sum_alpha; 
        S(i,idxa0) = EProjSimplex_new(ad);
    end;

    Sm = zeros(m,n);
    for i=1:m
        idxa0 = idxam{i};
        ai = A(idxa0,i);
        di = dist(idxa0,i);
        ad = (ai-0.5*lambda*di)/sum_alpha; 
        Sm(i,idxa0) = EProjSimplex_new(ad);
    end;

    S = sparse(S);
    Sm = sparse(Sm);
    SS = (S+Sm')/2;
    d1 = sum(SS,2);
    D1 = spdiags(1./sqrt(d1),0,n,n);
    d2 = sum(SS,1);
    D2 = spdiags(1./sqrt(d2'),0,m,m);
    SS1 = D1*SS*D2;

    SS2 = SS1'*SS1; 
    SS2 = full(SS2); 
    [V, ev0, ev]=eig1(SS2,c);
    U=(SS1*V)./(ones(n,1)*sqrt(ev0'));
    U = sqrt(2)/2*U; V = sqrt(2)/2*V;
    U_old = U;
    V_old = V;
    
    if length(ev) > c
        fn1 = sum(ev(1:c));
        fn2 = sum(ev(1:c+1));
        if fn1 < c-0.0000001
            lambda = 2*lambda; 
        elseif fn2 > c+1-0.0000001
            lambda = lambda/2; U = U_old;V = V_old;
        else
            break;
        end;
    else
        fn1 = sum(ev(1:c));
        if fn1 < c-0.0000001
            lambda = 2*lambda; 
        else
            break;
        end;
    end
end;
term2 = 0.1*(trace(U'*U)+trace(V'*V)-trace(U'*SS1*V));
SS0=sparse(n+m,n+m); SS0(1:n,n+1:end)=SS; SS0(n+1:end,1:n)=SS';
[~, y]=graphconncomp(SS0);
y1=y(1:n)';
y2=y(n+1:end)';



