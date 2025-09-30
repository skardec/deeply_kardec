%Para fase 1, vamos usar a mesma funcao da fase 2, so que com base = b e vetor c = zeros (zeramos todas as variáveis x_i) e deixamos w = b, logo temo uma solucao para Ax +w = b, com x=0 e w=b:
function [x,y,z,itera,base]=simplex2(A,b,c,base)
%[m,n]=size(A);
tol=1.0e-5;
iteramax=1e5;
%A=sparse(A);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%% Sem degenerecencia     %%%%%%%%%%%%%%%%%%%%
%b0=b;
%b=b+1.0e-5*rand(m,1);
itera=0;
while itera < iteramax
    % %%%%%%%%%%   %%%%%%%%%%   Achando as matrizes B e N   %%%%  %%%%%%%%
    B=A(:, base);
    M=[1:n]';
    M(base)=[];
    N=A(:,M);
    %%%%%%%%%%%%%%%%%%%    Cálculo de xb, y e u( costo relativo) %%%%%%%%%%%%%%
    [L,U]=lu(B);
    xb=U\(L\b);
    y=L'\(U'\c(base));
    u=c(M)-(N'*y);
    if min(u)>=-tol
        disp('solução ótima')
        break
    end
    %%%%%   critério de entrada na base (Encontrar o índice mínimo)%%%%%%%%%%%%%
    [r,j]=min(u);
    %%%%%%%%%%%%%%%%%%%%%%%%%% Cálculo de s   %%%%%%%%%%%%%%%%%%%%%%%%     %
    d=N(:,j);
    s=U\(L\d);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                verificando se a solução é ilimitada                     %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    p=find(s>tol);
    %p=find(s>0);
    if isempty(p)
        disp(' Solucao ilimitada!!!!')
        break
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                    Teste para saida da base                             %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    w=xb(p);    
    %[~,pmin]=min(w./s(p));         
    [nn,pmin]=min(w./s(p)); 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                      Atualizar a base                                   %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    aux=base(p(pmin));
    base(p(pmin))=M(j);
    M(j)=aux;
    itera=itera+1;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                      Função objetivo                                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x = zeros(n,1);
xb=U\(L\b0);
x(base)=xb;
x=sparse(x);
z=full(c'*x);
end

  e após rodar essa funcao nessas condicoes, teremos uma solucao basica, que vamos usar para para começar um novo problema, agora com c= c mesmo, e base = a solucao encontrada, logo em seguida temos a solucao otima.

variáveis a serem preencidas: 



function x = simplex_method_with_u(c, A, b, u)
    [m, n] = size(A);
    
    % Inicialização
    x = zeros(n, 1);
    B = eye(m);
    c_B = zeros(m, 1);
    B_inv = eye(m);
    
    % Encontre a solução básica inicial (viável)
    while true
        % Encontre a direção simplex com consideração de x_i <= u_i
        min_ratio = Inf;
        leaving_row = 0;
        for i = 1:m
            if B_inv(i, :) * A <= u && B_inv(i, :) * A >= 0
                d = B_inv(i, :) * A;
                ratio = (u - x(B' * (1:n))) ./ d; %esse é teuu tal de N_j
                [min_val, min_idx] = min(ratio);
                if min_val < min_ratio
                    min_ratio = min_val;
                    leaving_row = i;
                    entering_col = min_idx;
                end
            end
        end
        
        if leaving_row == 0
            break;  % A solução é ótima
        else
            d = B_inv(leaving_row, :) * A(:, entering_col);
            
            % Atualize a base
            B(:, leaving_row) = A(:, entering_col);
            
            % Atualize B_inv usando a decomposição LU
            [L, U] = lu(B);
            B_inv = U\(L);
            x=U\(L\b);
        end
    end
    
    %x(B' * (1:n)) = u;  % Garanta que x_i <= u_i para todo i
end