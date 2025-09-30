AI = [A  eye(m)];
cI=[zeros(n,1); ones(m,1)]; % custo da fase I
base=n+1:n+m;
itera = 0;
tol=1.0e-4;
iteramax = 1e+4;
b=-b; 

%%%%%%%%%%%%%%%%%%%%%%%%FASE 1%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
while itera < iteramax
   
    % %%%%%%%%%%   %%%%%%%%%%   Achando as matrizes B e N   %%%%  %%%%%%%%
    B=AI(:, base); %base é o vetor que diz as posições de A não-nulas
    M=(1:n+m)';
    M(base)=[]; % Pré- Processamento (Elimina as linhas e colunas nulas de A)
    N=AI(:,M);
    %%%%%%%%%%%%%%%%%%%    Cálculo de xb, y e u (custo relativo) %%%%%%%%%%%%%%
    [L,U]=lu(B);
    xb=U\(L\b);
    cbase = cI(base);
    y=L'\(U'\cbase);
    u=cI(M) -(N'*y);
    if min(u)>=-tol
        disp('solução ótima')
        quantas_iteracoes = itera;
        disp("quantas iteracoes:")
        disp(quantas_iteracoes)
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
    base(p(pmin))= M(j);
    M=aux;
    itera=itera+1;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                      Função objetivo                                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%dot(xb,cI(base))

itera = 0;
%%%%%%%%%%%%%%%%%%%%%%%% FASE 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
while itera < iteramax
   
    % %%%%%%%%%%   %%%%%%%%%%   Achando as matrizes B e N   %%%%  %%%%%%%%
    B=A(:, base); %base é o vetor que diz as posições de A não-nulas
    M=(1:n)';
    M(base)=[];
    N=A(:,M);
    %%%%%%%%%%%%%%%%%%%    Cálculo de xb, y e u( costo relativo) %%%%%%%%%%%%%%
    [L,U]=lu(B);
    base_old=base;
    x_reserva=xb;
    xb=U\(L\b);
    cbase = c(base);
    y=L'\(U'\cbase);
    u=c(M)-(N'*y);
    if min(u)>=-tol && abs(dot(xb,c(base))-dot(x_reserva,c(base_old)))<=tol
        disp('solução ótima');
        quantas_iteracoes = itera;
        disp("passei por aqui");
        break
    end

    %%%%%   critério de entrada na base (Encontrar o índice mínimo)%%%%%%%%%%%%%


    [r,j]=min(u);  %(posição j onde ocorre o menor custo relativo negativo)

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
    %[~,pmin]=min(w./s(p));  % s(p) é o equivalente a N~ij da teoria, onde p é a posição da coluna j onde occorreu uma saída de variável)       
    [nn,pmin]=min(w./s(p)); 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                      Atualizar a base                                   %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    aux=base(p(pmin));
    %base_old=base;
    base(p(pmin))= M(j);
    M=aux;
    itera=itera+1;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                      Função objetivo                                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dot(xb,c(base))
