function [robScoreNoUpdate, robScoreUpdate50] = MIOBIBreakNode(fname,k,topT)
% filename: in format: source dest (integer)
% k: budget
% top t: top t eigenvalue and eigenvector to compute the rob. score
% robScoreUpdate50 - recompute at every 50 operations and update using
% perturbation theory
% robScoreNoUpdate - no recompute at all 



% set the number of eigenpairs use in the perturbatin update 
t = 50;

format longG
% [fname, num2str(N), num2str(k), num2str(t)]
%
% data = textread(strcat('',fname));

[a,b,~] = find(fname);

data = [a-1 b-1];

N = max(max(data(:,1)),max(data(:,2)))+1;

%count the number of edges removed
edgesRM = 0;
edgesRMP = 0;
updateCount = 0;

% make the graph symmetric
A = sparse(data(:,1)+1, data(:,2)+1, ones(size(data,1),1), N, N);
A = A+A';

Aorg = sparse(data(:,1)+1, data(:,2)+1, ones(size(data,1),1), N, N);
Aorg = Aorg+Aorg';

ind = find(A);
Aorg(ind) = 1;
A(ind) = 1;


ep=1;

R = zeros(1+1,1);
Rtrue = zeros(1+1,1);

en = eigs(A,topT);

R(1) = log((1/topT)*sum(exp(en)));
Rtrue(1) = log((1/topT)*sum(exp(en)));




[V,D] = eigs(A, t);
V = [abs(V(:,1)) V(:,2:t)];

% returning the eignvalue
[e] = diag(D); %D(D~=0);


%True value
for j = 1:k

    score = zeros(N,1);
    
    [p,r,~] = find(Aorg);


    for n = 1:N
        c = exp(e)';
        indexP = find((p == n));
        indexR = find((r == n));
        Up = V(p([indexP;indexR]),:);
        Ur = V(r([indexP;indexR]),:);
        tnig = -1.*times(Up,Ur);
        colSum = sum(tnig,1);
        score(n) = sum(times(c,exp(colSum)));
    end


    %returning the index of all the lowest score edges
    [~,minIndex] = min(score);
    %minIndex = find(score==min(score));

    %pick an arbitrary one
    %remove = [p(minIndex(1)) r(minIndex(1))];

    remove = minIndex(1);

    Aorg(remove,:) = 1-ep;
    Aorg(:,remove) = 1-ep;
    
    
    
    deltaA = sparse(N,N);
    %update top t eigenvalues of A equation (4)
    deltaA(remove,:) = -ep;
    deltaA(:,remove) = -ep;
    
    deltaH = (V'*deltaA*V); %
    
    deltaE = diag(deltaH); %diagonal is the values we want
    
    tildaE = e + deltaE;
    
    deltaH(logical(eye(size(deltaH)))) = 0;
    
    %update top t eigenvectors of A equation (9)
    diffE = e';
    diffE = diffE(ones(1,t),:) - e(:,ones(1,t)); % compute lambda_j - lambda_i
    diffE = rdivide(1,diffE); % compute 1/(lambda_j-lambda_i)
    %diffE(logical(eye(size(diffE)))) = 0; %set all the diagonal to zero
    [naR, naC] = find(isfinite(diffE));
    diffE(naR, naC) = 0;
    
    
    
    %diffE = abs(diffE);
    
    innerSum = times(deltaH, diffE); %the corresponding entry multipication gives the constant terms inside the sum
    
    % second way of updating the eigenvectors using loops
    %    deltaV = [];
    
    %     for i=1:t
    %         deltaVTemp = 1/2 * V(:,i);
    %         %pair = [];
    %         for k=1:t
    %             if i ~= k
    %                 %pair = [pair; i k (V(:,k)'*deltaA*V(:,i))  e(i)-e(k) (V(:,k)'*deltaA*V(:,i))/(e(i)-e(k))];
    %                 %((V(:,k)'*deltaA*V(:,i))/((e(i)-e(k))))*V(:,k)
    %                 deltaVTemp = deltaVTemp + ((V(:,k)'*deltaA*V(:,i))/((e(i)-e(k))))*V(:,k);
    %             end
    %         end
    %         %pair
    %         deltaV = [deltaV deltaVTemp];
    %     end
    
    %    V + deltaV
    
    % another way of updating the eigenvectors
    deltaV = zeros(N,t);
    
    for i=1:t
        deltaVTemp = innerSum(:,i)';
        deltaVTemp = deltaVTemp(ones(length(V),1),:);
        %times(deltaVTemp, V)
        temp = sum(times(deltaVTemp, V),2);
        deltaV(:,i) = temp;
    end
    %V - deltaV
    
    tildaV = V + deltaV;
    
    for i=1:size(tildaV,2)
        tildaV(:,i) = tildaV(:,i)/norm(tildaV(:,i));
    end 
    
    e = tildaE;
    V = tildaV;
    
    V = [abs(V(:,1)) V(:,2:t)];
end

enTrue = eigs(Aorg,topT);
Rtrue(2) = log((1/topT)*sum(exp(enTrue)));
robScoreNoUpdate = ((Rtrue(1)-Rtrue(2))*100/Rtrue(1));







[V,D] = eigs(A, t);
V = [abs(V(:,1)) V(:,2:t)];

% returning the eignvalue
[e] = diag(D); %D(D~=0);

for j = 1:k
    
    score = zeros(N,1);
    
    [p,r,~] = find(A);
    
    
    for n = 1:N
        c = exp(e)';
        indexP = find((p == n));
        indexR = find((r == n));
        if ~isempty(indexP) && ~isempty(indexR)
            Up = V(p([indexP;indexR]),:);
            Ur = V(r([indexP;indexR]),:);
            tnig = -1.*times(Up,Ur);
            colSum = sum(tnig,1);
            score(n) = sum(times(c,exp(colSum)));
        else
            score(n) = Inf('double');
        end
    end
    
    [~,minIndex] = min(score);
    
    remove = minIndex;
    
    edgesRMP = max(edgesRMP, sum(A(remove,:))+sum(A(:,remove)));

    
    A(remove,:) = 1-ep;
    A(:,remove) = 1-ep;
    
    deltaA = sparse(N,N);
    %update top t eigenvalues of A equation (4)
    deltaA(remove,:) = -ep;
    deltaA(:,remove) = -ep;
    
    deltaH = (V'*deltaA*V); %
    
    deltaE = diag(deltaH); %diagonal is the values we want
    
    tildaE = e + deltaE;
    
    deltaH(logical(eye(size(deltaH)))) = 0;
    
    %update top t eigenvectors of A equation (9)
    diffE = e';
    diffE = diffE(ones(1,t),:) - e(:,ones(1,t)); % compute lambda_j - lambda_i
    diffE = rdivide(1,diffE); % compute 1/(lambda_j-lambda_i)
    %diffE(logical(eye(size(diffE)))) = 0; %set all the diagonal to zero
    [naR, naC] = find(isfinite(diffE));
    diffE(naR, naC) = 0;
    
    
    
    %diffE = abs(diffE);
    
    innerSum = times(deltaH, diffE); %the corresponding entry multipication gives the constant terms inside the sum
    
    % second way of updating the eigenvectors using loops
    %    deltaV = [];
    
    %     for i=1:t
    %         deltaVTemp = 1/2 * V(:,i);
    %         %pair = [];
    %         for k=1:t
    %             if i ~= k
    %                 %pair = [pair; i k (V(:,k)'*deltaA*V(:,i))  e(i)-e(k) (V(:,k)'*deltaA*V(:,i))/(e(i)-e(k))];
    %                 %((V(:,k)'*deltaA*V(:,i))/((e(i)-e(k))))*V(:,k)
    %                 deltaVTemp = deltaVTemp + ((V(:,k)'*deltaA*V(:,i))/((e(i)-e(k))))*V(:,k);
    %             end
    %         end
    %         %pair
    %         deltaV = [deltaV deltaVTemp];
    %     end
    
    %    V + deltaV
    
    % another way of updating the eigenvectors
    deltaV = zeros(N,t);
    
    for i=1:t
        deltaVTemp = innerSum(:,i)';
        deltaVTemp = deltaVTemp(ones(length(V),1),:);
        %times(deltaVTemp, V)
        temp = sum(times(deltaVTemp, V),2);
        deltaV(:,i) = temp;
    end
    %V - deltaV
    
    tildaV = V + deltaV;
    
    for i=1:size(tildaV,2)
        tildaV(:,i) = tildaV(:,i)/norm(tildaV(:,i));
    end 
    
    e = tildaE;
    V = tildaV;
    
    V = [abs(V(:,1)) V(:,2:t)];
    
    edgesRM = edgesRM + edgesRMP;
    edgesRMP = 0;
    
    if edgesRM >= 50
       [Vt,Dt] = eigs(A, t);
       [et] = diag(Dt);
       e=et;
       V=Vt;
        
       edgesRM = mod(edgesRM,50);
       updateCount = updateCount + 1;
       V = [abs(V(:,1)) V(:,2:t)];
   end
    
end

en = eigs(A,topT);
R(2) = log((1/topT)*sum(exp(en)));
robScoreUpdate50 = ((R(1)-R(2))*100/R(1));



end
