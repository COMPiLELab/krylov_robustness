function [robScoreNoUpdate, robScoreUpdate50, Aorg] = MIOBIBreakEdge2_weighted(fname,k,topT)
% filename: in format: source dest (integer)
% k: budget
% top t: top t eigenvalue and eigenvector to compute the rob. score
% robScoreUpdate50 - recompute at every 50 operations and update using
% perturbation theory
% robScoreNoUpdate - no recompute at all 

% FT: Aorg should be the modified matrix

% set the number of eigenpairs use in the perturbatin update 


t = topT; %50 <----- FT

format longG

[a,b,c] = find(fname);

data = [a b c];

N = max(max(data(:,1)),max(data(:,2)));

% build the matrix
A = sparse(data(:,1), data(:,2), data(:,3), N, N);
Aorg = A;

if ~ishermitian(A)
	error('MIOBIBreakEdge2_weighted:: adjacency matrix has to be symmetric');
end

R = zeros(1+1,1);
Rtrue = zeros(1+1,1);

en = eigs(A,topT);

R(1) = log((1/topT)*sum(exp(en)));
Rtrue(1) = log((1/topT)*sum(exp(en)));



[V,D] = eigs(A, t);
V = [abs(V(:,1)) V(:,2:t)];

% returning the eignvalue
e = diag(D); %D(D~=0);

%True value
for j = 1:k

    % returning the edges
    [p,r,~] = find(triu(Aorg,1));
    
    %procedure starts here
    Up = V(p,:); %getting the eigenvector value at p
    Ur = V(r,:); %getting the eignevector value at r
    c = exp(e)'; %-e(1))*exp(e(1)); %computing the constant c
    
    c = c(ones(1,length(p)),:);
    
    compScore = times(c,exp(-2.*times(Up,Ur))); %computing score for each (sum) component in Equ (8)
    
    score = sum(compScore,2);
    
    %returning the index of all the lowest score edges
    [~, minIndex] = min(score);
    %minIndex = find(score==min(score));
    
    %pick an arbitrary one
    remove = [p(minIndex(1)) r(minIndex(1))];
    
    
    Aorg(remove(1), remove(2)) = 0;
    Aorg(remove(2), remove(1)) = 0;
    
    deltaA = sparse(N,N);
    %update top t eigenvalues of A equation (4)
    deltaA(remove(1), remove(2)) = -A(remove(1), remove(2));
    deltaA(remove(2), remove(1)) = -A(remove(2), remove(1));
    
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
    
    
    deltaV = zeros(N,t);
    
    for i=1:t
        deltaVTemp = innerSum(:,i)';
        deltaVTemp = deltaVTemp(ones(length(V),1),:);
        %times(deltaVTemp, V)
        temp = sum(times(deltaVTemp, V),2);
        deltaV(:,i) = temp;
    end
    
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
%%%% FT modifies
robScoreNoUpdate = ((Rtrue(1)-Rtrue(2))*100/Rtrue(1));







[V,D] = eigs(A, t);
V = [abs(V(:,1)) V(:,2:t)];

% returning the eignvalue
e = diag(D); %D(D~=0);

for j = 1:k
    
    % returning the edges
    [p,r,~] = find(triu(A,1));
    
    %procedure starts here
    Up = V(p,:); %getting the eigenvector value at p
    Ur = V(r,:); %getting the eignevector value at r
    c = exp(e)'; %-e(1))*exp(e(1)); %computing the constant c
    
    c = c(ones(1,length(p)),:);
    
    compScore = times(c,exp(-2.*times(Up,Ur))); %computing score for each (sum) component in Equ (8)
    
    score = sum(compScore,2);
    
    %returning the index of all the lowest score edges
    [~, minIndex] = min(score);
    %minIndex = find(score==min(score));
    
    %pick an arbitrary one
    remove = [p(minIndex(1)) r(minIndex(1))];
    
    
    A(remove(1), remove(2)) = 0;
    A(remove(2), remove(1)) = 0;
    
    deltaA = sparse(N,N);
    %update top t eigenvalues of A equation (4)
    deltaA(remove(1), remove(2)) = -A(remove(1), remove(2));
    deltaA(remove(2), remove(1)) = -A(remove(2), remove(1));
    
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
    
    
    deltaV = zeros(N,t);
    
    for i=1:t
        deltaVTemp = innerSum(:,i)';
        deltaVTemp = deltaVTemp(ones(length(V),1),:);
        %times(deltaVTemp, V)
        temp = sum(times(deltaVTemp, V),2);
        deltaV(:,i) = temp;
    end
    
    tildaV = V + deltaV;
    
    for i=1:size(tildaV,2)
        tildaV(:,i) = tildaV(:,i)/norm(tildaV(:,i));
    end
    
    e = tildaE;
    V = tildaV;
    V = [abs(V(:,1)) V(:,2:t)];
    
    if mod(j,50) == 0
        [Vt,Dt] = eigs(A, t);
        [et] = diag(Dt);
        e=et;
        V=Vt;
        V = [abs(V(:,1)) V(:,2:t)];
    end
   
end


en = eigs(A,topT);
R(2) = log((1/topT)*sum(exp(en)));
robScoreUpdate50 = ((R(1)-R(2))*100/R(1));







end
