close all
clear
clc

load('A_mat.mat');

% % Calculate the rank of A.
% RR = rank(A);
% RRt = rank(A, 100);
% save('rankA.mat', 'RR', 'RRt');
% 
% % Calculate the condition number of A.
% N = cond(A);
% save('Acond.mat', 'N');
% 
% % Calculate the SVD.
% [U, S, V] = svd(A);
% sVec = diag(S);
% save('SVDval.mat', 'U', 'S', 'V');


% Moore-Penrose pseudo-inverse.
% [U, S, V] = svd(A, 'econ');
% save('svdComp.mat', 'U', 'S', 'V');

% % Alternative implementation.
% A_pinv = pinv(A);
% save('A_pinv_mat.mat', 'A_pinv');


% % Truncated Moore-Penrose pseudo-inverse. (tol = 200 good)
% [U, S, V] = svd(A, 'econ');
% save('svdComp.mat', 'U', 'S', 'V');

% % Alternative implementation.
% A_pinv = pinv(A,  200);
% save('A_pinv_trunc_mat.mat', 'A_pinv');


% % Moore-Penrose pseudo-inverse with zeroth-order Tikhonov
% % regularization.
% [U, S, V] = svd(A, 'econ');
% save('svdComp.mat', 'U', 'S', 'V');

% % Alternative implementation.
% lambda = 8e-8;
% A_new = A'*A + lambda*eye(size(A, 2));
% A_ZeroTik = pinv(A_new);
% save('A_ZeroTik_mat.mat', 'A_ZeroTik');


% % Truncated Moore-Penrose pseudo-inverse with zeroth-order Tikhonov
% % regularization.
% [U, S, V] = svd(A, 'econ');
% save('svdComp.mat', 'U', 'S', 'V');

% % Alternative implementation.
% lambda = 8e-8;
% A_new = A'*A + lambda*eye(size(A, 2));
% A_trunc_ZeroTik = pinv(A_new, 20000);
% save('A_trunc_ZeroTik_mat.mat', 'A_trunc_ZeroTik');

% % Randomized SVD.
% [U,S,V] = rsvd(A,size(A, 1));
% save('svdComp.mat', 'U', 'S', 'V');



% Find a good truncation parameter.
[U, S, V] = svd(A, 'econ');
sVec = diag(S);
CondMat = zeros(size(sVec));
for i = 1:length(sVec)
    sVecTemp = sVec(1:i);
    CondMat(i) = max(sVecTemp)/min(sVecTemp);
end

figure;
semilogy(1:length(sVec), CondMat, 'linewidth', 3);
set(gca, 'fontweight', 'bold', 'fontname', 'times', 'fontsize', 12);
xlabel('$\bf{truncation} \quad \bf{point} \quad \bf{\it{k}}$', 'Interpreter','Latex');
ylabel('$\bf{cond(G)}$', 'Interpreter','Latex');

figure;
Kstop = 500;
plot(1:Kstop, CondMat(1:Kstop), 'linewidth', 3);
hold on
yline(5, 'r--', 'linewidth', 2);
xlabel('$\bf{truncation} \quad \bf{point} \quad \bf{\it{k}}$', 'Interpreter','Latex');
ylabel('$\bf{cond(G)}$', 'Interpreter','Latex');
set(gca, 'fontweight', 'bold', 'fontname', 'times', 'fontsize', 16);
%}


%{
% Find a good lambda. (L-curve)
load('ReconsData.mat');
[U, S, V] = svd(A, 'econ');
sVec = diag(S);
lambda = logspace(-15, 15, 5000);
K = 1920;
MinTerm = zeros(length(lambda), 1);
Term1 = zeros(length(lambda), 1);
Term2 = zeros(length(lambda), 1);
for i = 1:length(lambda)
    x = zeros(length(size(A, 2)));
    for j = 1:K
        uVec = U(:, j);
        vVec = V(:, j);
        x = x + (sVec(j)/(sVec(j)^2 + lambda(i)^2))*(uVec'*b)*vVec;
    end
    MinTerm(i) = norm(b-A*x)^2 + lambda(i)*norm(x)^2;
    Term1(i) = norm(b-A*x);
    Term2(i) = norm(x);
end


lambdalin = log10(lambda);
NBin = 7;
xt = linspace(lambdalin(1), lambdalin(length(lambdalin)), NBin);

figure;
semilogx(lambda, MinTerm, 'linewidth', 3);
set(gca, 'xtick', 10.^xt, 'xticklabel', xt);
xlabel('$\bf{log_{10}(\eta)}$', 'Interpreter','Latex');
ylabel('$\bf{||G\rho-d||_2^2+\eta||x||_2^2}$', 'Interpreter','Latex');
set(gca, 'fontweight', 'bold', 'fontname', 'times', 'fontsize', 16);

figure;
semilogx(lambda, Term1, 'linewidth', 3);
set(gca, 'xtick', 10.^xt, 'xticklabel', xt);
xlabel('$\bf{log_{10}(\eta)}$', 'Interpreter','Latex');
ylabel('$\bf{||G\rho-d||_2}$', 'Interpreter','Latex');
set(gca, 'fontweight', 'bold', 'fontname', 'times', 'fontsize', 16);

figure;
semilogx(lambda, Term2, 'linewidth', 3);
set(gca, 'xtick', 10.^xt, 'xticklabel', xt);
xlabel('$\bf{log_{10}(\eta)}$', 'Interpreter','Latex');
ylabel('$\bf{||\rho||_2}$', 'Interpreter','Latex');
set(gca, 'fontweight', 'bold', 'fontname', 'times', 'fontsize', 16);

figure;
plot(Term1, Term2, 'linewidth', 3);
xlabel('$\bf{||G\rho-d||_2}$', 'Interpreter','Latex');
ylabel('$\bf{||\rho||_2}$', 'Interpreter','Latex');
set(gca, 'fontweight', 'bold', 'fontname', 'times', 'fontsize', 16);
%}


%{
% Find a good lambda. (generalized cross-validation)
load('ReconsData.mat');
[U, S, V] = svd(A);
sVec = diag(S);
lambda = logspace(-15, 15, 1000);
K = 162;
GCDTerm = zeros(length(lambda), 1);
tic;
for i = 1:length(lambda)
    x = zeros(length(size(A, 2)));
    for j = 1:K
        uVec = U(:, j);
        vVec = V(:, j);
        x = x + (sVec(j)/(sVec(j)^2 + lambda(i)^2))*(uVec'*b)*vVec;
    end
    bCal = U*S*V'*x;
    SDiag = S'*S+(lambda(i)^2)*eye(size(A, 2));
    AASharp = U*S*diag(1./diag(SDiag))*S'*U';
%     AASharp = (U*S/SDiag)*S'*U';
    GCDTermTemp = 0;
%     for k = 1:length(size(A, 1))
%         GCDTermTemp = GCDTermTemp + abs((bCal(k)-b(k))/(1-AASharp(k, k)))^2;
%     end
%     GCDTerm(i) = GCDTermTemp/size(A, 1);
    GCDTerm(i) = size(A, 1)*sum(abs(bCal-b).^2)/(abs(trace(eye(size(A, 1))-AASharp))^2);
end
tGCD = toc;

lambdalin = log10(lambda);
NBin = 7;
xt = linspace(lambdalin(1), lambdalin(length(lambdalin)), NBin);

figure;
semilogx(lambda, GCDTerm, 'linewidth', 2);
set(gca, 'xtick', 10.^xt, 'xticklabel', xt);
xlabel('$\bf{log_{10}(\alpha)}$', 'Interpreter','Latex');
ylabel('$\bf{g(\alpha)}$', 'Interpreter','Latex');
set(gca, 'fontweight', 'bold', 'fontname', 'times', 'fontsize', 16);

save('GCDTermPhantomApprox.mat', 'GCDTerm', 'lambda', 'lambdalin', 'NBin', 'xt', 'tGCD');
%}