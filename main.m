clc, clear, close all;

points_2D = load('data/points2D.txt');
points_3D = load('data/points3D.txt');
num = size(points_2D, 1);
points_2D_homo = [points_2D, ones(50, 1)]';
points_3D_homo = [points_3D, ones(50, 1)]';



% Part I: Normalize the data
% For 2D points
T = zeros(3, 3);
u_x = mean(points_2D(:, 1));
u_y = mean(points_2D(:, 2));
sigma_x = var(points_2D(:, 1));
sigma_y = var(points_2D(:, 2));
sigma = sqrt(sigma_x + sigma_y);
s = sqrt(2)/sigma;
T(1, 1) = s;
T(2, 2) = s;
T(3, 3) = 1;
T(1, 3) = -u_x * s;
T(2, 3) = -u_y * s;
points_2D_norm = T*points_2D_homo;

% For 3D points
U = zeros(4, 4);
u_x = mean(points_3D(:, 1));
u_y = mean(points_3D(:, 2));
u_z = mean(points_3D(:, 3));
sigma_x = var(points_3D(:, 1));
sigma_y = var(points_3D(:, 2));
sigma_z = var(points_3D(:, 3));
sigma = sqrt(sigma_x + sigma_y + sigma_z);
s = sqrt(3)/sigma;
U(1, 1) = s;
U(2, 2) = s;
U(3, 3) = s;
U(4, 4) = 1;
U(1, 4) = -u_x * s;
U(2, 4) = -u_y * s;
U(3, 4) = -u_z * s;
points_3D_norm = U*points_3D_homo;

% Part II: Create design matrix
A = [];
for i = 1:num
    cur_point = points_2D_norm(:, i);
    V = cur_point;
    V(1) = cur_point(1) + sign(cur_point(1)) * norm(cur_point);
    Hv = eye(3, 3) - 2 * (V*V')/(V'*V);
    m = Hv(2:3, :);
    left_null = m;
    A_i = kron(left_null, points_3D_norm(:, i)');
    A = [A; A_i];
end

% Part III: Use SVD to solve the right null space of design matrix
[~, ~, V] = svd(A);
pVec = V(:, 12);
P = reshape(pVec, 4, 3);
P = P';

% Part IV: Back to the original coordinate
P = inv(T) * P * U;
P = P / norm(P, 'fro');

% Part V: Show the result
format shortg
P = -P
 
% -------------- Problem II ----------------
% Config
lambda = 1000000;
iterNum = 10;
cost = [];


pVec_param = vector_parameterize(pVec);
P = reshape(pVec, 4, 3);
P = P';
delta = 0;
points_2D_proj = P * points_3D_norm;
covariance_matrix = inv(T(1, 1)^2 * eye(2*num, 2*num));
% Normalization
for i = 1:num
    points_2D_proj(:, i) = points_2D_proj(:, i) / points_2D_proj(3, i);
end
pre_error = points_2D_norm - points_2D_proj;
pre_error = pre_error(1:2, :);
iter = 1;
pre_error_norm = norm(pre_error(:)'*covariance_matrix *pre_error(:));
while iter <= iterNum
% Step I: Calculate Jacobian Matrix
J = zeros(2 * num, 11); % initial

for i = 1 : num
    Ai = zeros(2, 11);
    X_deri_P = zeros(2, 12);
    P_deri_P = zeros(12, 11);
    point_3D = points_3D_norm(:, i);
    w = P(3, :) * point_3D;
    % Calc X_deri_P
    X_deri_P(1, 1:4) = point_3D';
    X_deri_P(2, 5:8) = point_3D';
    X_deri_P(1, 9:12) = -points_2D_proj(1, i) * point_3D';
    X_deri_P(1, 9:12) = -points_2D_proj(2, i) * point_3D';
    X_deri_P = X_deri_P / w;
    temp = norm(pVec_param)/2;
    P_deri_P(1,:) = -1/4* (sin(temp))/temp * pVec_param';
    
    if temp == 0
        P_deri_P(2:12, :) = eye(11)/2;
    
    else
        P_deri_P(2:12, :) = mysinc(temp)*eye(11)/2 + ...
            (1/(4*temp)) * desinc(temp) * (pVec_param * pVec_param');
    end
    
    Ai = X_deri_P * P_deri_P;
    J(i:i+1, :) = Ai;
end
delta = linsolve(J'*covariance_matrix*J + lambda*eye(11), J'*covariance_matrix*pre_error(:));
pVec_param_new = pVec_param + delta;
pVec_new = vector_deparameterize(pVec_param_new);
P_new = reshape(pVec_new, 4, 3);
P_new = P_new';
points_2D_proj_new = P_new * points_3D_norm;

% Normalization
for i = 1:num
    points_2D_proj_new(:, i) = points_2D_proj_new(:, i) / points_2D_proj_new(3, i);
end
cur_error = points_2D_norm - points_2D_proj_new;
cur_error = cur_error(1:2, :);
cur_error_norm = norm(cur_error(:)'*covariance_matrix*cur_error(:));

if(pre_error_norm < cur_error_norm)
    lambda = 10 * lambda;
    cost = [cost; cur_error_norm];
else
    lambda = lambda/10;
    % save the cost
    cost = [cost; cur_error_norm];
    pVec = pVec_new;
    pVec_param = pVec_param_new;
    points_2D_proj = points_2D_proj_new;
    pre_error = cur_error;
    pre_error_norm = cur_error_norm;
    P = P_new;
    iter = iter + 1;
end
end

P_original = inv(T)*P*U;
P_original = P_original/norm(P_original, 'fro');

format shortg
P_original
cost