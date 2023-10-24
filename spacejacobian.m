S = [0 0 1 0 0 0;
     -1 0 0 0 -360 0;
     0 0 1 0 0 0;
     1 0 0 0 780 0;
     0 0 1 0 0 0;
     -1 0 0 0 -1180 0;
     0 0 1 0 0 0]';

q = [pi*0.25, pi*0.25, 0, -pi*0.25, 0, pi*0.25, 0];
R_home = [1 0 0; 0 1 0; 0 0 1]';
t_home = [0 0 1306]';
M = [R_home t_home; 0 0 0 1];
Js = jacob0(S,q);
Tsb = fkine(S,M,q);
fprintf('The forward kinematics transformation:\n');
disp(Tsb);
fprintf('The jacobian in space frame is:\n');
disp(Js);
rsb = Tsb(1:3,1:3);
p = Tsb(1:3,4);
Tbs = inverse(Tsb);
fprintf('The space frame w.r.t body frame is:\n');
disp(Tbs);
ADTbs = Adjoint(Tbs);
Jb = ADTbs*Js;
fprintf('The jacobian in body frame is:\n');
disp(Jb);
Jbw = Jb(1:3,1:7);
Jbv = Jb(4:6,1:7);
Jav = rsb*Jbv;
fprintf('The velocity component of analytic jacobian is:\n');
disp(Jav);
CJbw = condition(Jbw);
CJbv = condition(Jbv);
fprintf('The condition number for angular velocity component is:\n');
disp(CJbw);
fprintf('The condition number for linear velocity component is:\n');
disp(CJbv);
fprintf('The rank of the jacobian is:\n');
disp(rank(Jb));
F = [0, 0, -100]';
Tor = JointTorque(F,p,Js);
fprintf('The joint torques for the given force is:\n');
disp(Tor);
Tsd = [0.707, 0, -0.707, -500;
       0.707, 0, 0.707, 600;
       0, -1, 0, 400;
       0, 0, 0, 1];
Twbd = logm(Tbs*Tsd);
Vbd = [Twbd(3,2), Twbd(1,3), Twbd(2,1), Twbd(1,4), Twbd(2,4), Twbd(3,4)]';
fprintf('The twist required to move end effector from b to d in one second is Vb:\n');
disp(Vbd);
for i = 1:6
    Tbs = inverse(fkine(S, M, q));
    Twbd = logm(Tbs*Tsd);
    Vbd = [Twbd(3,2), Twbd(1,3), Twbd(2,1), Twbd(1,4), Twbd(2,4), Twbd(3,4)]';
    Jb = Adjoint(Tbs)*jacob0(S,q);
    Jinv = Jb.'*inv(Jb*Jb');
    Theta = q';
    Theta = Theta+Jinv*Vbd;
    q = Theta';
    fprintf('iteration/n');
    disp(i);
    disp(Theta);
end
fprintf('The final Joint variables list is: \n')
disp(Theta);
fprintf('The End effector position for this joint variable:\n')
disp(fkine(S, M, q));
function T = JointTorque(F,p,Js)
    Fs = [cross(p,F);
          F];
    T = Js.'*Fs;
end
function T = condition(J)
    A = J*J.';
    [V,D] = eig(A);
    Diag = diag(D);
    Lmax = max(Diag);
    Lmin = min(Diag);
    T = Lmax/Lmin;
end
function T = fkine(S, M, q)
    T = eye(4);
    for i = 1:size(S,2)
        twist = S(:,i);
        v = twist(4:6);
        omega = twist(1:3);
        H = POE(v,omega,q(i));
        T = T*H;
    end
    function T = POE(v,omega,q)
        I = eye(3);
        Om1_hat = [0, -omega(3), omega(2); omega(3), 0, -omega(1); -omega(2), omega(1), 0];
        exp1 = eye(3)+Om1_hat*sin(q)+(Om1_hat^2)*(1-cos(q));
        T = [exp1,(I*q+(1-cos(q))*Om1_hat+(q-sin(q))*Om1_hat^2)*v;0,0,0,1];
    end 
    T = T*M;
end

function Js = jacob0(Slist, thetalist)

    Js = Slist;
    T = eye(4);
    for i = 2: length(thetalist)
        T = T * MatrixExp6(VecTose3(Slist(:, i - 1) * thetalist(i - 1)));
	    Js(:, i) = Adjoint(T) * Slist(:, i);
    end
end

function T = MatrixExp6(se3mat_1)

    omgtheta = so3ToVec(se3mat_1(1: 3, 1: 3));
    if NearZero(norm(omgtheta))
        T = [eye(3), se3mat_1(1: 3, 4); 0, 0, 0, 1];
    else
        [omghat, theta] = AxisAng3(omgtheta);
        omgmat = se3mat_1(1: 3, 1: 3) / theta; 
        T = [MatrixExp3(se3mat_1(1: 3, 1: 3)), ...
             (eye(3) * theta + (1 - cos(theta)) * omgmat ...
              + (theta - sin(theta)) * omgmat * omgmat) ...
                * se3mat_1(1: 3, 4) / theta;
             0, 0, 0, 1];
    end
end

function  R = MatrixExp3(so3mat)
    omgtheta = so3ToVec(so3mat);
    if NearZero(norm(omgtheta))
        R = eye(3);
    else
        [omghat, theta] = AxisAng3(omgtheta);
        omgmat = so3mat / theta;
        R = eye(3) + sin(theta) * omgmat + (1 - cos(theta)) * omgmat * omgmat;
    end
end

function se3mat = VecTose3(V)
    se3mat = [VecToso3(V(1: 3)), V(4: 6); 0, 0, 0, 0];
end

function so3mat = VecToso3(omg)
    so3mat = [0, -omg(3), omg(2); omg(3), 0, -omg(1); -omg(2), omg(1), 0];
end

function omg = so3ToVec(so3mat)
    omg = [so3mat(3, 2); so3mat(1, 3); so3mat(2, 1)];
end

function V = se3ToVec(se3mat)
    V = [se3mat(3, 2); se3mat(1, 3); se3mat(2, 1); se3mat(1: 3, 4)];
end

function judge = NearZero(near)
    judge = norm(near) < 1e-6;
end

function [omghat, theta] = AxisAng3(expc3)
    theta = norm(expc3);
    omghat = expc3 / theta;
end

function AdT = Adjoint(T)
    [R, p] = TransToRp(T);
    AdT = [R, zeros(3); VecToso3(p) * R, R];
end

function [R, p] = TransToRp(T)
    R = T(1: 3, 1: 3);
    p = T(1: 3, 4);
end

function T = inverse(H)
    R = H(1: 3, 1: 3);
    p = H(1: 3, 4);
    T = [R.', -R.'*p;
        0, 0, 0, 1];
end