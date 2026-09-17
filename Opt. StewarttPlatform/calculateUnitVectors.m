
function [unitVectors, forceValues, cg] = calculateUnitVectors(fIDj, fixedPoints,eul)
    import casadi.*

    if length(fIDj) ~= 15
        error('fIDj debe tener 15 componentes.');
    end
    if size(fixedPoints, 1) ~= 6 || size(fixedPoints, 2) ~= 3
        error('fixedPoints debe ser una matriz 6x3 con las coordenadas de los puntos fijos.');
    end

    Rz = fIDj(1); 
    Rx = fIDj(2);
    Ry = fIDj(3);
    Fx = fIDj(4);
    Fy = fIDj(5);
    Fz = fIDj(6);
    Q1 = fIDj(7:9);   %  [Q1x, Q1y, Q1z]
    Q2 = fIDj(10:12); %  [Q2x, Q2y, Q2z]
    Q3 = fIDj(13:15); %  [Q3x, Q3y, Q3z]
    cg= (Q1+Q2+Q3)/3
   
   % Extraer ángulos de Euler: psi, theta, phi
    psi = eul(1);    % Rotación alrededor del eje Z
    theta = eul(2);  % Rotación alrededor del eje X
    phi = eul(3);    % Rotación alrededor del eje Y
    
% Definir las matrices de rotación
    Rz1 = [cos(psi), -sin(psi), 0; 
      sin(psi),  cos(psi), 0; 
      0,        0,        1];

    Rx1 = [1, 0,         0;
      0, cos(theta), -sin(theta);
      0, sin(theta),  cos(theta)];

    Ry1 = [cos(phi),  0, sin(phi);
      0,         1, 0;
      -sin(phi), 0, cos(phi)];

    Mx=Rz1*[Rx;0;0];
    My=Rz1*Rx1*[0;Ry;0];
    Mz=[0;0;Rz];

    M=Mx+My+Mz;

    % Inicializar matriz de vectores unitarios
    unitVectors = MX.zeros(6, 3);

    % O1 y O2 hacia Q1
    unitVectors(1, :) = (Q1 - fixedPoints(1, :)') / sqrt(sum((Q1 - fixedPoints(1, :)').^2)+1e-6);
    unitVectors(2, :) = (Q1 - fixedPoints(2, :)') / sqrt(sum((Q1 - fixedPoints(2, :)').^2)+1e-6);

    % O3 y O4 hacia Q2
    unitVectors(3, :) = (Q2 - fixedPoints(3, :)') / sqrt(sum((Q2 - fixedPoints(3, :)').^2)+1e-6);
    unitVectors(4, :) = (Q2 - fixedPoints(4, :)') / sqrt(sum((Q2 - fixedPoints(4, :)').^2)+1e-6);

    % O5 y O6 hacia Q3
    unitVectors(5, :) = (Q3 - fixedPoints(5, :)') / sqrt(sum((Q3 - fixedPoints(5, :)').^2)+1e-6);
    unitVectors(6, :) = (Q3 - fixedPoints(6, :)') / sqrt(sum((Q3 - fixedPoints(6, :)').^2)+1e-6);

    r1 = Q1 - cg; % Vector de cg a Q1
    r2 = Q2 - cg; % Vector de cg a Q2
    r3 = Q3 - cg; % Vector de cg a Q3


     % 1. Ecuaciones de fuerzas
     A_forces=unitVectors';

    % 2. Ecuaciones de momentos (producto cruzado con r)
    A_moments = [
        cross(r1, unitVectors(1, :)') , cross(r1, unitVectors(2, :)') , ...
        cross(r2, unitVectors(3, :)') , cross(r2, unitVectors(4, :)') , ...
        cross(r3, unitVectors(5, :)') , cross(r3, unitVectors(6, :)');
    ];

    % Matriz de ecuaciones combinadas
    A = [A_forces; A_moments];

    % Vector de términos independientes
    b = [Fx; Fy; Fz; M];

    % Resolver el sistema 
    forceValues = inv(A)* b;
end