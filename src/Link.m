classdef Link
    properties
        L %Compriment
        E %Módulo de Young
        A %Área da seção
        rho %Densidade
        a %Ângulo
        K %Matriz de rigidez
        M %Matriz de massa
    end

    methods
        function obj = Link(L, E, A, rho, a)
            obj.L = L;
            obj.E = E;
            obj.A = A;
            obj.rho = rho;
            obj.a = a;
            %Definição da matriz de rigidez conforme base teórica
            obj.K = ((E * A) / L) * ...
              [cosd(a)^2, cosd(a) * sind(a), 0, -cosd(a)^2, -cosd(a) * sind(a), 0;
               cosd(a) * sind(a), sind(a)^2, 0, -cosd(a) * sind(a), -sind(a)^2, 0;
               0, 0, 0, 0, 0, 0;
               -cosd(a)^2, -cosd(a) * sind(a), 0, cosd(a)^2, cosd(a) * sind(a), 0;
               -cosd(a) * sind(a), -sind(a)^2, 0, cosd(a) * sind(a), sind(a)^2, 0;
               0, 0, 0, 0, 0, 0];

            %Definição da matriz de massa conforme base teórica
            obj.M = ((rho * A * L) / 6) * ...
              [2 * cosd(a)^2, 2 * cosd(a) * sind(a), 0, cosd(a)^2, cosd(a) * sind(a), 0;
               2 * cosd(a) * sind(a), 2 * sind(a)^2, 0, cosd(a) * sind(a), sind(a)^2, 0;
               0, 0, 0, 0, 0, 0;
               cosd(a)^2, cosd(a) * sind(a), 0, 2 * cosd(a)^2, 2 * cosd(a) * sind(a), 0;
               cosd(a) * sind(a), sind(a)^2, 0, 2 * cosd(a) * sind(a), 2 * sind(a)^2, 0;
               0, 0, 0, 0, 0, 0];
        end
    end
end
