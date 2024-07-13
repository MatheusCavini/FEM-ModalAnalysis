classdef Portic
    properties
        L %Compriment
        E %Módulo de Young
        A %Área da seção
        I %Momento de inércia da seção
        r %Raio de giração
        rho %Densidade
        a %Ângulo
        K %Matriz de rigidez
        M %Matriz de massa
    end

    methods
          function obj = Portic(L, E, A, I, rho, a)
            obj.L = L;
            obj.E = E;
            obj.A = A;
            obj.I = I;
            obj.r = sqrt(I/A);
            obj.rho = rho;
            obj.a = a;

            %Matriz de rotação
            r = obj.r;
            c = cosd(a);
            s = sind(a);
            T = [ c, s, 0, 0, 0, 0;
                 s, -c, 0, 0, 0, 0;
                  0, 0, 1, 0, 0, 0;
                  0, 0, 0, c, s, 0;
                  0, 0, 0, s, -c, 0;
                  0, 0, 0, 0, 0, 1];

            %Definição da matriz de rigidez sem rotação
            k_local = E / L * ...
                [ A, 0, 0, -A, 0, 0;
                  0, 12*I/L^2, 6*I/L, 0, -12*I/L^2, 6*I/L;
                  0, 6*I/L, 4*I, 0, -6*I/L, 2*I;
                 -A, 0, 0, A, 0, 0;
                  0, -12*I/L^2, -6*I/L, 0, 12*I/L^2, -6*I/L;
                  0, 6*I/L, 2*I, 0, -6*I/L, 4*I];

            %Definição da matriz de massa sem rotação
            m_local = (rho * A * L) * ...
                [ 140/420 , 0, 0, 70/420, 0, 0;
                  0, 13/35 + 6*r^2/(5*L^2), 11/210 + r^2/(10*L), 0, 9/70 - 6*r^2/(5*L^2), -13*L/420 + r^2/(10*L);
                  0, 11/210 + r^2/(10*L), L^2/105 + 2*r^2/15, 0, 13*L/420 - r^2/(10*L), -L^2/140 - r^2/30;
                  70/420, 0, 0, 140/420, 0, 0;
                  0, 9/70 - 6*r^2/(5*L^2), 13*L/420 - r^2/(10*L), 0, 13/35 + 6*r^2/(5*L^2), -11*L/210 - r^2/(10*L);
                  0, -13*L/420 + r^2/(10*L), -L^2/140 - r^2/30, 0, -11*L/210 - r^2/(10*L), L^2/105 + 2*r^2/15];

            %Aplicação das rotações
            obj.M = T' * m_local * T;
            obj.K = T' * k_local * T;
        end
    end
end
