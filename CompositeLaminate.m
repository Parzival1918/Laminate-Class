classdef CompositeLaminate < handle
    %COMPOSITELAMINATE Summary: Create a laminate from laminae of class
    %CompositeLamina.
    %   Detailed explanation goes here
    
    properties (SetAccess = private)
        S (3,3) double {mustBeFloat} = zeros(3)
        C (3,3) double {mustBeFloat} = zeros(3)

        laminateE1 (1,1) double {mustBeFloat} = 0
        laminateE2 (1,1) double {mustBeFloat} = 0
        laminatev12 (1,1) double {mustBeFloat} = 0
        laminatev21 (1,1) double {mustBeFloat} = 0
        laminateG12 (1,1) double {mustBeFloat} = 0
        laminateIR121 (1,1) double {mustBeFloat} = 0
        laminateIR122 (1,1) double {mustBeFloat} = 0

        is_symmetric (1,1) logical = false
        is_balanced (1,1) logical = false
    end

    properties
        laminae (1,:) CompositeLamina 
        lamina_angles (1,:) double {mustBeFloat,mustBe90negpos(lamina_angles)} 
        thicknesses (1,:) double {mustBeFloat}
    end
    
    methods
        function obj = CompositeLaminate(laminae_in, angles, t)
            %COMPOSITELAMINATE Construct an instance of this class
            %   Detailed explanation goes here
            arguments
                laminae_in (1,:) CompositeLamina
                angles (1,:) double {mustBeFloat, ...
                       mustBeSameSize(angles,laminae_in)} = zeros([1 length(laminae_in)])
                t (1,:) double {mustBeFloat, ...
                  mustBeSameSize(angles,t)} = ones([1 length(laminae_in)])
            end
            
            obj.laminae = laminae_in;
            obj.lamina_angles = angles;
            obj.thicknesses = t;
        end

        function is_symmetric = get.is_symmetric(obj)
            angles = obj.lamina_angles;
            angles_fliped = flip(obj.lamina_angles);

            if sum(angles == angles_fliped)/length(angles) == 1
                is_symmetric = true;
            else
                is_symmetric = false;
            end
        end

        function is_balanced = get.is_balanced(obj)
            angles = obj.lamina_angles;
            angles(angles == 0) = [];
            angles(angles == 90) = [];

            if isempty(angles)
                is_balanced = false;
            else
                sumation = sum(angles);
                if sumation == 0
                    is_balanced = true;
                else
                    is_balanced = false;
                end
            end
        end
        
        function S = get.S(obj)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            S = ones(3);
            
            sumt = sum(obj.thicknesses);
            for i = 1:3
                for j = 1:3
                    sumSij = 0;
                    for k = 1:length(obj.laminae)
                        laminak = obj.laminae(k);
                        tk = obj.thicknesses(k);
                        Sk = smaller_mat(laminak.S);
                        anglek = obj.lamina_angles(k);

                        Sij_bar = inv(U(anglek))*Sk*T(anglek);
                        sumSij = sumSij + Sij_bar(i,j)*tk;
                    end
                    S(i,j) = sumSij/sumt;
                end
            end
        end

        function laminateE1 = get.laminateE1(obj)
            if ne(obj.S(1,1),0)
                laminateE1 = 1/obj.S(1,1);
            else
                laminateE1 = 0;
            end
        end

        function laminateE2 = get.laminateE2(obj)
            if ne(obj.S(2,2),0)
                laminateE2 = 1/obj.S(2,2);
            else
                laminateE2 = 0;
            end
        end

        function laminateG12 = get.laminateG12(obj)
            if ne(obj.S(3,3),0)
                laminateG12 = 1/obj.S(3,3);
            else
                laminateG12 = 0;
            end
        end

        function laminatev12 = get.laminatev12(obj)
            if ne(obj.S(1,2),0)
                laminatev12 = -obj.laminateE1*obj.S(1,2);
            else
                laminatev12 = 0;
            end
        end

        function laminatev21 = get.laminatev21(obj)
            if ne(obj.S(1,2),0)
                laminatev21 = -obj.laminateE2*obj.S(1,2);
            else
                laminatev21 = 0;
            end
        end

        function laminateIR121 = get.laminateIR121(obj)
            if ne(obj.S(1,3),0)
                laminateIR121 = obj.laminateE1*obj.S(1,3);
            else
                laminateIR121 = 0;
            end
        end

        function laminateIR122 = get.laminateIR122(obj)
            if ne(obj.S(2,3),0)
                laminateIR122 = obj.laminateE2*obj.S(2,3);
            else
                laminateIR122 = 0;
            end
        end

        function Sbar = S_angle(obj, angles)
            arguments
                obj CompositeLaminate
                angles (1,:) double {mustBeFloat} = 0
            end
            Sbar = cell([1 length(angles)]);

            for i = 1:length(angles)
                S_ = obj.S;
                U_ = U(angles(i));
                T_ = T(angles(i));
    
                Sbar{i} = inv(U_)*S_*T_;
            end
        end

        function E = E_angle(obj, angles, is_plot)
            arguments
                obj CompositeLaminate
                angles (1,:) double {mustBeFloat} = linspace(0,90,30)
                is_plot (1,1) logical = true
            end
            
            Sbar = obj.S_angle(angles);
            E = zeros([2 length(Sbar)]);
            
            for i = 1:length(Sbar)
                S_ = Sbar{i};
                S11 = S_(1,1);
                S22 = S_(2,2);
                E1 = 1/S11;
                E2 = 1/S22;
                
                E(:,i) = [E1;E2];
            end

            if is_plot
                figure;
                plot(angles,E(1,:),angles,E(2,:));
                title("Composite YMs gainst loading angle");
                xlabel("Angle (º)");
                ylabel("Major and Minor Modulus (" + obj.laminae(1).compositeComponents.E_units + ")");
                ylim([min(E,[],"all")-1 max(E,[],"all")+1])
                legend("E1'","E2'");
            end
        end

        function G = G_angle(obj, angles, is_plot)
            arguments
                obj CompositeLaminate
                angles (1,:) double {mustBeFloat} = linspace(0,90,30)
                is_plot (1,1) logical = true
            end
            
            Sbar = obj.S_angle(angles);
            G = zeros([1 length(Sbar)]);
            
            for i = 1:length(Sbar)
                S_ = Sbar{i};
                S66 = S_(3,3);
                G12 = 1/S66;
                
                G(1,i) = G12;
            end

            if is_plot
                figure;
                plot(angles,G);
                title("Composite G gainst loading angle");
                xlabel("Angle (º)");
                ylabel("Shear Modulus (" + obj.laminae(1).compositeComponents.G_units + ")");
                ylim([min(G,[],"all")-1 max(G,[],"all")+1])
                legend("G12'");
            end
        end

        function v = v_angle(obj, angles, is_plot)
            arguments
                obj CompositeLaminate
                angles (1,:) double {mustBeFloat} = linspace(0,90,30)
                is_plot (1,1) logical = true
            end
            
            Sbar = obj.S_angle(angles);
            v = zeros([2 length(Sbar)]);
            
            for i = 1:length(Sbar)
                S_ = Sbar{i};
                S11 = S_(1,1);
                S22 = S_(2,2);
                S12 = S_(1,2);
                E1 = 1/S11;
                E2 = 1/S22;
                v12 = -S12*E1;
                v21 = -S12*E2;
                
                v(:,i) = [v12;v21];
            end

            if is_plot
                figure;
                plot(angles,v(1,:),angles,v(2,:));
                title("Composite v gainst loading angle");
                xlabel("Angle (º)");
                ylabel("Major and Minor Poissons Ratio (" + obj.laminae(1).compositeComponents.v_units + ")");
                ylim([min(v,[],"all")-0.2 max(v,[],"all")+0.2])
                legend("v12'","v21'");
            end
        end

        function IR = IR_angle(obj, angles, is_plot)
            arguments
                obj CompositeLaminate
                angles (1,:) double {mustBeFloat} = linspace(0,90,30)
                is_plot (1,1) logical = true
            end
            
            Sbar = obj.S_angle(angles);
            IR = zeros([2 length(Sbar)]);
            
            for i = 1:length(Sbar)
                S_ = Sbar{i};
                S16 = S_(1,3);
                S26 = S_(2,3);
                S11 = S_(1,1);
                S22 = S_(2,2);
                E1 = 1/S11;
                E2 = 1/S22;
                nu121 = S16*E1;
                nu122 = S26*E2;
                
                IR(:,i) = [nu121;nu122];
            end

            if is_plot
                figure;
                plot(angles,IR(1,:),angles,IR(2,:));
                title("Composite interaction ratios gainst loading angle");
                xlabel("Angle (º)");
                ylabel("Interaction Ratios");
                ylim([min(IR,[],"all")-0.2 max(IR,[],"all")+0.2])
                legend("n121'","n122'");
            end
        end

        function stress_out = stress_lamina_angle(obj,angle,k,stress)
            arguments
                obj CompositeLaminate
                angle (1,1) double {mustBeFloat} = 0
                k (1,1) double {mustBeFloat,mustBeMaxLaminaeAmount(obj,k)} = 1
                stress (3,1) double {mustBeFloat} = [10;0;0]
            end
            
            stress_out = zeros([3 1]);
            S_c = obj.S_angle(angle);
            S_c = S_c{1};

            for i = 1:length(k)
                k_i = k(i);
                lamina_k = obj.laminae(k_i);
                C_k = smaller_mat(lamina_k.C);
                U_k = U(obj.lamina_angles(k_i) + angle);

                stress_out = C_k*U_k*S_c*stress;
            end
        end

        function stress_out = stress_lamina_angles(obj,k,angles,stress,is_plot)
            arguments
                obj CompositeLaminate
                k (1,1) double {mustBeFloat,mustBeMaxLaminaeAmount(obj,k)} = 1
                angles (1,:) double {mustBeFloat} = linspace(0,90,30)
                stress (3,1) double {mustBeFloat} = [10;0;0]
                is_plot (1,1) logical = true

            end

            stress_out = zeros([3 length(angles)]);

            for i = 1:length(angles)
                S_c = obj.S_angle(angles(i));
                S_c = S_c{1};
                for j = 1:length(k)
                    k_j = k(j);
                    lamina_k = obj.laminae(k_j);
                    C_k = smaller_mat(lamina_k.C);
                    U_k = U(obj.lamina_angles(k_j) + angles(i));

                    stress_out(:,i) = C_k*U_k*S_c*stress;
                end
            end

            if is_plot
                figure;
                plot(angles,stress_out(1,:));
                xlabel("Angle (º)");
                ylabel("In-plane Stress in the Fibre Direction");
                legend("stress");
            end
        end
    end
end

function mustBeSameSize(A,B)
    if ne(size(A,2),size(B,2))
        error("Inputs are not same sizes.")
    end
end

function mustBeMaxLaminaeAmount(obj,k)
    if k > length(obj.laminae)
        error("Trying to access properties of a Laminae that is not in the Laminate.")
    end
end

function mustBe90negpos(angles)
    for i = 1:length(angles) 
        angle = angles(i);
        if angle < -90 || angle > 90
            error("Angle must be between -90º and 90º.")
        end
    end
end

function val = U(angle)
    val = [cosd(angle)^2 sind(angle)^2 cosd(angle)*sind(angle); ...
           sind(angle)^2 cosd(angle)^2 -cosd(angle)*sind(angle); ...
           -2*cosd(angle)*sind(angle) 2*cosd(angle)*sind(angle) cosd(angle)^2 - sind(angle)^2];
end

function val = T(angle)
    val = [cosd(angle)^2 sind(angle)^2 2*cosd(angle)*sind(angle); ...
           sind(angle)^2 cosd(angle)^2 -2*cosd(angle)*sind(angle); ...
           -cosd(angle)*sind(angle) cosd(angle)*sind(angle) cosd(angle)^2 - sind(angle)^2];
end

function mat = smaller_mat(m)
    mat = [m(1,1) m(1,2) m(1,3); ...
           m(2,1) m(2,2) m(2,3); ...
           m(6,1) m(6,2) m(6,6)];
end