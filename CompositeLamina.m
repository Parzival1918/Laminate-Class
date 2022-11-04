classdef CompositeLamina < handle
    %COMPOSITELAMINA Summary: Use this class to create a lamina of aligned
    %fibre reinforcement in a matrix.
    %   Detailed explanation goes here
    
    properties (SetAccess = private)
        S (6,6) double {mustBeFloat} = zeros(6)
        C (6,6) double {mustBeFloat} = zeros(6)

        laminaE1 (1,1) double {mustBeFloat} = 0
        laminaE2 (1,1) double {mustBeFloat} = 0
        laminav12 (1,1) double {mustBeFloat} = 0
        laminav21 (1,1) double {mustBeFloat} = 0
        laminav23 (1,1) double {mustBeFloat} = 0
        laminaG12 (1,1) double {mustBeFloat} = 0
        laminaG13 (1,1) double {mustBeFloat} = 0
        laminaK (1,1) double {mustBeFloat} = 0
        laminaIR121 (1,1) double {mustBeFloat} = 0
        laminaIR122 (1,1) double {mustBeFloat} = 0

        n (1,1) double {mustBeFloat} = 0 %Dimensionless constant used in shear lag theory
        s (1,1) double {mustBeFloat} = 0 %Reinforcement aspect ratio
        st (1,1) double {mustBeFloat} = 0 %Stress transfer ratio, applicable for short fibres only
        modifiedE (1,1) double {mustBeFloat} = 0 %E'm value used in the modified shear lag theory
    end

    properties
        compositeComponents CompositeComponents = CompositeComponents('E-glass fibres');
        fibre_fraction (1,1) double {mustBeFloat} = 0.5
        Halpin_Tsai_EqnParameter (1,1) double {mustBeFloat} = 1

        lamina_type char {mustBeMember(lamina_type,{'isotropic',...
                    'orthotropic','transversely isotropic'})} = 'transversely isotropic'
        loading_type char {mustBeMember(loading_type,{'plane stress',...
                    })} = 'plane stress'
        short_fibre_theory char {mustBeMember(short_fibre_theory,{'shear lag theory',...
                    'modified shear lag theory'})} = 'shear lag theory'
    end
    
    methods
        function obj = CompositeLamina(CC)
            %COMPOSITELAMINA Construct an instance of this class
            %   Detailed explanation goes here
            arguments
                CC CompositeComponents = CompositeComponents('E-glass fibres');
            end

            obj.compositeComponents = CC;
        end

        function laminaE1 = get.laminaE1(obj)
            Ef = obj.compositeComponents.reinforcementE1;
            Em = obj.compositeComponents.matrixE1;
            ff = obj.fibre_fraction;

            laminaE1 = Ef*ff + (1-ff)*Em;
        end

        function laminaE2 = get.laminaE2(obj)
            Ef = obj.compositeComponents.reinforcementE1;
            Em = obj.compositeComponents.matrixE1;
            e = obj.Halpin_Tsai_EqnParameter;
            nu = (Ef-Em)/(Ef+e*Em);
            ff = obj.fibre_fraction;

            laminaE2 = Em*((1+e*nu*ff)/(1-nu*ff));
        end

        function laminaG12 = get.laminaG12(obj)
            Gf = obj.compositeComponents.reinforcementG12;
            Gm = obj.compositeComponents.matrixG12;
            e = obj.Halpin_Tsai_EqnParameter;
            nu = (Gf-Gm)/(Gf+e*Gm);
            ff = obj.fibre_fraction;

            laminaG12 = Gm*((1+e*nu*ff)/(1-nu*ff));
        end

        function laminav12 = get.laminav12(obj)
            vf = obj.compositeComponents.reinforcementv12;
            vm = obj.compositeComponents.matrixv12;
            ff = obj.fibre_fraction;

            laminav12 = vf*ff + (1-ff)*vm;
        end

        function laminav21 = get.laminav21(obj)
            laminav21 = obj.laminaE2*(obj.laminav12/obj.laminaE1);
        end

        function laminaK = get.laminaK(obj)
            Kf = obj.compositeComponents.reinforcementK;
            Km = obj.compositeComponents.matrixK;
            ff = obj.fibre_fraction;

            laminaK = 1/(ff/Kf + (1-ff)/Km);
        end

        function laminav23 = get.laminav23(obj)
            v21 = obj.laminav21;
            E2 = obj.laminaE2;
            K = obj.laminaK;

            laminav23 = 1 - v21 - E2/(3*K);
        end

        function n = get.n(obj)
            n = sqrt((2*obj.compositeComponents.matrixE1)/...
                (obj.compositeComponents.reinforcementE1*...
                (1+obj.compositeComponents.matrixv12)*log(1/obj.fibre_fraction)));
        end

        function s = get.s(obj)
            s = obj.compositeComponents.L/obj.compositeComponents.r;
        end

        function st = get.st(obj)
            st = 3/obj.n;
        end

        function modifiedE = get.modifiedE(obj)
            Ef = obj.compositeComponents.reinforcementE1;
            Em = obj.compositeComponents.matrixE1;
            modifiedE = (Ef*(1-sech(obj.n*obj.s))+Em)/2;
        end

        function S = get.S(obj)
            switch obj.lamina_type
                case 'isotropic'
                    E1 = obj.laminaE1;
                    v12 = obj.laminav12;
                    S11 = 1/E1;
                    S12 = -v12/E1;
                    S11S12 = 2*(S11-S12);
                    
                    S = [S11 S12 S12 0 0 0; ...
                         S12 S11 S12 0 0 0; ...
                         S12 S12 S11 0 0 0; ...
                         0 0 0 S11S12 0 0; ...
                         0 0 0 0 S11S12 0; ...
                         0 0 0 0 0 S11S12];
                case 'orthotropic'
                    S = zeros(6);
                case 'transversely isotropic'
                    E1 = obj.laminaE1;
                    E2 = obj.laminaE2;
                    v12 = obj.laminav12;
                    G12 = obj.laminaG12;
                    S11 = 1/E1;
                    S22 = 1/E2;
                    S12 = -v12/E1;
                    S66 = 1/G12;

                    S = [S11 S12 0 0 0 0; ...
                         S12 S22 0 0 0 0; ...
                         0 0 0 0 0 0; ...
                         0 0 0 0 0 0; ...
                         0 0 0 0 0 0; ...
                         0 0 0 0 0 S66];
            end
        end

        function C = get.C(obj)
            S_ = smaller_mat(obj.S);
            C_ = inv(S_);
            C =  bigger_mat(C_);
        end
        
        function strain = apply_stress(obj,stress,angle)
            arguments
                obj CompositeLamina
                stress (3,1) double {mustBeFloat} = [10;0;0]
                angle (1,1) double {mustBeFloat} = 0
            end
            
            S_ = smaller_mat(obj.S);
            U_ = U(angle);
            T_ = T(angle);

            Sbar = inv(U_)*S_*T_;
            
            if strcmp(obj.loading_type,'plane stress')
                strain = Sbar*stress;
            else
                error("Calculation of S with current Lamina properties not possible");
            end
        end

        function Sbar = S_angle(obj, angles)
            arguments
                obj CompositeLamina
                angles (1,:) double {mustBeFloat} = 0
            end
            Sbar = cell([1 length(angles)]);

            for i = 1:length(angles)
                S_ = smaller_mat(obj.S);
                U_ = U(angles(i));
                T_ = T(angles(i));
    
                Sbar{i} = inv(U_)*S_*T_;
            end
        end

        function E = E_angle(obj, angles, is_plot)
            arguments
                obj CompositeLamina
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
                ylabel("Major and Minor Modulus (" + obj.compositeComponents.E_units + ")");
                ylim([min(E,[],"all")-1 max(E,[],"all")+1])
                legend("E1'","E2'");
            end
        end

        function G = G_angle(obj, angles, is_plot)
            arguments
                obj CompositeLamina
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
                ylabel("Shear Modulus (" + obj.compositeComponents.G_units + ")");
                ylim([min(G,[],"all")-1 max(G,[],"all")+1])
                legend("G12'");
            end
        end

        function v = v_angle(obj, angles, is_plot)
            arguments
                obj CompositeLamina
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
                ylabel("Major and Minor Poissons Ratio (" + obj.compositeComponents.v_units + ")");
                ylim([min(v,[],"all")-0.2 max(v,[],"all")+0.2])
                legend("v12'","v21'");
            end
        end

        function IR = IR_angle(obj, angles, is_plot)
            arguments
                obj CompositeLamina
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

        function stress = short_fibre_stress(obj, position, strain)
            %SHORT_FIBRE_STRESS Summary: Calculate the stress at a point x
            %away from the fibre centre when the composite is under certain
            %strain. Only used for short fibres, as for long fibre
            %composites it is assumed stress is constant along the whole
            %fibre.
            %   Explain
            arguments
                obj CompositeLamina
                position (1,1) double {mustBeFloat,mustBeLessThanOrEqualPositive(obj,position)} = 0
                strain (1,1) double {mustBeFloat} = 0.01
            end

            if strcmp(obj.compositeComponents.reinforcement_type,'short fibres')
                if strcmp(obj.short_fibre_theory,'shear lag theory')
                    stress = obj.compositeComponents.reinforcementE1*strain*(1 - ...
                    cosh((obj.n*position)/obj.compositeComponents.r)*sech(obj.n*obj.s));
                else
                    Ef = obj.compositeComponents.reinforcementE1;
                    stress = strain*(Ef - (Ef - obj.modifiedE)* ...
                    cosh((obj.n*position)/obj.compositeComponents.r)*sech(obj.n*obj.s));
                end
            else
                error("Method only available for short fibre composites.")
            end
        end

        function shear = short_fibre_interfacial_shear(obj, position, strain)
            %INTERFACIAL_SHEAR Summary: Calculate the shear stress at a point x
            %away from the fibre centre when the composite is under certain
            %strain. Only used for short fibres, as for long fibre
            %composites it is assumed shear stress is 0.
            %   Explain
            arguments
                obj CompositeLamina
                position (1,1) double {mustBeFloat,mustBeLessThanOrEqualPositive(obj,position)} = 0
                strain (1,1) double {mustBeFloat} = 0.01
            end

            if strcmp(obj.compositeComponents.reinforcement_type,'short fibres')
                if strcmp(obj.short_fibre_theory,'shear lag theory')
                    shear = (obj.n/2)*obj.compositeComponents.reinforcementE1*strain* ...
                    sinh((obj.n*position)/obj.compositeComponents.r)*sech(obj.n*obj.s);
                else
                    Ef = obj.compositeComponents.reinforcementE1;
                    shear = (obj.n/2)*(Ef - obj.modifiedE)*strain* ...
                    sinh((obj.n*position)/obj.compositeComponents.r)*sech(obj.n*obj.s);
                end
            else
                error("Method only available for short fibre composites.")
            end
        end
    end
end

function mustBeLessThanOrEqualPositive(obj,B)
    A = obj.compositeComponents.L;
    if B < -A || B > A
        error("Argument value out of bounds. Position must be in range [-L, L].")
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
    mat = [m(1,1) m(1,2) 0; ...
           m(2,1) m(2,2) 0; ...
           0 0 m(6,6)];
end

function mat = bigger_mat(m)
    mat = [m(1,1) m(1,2) 0 0 0 0; ...
           m(2,1) m(2,2) 0 0 0 0; ...
           0 0 0 0 0 0; ...
           0 0 0 0 0 0; ...
           0 0 0 0 0 0; ...
           0 0 m(3,3) 0 0 0];
end