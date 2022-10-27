classdef Laminate < handle
    %LAMINATE Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (SetAccess = private)
        laminae (1,:) Lamina
        angles (1,:) double
        thickness (1,:) double
        stress (3,1) {mustBeFloat} = zeros([3 1])
        stressAngle (1,1) {mustBeFloat} = 0
        strain (3,1) {mustBeFloat} = zeros([3 1])
        C (3,3) {mustBeFloat} = zeros(3)
        S (3,3) {mustBeFloat} = zeros(3)

        compositeE1 (1,1) {mustBeFloat}
        compositeE2 (1,1) {mustBeFloat}
        compositev12 (1,1) {mustBeFloat}
        compositev21 (1,1) {mustBeFloat}
        compositeG12 (1,1) {mustBeFloat}
        interactionRatio121 (1,1) {mustBeFloat}
    end

    properties
        stress_units char {mustBeMember(stress_units,{'GPa','MPa','Pa'})} ...
            = 'MPa'
        stiffness_units char {mustBeMember(stiffness_units,{'GPa','MPa','Pa'})} ...
            = 'GPa'
    end
    
    methods
        function obj = Laminate(lamina, angle, t)
            arguments
                lamina = NaN
                angle = NaN
                t = 1
            end

            if ne(lamina,NaN) && ne(angle,NaN)
                obj.add_lamina(lamina,angle, t);
            end
        end

        function C = get.C(obj)
            C = zeros(3);
            sumtk = sum(obj.thickness);
            for i = 1:3
                for j = 1:3
                    sumCij = 0;
                    for k = 1:length(obj.laminae)
                        laminak = obj.laminae(k);
                        tk = obj.thickness(k);
                        Ck = [laminak.C(1,1), laminak.C(1,2), 0; ...
                              laminak.C(2,1), laminak.C(2,2), 0; ...
                              0, 0, laminak.C(6,6)];
                        angleZ = obj.angles(k) - obj.stressAngle;
                        cossin = [cosd(angleZ)^2 sind(angleZ)^2 cosd(angleZ)*sind(angleZ); ...
                          sind(angleZ)^2 cosd(angleZ)^2 -cosd(angleZ)*sind(angleZ); ...
                          -2*cosd(angleZ)*sind(angleZ) 2*cosd(angleZ)*sind(angleZ) cosd(angleZ)^2 - sind(angleZ)^2];
                        Cij = cossin*Ck*transpose(cossin);
                        
                        sumCij = sumCij + (Cij(i,j)*tk);
                    end
                    tempC = sumCij/sumtk;
                    C(i,j) = tempC;
                end
            end
        end

        function S = get.S(obj)
            S = zeros(3);
            sumtk = sum(obj.thickness);
            for i = 1:3
                for j = 1:3
                    sumSij = 0;
                    for k = 1:length(obj.laminae)
                        laminak = obj.laminae(k);
                        tk = obj.thickness(k);
                        Ck = [laminak.S(1,1), laminak.S(1,2), 0; ...
                              laminak.S(2,1), laminak.S(2,2), 0; ...
                              0, 0, laminak.S(6,6)];
                        angleZ = obj.angles(k) - obj.stressAngle;
                        cossin = [cosd(angleZ)^2 sind(angleZ)^2 cosd(angleZ)*sind(angleZ); ...
                          sind(angleZ)^2 cosd(angleZ)^2 -cosd(angleZ)*sind(angleZ); ...
                          -2*cosd(angleZ)*sind(angleZ) 2*cosd(angleZ)*sind(angleZ) cosd(angleZ)^2 - sind(angleZ)^2];
                        Sij = cossin*Ck*transpose(cossin);
                        
                        sumSij = sumSij + (Sij(i,j)*tk);
                    end
                    tempS = sumSij/sumtk;
                    S(i,j) = tempS;
                end
            end
        end

        function compositeE1 = get.compositeE1(obj)
            if ne(obj.S(1,1),0)
                compositeE1 = 1/obj.S(1,1);
            else
                compositeE1 = 0;
            end
        end

        function compositeE2 = get.compositeE2(obj)
            if ne(obj.S(2,2),0)
                compositeE2 = 1/obj.S(2,2);
            else
                compositeE2 = 0;
            end
        end

        function compositeG12 = get.compositeG12(obj)
            if ne(obj.S(3,3),0)
                compositeG12 = 1/obj.S(3,3);
            else
                compositeG12 = 0;
            end
        end

        function compositev12 = get.compositev12(obj)
            if ne(obj.S(1,2),0)
                compositev12 = -obj.compositeE1*obj.S(1,2);
            else
                compositev12 = 0;
            end
        end

        function interactionRatio121 = get.interactionRatio121(obj)
            if ne(obj.S(1,3),0)
                interactionRatio121 = obj.compositeE1*obj.S(1,3);
            else
                interactionRatio121 = 0;
            end
        end

        function obj = add_lamina(obj, lamina, angle, t)
            arguments
                obj
                lamina Lamina
                angle double = 0
                t = 1
            end

            if class(lamina) == "Lamina"
                obj.laminae(end+1) = lamina;
                if length(obj.laminae) == 1
                    obj.angles(end+1) = 0;
                else
                    obj.angles(end+1) = angle;
                end
                obj.thickness(end+1) = t;
            else
                error("Wrong data type. Input arguments must be of class Lamina.");
            end
        end

        function plot_YMs_angle(obj,angles)
            arguments
                obj
                angles (1,:) {mustBeFloat} = linspace(0,90,30)
            end
            
            originalAngle = obj.stressAngle;
            YMsAngle = zeros([length(angles), 2]);
            pos = 1;
            for angle = angles
                obj.stressAngle = angle;
                YMsAngle(pos,:) = [obj.compositeE1, obj.compositeE2];
                pos = pos + 1;
            end
            
            figure;
            plot(angles,YMsAngle);
            xlabel("Angle (º)");
            ylabel("Major and Minor Modulus (" + obj.stiffness_units + ")");
            ylim([min(YMsAngle,[],"all")-1 max(YMsAngle,[],"all")+1])
            legend("E1","E2");

            obj.stressAngle = originalAngle;
        end

        function plot_G_angle(obj,angles)
            arguments
                obj
                angles (1,:) {mustBeFloat} = linspace(0,90,30)
            end
            
            originalAngle = obj.stressAngle;
            GAngle = zeros([length(angles), 1]);
            pos = 1;
            for angle = angles
                obj.stressAngle = angle;
                GAngle(pos,:) = obj.compositeG12;
                pos = pos + 1;
            end
            
            figure;
            plot(angles,GAngle);
            xlabel("Angle (º)");
            ylabel("G12 (" + obj.stiffness_units + ")");
            ylim([min(GAngle,[],"all")-1 max(GAngle,[],"all")+1])

            obj.stressAngle = originalAngle;
        end

        function plot_v_angle(obj,angles)
            arguments
                obj
                angles (1,:) {mustBeFloat} = linspace(0,90,30)
            end
            
            originalAngle = obj.stressAngle;
            vAngle = zeros([length(angles), 1]);
            pos = 1;
            for angle = angles
                obj.stressAngle = angle;
                vAngle(pos,:) = obj.compositev12;
                pos = pos + 1;
            end
            
            figure;
            plot(angles,vAngle);
            xlabel("Angle (º)");
            ylabel("Major Poisson's Ratio");
            legend("v12");
            ylim([min(vAngle,[],"all")-1 max(vAngle,[],"all")+1])

            obj.stressAngle = originalAngle;
        end

        function plot_interaction_ratio_angle(obj,angles)
            arguments
                obj
                angles (1,:) {mustBeFloat} = linspace(0,90,30)
            end
            
            originalAngle = obj.stressAngle;
            IRAngle = zeros([length(angles), 1]);
            pos = 1;
            for angle = angles
                obj.stressAngle = angle;
                IRAngle(pos,:) = obj.interactionRatio121;
                pos = pos + 1;
            end
            
            figure;
            plot(angles,IRAngle);
            xlabel("Angle (º)");
            ylabel("Interaction Ratios");
            legend("n121");
            ylim([min(IRAngle,[],"all")-1 max(IRAngle,[],"all")+1])

            obj.stressAngle = originalAngle;
        end
    end
end

