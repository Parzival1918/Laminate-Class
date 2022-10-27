classdef CompositeLamina < handle
    %COMPOSITELAMINA Summary: Use this class to create a lamina of aligned
    %fibre reinforcement in a matrix.
    %   Detailed explanation goes here
    
    properties (SetAccess = private)
        S
        C
    end

    properties
        compositeComponents CompositeComponents = CompositeComponents(100,10);
    end
    
    methods
        function obj = CompositeLamina(CC)
            %COMPOSITELAMINA Construct an instance of this class
            %   Detailed explanation goes here
            arguments
                CC CompositeComponents = CompositeComponents(100,10);
            end

            obj.compositeComponents = CC;
        end
        
        function outputArg = method1(obj,inputArg)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            outputArg = obj.Property1 + inputArg;
        end
    end
end

function val = U(angle)
    val = [cosd(angle)^2 sind(angle)^2 cosd(angle)*sind(angle); ...
           sind(angle)^2 cosd(angle)^2 -cosd(angle)*sind(angle); ...
           -2*cosd(angle)*sind(angle) 2*cosd(angle)*sind(angle) cosd(angle)^2 - sind(angle)^2];
end