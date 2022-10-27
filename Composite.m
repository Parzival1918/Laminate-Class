classdef Composite < handle
    %COMPOSITE Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (SetAccess = private)
        ratio (1,1) {mustBeFloat} = 1e-3
        Halpin_Tsai_EqnParameter (1,1) {mustBeFloat} = 1
        compositeE1 (1,1) {mustBeFloat}
        compositeE2 (1,1) {mustBeFloat}
        compositev12 (1,1) {mustBeFloat}
        compositev21 (1,1) {mustBeFloat}
        compositev23 (1,1) {mustBeFloat}
        compositeG12 (1,1) {mustBeFloat}
        compositeG23 (1,1) {mustBeFloat}
        compositeB (1,1) {mustBeFloat}
    end
    
    properties
        FibreE (1,1) {mustBeFloat}
        Fibrev (1,1) {mustBeFloat}
        FibreG (1,1) {mustBeFloat}
        FibreFraction (1,1) {mustBeFloat,mustBe01} = 0.5
        MatrixE (1,1) {mustBeFloat}
        Matrixv (1,1) {mustBeFloat}
        MatrixG (1,1) {mustBeFloat}

        stress_units char {mustBeMember(stress_units,{'GPa','MPa','Pa'})} ...
            = 'MPa'
        stiffness_units char {mustBeMember(stiffness_units,{'GPa','MPa','Pa'})} ...
            = 'GPa'
    end
    
    methods
        function ratio = get.ratio(obj)
            ratio = units_difference(obj);
        end

        function compositeE1 = get.compositeE1(obj)
            if ne(obj.FibreE,0) && ne(obj.Fibrev,0) && ne(obj.FibreG,0) && ...
                    ne(obj.MatrixE,0) && ne(obj.Matrixv,0) && ne(obj.MatrixG,0)
                compositeE1 = obj.FibreFraction*obj.FibreE + ...
                    (1-obj.FibreFraction)*obj.MatrixE;
            end
        end

        function compositeE2 = get.compositeE2(obj)
            if ne(obj.FibreE,0) && ne(obj.Fibrev,0) && ne(obj.FibreG,0) && ...
                    ne(obj.MatrixE,0) && ne(obj.Matrixv,0) && ne(obj.MatrixG,0)
                nu = (obj.FibreE - obj.MatrixE)/(obj.FibreE - obj.Halpin_Tsai_EqnParameter*obj.MatrixE);
                compositeE2 = obj.MatrixE*((1 + obj.Halpin_Tsai_EqnParameter*nu*obj.FibreFraction)/(1-nu*obj.FibreFraction));
            end
        end

        function compositeG12 = get.compositeG12(obj)
            if ne(obj.FibreE,0) && ne(obj.Fibrev,0) && ne(obj.FibreG,0) && ...
                    ne(obj.MatrixE,0) && ne(obj.Matrixv,0) && ne(obj.MatrixG,0)
                nu = (obj.FibreG - obj.MatrixG)/(obj.FibreG - obj.Halpin_Tsai_EqnParameter*obj.MatrixG);
                compositeG12 = obj.MatrixG*((1 + obj.Halpin_Tsai_EqnParameter*nu*obj.FibreFraction)/(1-nu*obj.FibreFraction));
            end
        end

        function compositev12 = get.compositev12(obj)
            if ne(obj.FibreE,0) && ne(obj.Fibrev,0) && ne(obj.FibreG,0) && ...
                    ne(obj.MatrixE,0) && ne(obj.Matrixv,0) && ne(obj.MatrixG,0)
                compositev12 = obj.FibreFraction*obj.Fibrev + ...
                    (1-obj.FibreFraction)*obj.Matrixv;
            end
        end

        function compositev21 = get.compositev21(obj)
            if ne(obj.FibreE,0) && ne(obj.Fibrev,0) && ne(obj.FibreG,0) && ...
                    ne(obj.MatrixE,0) && ne(obj.Matrixv,0) && ne(obj.MatrixG,0)
                compositev21 = obj.compositeE2*(obj.compositev12/obj.compositeE1);
            end
        end

        function compositeB = get.compositeB(obj)
            if ne(obj.FibreE,0) && ne(obj.Fibrev,0) && ne(obj.FibreG,0) && ...
                    ne(obj.MatrixE,0) && ne(obj.Matrixv,0) && ne(obj.MatrixG,0)
                Bf = obj.FibreE/(3*(1-2*obj.Fibrev));
                Bm = obj.MatrixE/(3*(1-2*obj.Matrixv));
                compositeB = 1/(obj.FibreFraction/Bf + (1-obj.FibreFraction)/Bm);
            end
        end

        function compositev23 = get.compositev23(obj)
            if ne(obj.FibreE,0) && ne(obj.Fibrev,0) && ne(obj.FibreG,0) && ...
                    ne(obj.MatrixE,0) && ne(obj.Matrixv,0) && ne(obj.MatrixG,0)
                compositev23 = 1 - obj.compositev21 - obj.compositeE2/(3*obj.compositeB);
            end
        end

        function compositeG23 = get.compositeG23(obj)
            if ne(obj.FibreE,0) && ne(obj.Fibrev,0) && ne(obj.FibreG,0) && ...
                    ne(obj.MatrixE,0) && ne(obj.Matrixv,0) && ne(obj.MatrixG,0)
                compositeG23 = obj.compositeE2/(2*(1 + obj.compositev23));
            end
        end

        function obj = update_Fibres(obj)
            if ne(obj.FibreE,0) && ne(obj.Fibrev,0)
                obj.FibreG = obj.FibreE/(2*(1+obj.Fibrev));
            elseif ne(obj.FibreG,0) && ne(obj.Fibrev,0)
                obj.FibreE = obj.FibreG*(2*(1+obj.Fibrev));
            elseif ne(obj.FibreG,0) && ne(obj.FibreE,0)
                obj.Fibrev = (obj.FibreE/(2*obj.FibreG)) - 1;
            end
        end

        function obj = update_Matrix(obj)
            if ne(obj.MatrixE,0) && ne(obj.Matrixv,0)
                obj.MatrixG = obj.MatrixE/(2*(1+obj.Matrixv));
            elseif ne(obj.MatrixG,0) && ne(obj.Matrixv,0)
                obj.MatrixE = obj.MatrixG*(2*(1+obj.Matrixv));
            elseif ne(obj.MatrixG,0) && ne(obj.MatrixE,0)
                obj.Matrixv = (obj.MatrixE/(2*obj.MatrixG)) - 1;
            end
        end

        function obj = update_Composite(obj)
            if ne(obj.FibreE,0) && ne(obj.Fibrev,0)
                obj.FibreG = obj.FibreE/(2*(1+obj.Fibrev));
            elseif ne(obj.FibreG,0) && ne(obj.Fibrev,0)
                obj.FibreE = obj.FibreG*(2*(1+obj.Fibrev));
            elseif ne(obj.FibreG,0) && ne(obj.FibreE,0)
                obj.Fibrev = (obj.FibreE/(2*obj.FibreG)) - 1;
            end

            if ne(obj.MatrixE,0) && ne(obj.Matrixv,0)
                obj.MatrixG = obj.MatrixE/(2*(1+obj.Matrixv));
            elseif ne(obj.MatrixG,0) && ne(obj.Matrixv,0)
                obj.MatrixE = obj.MatrixG*(2*(1+obj.Matrixv));
            elseif ne(obj.MatrixG,0) && ne(obj.MatrixE,0)
                obj.Matrixv = (obj.MatrixE/(2*obj.MatrixG)) - 1;
            end
        end

        function obj = Composite()
        end
    end
end

function mustBe01(val)
    if val < 0 || val > 1
        error("Invalid FibreFraction value. Must be in [0,1] range.")
    end
end

function multiplyer = units_difference(obj)
    switch obj.stress_units
        case 'GPa'
            mult1 = 10^9;
        case 'MPa'
            mult1 = 10^6;
        otherwise
            mult1 = 1;
    end

    switch obj.stiffness_units
        case 'GPa'
            mult2 = 10^9;
        case 'MPa'
            mult2 = 10^6;
        otherwise
            mult2 = 1;
    end

    multiplyer = mult1/mult2;
end