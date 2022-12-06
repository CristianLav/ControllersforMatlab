classdef ServoController < Controller
    %ServoController Calcula un servosistema para entrada escalon y
    %retroalimentacino de estados

    properties
        K_acker;
        L_acker;
        G,H;
    end

    properties(Access = private)
        Ahat;
        Bhat;
        AOb;
        desireObservatorPolinomial;
        desireDiscreteObservatorPolinomial;
    end

    methods
        %Constructor of ServoController class-------------------
        function obj =ServoController(input_,DESIGN_CRITERIA_,samplingTime_)
            obj  = obj@Controller(input_,DESIGN_CRITERIA_,samplingTime_);

            %obj.sizeKvector();
        end

        %Calculate Everything-------------------
        function computeServo(obj)
            obj.defineDesirePolinomialSize();
            obj.stablishAhatBhat();
            if obj.samplingTime == 0
                obj.getDesirePolinomial();
                obj.getDesirePolinomialforObservator();
            else
                obj.getDesireDiscretePolinomial();
                obj.getDesireDiscreteObservatorPolinomial();
            end

            obj.getKforAckerMethod();
            obj.getLforObserverAckerMan();
        end

    end

    methods(Access = private)

        %Calculate the hat matriz-------------------
        function stablishAhatBhat(obj)
            if obj.samplingTime == 0
                obj.AOb=obj.A;
                if obj.input ==0
                    obj.Ahat=obj.A;
                    obj.Bhat=obj.B;
                    return;
                elseif obj.input == 1
                    obj.Ahat=[obj.A zeros(size(obj.A,1) , 1); -obj.C 0];
                    obj.Bhat=[obj.B;0];

                    return;
                end
                disp("Otros casos siguen sin ser programados")
            else
                [G1,H1]=c2d(obj.A,obj.B,obj.samplingTime);
                obj.G=G1;
                obj.H=H1;
                obj.AOb=obj.G;
                if obj.input == 0

                    obj.Ahat=obj.G;
                    obj.Bhat=obj.H;
                    return;
                elseif obj.input == 1

                    obj.Ahat=[obj.G zeros(size(obj.G,2),1);obj.C*obj.G 1];
                    obj.Bhat=[obj.H;obj.C*obj.H];

                    return;
                end
                disp("Otros casos siguen sin ser programados")
            end
        end

        %Size of the diserePolinomial-------------------
        function defineDesirePolinomialSize(obj)
            obj.sizeDiserePol = size(obj.A,1) + obj.input + 1;
        end

        %Calculate  the Discrete diserePolinomial for observator-------------------
        function getDesireDiscreteObservatorPolinomial(obj)

            m=exp(-obj.samplingTime *obj.zeta*(obj.Wn*10));
            phi= obj.samplingTime*(obj.Wn*10);

            obj.desireDiscreteObservatorPolinomial=[1 -2*m*cos(phi) m^2];

            MatrizSize=size(obj.G,2) ;
            PolSize=size(obj.desireDiscreteObservatorPolinomial,2)-1;
            if MatrizSize > PolSize
                for i=PolSize:MatrizSize-1
                    obj.desireDiscreteObservatorPolinomial=conv(obj.desireDiscreteObservatorPolinomial,[1 -0.05]);
                end
            end
        end

        %Calculate  the diserePolinomial for observator-------------------
        function getDesirePolinomialforObservator(obj)

            %Calcular Polinomio Deseado
            Beta=5;
            obj.desireObservatorPolinomial=[1 2*obj.zeta*(10*obj.Wn) (10*obj.Wn)^2];
            %calcular si es requerido agregar polos no dominantes para igualar el orden
            MatrizSize=size(obj.A,1);
            PolSize=size(obj.desireObservatorPolinomial,2)-1;
            if MatrizSize > PolSize
                for i=PolSize:MatrizSize-1
                    obj.desireObservatorPolinomial=conv(obj.desireObservatorPolinomial,[1 Beta*obj.zeta*(10*obj.Wn)]);
                end
            end
        end

        %Calculate Controlator vector-------------------
        function getKforAckerMethod(obj)
            MatrixRange=size(obj.Ahat,1);
            S=[];
            for i=0:MatrixRange - 1
                S=[ S ,obj.Ahat^i*obj.Bhat];
            end

            RangeA=size(obj.Ahat,2);
            ZerosAckerMatriz=zeros(1,RangeA);
            ZerosAckerMatriz(RangeA)=1;
            %Calcular Polinomio Deseado

            %Calcular los Valores de K
            if obj.samplingTime == 0
                PhiA=0;
                for i=1:RangeA
                    PhiA= PhiA+ obj.desirePolinomial(i)*obj.Ahat^(RangeA-(i-1)) ;
                end
                I=eye(RangeA);
                PhiA=PhiA+obj.desirePolinomial(RangeA+1)*I;
                obj.K_acker=ZerosAckerMatriz*inv(S)*PhiA;
                return;
            end

            PhiA=0;
            for i=1:RangeA
                PhiA= PhiA+ obj.desirePolinomialDiscreto(i)*obj.Ahat^(RangeA-(i-1)) ;
            end
            I=eye(RangeA);
            PhiA=PhiA+obj.desirePolinomialDiscreto(RangeA+1)*I;
            obj.K_acker=ZerosAckerMatriz*inv(S)*PhiA;
            obj.K_acker=ZerosAckerMatriz*inv(S)*PhiA;
        end

        %Calculate L for observator-------------------
        function  getLforObserverAckerMan(obj)
            MatrixRange=size(obj.AOb,1);
            O=[];
            for i=0:MatrixRange - 1
                O=[ O ;obj.C*obj.AOb^i];
            end

            RangeA=size(obj.AOb,1);
            ZerosAckerMatriz=zeros(1,RangeA);
            ZerosAckerMatriz(RangeA)=1;
            ZerosAckerMatriz=ZerosAckerMatriz';

            %Calcular los Valores de L
            PhiA=0;
            if obj.samplingTime == 0
                for i=1:RangeA
                    PhiA= PhiA+ obj.desireObservatorPolinomial(i)*obj.AOb^(RangeA-(i-1)) ;
                end
                I=eye(RangeA);
                PhiA=PhiA+obj.desireObservatorPolinomial(RangeA+1)*I;
            else
                for i=1:RangeA
                    PhiA= PhiA+ obj.desireDiscreteObservatorPolinomial(i)*obj.AOb^(RangeA-(i-1)) ;
                end
                I=eye(RangeA);
                PhiA=PhiA+obj.desireDiscreteObservatorPolinomial(RangeA+1)*I;
            end
           
            obj.L_acker = PhiA*inv(O)*ZerosAckerMatriz;
            
        end

    end
end