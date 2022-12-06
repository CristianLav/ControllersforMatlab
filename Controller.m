classdef Controller <handle

    %Controller this SuperClass defines the core part of every single
    %Controller seen in the TAC class, where all the plant and modeling
    %part is define.
    properties(Access=public)

        input= 1; %if the type of signal it have to follow is an step, ramp or an parabolla choose 1,2 or 3 respectibly
        transferFunction ;%this eather way will be calculated by the methods or be an input from the user
%         discreteTransferFunction;
        samplingTime=0 ;%Usefull for using an discrete type of controller let in 0 for continuos controller
        controlability=0; % shows the controlability of a Tf
        observability=0; % shows at the obserbavility of a Tf
        A,B,C,D,E,F; %Spaces stades
        num,den; %num and den
        plantRoots; %roots
        plantPoles; %Poles
        Mp,zeta,Wn,ts; %desing criteria
    end
    properties(Hidden=false)
        desirePolinomial;
        desirePolinomialDiscreto;
        M,Q,T;
        sizeDiserePol;
        systemType;
        S,O; %Controlability and Observability
        As,Bs,Cs,Ds,Es,Fs;%Spaces stades Symbolic
        design_criteria; %this paramaeter must be an containers.map type for the program to recognize the value and the inputs are Mp,ts,zeta and wn
        %       Mp,zeta,Wn,ts;
    end
    methods
        %Constructor-------------------------------------------------------
        function obj = Controller(input_,design_criteria_,samplingTime_)
            obj.input=input_;
            obj.design_criteria=design_criteria_;
            obj.samplingTime=samplingTime_;
            obj.Obtain_zeta_and_wn_from();

        end

        %ObtenerMatricez A,B,C,D------------------------publico
        function getSpaceStadeSymsModel(obj,Estados,derivadas_estados,Entrada,Salida,Perturbaciones)
            %Sacar las Funciones haciendo uso del comando Jacobian
            obj.As=jacobian(derivadas_estados,Estados);
            obj.Bs=jacobian(derivadas_estados,Entrada);
            obj.Cs=jacobian(Salida,Estados);
            obj.Ds=jacobian(Salida,Entrada);
            obj.Es=jacobian(derivadas_estados,Perturbaciones(1));
            obj.Fs=jacobian(derivadas_estados,Perturbaciones(2));
        end

        %Complete all stats fi the trasfer function is an input
        function initializeWithTf(obj)
            obj.getTransferFuncionNumDen();
            obj.setABCDMatrixwithtf();
            obj.getControlability();
            obj.getObservability();
        end

        %Complete all stats fi the Equations as is an input
        function initializeWithSpaceStades(obj)
            space_stade=ss(obj.A,obj.B,obj.C,obj.D);
            obj.transferFunction = tf(space_stade);
            obj.getTransferFuncionNumDen();
            obj.getControlability();
            obj.getObservability();
        end

    end
    methods
        %Set and Get Transfer function----------------------------------
        %set Trasnfer function
        function  setTransferfunction(obj,transferFunction_)
            obj.transferFunction=transferFunction_;
        end
        %get Trasnfer function
        function transferFunction = getTransferfunction(obj)
            transferFunction=obj.transferFunction;
        end
    end
    methods(Access = private)
        %calculate zeta if overshoot is a control design parameter
        function  get_zeta_by_overshoot(obj)
            obj.zeta = abs(log(obj.Mp/100)) / sqrt(pi^2 + log(obj.Mp/100)^2);
        end
        %calculate zeta if ts and Wn is a control design parameter
        function  get_zeta_by_ts_and_wn(obj)
            obj.zeta = 4 / (obj.ts * obj.Wn);
        end
        %calculate Wn if ts and zeta is a control design parameter
        function get_wn_by_ts_and_zeta(obj)
            obj.Wn = 4 / (obj.ts * obj.zeta);
        end       
        %Get zeta and Wn for the control calculous------------
        function Obtain_zeta_and_wn_from(obj)

            if isKey(obj.design_criteria, 'Mp')

                obj.Mp = obj.design_criteria('Mp');
                obj.get_zeta_by_overshoot();

                if isKey(obj.design_criteria, 'ts')

                    obj.ts = obj.design_criteria('ts');
                    obj.get_wn_by_ts_and_zeta();

                elseif isKey(obj.design_criteria, 'Wn')

                    obj.Wn = obj.design_criteria('Wn');
                end

            elseif isKey(obj.design_criteria, 'ts')

                obj.ts = obj.design_criteria('ts');

                if isKey(obj.design_criteria, 'Wn')

                    obj.Wn = obj.design_criteria('Wn');
                    obj.get_zeta_by_ts_and_wn();

                elseif isKey(obj.design_criteria, 'zeta')

                    obj.zeta = obj.design_criteria('zeta');
                    obj.get_wn_by_ts_and_zeta();
                end

            elseif isKey(obj.design_criteria, 'Wn') && isKey(obj.design_criteria, 'zeta')
                obj.Wn = obj.design_criteria('Wn');
                obj.zeta = obj.design_criteria('zeta');
            end

        end
    end
    methods(Access=protected)
         %set A B C D Matrix from TrasnferFuncion------------
        function setABCDMatrixwithtf(obj)
            [A1,B1,C1,D1]=tf2ss(obj.num,obj.den);
            obj.A=A1;
            obj.B=B1;
            obj.C=C1;
            obj.D=D1;
        end

        %get numerator and denominator of trasnfer funcion --------------
        function getTransferFuncionNumDen(obj)
            [ numft , denft ] = tfdata(obj.transferFunction);
            obj.num =cell2mat(numft);
            obj.den =cell2mat(denft);
            obj.plantRoots=roots(obj.den);
            obj.plantPoles=roots(obj.num);
        end

        %Generation of  continuous Desire Polinomial--------
        function getDesirePolinomial(obj)
            %Calcular Polinomio Deseado
            Beta=5;
            obj.desirePolinomial=[1 2*obj.zeta*obj.Wn obj.Wn^2];
            %calcular si es requerido agregar polos no dominantes para igualar el orden
            PolSize=size(obj.desirePolinomial,2)-1;
            if obj.sizeDiserePol > PolSize
                for i=PolSize:obj.sizeDiserePol-1
                    obj.desirePolinomial=conv(obj.desirePolinomial,[1 Beta*obj.zeta*obj.Wn]);
                end
            end
        end
       
        %Generation of  Discrete Desire Polinomial--------
        function getDesireDiscretePolinomial(obj)
            m=exp(- obj.samplingTime*obj.zeta*obj.Wn);
            phi=obj.samplingTime*obj.Wn;

            obj.desirePolinomialDiscreto=[1 -2*m*cos(phi) m^2];

            MatrizSize=obj.sizeDiserePol ;
            PolSize=size(obj.desirePolinomialDiscreto,2)-1;
            if MatrizSize > PolSize
                for i=PolSize:MatrizSize-1
                    obj.desirePolinomialDiscreto=conv(obj.desirePolinomialDiscreto,[1 -0.05]);
                end
            end
        end
        
        %Get the controlability-------------------------------------
        
        function getControlability(obj)
            syms s;
            MatrixRange=size(obj.A,1);
            obj.S=zeros(MatrixRange);

            for i=0:MatrixRange - 1
                obj.S(:, i+1)=[obj.A^i*obj.B];
            end

            %S=[B A*B A^2*B]
            %Sacar los Subindices de la Determinante  SI-A
            SI=eye(size(obj.A))*s;
            determinateM=det(SI-obj.A);
            determinateM=collect(determinateM,s);
            determinateM=eval(determinateM);
            subexpresionDetermianteM =double(coeffs(determinateM,s));

            %En Caso de inconsistencia
            PolSize=size(subexpresionDetermianteM,2);
            if size(subexpresionDetermianteM) < MatrixRange+1
                for i=PolSize:MatrixRange
                    subexpresionDetermianteM=[0,subexpresionDetermianteM];
                end
            end
            %Crear la Matriz M
            obj.M=zeros(MatrixRange);
            for i=1:MatrixRange
                subexpresionDetermianteM(1)=[]  ;
                obj.M(i,:)=[(subexpresionDetermianteM)];
                subexpresionDetermianteM(MatrixRange+1)=[0];
            end
            %Matriz T de trasformacion
            obj.T=obj.S*obj.M;
            obj.controlability=det(obj.S);
        end
        
        %Get the Observability-------------------------------
        function getObservability(obj)
            syms s;
            MatrixRange=size( obj.A,1);
            obj.O=zeros(MatrixRange);
            for i=0:MatrixRange - 1
                obj.O(i+1,:)=[obj.C * obj.A^i];
            end
            obj.Q=inv( obj.M* obj.O);
            obj.observability=det(obj.O);
        end

    end
end