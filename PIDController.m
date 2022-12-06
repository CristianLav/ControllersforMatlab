classdef PIDController <Controller
    %PIDController is the class that take care of creating the PID constands
    %   Detailed explanation goes here

    properties       
        numPID;
        denPID;
        PIDValues;
        PIDType;
        PIDTransferFunction;
        numOfIntegrators;
    end
   

    methods
        %Constructor of PIDController class-------------------
        function obj =PIDController(input_,DESIGN_CRITERIA_,samplingTime_)
            obj  = obj@Controller(input_,DESIGN_CRITERIA_,samplingTime_);
        end
        
        %literally makes everything for you, thanks ;)
        function startComputesPID(obj)
            if obj.samplingTime == 0
                obj.getPIDtype();
                obj.getDesirePolinomial();
                obj.getPIDvalues();
                return;
            end
            disp('No programo todavia la parte discreta gg')
        end

        %Generation PID type---------------------------------
        function getPIDtype(obj)
            syms s;
            obj.systemType=size(find(obj.plantRoots==0),1);
            %Definir el numero de integradores y el orden de la funcion de
            %trasnferencia

            obj.numOfIntegrators= obj.input-obj.systemType;
            ordenDenominador=size(obj.den,2)-1;

            %Definir la variabes del PID y calcular la funcion de trasferencia que se
            %requiere
            syms kd4 kd3 kd2 kd1 kp ki1 ki2 ki3 ki4
            definirtipoPID=[kd4 kd3 kd2 kd1 kp ki1 ki2 ki3 ki4];
            posicionfinal=5+obj.numOfIntegrators;
            posicionfinalInicial=posicionfinal-(ordenDenominador+obj.numOfIntegrators)+1;
            definirtipoPID=definirtipoPID(posicionfinalInicial:posicionfinal);
            obj.sizeDiserePol=size(definirtipoPID,2);

            %arma la funcion de trasferencia de manera simbolica
            numeradorPID=0;
            for i=1 :obj.sizeDiserePol
                numeradorPID=numeradorPID+ definirtipoPID(i)*s^(obj.sizeDiserePol-i);
            end
            obj.PIDTransferFunction=numeradorPID/s^(obj.numOfIntegrators);

        end
            
        %Computes the values of a PID
        function getPIDvalues(obj)
            
            %Obtener el la funcion de trasferencia de manera simbolica
            syms s;
            [numeradorPID, denominadorPID] = numden(obj.PIDTransferFunction);
            symbolicNumerator=0;
            sizenumerator=size(obj.num,2)-1;
            for i=0:sizenumerator
                symbolicNumerator=symbolicNumerator+obj.num(i+1)*s^(sizenumerator-i);
            end
            symbolicDenominator=0;
            sizedenominator=size(obj.den,2)-1;
            for i=0:sizedenominator
                symbolicDenominator=symbolicDenominator+obj.den(i+1)*s^(sizedenominator-i);
            end
            %Sacar la funcion de trasnferencia en Lazo cerrado para hacer igualacion de
            %polinomio
            LazoCerradoPID=symbolicDenominator*denominadorPID + numeradorPID*symbolicNumerator; %denominador de la funcion
            subexpresionsofPID=collect(LazoCerradoPID,s);
            subexpresionsofPID=children(subexpresionsofPID);
            subexpresionsofPID=[subexpresionsofPID{:}];
            OrdenPID=size(subexpresionsofPID,2)-1;
            subexpresionsofPID=subexpresionsofPID/(subexpresionsofPID(1)/s^OrdenPID);

            for i=0:OrdenPID
                subexpresionsofPID(i+1)=subexpresionsofPID(i+1)/s^(OrdenPID-i);
            end
            %Solociona los valores del PID
            obj.PIDValues = vpasolve(subexpresionsofPID== obj.desirePolinomial);
        end

    end
end