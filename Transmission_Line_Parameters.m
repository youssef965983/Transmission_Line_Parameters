

function Transmission_Line_Parameters
    fprintf('\nThis program calculate Transmission Line Parameters (R,L,C)/phase \n\n');
    
    %%%%%%%%%%%%%%%%% get parameters %%%%%%%%%%%%%%%%%%%%%%%%%
    conductor_resistivity = get_conductor_resistivity_func();
    conductor_length = get_conductor_length_func();
    conductor_diameter = get_conductor_diameter_func();
    freq = input('Enter The operating frequency: ');
  
    %%%%%%%%%%%%%%%%%% Resistance calculation %%%%%%%%%%%%%%%%

    dc_resistance = calc_dc_resistance_func(conductor_resistivity,conductor_length,conductor_diameter);
    ac_resistance = calc_ac_resistance_func(conductor_resistivity,conductor_length,conductor_diameter);
    
    %%%%%%%%%%%%%%%%%% get D_equivalent   %%%%%%%%%%%%%%%%%%%%
    D_eq = get_equivalent_distance();

    %%%%%%%%%%%%%%%%%% Inductance calculation %%%%%%%%%%%%%%%%
    inductance_per_phase = calc_inductance_per_phase_func(conductor_diameter,conductor_length,D_eq);

    %%%%%%%%%%%%%%%%%% Capacitance calculation %%%%%%%%%%%%%%%%
    capacitance_per_phase = calc_capacitance_per_phase_func(conductor_diameter,conductor_length,D_eq);

    %%%%%%%%%%%%%%%%%% show results %%%%%%%%%%%%%%%%%%%%%%%%%%%
    fprintf('DC_Resistance/phase is %d ohm\n',dc_resistance);
    fprintf('AC_Resistance/phase is %d ohm\n',ac_resistance);
    fprintf('Inductance/phase is %d Henry\n',inductance_per_phase);
    fprintf('Capacitance/phase is %d Farad\n',capacitance_per_phase);
 
    %%%%%%%%%%%%%%%%%%% ABCD parameters calculations%%%%%%%%%%%%%%%%%%%%%%%
    Z=complex(ac_resistance,2*pi*freq*inductance_per_phase);
    Y=complex(0,2*pi*freq*capacitance_per_phase);
    [A,B,C,D]=TL_Model(conductor_length,Z,Y);
    
    %%%%%%%%%%%%%%%%%% show ABCD results  %%%%%%%%%%%%%%%%%%%%%%
    fprintf('A parameter = %d +j %d \n',real(A),imag(A));
    fprintf('B parameter = %d +j %d \n',real(B),imag(B));
    fprintf('C parameter = %d +j %d \n',real(C),imag(C));
    fprintf('D parameter = %d +j %d \n',real(D),imag(D));
    
    %%%%%%%%%%%%%%%%%% Transmission Line Performance %%%%%%%%%%%
    Receiving_Voltage = input('Enter The Receiving End Voltage in KV: ') /sqrt(3);
    CaseChoice = 0;
    
    while(CaseChoice ~= 3)
        fprintf('please select: \n');
        fprintf('1- efficiency and voltage regulation to the active power curve . \n');
        fprintf('2- efficiency and voltage regulation to power factor curve.\n');
        fprintf('3- Quit.\n');
        CaseChoice = input('\nYour choice: ');
         
        eff= zeros (1000);
        VR= zeros (1000); 
        Receiving_Current= zeros (1000);
        Sending_Voltage= zeros (1000); 
        Sending_Current= zeros (1000); 
        Sending_PowerFactor= zeros (1000); 
        Sending_Power= zeros (1000); 
    
        switch CaseChoice
            case 1
                Receiving_PowerFactor = 0.8;
                sign =-1;
                Receiving_Power = linspace(0, 100000, 1000); 
    
                for i = 1:1000
                    Receiving_Current(i) =( Receiving_Power(i) / (3 * Receiving_Voltage * Receiving_PowerFactor))*exp(1i * acos(Receiving_PowerFactor)* (pi/180) * sign);
                    Sending_Voltage(i)= A* Receiving_Voltage + B*Receiving_Current(i);
                    Sending_Current(i)= C* Receiving_Voltage + D*Receiving_Current(i);
                    Sending_PowerFactor(i)= cos(angle(Sending_Voltage(i)) - angle(Sending_Current(i))); 
                    Sending_Power(i) = 3 * abs(Sending_Voltage(i)) * abs(Sending_Current(i)) * Sending_PowerFactor(i);
                    eff(i) = Receiving_Power(i)*100 / Sending_Power(i);
                    VR(i) = (abs(Sending_Voltage(i))/abs(A) - abs( Receiving_Voltage))*100 / abs(Receiving_Voltage) ;
                end
    
                figure
                title('Efficiency')
                plot(Receiving_Power, eff);
                ylim([0 100]);
                ylabel('Efficiency')
                xlabel('Active Power (W)')
    
                figure
                title('Voltage Regulation')
                plot(Receiving_Power, VR);
                ylabel('Voltage Regulation')
                xlabel('Active Power (W)')
    
            case 2  
                Receiving_PowerFactor = linspace(0.3, 1, 1000);
                sign=-1;
                Receiving_Power  = 100000;

                for i = 1:1000
                    Receiving_Current(i) =( Receiving_Power / (3 * Receiving_Voltage * Receiving_PowerFactor(i)))*exp(1i * acos(Receiving_PowerFactor(i))* (pi/180) * sign);
                    Sending_Voltage(i)= A* Receiving_Voltage + B*Receiving_Current(i);
                    Sending_Current(i)= C* Receiving_Voltage + D*Receiving_Current(i);
                    Sending_PowerFactor(i)= cos(angle(Sending_Voltage(i)) - angle(Sending_Current(i))); 
                    Sending_Power(i) = 3 * abs(Sending_Voltage(i)) * abs(Sending_Current(i)) * Sending_PowerFactor(i);
                    eff(i) = Receiving_Power*100 / Sending_Power(i);
                    VR(i) = (abs(Sending_Voltage(i))/abs(A) - abs( Receiving_Voltage))*100 / abs(Receiving_Voltage) ;
                end
            
                figure
                title('Efficiency')
                plot(abs(Receiving_PowerFactor), eff);
                ylabel('Efficiency')
                xlabel('Power Factor - Lagging')
    
                figure
                title('Voltage Regulation')
                plot(abs(Receiving_PowerFactor), VR);
                ylabel('Voltage Regulation')
                xlabel('Power Factor - Lagging')
 
                sign=1; 
  

                for i = 1:1000
                    Receiving_Current(i) =( Receiving_Power / (3 * Receiving_Voltage * Receiving_PowerFactor(i)))*exp(1i * acos(Receiving_PowerFactor(i))* (pi/180) * sign);
                    Sending_Voltage(i)= A* Receiving_Voltage + B*Receiving_Current(i);
                    Sending_Current(i)= C* Receiving_Voltage + D*Receiving_Current(i);
                    Sending_PowerFactor(i)= cos(angle(Sending_Voltage(i)) - angle(Sending_Current(i))); 
                    Sending_Power(i) = 3 * abs(Sending_Voltage(i)) * abs(Sending_Current(i)) * Sending_PowerFactor(i);
                    eff(i) = Receiving_Power*100 / Sending_Power(i);
                    VR(i) = (abs(Sending_Voltage(i))/abs(A) - abs( Receiving_Voltage))*100 / abs(Receiving_Voltage) ;
                end
  
                figure
                title('Efficiency')
                plot(abs(Receiving_PowerFactor), eff);
                ylabel('Efficiency ')
                xlabel('Power Factor - Leading')

                figure
                title('Voltage Regulation')
                plot(abs(Receiving_PowerFactor), VR);
                ylabel('Voltage Regulation')
                xlabel('Power Factor - Leading')
            case 3
            
            otherwise
                fprintf('Enter a valid Choice ..\n');
        end
    end
end

function conductor_resistivity = get_conductor_resistivity_func()
    conductor_resistivity = input('please enter Conductor Resistivity in ?.m : ');
    while (conductor_resistivity <= 0)
        conductor_resistivity = input('please enter Conductor Resistivity in ?.m above 0: ');
    end   
end

function conductor_length = get_conductor_length_func()
    conductor_length = input('please enter Conductor Length in Km : ')*1000;
    while (conductor_length <= 0)
        conductor_length = input('please enter Conductor Length in Km above 0: ')*1000;
    end   
end
function conductor_diameter = get_conductor_diameter_func()
  fprintf('Enter the Conductor Diameter unit: \n');
  fprintf('1- mm\n');
  fprintf('2- cm\n');

    conductor_diameter_unit = input('');
    
    switch conductor_diameter_unit
            case 1
                conductor_diameter = input('please enter Conductor Diameter in mm : ')/1000;
            case 2
                conductor_diameter = input('please enter Conductor Diameter in cm : ')/100;
             otherwise
                fprintf('wrong choice again!! \n');
                conductor_diameter = get_conductor_diameter_func();
    end
  fprintf('\n');
 
end

function dc_resistance = calc_dc_resistance_func(conductor_resistivity, conductor_length, conductor_diameter)
    dc_resistance = (conductor_resistivity*conductor_length)/( pi.*(conductor_diameter/2).*(conductor_diameter/2));
end

function ac_resistance = calc_ac_resistance_func(conductor_resistivity, conductor_length, conductor_diameter)
    ac_resistance = 1.1.*calc_dc_resistance_func(conductor_resistivity,conductor_length,conductor_diameter);
end


function D_eq = get_equivalent_distance()
    fprintf('please select: \n');
    fprintf('1- symmetrical configuration\n');
  fprintf('2- unsymmetrical configuration\n');
    symmetrical_choice = input('');
    
    switch symmetrical_choice
        case 1
            D_eq = get_equivalent_distance_for_smmetrical_func();
        case 2
            D_eq = get_equivalent_distance_for_unsmmetrical_func();
        otherwise
            fprintf('wrong choice ,try again!! \n');
            D_eq = get_equivalent_distance();
    end 
end

function D_eq = get_equivalent_distance_for_smmetrical_func()
    D_eq = input('Enter Equivalent Distance : ');
end

function D_eq = get_equivalent_distance_for_unsmmetrical_func()
    D1 = input('Enter distance between line 1 and line 2 in m: ');
    D2 = input('Enter distance between line 1 and line 3 in m: ');
    D3 = input('Enter distance between line 2 and line 2 in m: ');
    
    D_eq = (D1*D2*D3).^(1/3);
end

function inductance_per_phase = calc_inductance_per_phase_func(conductor_diameter,conductor_length,D_eq)
    inductance_per_phase = 2*(10)^(-7).*log( D_eq / ( 0.7788.*(conductor_diameter/2) ) )*conductor_length;
end

function capacitance_per_phase = calc_capacitance_per_phase_func(conductor_diameter,conductor_length,D_eq)
    capacitance_per_phase = ( ( 2.*pi.*(8.854187817.*10.^-12) ) ./ log(D_eq ./ (conductor_diameter./2) ) ) * conductor_length;
end

function [A,B,C,D] = TL_Model(Length, Z, Y)
	if Length < 80000
        fprintf("Using Short model\n");
        A =1;
        B = Z;
        C = 0;
        D = 1;
    elseif  Length >= 80000
        fprintf("Using Medium model\n");
        t1 = 1 + Z*Y/2;
        t2 = 1 + Z*Y/4;
        fprintf('Choose the Line model : \n 1- PI model \n 2- T model');
        Line_Model = input('\nYour choice: ');
        
        switch Line_Model
            case 1
                fprintf("Using Pi model\n");
                A = t1;
                B = Z;
                C = Y.*t2;
                D = t1;
            case 2  
                fprintf("Using T model\n");
                A = t1;
                B = Z.*t2;
                C = Y;
                D = t1;
            otherwise
                fprintf('Enter a valid Choice .\n');
                [A,B,C,D] = TL_Model(Length, Z, Y);
        end   
    end
end
