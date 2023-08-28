% This code is for model calibration of perfusion bioreactor via the single shooting approach
% CAUTION: Model calibration here is based on our based on our first attempt on perfusion bioreactor mechanistic model

clear all; close all; clc

% % How To:
% Substitute perfusionMeasuredData11to6 here with your data fileName
% Substitute sheet1 here with your data sheetName

datafileName = "perfusionMeasuredDataDay1to6";
sheetName = "Sheet1";

dataCreatedTable = readtable(datafileName,"Sheet",sheetName);

% Extract variable name from the table. If a column deos not have a name, corresponding column number is assigned!
% Get the number of columns in the table
numOfColumns = size(dataCreatedTable, 2);
% Assign names to columns without names and use existing names for others
for i = 1:numOfColumns
    if isempty(dataCreatedTable.Properties.VariableNames{i})
        dataCreatedTable.Properties.VariableNames{i} = ['column' num2str(i)];
    end
end
columnNames = dataCreatedTable.Properties.VariableNames;

% Extract those for variables (that is remove time)
variablecolumnNames = columnNames{2:end};

% State variables considered here are:
% Xv = Viable Cells Density, E6 Cells/mL
% Xd = Dead Cells Density, E6 Cells/mL
% Xt = Total Cells Density, E6 Cells/mL
% Glc = Glucose concentration, g/L
% Gln = Glutamine concentration, mM
% Ami = Other amino acids (lumped minus glutamine and glutamate) concentration, mM
% Lac = Lactate concentration, g/L
% Amm = Ammonia concentration, mM
% Titer = mAb concnetration, E6 Cells/mL
% V = Bioreactor Volume, mL
% Glu = Glutamate concentration, mM

% Store the state variable units into a cell array:
stateVariableUnits = {'E6 Cells/mL','E6 Cells/mL','E6 Cells/mL','g/L','mM','mM','g/L','mM','g/L','mL','mM'};

% % Concatenate columnNames with stateVariableUnits:
% concatVarNamesandUnits = cell(1, numel(variablecolumnNames));  % Initialize A3 as a cell array with the same size as A1
% 
% for i = 1:numel(variablecolumnNames)
%     concatVarNamesandUnits(i) = [variablecolumnNames{i}, ', ', stateVariableUnits{i}];
% end

t = dataCreatedTable.time; % XDATA
c1 = dataCreatedTable.Xv;
c2 = dataCreatedTable.Xd;
c3 = dataCreatedTable.Xt;
c4 = dataCreatedTable.Glc;
c5 = dataCreatedTable.Gln;
c6 = dataCreatedTable.Ami;
c7 = dataCreatedTable.Lac;
c8 = dataCreatedTable.Amm;
c9 = dataCreatedTable.Titre;
c10 = dataCreatedTable.Vol;
c11 = dataCreatedTable.Glu;

c = [c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11]; % YDATA

% Define your state variable initial values
c1initial = 1.75; % Viable cells density, [E6 Cells/mL]
c2initial = 1e-2; % Dead cells density, [E6 Cells/mL]
c3initial = c1initial+c2initial; % Total cells density, [E6 Cells/mL]
c4initial = 4.94; % Glucose conc, [g/L]
c5initial = 0.00; %0.50*c4initial; % Glutamine conc, [mM]
c6initial = 0.25*c4initial; % Amino acid conc, [mM]
c7initial = 0.51; % Lactate conc, [g/L]
c8initial = 0.96; % Ammonia conc, [mM]
c9initial = 0.00; % Titre, [g/L]
c10initial = 2500; % Culture volume, [mL]
c11initial = 2.46; % Glutamate conc, [mM]

% Lumped your initial values into initial condition vector
cinitial = [c1initial c2initial c3initial c4initial c5initial c6initial c7initial c8initial c9initial c10initial c11initial];

% Initialize your model parameters, parameters that would want to estimate their value!
theta0 = [1.00;1.00;1.00;1.00;1.00;1.00;1.00;1.00]; % X0, you can individually change these values!

LB = rand(length(theta0),1);
for i = 1:length(theta0)
    LB(i) = 1e-3; % You could play around with this value!
end
UB = []; % No upper bound, but you could play around with this value, say change to 1e5!

% options = optimoptions('lsqcurvefit', 'OptimalityTolerance', 1e-10, 'FunctionTolerance', 1e-10);
tic
[theta,Rsdnrm,Rsd,ExFlg,OptmInfo,Lmda,Jmat] = lsqcurvefit(@paramEstCperfusionDesign,theta0,t,c,LB,UB);
toc

fprintf(1,'\tRate Constants:\n')
for k1 = 1:length(theta)
    fprintf(1, '\t\tTheta(%d) = %8.5f\n', k1, theta(k1))
end
tfit = (linspace(min(t), max(t)))';
Cfit = paramEstCperfusionDesign(theta, tfit);

for i = 1:length(c)
    subplot(3,4,i)
    plot(t, c(:,i), 'p')
    hold on
    hlp = plot(tfit, Cfit(:,i));
    hold off
    grid
    xlabel('Time, [day]')
    ylabel({columnNames{i+1},stateVariableUnits{i}})
    %ylabel(concatVarNamesandUnits)
    legend(hlp, columnNames{i}, 'Location','NE')
end

function C = paramEstCperfusionDesign(theta,t)
% Define your state variable initial values
c1initial = 1.75; % Viable cells density, [E6 Cells/mL]
c2initial = 1e-2; % Dead cells density, [E6 Cells/mL]
c3initial = c1initial+c2initial; % Total cells density, [E6 Cells/mL]
c4initial = 4.94; % Glucose conc, [g/L]
c5initial = 0.00; %0.50*c4initial; % Glutamine conc, [mM]
c6initial = 0.25*c4initial; % Amino acid conc, [mM]
c7initial = 0.51; % Lactate conc, [g/L]
c8initial = 0.96; % Ammonia conc, [mM]
c9initial = 0.00; % Titre, [g/L]
c10initial = 2500; % Culture volume, [mL]
c11initial = 2.46; % Glutamate conc, [mM]

% Lumped your initial values into initial condition vector
cinitial = [c1initial c2initial c3initial c4initial c5initial c6initial c7initial c8initial c9initial c10initial c11initial];

% simulation here:
[t,C] = ode45(@perfusionDesign,t,cinitial);

% ODE function here:
    function dcdt = perfusionDesign(t,c)

        % Define all these model contants. These values usually come from experiments or literature.
        % Anyone that is unkmown can be included as part of model parameter to be estimated
        mugmax = 2*0.481360322;
        mudmax = 0.05;
        mudmin = 0.1*0.05;
        Ks4 = 0.10;
        %Ks5 = 0.5*0.10;
        Ks6 = 0.25*0.10;
        Ks11 = 0.5*0.10;
        Ki7 = 50;
        Ki8 = 10;
        Kd7 = 25;
        Kd8 = 5;
        Km4 = 0.005;
        Km5 = 0.00025*0.005;
        Km6 = 0.025*0.005;
        Km11 = 0.5*0.005;
        Kdi4 = 0.00;
        Kdi5 = 0.000025;
        Kdi6 = 0.025;
        Kdi11 = 0.025;
        Ku7 = 0.05;

        % Define your feed concentration values:
        c1feed = 0.00;
        c2feed = 0.00;
        c3feed = c1feed+c2feed;
        c4feed = 9.030312;
        c5feed = 0.00;
        c6feed = 0.25*c4feed;
        c7feed = 0.00;
        c8feed = 0.00;
        c9feed = 0.00;
        c10feed = 0.00; % Keep this at zero!
        c11feed = 2.5269;

        % Lumped your cfeed values into a vector
        cfeed = [c1feed,c2feed,c3feed,c4feed,c5feed,c6feed,c7feed,c8feed,c9feed,c10feed,c11feed];

        Fmedia = 0.05e-6*c(1)*c(10); % [mL/d]
        Foutlet = 0.05e-6*c(1)*c(10); % [mL/d]
        Fbleed = Fmedia-Foutlet;

        r1 = mugmax*(c(4)/(Ks4+c(4)))*(c(11)/(Ks11+c(11)))*(Ki7/(Ki7+c(7)))*(Ki8/(Ki8+c(8)))*c(1)-(mudmin*(Kd7/(Kd7+c(7)))*(Kd8/(Kd8+c(8)))+mudmax*(c(7)/(Kd7+c(7)))*(c(8)/(Kd8+c(8)))*c(1));
        r2 = mudmin*(Kd7/(Kd7+c(7)))*(Kd8/(Kd8+c(8)))+mudmax*(c(7)/(Kd7+c(7)))*(c(8)/(Kd8+c(8)))*c(1);
        r3 = r1+r2;
        r4 = -theta(1)*(c(4)/(Ks4+c(4)))*c(1)-Km4*c(1)-Kdi4*c(4);
        r5 = theta(2)*c(1)-Km5*c(1)-Kdi5*c(5);
        r6 = -theta(3)*(c(6)/(Ks6+c(6)))*c(1)-Km6*c(1)-Kdi6*c(6);
        r7 = theta(4)*(Ki7/(Ki7+c(7)))*c(1)-Ku7*c(1);
        r8 = theta(5)*(Ki8/(Ki8+c(8)))*c(1)+Kdi5*c(5)+Kdi6*c(6)+Kdi11*c(11);
        r9 = theta(6)*c(1);
        r10 = theta(7)*(Fmedia-Foutlet-Fbleed);
        r11 = -theta(8)*(c(11)/(Ks11+c(11)))*c(1)-Km11*c(1)-Kdi11*c(11);
        dcdt = zeros(11,1);
        dcdt(1) = (Fmedia/c(10))*cfeed(1)-(Foutlet/c(10))*c(1)*0-(Fbleed/c(10))*c(1)+r1-(c(1)/c(10))*dcdt(10); % Viable cells density
        dcdt(2) = (Fmedia/c(10))*cfeed(2)-(Foutlet/c(10))*c(2)*0-(Fbleed/c(10))*c(2)+r2-(c(2)/c(10))*dcdt(10); % Dead cells density
        dcdt(3) = (Fmedia/c(10))*cfeed(3)-(Foutlet/c(10))*c(3)*0-(Fbleed/c(10))*c(3)+r3-(c(3)/c(10))*dcdt(10); % Total cells density
        dcdt(4) = (Fmedia/c(10))*cfeed(4)-(Foutlet/c(10))*c(4)-(Fbleed/c(10))*c(4)+r4-(c(4)/c(10))*dcdt(10); % Glucose concentration
        dcdt(5) = (Fmedia/c(10))*cfeed(5)-(Foutlet/c(10))*c(5)-(Fbleed/c(10))*c(5)+r5-(c(5)/c(10))*dcdt(10); % Glutamine concentration
        dcdt(6) = (Fmedia/c(10))*cfeed(6)-(Foutlet/c(10))*c(6)-(Fbleed/c(10))*c(6)+r6-(c(6)/c(10))*dcdt(10); % Other amino acids concentration
        dcdt(7) = (Fmedia/c(10))*cfeed(7)-(Foutlet/c(10))*c(7)-(Fbleed/c(10))*c(7)+r7-(c(7)/c(10))*dcdt(10); % Lactate concentration
        dcdt(8) = (Fmedia/c(10))*cfeed(8)-(Foutlet/c(10))*c(8)-(Fbleed/c(10))*c(8)+r8-(c(8)/c(10))*dcdt(10); % Ammonia concentration
        dcdt(9) = (Fmedia/c(10))*cfeed(9)-(Foutlet/c(10))*c(9)-(Fbleed/c(10))*c(9)+r9-(c(9)/c(10))*dcdt(10); % Titre (mAb concentration)
        dcdt(10) = cfeed(10)+r10;                                                                            % Reactor volume
        dcdt(11) = (Fmedia/c(10))*cfeed(11)-(Foutlet/c(10))*c(11)-(Fbleed/c(10))*c(11)+r11-(c(11)/c(10))*dcdt(10); % Glutamate concentration
    end
end