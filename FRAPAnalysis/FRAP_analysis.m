function FRAP_analysis()
    close all
    clear all

% ======= FRAP_analysis.m ===============
% Script written by Bob Rawle, Williams College, 2025
% Model used to fit FRAP curve is the classic equation by Soumpasis (1983), 
% "Theoretical analysis of fluorescence photobleaching recovery 
% experiments", Biophysical Journal.
%
% Input: This script uses the input of an excel sheet (see the example in this
% same folder for how it should be formatted).
% =======================================

% Options/Inputs/Guesses
    radius = 20.1; %r of spot being photobleached, in um (measure in FIJI)
    diffCoeff_guess = 1; %diffision coefficient, in um^2/sec
    fractImmobile_guess = 0.01; %fraction of molecules that are immobile. 0 to 1.
    correctForPhotofading = 'y'; 

    diffCoeff_bounds = [0, 20];
    fractImmobile_bounds = [0, 1];

% Generate test data (for debugging only)
    % [normFluor_backcorr, time_vector]= Generate_Test_Data();

% Import Data and Process
    
    data = readtable("Example input data.xlsx");
    data = data{:,:}; %convert to matrix format
        time_fromdata = data(:,2); %time vector should have t=0 as the frame after photobleaching
        fluorint_photobleachspot = data(:,3); %fluor intensity of photobleached spot
            fluorint_photobleachspot = fluorint_photobleachspot - 400; %subtract 400 b/c of offset at 2x2 binning
        fileInfo = num2str(data(1,5))+ "-" + num2str(data(2,5));


    idxoftzero = find(time_fromdata(:)==0);
    fluorint_beforebleach = mean(fluorint_photobleachspot(1:idxoftzero-1));
    fluorint_rightafterbleach = fluorint_photobleachspot(idxoftzero); 

    normFluor_photobleachspot = (fluorint_photobleachspot - fluorint_rightafterbleach)./(fluorint_beforebleach-fluorint_rightafterbleach);
    plot(time_fromdata,normFluor_photobleachspot,'o')
    hold on


    idxtostartfit = idxoftzero+1;
    time_vector = time_fromdata(idxtostartfit:end); %exclude data before photobleach

    % Correct for photofading if chosen
    if strcmp(correctForPhotofading,"y")
        fluorint_backgroundspot = data(:,4); %fluor intensity of background spot somewhere else (corner probably)
                fluorint_backgroundspot = fluorint_backgroundspot - 400; %subtract 400 b/c of offset at 2x2 binning
        
        normFluor_backgroundspot = fluorint_backgroundspot./fluorint_backgroundspot(1);
            plot(time_fromdata,normFluor_backgroundspot,'x')
            hold on
    
        LinEq = fittype('m*x+b');
        [LinFit] = fit(time_vector,normFluor_backgroundspot(idxtostartfit:end),LinEq, 'Lower', [0.999 -inf], 'Upper', [1.0001, inf]);
        plot(LinFit,'-')
        Slope = LinFit.m;
    
        correctionFactors = time_vector*Slope+1;
        normFluor_backcorr = normFluor_photobleachspot(idxtostartfit:end)./correctionFactors;
    
        plot(time_vector,normFluor_backcorr,'s')
        hold on
        legend("PhotoBleach Spot, Full Time", "Back Spot", "Back Fit", "PhotoBleach Spot Corrected")

    else
        normFluor_backcorr = normFluor_photobleachspot(idxtostartfit:end);
        plot(time_vector,normFluor_backcorr,'s')
        hold on
        legend("PhotoBleach Spot, Full Time", "PhotoBleach Spot UnCorrected")
    end

    title("Diagnostic Plot")

%Plot the data on the clean figure for the fit
        figure(2)
        plot(time_vector,normFluor_backcorr,'bs')
        hold on


% Fit the Data
    Init_Guess = [diffCoeff_guess fractImmobile_guess];
    Low_Bounds = [diffCoeff_bounds(1) fractImmobile_bounds(1)];
    Up_Bounds = [diffCoeff_bounds(2) fractImmobile_bounds(2)];
    OptimizationOptions = optimset('Display','off','TolFun',1e-9);

    [FitValues, ResNorm, Resid_sglwlag, ExFlag]...
           = lsqnonlin(@Calc_Error,Init_Guess,Low_Bounds,...
           Up_Bounds,...
           OptimizationOptions,time_vector,normFluor_backcorr,radius); %this line is ex params

% Plot and display
    diffCoeff_fromfit = FitValues(1);
    fractImmobile_fromfit = FitValues(2);
    time_for_plot = linspace(.1,max(time_vector),1000);
        %add lots of time points to make it look smoother

    YVals_FromFit = FRAP_Equation(time_for_plot, radius, diffCoeff_fromfit, fractImmobile_fromfit);

    disp("Diff Coeff = " + num2str(diffCoeff_fromfit) + " um^2/sec")
    disp("Fraction Immobile = " + num2str(fractImmobile_fromfit))
    
    hold on
    plot(time_for_plot, YVals_FromFit, "-")
    plot(0,0,"bo")
    ylim([0 1])
    xlabel("t (sec)")
    ylabel("f(t)")
    title("Data and fit; File info = " + fileInfo)

end

function [normFluor_backcorr, time_vector]= Generate_Test_Data()
% Test Data
    radius = 10; %r of spot being photobleached, in um
    diffCoeff = 1; %diffision coefficient, in um^2/sec
    time_vector = 0.1:1:20; %time vector, in sec
    fractImmobile = 0; %fraction of molecules that are immobile. 0 to 1.

    YVals = FRAP_Equation(time_vector, radius, diffCoeff, fractImmobile);

    %add noise
    noise = 0.1;
    normFluor_backcorr = YVals + (rand(size(YVals))-0.5)*noise;

    plot(time_vector, normFluor_backcorr, "bs")
    ylim([0 1])
    xlabel("t (sec)")
    ylabel("f(t)")

end

function Error = Calc_Error(p,x,Data,radius)

    diffCoeff = p(1);
    fractImmobile = p(2);

    YVals_model = FRAP_Equation(x, radius, diffCoeff, fractImmobile);

    Error = Data - YVals_model;

end

function YVals = FRAP_Equation(t, radius, diffCoeff, fractImmobile)
    %To Bob: List reference

    tau = radius.^2/(4*diffCoeff);

    YVals = (1-fractImmobile).*exp(-2*tau./t).*(besseli(0,2*tau./t) + besseli(1,2*tau./t));

    
end
