% Description:
% This script will extract the track data and metadata from the output .xml 
% file from the TrackMate plug-in of Image J. It will then use that extracted
% data to determine the 2D diffusion coefficient for each track, by fitting
% into a linear model. The best number of points to include in the fit is
% determined following the algorithm in Michalet, Phys Rev. E., 2010, 82: 041914
% (https://doi.org/10.1103/PhysRevE.82.041914). The outputs are a histogram
% showing all the measured diffusion coefficients, and a histogram showing
% only the diffusion coefficients of those above the immobile cutoff.
% All of the data is saved as a .MAT file.

% Code written by Bob Rawle, Williams College, 2021
% Published online in conjunction with the manuscript by Lam et al. (2021).

% Special thanks to example code at:
% https://measurebiology.org/wiki/Calculating_MSD_and_Diffusion_Coefficients

close all
clear all

Filename = '211030- 104- SPLAT TR mobility 5s timelapse- no rinse_Tracks.xml';
% Filename = '210702- 303- TR mobility 5s timelapse after 30 min-1_Tracks.xml';
% Filename = '210702- 403- TR mobility 1s timelapse_Tracks.xml';
Filename = 'test best (210122-406)_Tracks.xml';

ImmobileDiffCoeff = 0.000060496;
HistogramBins = 100;

[Trackdata, Metadata] = importTrackMateTracks(Filename);

NumTracks = numel(Trackdata);

TimeVector = [1:301]'; % skip deleted frames, go from back to front ; [1:9,12:40]
DefaultNumPtsToIncludeInFit = 10; %round(0.1*(length(TimeVector)-1));
FrameDuration = 0.1; %in sec

DiffusionFitEq = fittype('4*D*x+Yint');
    %MSD = 4D?t for 2D diffusion

DiffCoeffList = [];
LocErrorList = [];
ConvergedList = [];
NumExcluded = 0;

for b = 1:NumTracks
   
    CurrTrackLength = length(Trackdata{b}(:,1));
    NumMSDsToCalc = CurrTrackLength - 1;
    
    ParticleData(b).MSDVector = [];
    BigdX = [];
    BigdY = [];
    BigdT = [];
    deltaFrame = 1:NumMSDsToCalc;
    
    for dF = deltaFrame
       dx = Trackdata{b}((dF+1):1:end,2) - Trackdata{b}(1:1:(end-dF),2);
       dy = Trackdata{b}((dF+1):1:end,2) - Trackdata{b}(1:1:(end-dF),2);
       dT = TimeVector((dF+1):1:CurrTrackLength) - TimeVector(1:1:(CurrTrackLength-dF));
       
       BigdX = [BigdX; dx];
       BigdY = [BigdY; dy];
       BigdT = [BigdT; dT];
       
    end
    
    UniquedT = unique(BigdT);
    
    for i = 1:length(UniquedT)
        
        deltaT = UniquedT(i);
        
        Idx = BigdT == deltaT;
        
        dX_vector = BigdX(Idx);
        dY_vector = BigdY(Idx);
        
        drSquared = dX_vector.^2 + dY_vector.^2;
        MSD = mean(drSquared);
        
        ParticleData(b).MSDVector = [ParticleData(b).MSDVector, MSD];
        
    end
    
    % TimeVector = Trackdata{b}(:,1);
    figure(1)
    plot(UniquedT, ParticleData(b).MSDVector);
    
    %% Fit w/ best num of points
            MaxLoops = 20;
            KeepGoing = 'y';
            Status = 'Not Yet Converged';
            LoopCounter = 0;
            NumPtsToIncludeInFit = DefaultNumPtsToIncludeInFit;
            OldMinNumPoints = 0;
            KeepThisOne = 'n';

            while KeepGoing == 'y'
                LoopCounter = LoopCounter + 1;

                EndPtOfFit = min([length(ParticleData(b).MSDVector),NumPtsToIncludeInFit]);
                FitResults = fit(UniquedT(1:EndPtOfFit), ParticleData(b).MSDVector(1:EndPtOfFit)',DiffusionFitEq,'StartPoint',[0 0]);

                    CoeffOutputs = coeffvalues(FitResults);
                    DiffCoeff = CoeffOutputs(1);
                    LocError = CoeffOutputs(2);

                    ReducedLocError = LocError^2/(DiffCoeff*FrameDuration);

                    NewMinNumPoints = round(2+2.7*ReducedLocError^0.5);

                    if NewMinNumPoints == NumPtsToIncludeInFit &&...
                            NewMinNumPoints == OldMinNumPoints

                        %Fit has converged
                         KeepGoing = 'n';
                         Status = 'Converged';
                         
                         if isreal(NewMinNumPoints)
                             KeepThisOne = 'y';
                         else
                             %KeepThisOne = 'n';
                             %NumExcluded = NumExcluded + 1;
                             
                            Status = 'Converged To Imaginary Number; Use Default';
                            NumPtsToIncludeInFit = DefaultNumPtsToIncludeInFit;
                            EndPtOfFit = min([length(ParticleData(b).MSDVector),NumPtsToIncludeInFit]);
                            FitResults = fit(UniquedT(1:EndPtOfFit), ParticleData(b).MSDVector(1:EndPtOfFit)',DiffusionFitEq,'StartPoint',[0 0]);
                            KeepThisOne = 'y';
                            
                         end

                    elseif NewMinNumPoints >= length(ParticleData(b).MSDVector)
                        %Never converged
                        KeepGoing = 'n';
                        Status = 'Converged To Large Number; Use Default';
                        NumPtsToIncludeInFit = DefaultNumPtsToIncludeInFit;
                        EndPtOfFit = min([length(ParticleData(b).MSDVector),NumPtsToIncludeInFit]);
                        FitResults = fit(UniquedT(1:EndPtOfFit), ParticleData(b).MSDVector(1:EndPtOfFit)',DiffusionFitEq,'StartPoint',[0 0]);
                        KeepThisOne = 'y';
                        
                    else
                        %Try another round;
                        OldMinNumPoints = NumPtsToIncludeInFit;
                        NumPtsToIncludeInFit = NewMinNumPoints;
                    end

                if LoopCounter >= MaxLoops
                    %Never converged
                    KeepGoing = 'n';
                    Status = 'Never Converged, Use Default';
                    NumPtsToIncludeInFit = DefaultNumPtsToIncludeInFit;
                    EndPtOfFit = min([length(ParticleData(b).MSDVector),NumPtsToIncludeInFit]);
                    FitResults = fit(UniquedT(1:EndPtOfFit), ParticleData(b).MSDVector(1:EndPtOfFit)',DiffusionFitEq,'StartPoint',[0 0]);
                    KeepThisOne = 'y';
                end
            end
    
    %% Plot results
    
    CoeffOutputs = coeffvalues(FitResults);
    
    if KeepThisOne == 'y'
        DiffCoeffList = [DiffCoeffList; CoeffOutputs(1)];
        LocErrorList = [LocErrorList; CoeffOutputs(2)];
        if strcmp(Status,'Converged')
            ConvergedList = [ConvergedList; 1];
        else
            ConvergedList = [ConvergedList; 0];
        end
    end
    
    hold on
    plot(FitResults,'-')
    xlabel('deltat')
    ylabel('MSD')
    title(strcat(Status, "; Keep This One = ", KeepThisOne, "; Num Pts = ", num2str(NumPtsToIncludeInFit)))
    hold off
    
    figure(2)
        try
            plot(UniquedT(1:NumPtsToIncludeInFit), ParticleData(b).MSDVector(1:NumPtsToIncludeInFit));
        catch
            ErrorHasOccurred = 1;
        end
        hold on
        plot(FitResults,'-')
        xlabel('deltat')
        ylabel('MSD')
        hold off
        title(strcat("Diff Coeff = ", num2str(CoeffOutputs(1)), "; Loc Error = ", num2str(CoeffOutputs(2)),...
            "; Reduced Loc Error = ", num2str(ReducedLocError)));
    drawnow
    
end

figure(3)
% IndxToKeep = DiffCoeffList > ImmobileDiffCoeff; %0.03612
% DiffCoeffList = DiffCoeffList(IndxToKeep);
histogram(DiffCoeffList); %Look up specif, bin width ;
% histogram(DiffCoeffList,HistogramBins)
xlabel('Diffusion Coeff')
ylabel('Num Observed')
%end

figure(4) % for mobile particles only
IndxToKeep = DiffCoeffList > ImmobileDiffCoeff;
Mobile = DiffCoeffList(IndexToKeep);
histogram(Mobile,HistogramBins)
xlabel('Diffusion Coeff')
ylabel('Num Observed')

% immobile = DiffCoeffList(DiffCoeffList < ImmobileDiffCoeff)
% total = DiffCoeffList
% percentImmobile = (immobile/total)*100

save(strcat(Filename,'-AnalysisFile','.mat'))