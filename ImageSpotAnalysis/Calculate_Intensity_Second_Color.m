function [IsVirusGood,LineColor,RoughIntensity2,GaussianIntensity2] = Calculate_Intensity_Second_Color(...
    CurrentImage2,CurrentVirusBox,RoughBackground,...
    FigureHandles,BoxToPlot,LineColor,CurrFrameNum,IsVirusGood)

    OffsetTolerance = 0.1;
    
    set(0,'CurrentFigure',FigureHandles.Image2Window);
    plot(BoxToPlot(:,2),BoxToPlot(:,1),LineColor)
    hold on
    drawnow
    
    % Rough background intensity calculation
        CurrentVirusArea = CurrentImage2(...
            CurrentVirusBox.Top:CurrentVirusBox.Bottom,...
            CurrentVirusBox.Left:CurrentVirusBox.Right);

        CurrentRawIntensity = sum(sum((CurrentVirusArea)));

        RoughBackgroundValue = RoughBackground(CurrFrameNum+1);
        RoughIntensity2 = CurrentRawIntensity -...
            RoughBackgroundValue.*(CurrentVirusBox.Bottom - CurrentVirusBox.Top + 1)^2;

    % Gaussian intensity calculation (the background is calculated from a
    % Gaussian fit and then subtracted from the raw intensity values)
        ImageToFit = CurrentVirusArea;
        try
            [OffsetFrom2DFit, Noise] = Vesicle_Gaussian_Fit(ImageToFit);

            GaussianIntensity2 = CurrentRawIntensity -...
                OffsetFrom2DFit*(CurrentVirusBox.Bottom - CurrentVirusBox.Top + 1)^2;
            
            if OffsetFrom2DFit < RoughBackgroundValue - OffsetTolerance*RoughBackgroundValue
            
              IsVirusGood = 'n';
              disp('Failed Gaussian fit-second color')
              
              LineColor = 'r-';
              set(0,'CurrentFigure',FigureHandles.Image2Window);
                plot(BoxToPlot(:,2),BoxToPlot(:,1),LineColor)
                hold on
                drawnow
            end
        
        catch
            IsVirusGood = 'n';
            disp('Failed Gaussian fit-second color')
              
              LineColor = 'r-';
              set(0,'CurrentFigure',FigureHandles.Image2Window);
                plot(BoxToPlot(:,2),BoxToPlot(:,1),LineColor)
                hold on
                drawnow
            
            GaussianIntensity2 = 0;
        end
        
        
end