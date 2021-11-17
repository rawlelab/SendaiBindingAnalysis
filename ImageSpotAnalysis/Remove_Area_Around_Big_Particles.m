function [CurrentImage,BinaryCurrentImage] = Remove_Area_Around_Big_Particles(Options,...
    ImageWidth,ImageHeight,CurrentImage,BinaryCurrentImage,FigureHandles)

        %Plot the image
            set(0,'CurrentFigure',FigureHandles.ImageWindow);
            hold off
            imshow(CurrentImage, [Options.MinImageShow, Options.MaxImageShow], 'InitialMagnification', 'fit','Border','tight');
            hold on
            
            set(0,'CurrentFigure',FigureHandles.BinaryImageWindow);
            imshow(BinaryCurrentImage, 'InitialMagnification', 'fit','Border','tight');
            drawnow
            
        VirusComponentArray = bwconncomp(BinaryCurrentImage,8);

    %The properties associated with each virus in the binary image are
    %extracted.
        VirusProperties = regionprops(VirusComponentArray, CurrentImage, 'Centroid',...
            'Eccentricity', 'PixelValues', 'Area','PixelIdxList');
        NumberOfVirusesFound = length(VirusProperties);

    for j = 1:NumberOfVirusesFound
        if VirusProperties(j).Area > Options.MinAreaBig
            CurrentVirusArea = VirusProperties(j).Area;
            CurrentVirusCoordinates = VirusProperties(j).Centroid;
                CurrentVirusX = CurrentVirusCoordinates(1);
                CurrentVirusY = CurrentVirusCoordinates(2);

            SizeOfSquareAroundCurrVesicle = (sqrt(CurrentVirusArea)*2+20);

                %Coordinates are set up to define the area around the vesicle
                %of interest (i.e. the ROI)  
                CurrentBoxLeft = max([round(CurrentVirusX) - round(SizeOfSquareAroundCurrVesicle/2),...
                    1]);
                CurrentBoxRight = min([round(CurrentVirusX) + round(SizeOfSquareAroundCurrVesicle/2),...
                    ImageWidth]);
                CurrentBoxTop = max([round(CurrentVirusY) - round(SizeOfSquareAroundCurrVesicle/2),...
                    1]);
                CurrentBoxBottom = min([round(CurrentVirusY) + round(SizeOfSquareAroundCurrVesicle/2),...
                    ImageHeight]);


                CurrentImage(CurrentBoxTop:CurrentBoxBottom,...
                    CurrentBoxLeft:CurrentBoxRight) = 0;

                BinaryCurrentImage(CurrentBoxTop:CurrentBoxBottom,...
                    CurrentBoxLeft:CurrentBoxRight) = 0;
        end

    end

   
        %Plot the corrected image
            set(0,'CurrentFigure',FigureHandles.ImageWindow);
            hold off
            imshow(CurrentImage, [Options.MinImageShow, Options.MaxImageShow], 'InitialMagnification', 'fit','Border','tight');
            hold on
            
            set(0,'CurrentFigure',FigureHandles.BinaryImageWindow);
            imshow(BinaryCurrentImage, 'InitialMagnification', 'fit','Border','tight');
            drawnow
            
end