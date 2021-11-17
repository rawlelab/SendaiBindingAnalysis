function [FigureHandles] = Setup_Figures(Options)

        FigureHandles.BinaryImageWindow = figure(3);
        set(FigureHandles.BinaryImageWindow,'Position',[1 -50 450 341]);
        FigureHandles.BackgroundTraceWindow = figure(4);
        set(FigureHandles.BackgroundTraceWindow,'Position',[452 -130 450 341]);
        FigureHandles.ImageWindow = figure(1);
        set(FigureHandles.ImageWindow,'Position',[    1.6667  301.6667  477.3333  338.6667]); %6   479   451   338]);
        FigureHandles.CurrentTraceWindow = figure(2);
        set(FigureHandles.CurrentTraceWindow,'Position',[472   476   450   341]); 
        
        if strcmp(Options.RemoveBleedthroughAreas,'y')
            FigureHandles.BilayerWindow = figure(6);
            set(FigureHandles.BilayerWindow,'Position',[472   196   450   341]);
        end
        
        if strcmp(Options.DualColor,'y') 
            FigureHandles.Image2Window= figure(7);
            set(FigureHandles.Image2Window,'Position',[472.3333  298.3333  503.3333  338.0000]); %472   476    450   341
        end

end