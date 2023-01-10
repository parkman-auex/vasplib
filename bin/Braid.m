classdef Braid < matlab.mixin.CustomDisplay
    %UNTITLED7 Summary of this class goes here 
    %   Detailed explanation goes here

    properties
        BraidWord;
        Nstrings;
        Crossings;
    end
    properties
        Generator;
        BraidNum;
    end
    properties
        NaiveSeperateBands;
        PlotSeqL;
        PlotSeqLMirror;
    end
    properties
        Permutation;
        CycleDecompositionStr;
        CycleDecomposition;
    end
    %% Define which properties show
    methods (Access = protected)
        function propgrp = getPropertyGroups(~)
            proplist = {'BraidWord','Nstrings','Crossings','CycleDecompositionStr'};
            propgrp = matlab.mixin.util.PropertyGroup(proplist);
        end
    end
    methods
        function BraidObj = Braid(BraidWord)
            
            DoubleBraidWord_col = Braid.BraidWord2Num(BraidWord);
            Nstrings = Braid.NumBraidWord2Nstrings(DoubleBraidWord_col);
            Crossings = length(DoubleBraidWord_col);
            [Generator,GeneratorStr] = Braid.Nstrings2Generator(Nstrings);
            Permutation =  Braid.PermutationBraidNum(DoubleBraidWord_col);
            [CyclePresentation,CyclePresentationStr] = Braid.Cauchy2Cycle(Permutation);
            %
            BraidObj.BraidWord = BraidWord;
            BraidObj.BraidNum = DoubleBraidWord_col;
            BraidObj.Nstrings = Nstrings;
            BraidObj.Crossings = Crossings;
            BraidObj.Generator = Generator;
            BraidObj.Permutation = Permutation;
            BraidObj.CycleDecomposition = CyclePresentation;
            BraidObj.CycleDecompositionStr = CyclePresentationStr;
            [BraidObj.NaiveSeperateBands,BraidObj.PlotSeqL,BraidObj.PlotSeqLMirror] = NaiveBands(BraidObj);
        end

    end
    methods
        function [NaiveSeperateBands,PlotSeqL,PlotSeqLVertical] = NaiveBands(BraidObj)
            NaiveSeperateBands{BraidObj.Crossings} = (1:BraidObj.Nstrings).';
            PermutationList = abs(BraidObj.BraidNum);
            BraidList = BraidObj.BraidNum;
            Nbands = BraidObj.Nstrings;
            DefaultL = (1:BraidObj.Nstrings).';
            PlotSeqL = repmat(DefaultL,[1 BraidObj.Crossings]);
            PlotSeqLVertical = PlotSeqL;
            for i = 1:BraidObj.Crossings
                if i > 1
                    NaiveSeperateBands{i}(:,1) =  NaiveSeperateBands{i-1}(:,2);
                else
                    NaiveSeperateBands{i}(:,1) = DefaultL;
                end
                PlotSeqL(:,i) =  NaiveSeperateBands{i}(:,1);
                PlotSeqLVertical(:,i)  = NaiveSeperateBands{i}(:,1);
                if BraidList(i) > 0
                    PlotSeqL([PermutationList(i),PermutationList(i)+1],i) = ...
                        PlotSeqL([PermutationList(i)+1,PermutationList(i)],i);
                else
                   PlotSeqLVertical([PermutationList(i),PermutationList(i)+1],i) = ...
                        PlotSeqLVertical([PermutationList(i)+1,PermutationList(i)],i);
                end
                for n = 1:Nbands
                    if NaiveSeperateBands{i}(n,1) == PermutationList(i)
                        NaiveSeperateBands{i}(n,2) = NaiveSeperateBands{i}(n,1) + 1;
                    elseif NaiveSeperateBands{i}(n,1) == PermutationList(i) +1
                        NaiveSeperateBands{i}(n,2) = NaiveSeperateBands{i}(n,1) - 1;
                    else
                        NaiveSeperateBands{i}(n,2) = NaiveSeperateBands{i}(n,1);
                    end
                end
            end
        end
    end
    methods
        function ax = LinePlot(BraidObj,options,opts)
            arguments
                BraidObj Braid
                options.ax = gca;
                options.Color = @parula;
                options.LineSpec = '-';
                options.LineWidth = 10;
                options.MarkerSize = 3;
                options.MarkerEdgeColor = 'none';
                options.MarkerFaceColor = 'none';
                opts.vertical = false;
            end
            Nbands=BraidObj.Nstrings;
            ax = options.ax;
            hold(ax,'on');
            if isstring(options.Color)
                if isrow
                    colorL = options.Color.';
                end
            elseif isnumeric(options.Color)
                x = linspace(0,1,size(options.Color,1));
                xq = linspace(0,1,Nbands);
                colorL = [...
                    interp1(x,options.Color(:,1),xq).',...
                    interp1(x,options.Color(:,2),xq).',...
                    interp1(x,options.Color(:,3),xq).',...
                    ];
            elseif isa(options.Color,'function_handle')
                colorL = options.Color(Nbands);
            else
                for i = 1:Nbands
                    colorL(i,:) = [rand,rand,rand];
                end
            end
            HSV = rgb2hsv(colorL);
            HSV(:,2) = HSV(:,2);
            HSV(:,3) = HSV(:,3)-0.2;
            ModifycolorL = hsv2rgb(HSV);
            PlotSeqList = BraidObj.PlotSeqL;
            EIGENCARCELL = BraidObj.NaiveSeperateBands;
            Ncross = BraidObj.Crossings;       
            EIGENCAR = zeros(Nbands,Ncross+1);
            EIGENCAR(:,1) = EIGENCARCELL{1}(:,1);
            if opts.vertical
                PlotSeqList = BraidObj.PlotSeqL;
                for i = 1:Ncross
                    EIGENCAR(:,i+1) = EIGENCARCELL{i}(:,2);
                    for n = 1:BraidObj.Nstrings
                        
                        Ei = PlotSeqList(n,i);
                        plot(ax,EIGENCARCELL{i}(Ei,:),[i-1 i],options.LineSpec,...
                            'LineWidth',options.LineWidth,...
                            'Color',colorL(Ei,:),...
                            'MarkerSize',options.MarkerSize,...
                            'MarkerEdgeColor',options.MarkerEdgeColor,...
                            'MarkerFaceColor',options.MarkerFaceColor,...
                            'DisplayName',['E_',num2str(Ei)]);
                    end
                end
                for n = 1:Nbands
                    scatter(ax,EIGENCAR(n,:),0:Ncross,...
                        ones(1,Ncross+1)*options.LineWidth*20,ModifycolorL(n,:),"filled",'DisplayName',['String_',num2str(n)]);
                end
            else
                for i = 1:Ncross
                    EIGENCAR(:,i+1) = EIGENCARCELL{i}(:,2);
                    for n = 1:BraidObj.Nstrings
                        
                        Ei = PlotSeqList(n,i);
                        plot(ax,[i-1 i],-EIGENCARCELL{i}(Ei,:),options.LineSpec,...
                            'LineWidth',options.LineWidth,...
                            'Color',colorL(Ei,:),...
                            'MarkerSize',options.MarkerSize,...
                            'MarkerEdgeColor',options.MarkerEdgeColor,...
                            'MarkerFaceColor',options.MarkerFaceColor,...
                            'DisplayName',['E_',num2str(Ei)]);
                    end
                end
                for n = 1:Nbands
                    scatter(ax,0:Ncross,-EIGENCAR(n,:),...
                        ones(1,Ncross+1)*options.LineWidth*20,ModifycolorL(n,:),"filled",'DisplayName',['String_',num2str(n)]);
                end
            end
            ylabel(ax,'');
            yticks(ax,[]);
            xlabel(ax,'');
            xticks(ax,[]);
            axis(ax,'off');
        end
    end
    methods(Static)
        function DoubleBraidWord_col = BraidWord2Num(BraidWord)
            CharBraidWord = char(BraidWord);
            CharBraidWord_col = (CharBraidWord.');
            DoubleBraidWord_col = double(CharBraidWord_col);
            for i = 1:numel(DoubleBraidWord_col)
                if DoubleBraidWord_col(i) < 97 
                    DoubleBraidWord_col(i) =  DoubleBraidWord_col(i) - 64;
                else
                    DoubleBraidWord_col(i) = -DoubleBraidWord_col(i) + 96;
                end
                
            end
            
        end
        function Nstrings = NumBraidWord2Nstrings(DoubleBraidWord_col)
            Nstrings = max(abs(DoubleBraidWord_col)) + 1;
        end
        function [Generator,GeneratorStr] = Nstrings2Generator(Nstrings)
            TmpList =  string((1:Nstrings-1).');
            GeneratorStr = "sigma_" + TmpList;
            Generator = str2sym(GeneratorStr);
        end
        function Permutation = PermutationBraidNum(DoubleBraidWord_col)
            PermutationList = abs(DoubleBraidWord_col);
            Nstrings = max(abs(DoubleBraidWord_col)) + 1;
            Permutation(1,:) = 1:Nstrings;
            for iString = Permutation(1,:)
                Permutation(2,iString) = Braid.FinalPosition(iString,PermutationList);
            end
        end
        function oString = FinalPosition(iString,PermutationList)
            tmpString =  iString;
            for i = 1:numel(PermutationList)
                if tmpString == PermutationList(i)
                    tmpString = tmpString + 1;
                elseif tmpString == PermutationList(i) +1
                    tmpString = tmpString - 1;
                else

                end
            end
            oString = tmpString;
        end
        function [CyclePresentation,CyclePresentationStr] = Cauchy2Cycle(Permutation)
            %CyclePresentation{1} = 1;
            PermutationStore = Permutation(1,:);
            PermutationFunction = Permutation(2,:);
            Ndivide = 0;
            while ~isempty(PermutationStore)
                TriceElement = PermutationStore(1);
                Converge = false;
                RefTriceElement = TriceElement;
                Ndivide = Ndivide + 1;
                count = 0;
                while ~Converge
                    count = count + 1;
                    CyclePresentation{Ndivide}(count) = TriceElement;
                    PermutationStore(PermutationStore == TriceElement) = [];
                    TriceElement = PermutationFunction(TriceElement);
                    Converge = TriceElement == RefTriceElement;
                end            
            end
            CyclePresentationStr = ''; 
            for i = 1:numel(CyclePresentation)
                CyclePresentationStr = [CyclePresentationStr,'('];
                CyclePresentationStr = [CyclePresentationStr,num2str(CyclePresentation{i})];
                CyclePresentationStr = [CyclePresentationStr,')'];
            end
        end
    end
    methods(Static) % script
        function Permutation = PermutationBraidWord(BraidWord)
            DoubleBraidWord_col = Braid.BraidWord2Nstrings(BraidWord);
            Permutation = Braid.PermutationBraidNum(DoubleBraidWord_col);
        end
        function Nstrings = BraidWord2Nstrings(BraidWord)
            DoubleBraidWord_col = Braid.BraidWord2Nstrings(BraidWord);
            Nstrings = Braid.NumBraidWord2Nstrings(DoubleBraidWord_col);
        end
    end
end