classdef POINavigator < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        UIFigure      matlab.ui.Figure
        UIAxes        matlab.ui.control.UIAxes
        UITable       matlab.ui.control.Table
        UIAxes2       matlab.ui.control.UIAxes
        SetPOIButton  matlab.ui.control.Button
        SetORIButton  matlab.ui.control.Button
        LoadButton    matlab.ui.control.Button
        SaveButton    matlab.ui.control.Button
    end

    methods (Access = private)

        % Button pushed function: SetPOIButton
        function SetPOIButtonPushed(app, event)
            
        end

        % Button pushed function: SetORIButton
        function SetORIButtonPushed(app, event)
            
        end

        % Button pushed function: LoadButton
        function LoadButtonPushed(app, event)
            
        end

        % Button pushed function: SaveButton
        function SaveButtonPushed(app, event)
             
        end
    end

    % App initialization and construction
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Create UIFigure
            app.UIFigure = uifigure;
            app.UIFigure.Position = [100 100 823 680];
            app.UIFigure.Name = 'UI Figure';
            setAutoResize(app, app.UIFigure, true)

            % Create UIAxes
            app.UIAxes = uiaxes(app.UIFigure);
            title(app.UIAxes, 'Live');
            app.UIAxes.Box = 'on';
            app.UIAxes.XColor = [0.149 0.149 0.149];
            app.UIAxes.XColorMode = 'manual';
            app.UIAxes.XTick = [];
            app.UIAxes.XTickMode = 'manual';
            app.UIAxes.YTick = [];
            app.UIAxes.YTickMode = 'manual';
            app.UIAxes.Position = [10 385 400 269];

            % Create UITable
            app.UITable = uitable(app.UIFigure);
            app.UITable.ColumnName = {'X'; 'Y'; 'Z'; 'Comment'};
            app.UITable.RowName = {};
            app.UITable.Position = [58 144 513 185];

            % Create UIAxes2
            app.UIAxes2 = uiaxes(app.UIFigure);
            title(app.UIAxes2, 'Template');
            app.UIAxes2.Box = 'on';
            app.UIAxes2.XTick = [];
            app.UIAxes2.XTickMode = 'manual';
            app.UIAxes2.YTick = [];
            app.UIAxes2.YTickMode = 'manual';
            app.UIAxes2.Position = [409 385 400 269];

            % Create SetPOIButton
            app.SetPOIButton = uibutton(app.UIFigure, 'push');
            app.SetPOIButton.ButtonPushedFcn = createCallbackFcn(app, @SetPOIButtonPushed, true);
            app.SetPOIButton.Position = [140 348 100 22];
            app.SetPOIButton.Text = 'Set POI';

            % Create SetORIButton
            app.SetORIButton = uibutton(app.UIFigure, 'push');
            app.SetORIButton.ButtonPushedFcn = createCallbackFcn(app, @SetORIButtonPushed, true);
            app.SetORIButton.Position = [281 348 100 22];
            app.SetORIButton.Text = 'Set ORI';

            % Create LoadButton
            app.LoadButton = uibutton(app.UIFigure, 'push');
            app.LoadButton.ButtonPushedFcn = createCallbackFcn(app, @LoadButtonPushed, true);
            app.LoadButton.Position = [414 348 100 22];
            app.LoadButton.Text = 'Load';

            % Create SaveButton
            app.SaveButton = uibutton(app.UIFigure, 'push');
            app.SaveButton.ButtonPushedFcn = createCallbackFcn(app, @SaveButtonPushed, true);
            app.SaveButton.Position = [559 348 100 22];
            app.SaveButton.Text = 'Save';
        end
    end

    methods (Access = public)

        % Construct app
        function app = POINavigator()

            % Create and configure components
            createComponents(app)

            % Register the app with App Designer
            registerApp(app, app.UIFigure)

            if nargout == 0
                clear app
            end
        end

        % Code that executes before app deletion
        function delete(app)

            % Delete UIFigure when app is deleted
            delete(app.UIFigure)
        end
    end
end