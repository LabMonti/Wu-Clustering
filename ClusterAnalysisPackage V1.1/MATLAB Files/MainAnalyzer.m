%Hierarchical clustering program for quantum transport histogram data
%Copyright 2017 LabMonti, The University of Arizona

%Original Article: Uncovering Hierarchical Data Structure in Singe Molecule 
%Transport, B.H. Wu, J.A. Ivie, T.K. Johnson, O.L.A. Monti, J. Chem. Phys.,
%146, 092321 (2017) DOI: http://dx.doi.org/10.1063/1.4974937

%Created & Written by Ben Wu
%Edited by Nathan Bamberger

%Version 1.1 Changes:
%1) Added line for removing zero-count bins before processing
%2) Added wait bar for optics clustering
%3) Optimized memory usage in CDbw index calculation

function MainAnalyzer()

ss = get(0,'screensize');
ss = ss(ss>1);
dim = min(ss);
ui = figure('Visible','off','outerposition',[0 0 dim dim]);

pnl = uipanel('Title','Main Panel','FontSize',12,...
                'BackgroundColor','white',...
                'Position',[.05 .1 .4 .3]);
subpnl = uipanel('parent',pnl,'Title','Input Used','FontSize',10,...
                'Position',[0 0.2 .25 .6]);

subpnl1 = uipanel('parent',pnl,'Title','Operations','FontSize',10,...X1
                'Position',[.25 0 .75 .8]);

stop_button = uicontrol('parent',pnl,'Style','togglebutton',...
    'String','Exit',...
    'units','normalized',...
    'Position',[0 0 .25 .2]);

tbox = uicontrol('parent',ui,'style','edit','enable','inactive',...
    'units','normalized',...
    'horizontalalignment','left',...
    'fontsize',14,...
    'Position',[.05 .41 .4 .05]);
uicontrol('parent',ui,'Style','Text','String','Status',...
    'units','normalized',...
    'horizontalalignment','left',...
    'fontsize',12,...
    'fontweight','bold',...
    'Position',[.05 .46 .4 .03]); 

filepath = uicontrol('parent',pnl,'Style','Edit','enable',...
    'inactive','units','normalized','Position',[0 .8 1 .1]); 
uicontrol('parent',pnl,'Style','Text','String','File Path',...
    'units','normalized',...
    'fontsize',12,...
    'Position',[0 .9 1 .1]); 

min_pts = uicontrol('parent',subpnl,'Style','Edit','enable','inactive',...
    'units','normalized',...
    'Position',[0 .75 .3 .25]); 
uicontrol('parent',subpnl,'Style','Text','String','min pts',...
    'units','normalized',...
    'Position',[.3 .7 .7 .25]); 

dbscan_eps = uicontrol('parent',subpnl,'Style','Edit','enable','inactive',...
    'units','normalized',...
    'Position',[0 .5 .3 .25]); 
uicontrol('parent',subpnl,'Style','Text','String','eps',...
    'units','normalized',...
    'Position',[.3 .45 .7 .25]); 

optics_eps = uicontrol('parent',subpnl,'Style','Edit','enable','inactive',...
    'units','normalized',...
    'Position',[0 .25 .3 .25]); 
uicontrol('parent',subpnl,'Style','Text','String','geneps',...
    'units','normalized',...
    'Position',[.3 .2 .7 .25]); 

weight = uicontrol('parent',subpnl,'Style','Edit','enable','inactive',...
    'units','normalized',...
    'Position',[0 .0 .3 .25]); 
uicontrol('parent',subpnl,'Style','Text','String','weight',...
    'units','normalized',...
    'Position',[.3 -.05 .7 .25]); 

importdata = uicontrol('parent',subpnl1,'Style','pushbutton',...
    'String','Load Data Matrix','Callback',@ImportDataCallback,...
    'units','normalized',...
    'Position',[.05 .83 .45 .16]); 

datastand = uicontrol('parent',subpnl1,'Style','PushButton',...
    'String','Standardize Data','Callback',@DataStandardizationCallback,...
    'units','normalized',...
    'Position',[.05 .67 .45 .16]); 

opticsdata = uicontrol('parent',subpnl1,'Style','PushButton',...
    'String','OPTICS Clustering','Callback',@OpticsDataCallback,...
    'units','normalized',...
    'Position',[.05 .5 .45 .16]); 

opticsextract = uicontrol('parent',subpnl1,'Style','PushButton',...
    'String','Extract Optics Partition','Callback',@ExtractOpticsPartitionCallback,...
    'units','normalized',...
    'Position',[.05 .33 .45 .16]);

tableview = uicontrol('parent',subpnl1,'Style','PushButton',...
    'String','View Cluster Size Table','Callback',@ViewTableCallback,...
    'units','normalized',...
    'Position',[.05 0.16 .45 .16]); 

extracclust = uicontrol('parent',subpnl1,'Style','PushButton',...
    'String','Extract Clusters','Callback',@ExtractClustCallback,...
    'units','normalized',...
    'Position',[.05 0 .45 .16]); 

CDbwtest = uicontrol('parent',subpnl1,'Style','PushButton',...
    'String','CDbw Analysis','Callback',@CDbwTestCallback,...
    'units','normalized',...
    'Position',[.525 .83 .45 .16]); 

printCDbw = uicontrol('parent',subpnl1,'Style','pushbutton',...
    'Callback',@Print_CDbwCallback,...
    'String','Print CDbw Partition',...
    'units','normalized',...
    'Position',[.525 .67 .45 .16]); 

gammatest = uicontrol('parent',subpnl1,'Style','PushButton',...
    'String','Gamma Statistic','Callback',@GammaTestCallback,...
    'units','normalized',...
    'Position',[.525 .5 .45 .16]);

print2Dhist = uicontrol('parent',subpnl1,'Style','pushbutton',...
    'Callback',@Print_2D_HistCallback,...
    'String','Print 2D Histogram',...
    'units','normalized',...
    'Position',[.525 .33 .45 .16]); 

printoptics = uicontrol('parent',subpnl1,'Style','pushbutton',...
    'Callback',@Print_OpticsCallback,...
    'String','Print Reachability Plot',...
    'units','normalized',...
    'Position',[.525 0.17 .45 .16]); 

printclust = uicontrol('parent',subpnl1,'Style','pushbutton',...
    'Callback',@Print_ClustCallback,...
    'String','Print Clustering',...
    'units','normalized',...
    'Position',[.525 0 .45 .16]); 

ui.Visible = 'on';
setappdata(subpnl1,'val',[]);
setappdata(subpnl1,'partitions',[]);

pause(.1)

%This loop allows the code to continuously wait for buttons to be pressed
while get(stop_button,'Value') < 1    
    pause(.1)
end

close all

end

%Imports histogram data from the file specified in the filepath text box.
%Expected format of file is three columns of data: distance,
%log(conductance), count (i.e. number of points in that
%distance-conductance bin)
function ImportDataCallback(hObject,~)
    ClearAppData('X1',hObject)
    ClearAppData('X2',hObject)
    ClearAppData('RD',hObject)
    ClearAppData('CD',hObject)
    ClearAppData('order',hObject)
    ClearAppData('Y',hObject)
    ClearAppData('T',hObject)
    ClearAppData('eps_range',hObject)
    ClearAppData('partitions',hObject)
    ClearAppData('index',hObject)
    
    set(evalin('caller','tbox'),'string','Running File IO')
    drawnow
    
    prompt = {'Enter file path for input data file:'};
    dlg_title = 'File path input';
    defaultans = {'Cluster Analysis Package for GitHub\Sample Data and Solutions\Fig5.dat'};
    usrinput = inputdlg(prompt,dlg_title,[1 60],defaultans); 
    path = usrinput{1};
    set(evalin('caller','filepath'),'string',usrinput{1});
    drawnow
    
    try
        importdata(path);
    catch
        set(evalin('caller','tbox'),'string','File not Found')
    end
    
    %X1 is a matrix of this 3-column histogram data
    X1 = importdata(path);
    
    %Remove zero-count bins:
    X1 = X1(X1(:,3) ~= 0, :);
    
    setappdata(hObject.Parent,'X1',X1);
    set(evalin('caller','tbox'),'string','File IO Completed')
end

%Takes X1 histogram data, inverts counts, standardizes each column, and
%weights conductances by weighting factor w.  The new data is put into
%matrix X2 (with X1 still having the original data)
%From now on, each row of X2 is considered a single data point
function DataStandardizationCallback(hObject,~)
    ClearAppData('X2',hObject)
    ClearAppData('RD',hObject)
    ClearAppData('CD',hObject)
    ClearAppData('order',hObject)
    ClearAppData('Y',hObject)
    ClearAppData('T',hObject)
    ClearAppData('eps_range',hObject)
    ClearAppData('partitions',hObject)
    ClearAppData('index',hObject)

    set(evalin('caller','tbox'),'string','Running Data Standardization')
    drawnow
    
    %Get w from user input, then display and store value
    prompt = {'Enter scaling weight for conductance data (recommended value is 1.5):'};
    dlg_title = 'Scaling Weight Input';
    defaultans = {'1.5'};
    usrinput = inputdlg(prompt,dlg_title,1,defaultans);
    w = str2double(usrinput{1});
    setappdata(hObject.Parent,'w',w);
    set(evalin('caller','weight'),'string',usrinput{1});
    
    X1 = getappdata(hObject.Parent,'X1');
    X2 = DataStandardization(X1,w);
    setappdata(hObject.Parent,'X2',X2);
    set(evalin('caller','tbox'),'string','Data Standardization Completed')
end

%The heart of the program!  This function uses the OPTICS algorithm on the
%X2 data to find the cluster order, as well as the reachability and core
%distances of each point in X2.  
%Also plots the reachability plot on-screen.  
%order is a vector containing the order in which each point (row) of X2 was
%processed (e.g. [5, 2, 3, ...] => 1st points is 5th in cluster order,
%etc.)
%RD is the a vector containing the reachability distance for each point
%(row) in X2, CD is the same thing for the core distance
function OpticsDataCallback(hObject,~)
    ClearAppData('RD',hObject)
    ClearAppData('CD',hObject)
    ClearAppData('order',hObject)
    ClearAppData('Y',hObject)
    ClearAppData('T',hObject)
    ClearAppData('eps_range',hObject)
    ClearAppData('partitions',hObject)
    ClearAppData('index',hObject)

    set(evalin('caller','tbox'),'string','Running Optics Clustering')
    drawnow
    
    %Get k and geneps values from user input (then display and store
    %values)
    prompt = {'Enter value for min points (recommended value is 25):', 'Enter value for generating epsilon (recommended value is 1):'};
    dlg_title = 'Min Points and Generating Epsilon Input';
    defaultans = {'25', '1'};
    usrinput = inputdlg(prompt,dlg_title,1,defaultans);
    k = str2double(usrinput{1});
    geneps = str2double(usrinput{2});
    setappdata(hObject.Parent, 'k', k);
    setappdata(hObject.Parent, 'geneps', geneps);
    set(evalin('caller','min_pts'), 'string', usrinput{1});
    set(evalin('caller','optics_eps'), 'string', usrinput{2});
    drawnow
    
    X2 = getappdata(hObject.Parent,'X2');
    [RD,CD,order] = OpticsAnalysis(X2,geneps,k);
    setappdata(hObject.Parent,'order',order);
    setappdata(hObject.Parent,'RD',RD);
    setappdata(hObject.Parent,'CD',CD);
    fig = gcf;
    setappdata(hObject.Parent,'opt',fig.CurrentAxes);
    set(evalin('caller','tbox'),'string','Optics Clustering Complete')
end

%Uses the reachability plot to assign points to clusters for the
%user-defined value of epsilon (eps).  
%Also plots the clustering solution on-screen.  
%Y is a vector containing the cluster ID# of each point (row) in X2. Note
%that points in clusters representing <0.35% of total points were
%reassigned to noise cluster (Id# -1)
%T is a table with one row for each cluster ID, listing the ID# of the
%cluster, the number of points assigned to the cluster, and the percentage
%of all points that belong to that cluster.
%CDbw is the CDbw index value for the specific clustering solution
function ExtractOpticsPartitionCallback(hObject,~)
    ClearAppData('Y',hObject)
    ClearAppData('T',hObject)

    set(evalin('caller','tbox'),'string','Extracting OPTICS Partition')
    drawnow
    
    X1 = getappdata(hObject.Parent,'X1');
    X2 = getappdata(hObject.Parent,'X2');
    order = getappdata(hObject.Parent,'order');
    RD = getappdata(hObject.Parent,'RD');
    CD = getappdata(hObject.Parent,'CD');
    
    %Get epsilon value from user prompt (and display and store that value)
    prompt = {'Enter epsilon value for cluster extraction:'};
    dlg_title = 'Extraction Epsilon Input';
    defaultans = {'0.1'};
    usrinput = inputdlg(prompt,dlg_title,1,defaultans);
    eps = str2double(usrinput{1});
    setappdata(hObject.Parent,'eps',eps);
    set(evalin('caller','dbscan_eps'),'string',usrinput{1});
    drawnow
    
    geneps = getappdata(hObject.Parent,'geneps');
    [Y,T,CDbw] = ExtractOpticsPartition(X1,X2,RD,CD,order,eps,geneps);
    
    %Plot clustering solution to screen
    map = distinguishable_colors(max(Y(:))+2);
    map1 = distinguishable_colors(max(Y(:)));
    subplot(2,2,2);
    scatter3(X1(order,1),10.^(X1(order,2)),X1(order,3), 50, Y, 'filled')
    view(2)
    set(gca,'yscale','log')
    grid off
    if min(Y(:)) == -1
        caxis(caxis)
        colormap(map) 
        labels(1) = -1;
        for l = 2:(max(Y(:))+2)
            labels(l) = l-2;
        end
    else
        colormap(map1) 
        for l = 1:(max(Y(:)))
            labels(l) = l;
        end
    end
    lcolorbar(labels,'FontSize', 14,'fontweight','bold');
    title('OPTICS Extracted Clustering','FontWeight','bold','fontsize',18)
    ylabel('Conductance (G/Go)','Fontweight','bold','fontsize',14)
    xlabel('Distance (m)','FontWeight','bold','fontsize',14) % x-axis label
    set(gca,'fontweight','bold')
    
    val = getappdata(hObject.Parent,'val');
    setappdata(hObject.Parent,'val',[val,CDbw])
    setappdata(hObject.Parent,'Y',Y)
    setappdata(hObject.Parent,'T',T)
    set(evalin('caller','tbox'),'string',strcat('Partition Extracted:',{' '},'CDbw =',{' '},num2str(CDbw)))
end

%Displays the table T, containing number of points and percentage of points
%assigned to each cluster ID, on screen
function ViewTableCallback(hObject,~)
    set(evalin('caller','tbox'),'string','Creating Table')
    drawnow
    T = getappdata(hObject.Parent,'T');
    pos = get(subplot(2,2,4),'position');
    delete(subplot(2,2,4))
    setappdata(hObject.Parent,'table',uitable('data',table2array(T),'ColumnName',T.Properties.VariableNames,...
        'units','normalized',...
        'position',pos))
    
    set(evalin('caller','tbox'),'string','Table Created')  
end

%Takes only the data points belonging to user-specified cluster ID#s and
%makes a clustering plot and reachability plot from just those data (both
%in new windows)
function ExtractClustCallback(hObject,~)
    set(evalin('caller','tbox'),'string','Extracting Clusters')
    drawnow
    X1 = getappdata(hObject.Parent,'X1');
    Y = getappdata(hObject.Parent,'Y');
    order = getappdata(hObject.Parent,'order');
    RD = getappdata(hObject.Parent,'RD');
    ExtractClust(X1,Y,order,RD)
    set(evalin('caller','tbox'),'string','Clusters Extracted')
end

%Prompts user for a range of epsilon (eps) values, then re-runs
%"ExtractOpticsPartition" once for each of those epsilon values.  The
%cluster ID vector for each extraction is stored in the "parition" vector,
%while the CDbw value for each extraction is stored in the "index" vector.
%Creates two new windows when completed: one shows the clustering solution
%corresponding to the highest CDbw value, the other shows a plot of CDbw
%versus partition number.  
function CDbwTestCallback(hObject,~)
    ClearAppData('eps_range',hObject)
    ClearAppData('partitions',hObject)
    ClearAppData('index',hObject)

    set(evalin('caller','tbox'),'string','Performing CDbw Analysis')
    drawnow
    
    X1 = getappdata(hObject.Parent,'X1');
    X2 = getappdata(hObject.Parent,'X2');
    order = getappdata(hObject.Parent,'order');
    RD = getappdata(hObject.Parent,'RD');
    CD = getappdata(hObject.Parent,'CD');
    geneps = getappdata(hObject.Parent,'geneps');
    prompt = {'Enter Number of Epsilon Steps:','Enter Starting Epsilon Value:', 'Enter Ending Epsilon Value:', 'Enter Epsilon value to compare against:'};
    dlg_title = 'Input';
    defaultans = {'10','.05','.225','0.10'};
    usrinput = inputdlg(prompt,dlg_title,1,defaultans);
    n = str2double(usrinput{1});
    epsi = str2double(usrinput{2});
    epsf = str2double(usrinput{3});
    eps = linspace(epsi,epsf,n);
    
    eps_compare = str2double(usrinput{4});
    [~,~,CDbw_compare] = ExtractOpticsPartition(X1,X2,RD,CD,order,eps_compare,geneps);
    
    [partitions,index] = CDbwAnalysis(X1,X2,RD,CD,order,geneps,eps,n);
    setappdata(hObject.Parent,'eps_range',eps)
    setappdata(hObject.Parent,'partitions',partitions)
    setappdata(hObject.Parent,'index',index)
    [~,I] = max(index);
    eps_str = num2str(eps(I));
    set(evalin('caller','tbox'),'string',strcat('Optimal Partition was Extracted at epsilon =',{' '},eps_str))
    
    t = sum(index >= CDbw_compare);
    message = strcat(num2str(t),{' '}, 'out of', {' '}, num2str(n),...
    {' '},'epsilon values in range have CDbw values >= CDbw value for comparison epsilon (i.e. eps =',...
    {' '}, num2str(eps_compare),')');
    
    d = dialog('Position',[300 300 250 150],'Name','My Dialog');
    txt = uicontrol('Parent',d,'Style','text','Position',[20 80 210 40],...
               'String',message);
    btn = uicontrol('Parent',d,'Position',[85 20 70 25],'String','OK',...
               'Callback','delete(gcf)');
end

%Calculates the gamma statistic for the clustering solution displayed in the
%main GUI window.  In addition, ranomly permutes cluster membership to
%derive a Monte Carlo simulation of the sampling distribution for the gamma
%statistic, and from that distribution calcualtes a (one-sided) p-value for
%the gamma of the clustering solution shown in the mian GUI window.  
function GammaTestCallback(hObject,~)
    
    set(evalin('caller','tbox'),'string','Calculating Gamma Statistic')
    drawnow
    
    %Get # of MC samples from user
    prompt = {'Enter # of Monte Carlo samples to use for Gamma Statistic distribution (recommended value 100):'};
    dlg_title = 'Enter Monte Carlo Samples';
    defaultans = {'100'};
    usrinput = inputdlg(prompt,dlg_title,1,defaultans);
    n_MC = str2double(usrinput{1});
    
    Y = getappdata(hObject.Parent,'Y');
    X2 = getappdata(hObject.Parent,'X2');
    [gamma,gammadist] = gammatest(Y,X2,n_MC);
    setappdata(hObject.Parent,'gamma',gamma)
    setappdata(hObject.Parent,'gammadist',gammadist)
    t = sum(gammadist >= gamma);
    gammastr = num2str(gamma);
    set(evalin('caller','tbox'),'string',strcat('Gamma =', {' '},...
        gammastr,';',{' '},num2str(t),{' '},'out of',{' '},num2str(n_MC),...
        {' '},'MC samples have higher Gamma values'));
end

%Create new window displaying a 2D histogram of the
%distance-log(conductance) data
function Print_2D_HistCallback(hObject,~)

    prompt = {'Enter Image Title: ','Enter x-axis label: ','Enter y-axis label: '};
    dlg_title = 'Input';
    defaultans = {'2D Conductance-Displacement Histogram','Displacement (nm)','Conductance (G/G_o)'};
    usrinput = inputdlg(prompt,dlg_title,1,defaultans);
    usrtitle = usrinput{1};
    usrx = usrinput{2};
    usry = usrinput{3};
    X1 = getappdata(hObject.Parent,'X1');
    Print2DHist(X1,usrtitle,usrx,usry);
    
end

%Opens a new window displaying the clustering solution most recently
%extracted (not counting extractions performed within CDbw Analysis)
function Print_ClustCallback(hObject,~)
    prompt = {'Enter Title: ', 'Enter x-axis label: ','Enter y-axis label: '};
    dlg_title = 'Input';
    defaultans = {'Clustering','Displacement (nm)','Conductance (G/Go)'};
    usrinput = inputdlg(prompt,dlg_title,1,defaultans);
    usrtitle = usrinput{1};
    usrx = usrinput{2};
    usry = usrinput{3};
    X1 = getappdata(hObject.Parent,'X1');
    order = getappdata(hObject.Parent,'order');
    Y = getappdata(hObject.Parent,'Y');
    X1 = X1(order,:);
    PrintClust(X1,Y,usrtitle,usrx,usry)
end

%Creates a new window displaying the reachability plot
function Print_OpticsCallback(hObject,~)
    
    usrtitle = inputdlg('Enter Image Title');
    order = getappdata(hObject.Parent,'order');
    RD = getappdata(hObject.Parent,'RD');
    PrintRP(order,RD,usrtitle)
    
end

%Allows user to specify a partition ID corresponding to one of the
%extracted partitions calculated by CDbwAnalysis.  Then creates a new
%window displaying that clustering solution, and prints the CDbw and
%epsilon value associated with that solution to the status bar.  
function Print_CDbwCallback(hObject,~) 
    X1 = getappdata(hObject.Parent,'X1');
    order = getappdata(hObject.Parent,'order');
    eps = getappdata(hObject.Parent,'eps_range');
    partitions = getappdata(hObject.Parent,'partitions');
    [partition_index] = PrintCDbwPartition(X1,order,partitions);
    eps_str = num2str(eps(partition_index));    
    index = getappdata(hObject.Parent, 'index');
    cdbw_str = num2str(index(partition_index));
    set(evalin('caller','tbox'),'string',strcat('Partition was Extracted at epsilon =',{' '},eps_str, {', '}, 'CDbw =', {' '},cdbw_str))
end