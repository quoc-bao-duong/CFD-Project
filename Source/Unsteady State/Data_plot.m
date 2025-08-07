function Data_plot(PHI,delta_x,k)
%{

%}
% Create data
    x=0:delta_x:1;
    y=0:delta_x:1;
    n=length(x);
    IMax=n;
    JMax=n;
    [X, Y] = meshgrid(x, y);

    PHI = flipud(PHI); 

    % Output file
    if k==0
        filename ='Results_Initial_PHI.dat';
    end

    if k==2
        filename ='Results_CENTRAL.dat';
    end



    fileID = fopen(filename, 'w');

    % Write header
    fprintf(fileID, 'TITLE = "IJ-Ordered Data with Cell-Centered Variables"\n');
    fprintf(fileID, 'VARIABLES = "X", "Y", "Temperature" \n');
    fprintf(fileID, 'ZONE I=%d, J=%d, DATAPACKING=BLOCK, VARLOCATION=(3=CELLCENTERED)\n', IMax, JMax);

    % Write nodal data in BLOCK format
    % Write X-coordinates block
    for j = 1:JMax
        for i = 1:IMax
            fprintf(fileID, '%f ', X(j, i));
        end
        fprintf(fileID, '\n');
    end
    
    % Write Y-coordinates block
    for j = 1:JMax
        for i = 1:IMax
            fprintf(fileID, '%f ', Y(j, i));
        end
        fprintf(fileID, '\n');
    end

    % Write cell-centered data in BLOCK format
    % Cell-centered data has size (IMax-1) x (JMax-1)
    for j = 1:(JMax-1)
        for i = 1:(IMax-1)
            fprintf(fileID, '%f ', PHI(j, i));
        end
        fprintf(fileID, '\n');
    end
    
    % Close the file
    fclose(fileID);
end