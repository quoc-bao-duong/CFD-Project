function Mesh(delta_x,n)
%{
 INPUT
 delta_x,n
%}
    % Create data
    x=0:delta_x:1;
    y=0:delta_x:1;
    np=length(x); % number of points along x-axis or y-axis
    % Create meshgrid for X, Y coordinates
    [X, Y] = meshgrid(x, y);
    % Open a file for writing
    fileID = fopen('mesh.dat', 'w');
    % Write Tecplot header
    fprintf(fileID, 'TITLE = "Mesh"\n');
    fprintf(fileID, 'VARIABLES = "X", "Y"\n');
    fprintf(fileID, 'ZONE T="%d x %d Cells", I=%d, J=%d, F=POINT\n',n,n,np,np);
    % Write the data in POINT format
    for j = 1:np
        for i = 1:np
            fprintf(fileID, '%f %f\n', X(j, i), Y(j, i));
        end
    end

    % Close the file
    fclose(fileID);
end