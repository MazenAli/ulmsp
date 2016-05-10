clear all

load 'coordinates_fine.dat'
load 'elements_fine.dat'
load 'sol_fine.dat'

trisurf(elements_fine,coordinates_fine(:,1), ...
                      coordinates_fine(:,2), ...
                      sol_fine(:,1),'facecolor','interp');
view(3);
