function load_dm

% Copyright 2007 John L Baker. All rights reserved.
% This software is provided AS IS under the terms of the Open Source
% MIT License. See http://www.opensource.org/licenses/mit-license.php.

global SRCDIR NCOMP DM

% Load dendrite morphology table into a global
% DM columns are:
% 1 type
% 2 id (should be array index-1)
% 3 parent
% 4 branch id
% 5 radius (microns)
% 6 length (microns)
% 7 distance from soma (microns)
% 8 x coordinate
% 9 y coordinate
% 10 z coordinate

if size(DM,1)>0
   return
end
DM=[];

% Identify the appropriate morphology table
% based on the number of compartments in the model.

if NCOMP==26
   DM=dlmread([SRCDIR 'cell_pyramidal_ball_and_stick.csv'],',',1,0);
end
if NCOMP==368
   DM=dlmread([SRCDIR 'cell_l56a_50_micron.csv'],',',1,0);
end
if NCOMP==436
   DM=dlmread([SRCDIR 'cell_n123_50_micron.csv'],',',1,0);
end
if NCOMP==631
   DM=dlmread([SRCDIR 'cell_l56a_25_micron.csv'],',',1,0);
end
if NCOMP==2776
   DM=dlmread([SRCDIR 'cell_l56a_5_micron.csv'],',',1,0);
end
if NCOMP==504
   DM=dlmread([SRCDIR 'cell_l51_50_micron.csv'],',',1,0);
end
if NCOMP==842 
   DM=dlmread([SRCDIR 'cell_l51_25_micron.csv'],',',1,0);
end
if NCOMP==1918 
   DM=dlmread([SRCDIR 'cell_l51_10_micron.csv'],',',1,0);
end
if NCOMP==440
   DM=dlmread([SRCDIR 'cell_ml02_50_micron.csv'],',',1,0);
end

if size(DM,1)>0
   % make distance negative for basal dendrites
   bi=find(DM(:,1)==3);
   DM(bi,7)=-abs(DM(bi,7));
end
