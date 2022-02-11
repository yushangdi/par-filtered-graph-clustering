
These matlab codes are the implementation of Directed Bubble Hierarchical Tree (DBHT) technique as introduced in the: "Hierarchical information clustering by means of topologically embedded graphs", Plos One. 

The DBHT code, makes use of the graph library package 'Matlab BGL', which uses Boost Graph Library and Statistics Toolbox, Under Boost Software License (http://www.mathworks.com/matlabcentral/fileexchange/10922). The Matlab BGL package must be installed and its path should be added (for instance by using: addpath('matlab_bgl')).
Note that in some platforms (i.e. 64-bit Mac's with R2009b or higher) the .mex files must be re-compiled (see: http://dgleich.wordpress.com/2010/07/08/matlabbgl-osx-64-bit/)

----------------Clustering with DBHT-----------------------

Discrete clustering, requires two inputs: 

- N x N similarity matrix S, 
- N x N dissimilarity matrix D. 

The function call to perform clustering is: 

T=DBHT(D,R);

where T is N x 1 column vector containing cluster membership of vertices. T(n)=k indicates cluster membership of nth vertex to kth cluster.

-------Hierarchical structure by DBHT------------------

To obtain intra- and inter-cluster hierarchical structures, the function calls are:

[T8,Rpm,Adjv,Dpm,Mv]=DBHT(D,R);
Z=HierarchyConstruct4(Rpm,Dpm,T8,Adjv,Mv);

The hierarchical structure is contained in the output 'Z', a (N-1) x N linkage matrix, and follows the same format as the output from matlab built-in funcion 'linkage'. To view the hierarchy, call:

dendrogram(Z,0);


%-----------------------------------------------------------------------------------------
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should receive a copy of the GNU General Public License
% along with this program.  See also <http://www.gnu.org/licenses/>.
%-----------------------------------------------------------------------------------------


