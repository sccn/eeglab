function [g, c, tmp, m]=kMeansCluster(m,k,isRand)
% kMeansCluster - Simple k means clustering algorithm                                                              
% Author: Kardi Teknomo, Ph.D.                                                                  
%                                                                                                                    
% Purpose: classify the objects in data matrix based on the attributes    
% Criteria: minimize Euclidean distance between centroids and object points                    
% For more explanation of the algorithm, see http://people.revoledu.com/kardi/tutorial/kMean/index.html    
% Output: matrix data plus an additional column represent the group of each object               
%                                                                                                                
% Example: m = [ 1 1; 2 1; 4 3; 5 4]  or in a nice form                         
%          m = [ 1 1;                                                                                     
%                2 1;                                                                                         
%                4 3;                                                                                         
%                5 4]                                                                                         
%          k = 2                                                                                             
% kMeansCluster(m,k) produces m = [ 1 1 1;                                        
%                                   2 1 1;                                                                   
%                                   4 3 2;                                                                   
%                                   5 4 2]                                                                   
% Input:
%   m      - required, matrix data: objects in rows and attributes in columns                                                 
%   k      - optional, number of groups (default = 1)
%   isRand - optional, if using random initialization isRand=1, otherwise input any number (default)
%            it will assign the first k data as initial centroids
%
% Local Variables
%   f      - row number of data that belong to group i
%   c      - centroid coordinate size (1:k, 1:maxCol)
%   g      - current iteration group matrix size (1:maxRow)
%   i      - scalar iterator 
%   maxCol - scalar number of rows in the data matrix m = number of attributes
%   maxRow - scalar number of columns in the data matrix m = number of objects
%   temp   - previous iteration group matrix size (1:maxRow)
%   z      - minimum value (not needed)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin<3,        isRand=0;   end
if nargin<2,        k=1;        end
tmp = [];

[maxRow, maxCol]=size(m);
if maxRow<=k, 
    y=[m, 1:maxRow];
else
	
	% initial value of centroid
    if isRand,
        p = randperm(size(m,1));      % random initialization
        for i=1:k
            c(i,:)=m(p(i),:)  
    	end
    else
        for i=1:k
           c(i,:)=m(i,:);        % sequential initialization
    	end
    end
    
	temp=zeros(maxRow,1);   % initialize as zero vector
    
	while 1,
        d=DistMatrix(m,c);  % calculate objcets-centroid distances
        [z,g]=min(d,[],2);  % find group matrix g
        if g==temp,
            break;          % stop the iteration
        else
            temp=g;         % copy group matrix to temporary variable
        end
        for i=1:k
            f=find(g==i);
            if f            % only compute centroid if f is not empty
                c(i,:)=mean(m(find(g==i),:),1);
            end
        end
	end
    
	y=[m,g];
    
end

%The Matlab function kMeansCluster above call function DistMatrix as shown in the code below. It works for multi-dimensional Euclidean distance. Learn about other type of distance here.

    function d=DistMatrix(A,B)
                 %%%%%%%%%%%%%%%%%%%%%%%%%
                 % DISTMATRIX return distance matrix between points in A=[x1 y1 ... w1] and in B=[x2 y2 ... w2]
                 % Copyright (c) 2005 by Kardi Teknomo,  http://people.revoledu.com/kardi/
                 %
                 % Numbers of rows (represent points) in A and B are not necessarily the same.
                 % It can be use for distance-in-a-slice (Spacing) or distance-between-slice (Headway),
                 %
                 % A and B must contain the same number of columns (represent variables of n dimensions),
                 % first column is the X coordinates, second column is the Y coordinates, and so on.
                 % The distance matrix is distance between points in A as rows
                 % and points in B as columns.
                 % example: Spacing= dist(A,A)
                 % Headway = dist(A,B), with hA ~= hB or hA=hB
                 %          A=[1 2 3; 4 5 6; 2 4 6; 1 2 3]; B=[4 5 1; 6 2 0]
                 %          dist(A,B)= [ 4.69   5.83;
                 %                       5.00   7.00;
                 %                       5.48   7.48;
                 %                       4.69   5.83]
                 %
                 %          dist(B,A)= [ 4.69   5.00     5.48    4.69;
                 %                       5.83   7.00     7.48    5.83]
                 %%%%%%%%%%%%%%%%%%%%%%%%%%%
                 [hA,wA]=size(A);
                 [hB,wB]=size(B);
                 if wA ~= wB,  error(' second dimension of A and B must be the same'); end
                 for k=1:wA
                      C{k}= repmat(A(:,k),1,hB);
                      D{k}= repmat(B(:,k),1,hA);
                 end
                 S=zeros(hA,hB);
                 for k=1:wA
                      S=S+(C{k}-D{k}').^2;
                 end
                 d=sqrt(S);

