function IndependentEdgeVector = hl_MatrixToIndependentEdgeVector(Matrix)


Temp=[];
IndependentEdgeVector=[];

for i=1:size(Matrix,1)

    Temp=Matrix(i,i+1:end);
    IndependentEdgeVector=[IndependentEdgeVector,Temp];

end