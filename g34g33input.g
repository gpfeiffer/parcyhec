# the complex reflection group.
W:= ComplexReflectionGroup(34);

# the relations.
relations1:= [
            [ [ 1, 2, 1 ], [ 2, 1, 2 ] ], [ [ 3, 2, 3 ], [ 2, 3, 2 ] ], 
  [ [ 4, 2, 4 ], [ 2, 4, 2 ] ], [ [ 4, 3, 4 ], [ 3, 4, 3 ] ], 
  [ [ 4, 5, 4 ], [ 5, 4, 5 ] ], [ [ 1, 3 ], [ 3, 1 ] ], 
  [ [ 1, 4 ], [ 4, 1 ] ], [ [ 1, 5 ], [ 5, 1 ] ], [ [ 2, 5 ], [ 5, 2 ] ], 
  [ [ 3, 5 ], [ 5, 3 ] ], [ [ 5, 6, 5 ], [ 6, 5, 6 ] ], 
  [ [ 1, 6 ], [ 6, 1 ] ], [ [ 2, 6 ], [ 6, 2 ] ], [ [ 3, 6 ], [ 6, 3 ] ], 
  [ [ 4, 6 ], [ 6, 4 ] ], [ [ 4, 2, 3, 4, 2, 3 ], [ 3, 4, 2, 3, 4, 2 ] ], 
[  [ 3, 4, 2, 3, 4, 2 ],[ 2, 3, 4, 2, 3, 4 ] ], 
            ];
            
# the parabolic Coxeter subgroup, and how it relates to W.
subset1:= [1,2,3,4,5];
W1:= ComplexReflectionGroup(33);

# the relations.
relations2:= [
            [[1,2,1], [2,1,2]],
            [[1,3],[3,1]],
            [[1,4],[4,1]],
            [[1,5],[5,1]],
            [[2,3,2],[3,2,3]],
            [[2,4,2],[4,2,4]],
            [[2,5],[5,2]],
            [[3,4,3],[4,3,4]],
            [[3,5],[5,3]],
            [[4,5,4],[5,4,5]],
            [[4,2,3,4,2,3],[3,4,2,3,4,2]],
            [[3,4,2,3,4,2],[2,3,4,2,3,4]],
            ];
            
# the parabolic Coxeter subgroup, and how it relates to W.
subset2:= [1,2,4,5];
W2:= CoxeterGroup("A", 4);
q:= X(Rationals);  q.name:="q";
H2:= Hecke(W2, q);

H1:= HeckeMatrices(W1, relations2, H2, subset2, 1);;
H1.name:= "H1";
Print("Runtime: ", Runtime(), "ms\n");


H0:= HeckeMatrices(W, relations1, H1, subset1, 1);;
Print("Runtime: ", Runtime(), "ms\n");
