# get the complex reflection group.
W:= ComplexReflectionGroup(34);  ##  <- INPUT

# the real reflection subgroup
subset:= [1,2,4,5,6];   ##  <- INPUT

# the relations.
relations:= [     #  <- INPUT !
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
W1:= CoxeterGroup("A", 5);   ##  <- INPUT

q:= X(Rationals);  q.name:="q";
H1:= Hecke(W1, q);

mmm:= HeckeMatrices(W, relations, H1, subset, 2);;
Print("Runtime: ", Runtime(), "ms\n");
