# the complex reflection group.
W:= ComplexReflectionGroup(12);

# the relations.
relations:= [
                  [ [ 1, 2, 3, 1 ], [ 2, 3, 1, 2 ] ], 
                  [ [ 1, 2, 3, 1 ], [ 3, 1, 2, 3 ] ],
            ];
            
# the parabolic Coxeter subgroup, and how it relates to W.
subset:= [1];
W1:= CoxeterGroup("A", 1);
q:= X(Rationals);  q.name:="q";
H1:= Hecke(W1, q);

mmm:= HeckeMatrices(W, relations, H1, subset, 2);;
Print("Runtime: ", Runtime(), "ms\n");
