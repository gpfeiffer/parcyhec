# the complex reflection group.
W:= ComplexReflectionGroup(29);

# the relations.
relations:= [
            [[1,2,1], [2,1,2]],
            [[1,3], [3,1]],
            [[1,4], [4,1]],
            [[2,3,2,3], [3,2,3,2]],
            [[2,4,2], [4,2,4]],
            [[3,4,3], [4,3,4]],
            [[3,2,4,3,2,4], [4,3,2,4,3,2]],
            ];
            
# the parabolic Coxeter subgroup, and how it relates to W.
subset:= [3,2,1];
W1:= CoxeterGroup("B", 3);
q:= X(Rationals);  q.name:="q";
H1:= Hecke(W1, q);

mmm:= HeckeMatrices(W, relations, H1, subset, 2);;
Print("Runtime: ", Runtime(), "ms\n");
