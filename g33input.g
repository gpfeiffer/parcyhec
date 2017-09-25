# the complex reflection group.
W:= ComplexReflectionGroup(33);

# the relations.
relations:= [
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
subset:= [1,2,4,5];
W1:= CoxeterGroup("A", 4);
q:= X(Rationals);  q.name:="q";
H1:= Hecke(W1, q);

mmm:= HeckeMatrices(W, relations, H1, subset, 1);;
Print("Runtime: ", Runtime(), "ms\n");

