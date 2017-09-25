# the complex reflection group.
W:= ComplexReflectionGroup(31);

# the relations.
relations:= [
                  [[1,4,1],[4,1,4]],
                  [[1,5],[5,1]],
                  [[2,4,2],[4,2,4]],
                  [[2,5,2],[5,2,5]],
                  [[3,4],[4,3]],
                  [[3,5,3],[5,3,5]],
                  [[4,5],[5,4]],
                  [[1,2,3],[2,3,1]],
                  [[2,3,1],[3,1,2]],
                  [[1,2,3],[3,1,2]],
                  ];
            
# the parabolic Coxeter subgroup, and how it relates to W.
subset:= [4,2,5];
W1:= CoxeterGroup("A", 3);
q:= X(Rationals);  q.name:="q";
H1:= Hecke(W1, q);

mmm:= HeckeMatrices(W, relations, H1, subset, 1);;
Print("Runtime: ", Runtime(), "ms\n");

