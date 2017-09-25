# the complex reflection group.
W:= ComplexReflectionGroup(22);

# the relations.
relations:= [
                  [[1,2,3,1,2],[2,3,1,2,3]],
                  [[2,3,1,2,3],[3,1,2,3,1]],
                  [[1,2,3,1,2],[3,1,2,3,1]],
            ];
            
# the parabolic Coxeter subgroup, and how it relates to W.
subset:= [1];
W1:= CoxeterGroup("A", 1);
q:= X(Rationals);  q.name:="q";
H1:= Hecke(W1, q);

mmm:= HeckeMatrices(W, relations, H1, subset, 1);;
Print("Runtime: ", Runtime(), "ms\n");
