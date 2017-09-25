# the complex reflection group.
W:= CoxeterGroup("A", 2);

# the relations.
relations:= [
             [[1,2,1],[2,1,2]],
            ];
            
# the parabolic Coxeter subgroup, and how it relates to W.
subset:= [1];
W1:= CoxeterGroup("A", 1);
q:= X(Rationals);  q.name:="q";
H1:= Hecke(W1, q);

mmm:= HeckeMatrices(W, relations, H1, subset, 2);;
Print("Runtime: ", Runtime(), "ms\n");
