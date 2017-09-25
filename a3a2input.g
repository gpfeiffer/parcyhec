# the complex reflection group.
W:= CoxeterGroup("A", 3);

# the relations.
relations1:= [
             [[1,2,1],[2,1,2]],
             [[2,3,2],[3,2,3]],
             [[1,3],[3,1]],
            ];
            
# the parabolic Coxeter subgroup, and how it relates to W.
subset1:= [1,2];
W1:= CoxeterGroup("A", 2);

# the relations.
relations2:= [
             [[1,2,1],[2,1,2]],
            ];
            
# the parabolic Coxeter subgroup, and how it relates to W.
subset2:= [1];
W2:= CoxeterGroup("A", 1);
q:= X(Rationals);  q.name:="q";
H2:= Hecke(W2, q);

H1:= HeckeMatrices(W1, relations2, H2, subset2, 2);;
H1.name:= "H1";
Print("Runtime: ", Runtime(), "ms\n");


H0:= HeckeMatrices(W, relations1, H1, subset1, 2);;
Print("Runtime: ", Runtime(), "ms\n");
