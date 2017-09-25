# the complex reflection group.
W:= CoxeterGroup("A", 4);
relations1:= [
             [[1,2,1],[2,1,2]],
             [[2,3,2],[3,2,3]],
              [[1,3],[3,1]],
              [[3,4,3],[4,3,4]],
              [[2,4],[4,2]],
              [[1,4],[4,1]],
              ];

# the parabolic Coxeter subgroup, and how it relates to W.
subset1:= [1,2,3];
W1:= CoxeterGroup("A", 3);

# the relations.
relations2:= [
             [[1,2,1],[2,1,2]],
             [[2,3,2],[3,2,3]],
             [[1,3],[3,1]],
            ];
            
# the parabolic Coxeter subgroup, and how it relates to W.
subset2:= [1,2];
W2:= CoxeterGroup("A", 2);

# the relations.
relations3:= [
             [[1,2,1],[2,1,2]],
            ];
            
# the parabolic Coxeter subgroup, and how it relates to W.
subset3:= [1];
W3:= CoxeterGroup("A", 1);
q:= X(Rationals);  q.name:="q";
H3:= Hecke(W3, q);

H2:= HeckeMatrices(W2, relations3, H3, subset3, 1);;
H2.name:= "H2";
Print("Runtime: ", Runtime(), "ms\n");

H1:= HeckeMatrices(W1, relations2, H2, subset2, 1);;
H1.name:= "H1";
Print("Runtime: ", Runtime(), "ms\n");

H0:= HeckeMatrices(W, relations1, H1, subset1, 1);;
H0.name:= "H0";
Print("Runtime: ", Runtime(), "ms\n");
