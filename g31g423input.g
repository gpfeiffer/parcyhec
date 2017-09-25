# the complex reflection group.
W:= ComplexReflectionGroup(31);

# the relations.
relations1:= [
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
                  [[2,3,-2,4,1,2,4],[4,1,2,4,-2,3,2]], # <- essential
            ];

# the parabolic Coxeter subgroup, and how it relates to W.
#subset1:= [2,1,3,4];
subset1:= [1,2,3,5];

W1:= ComplexReflectionGroup(4,2,3);

# the relations.
relations2:= [
                  [[3,4,3],[4,3,4]],
                  [[2,4,2],[4,2,4]],
                  [[1,4],[4,1]],
                  [[1,2,3],[2,3,1]],
                  [[2,3,1],[3,1,2]],
                  [[1,2,3],[3,1,2]],
            ];

# the parabolic Coxeter subgroup, and how it relates to W.
subset2:= [2,4];
W2:= CoxeterGroup("A", 2);
q:= X(Rationals);  q.name:="q";
H2:= Hecke(W2, q);

H1:= HeckeMatrices(W1, relations2, H2, subset2, 1);;
H1.name:= "H1";
Print("Runtime: ", Runtime(), "ms\n");


H0:= HeckeMatrices(W, relations1, H1, subset1, 1);;
H0.name:= "H0";
Print("Runtime: ", Runtime(), "ms\n");
