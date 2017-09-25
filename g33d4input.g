#(Le but est de faire G33 a partir de D4)
# the complex reflection group.
W:= ComplexReflectionGroup(33);
pp := [W.1,W.2,W.4,W.3,W.5^((W.3*W.4)^-1)];
W:= Group(pp,());
W.generatingReflections:= [1..Length(pp)];

# the relations. (G33-D4)(celles de l'article de Bessis-Michel)
relations := [
	[[1,2,1],[2,1,2]],
	[[1,3],[3,1]],
	[[1,4],[4,1]],
	[[1,5],[5,1]],
	[[2,3,2],[3,2,3]],
	[[2,4,2],[4,2,4]],
	[[2,5,2],[5,2,5]],
	[[3,4,3],[4,3,4]],
	[[3,5],[5,3]],
	[[4,5,4],[5,4,5]],
	[[5,4,3,2,5,4],[4,3,2,5,4,3]], # 424
#	[[3,2,4,3,2,4],[4,3,2,4,3,2]],
#	[[3,2,4,3,2,4],[2,4,3,2,4,3]],
#	[[4,3,2,4,3,2],[2,4,3,2,4,3]], #410
#	[[ 4, 2, 1, 5, 4, 2, 1 ],[ 2, 5, -2, 4, 2, 1, 5, 4, 2 ]], # 336
        [[4,2,5,4,3,2,5],[3,2,5,4,3,2,4]],
	];
	
# the parabolic Coxeter subgroup, and how it relates to W.
subset:= [1,3,2,5];
W1:= CoxeterGroup("D", 4);
q:= X(Rationals);  q.name:="q";
H1:= Hecke(W1, q);

mmm:= HeckeMatrices(W, relations, H1, subset, 1);;
Print("Runtime: ", Runtime(), "ms\n");
