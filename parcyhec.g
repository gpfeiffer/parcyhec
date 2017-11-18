#############################################################################
##
#A  parcyhec.g                                    goetz.pfeiffer@nuigalway.ie
##
##  represent a cyclotomic Hecke algebra as matrices over a parabolic.
##

#############################################################################
##
##  the main procedure
##
##  Input:
##
##    W -- complex reflection group
##    relations -- a list of relations
##    W1 -- the parabolic subgroup
##    subset -- injection of W1 into W
##    strategy -- 1 or 2: cosets or double cosets
##
HeckeMatrices:= function(W, relations, H1, subset, strategy)
    local   simples,  H,  T,  q,  one,  zero,  cosets,
            words,  rules,  process,  i,  s,  complement,  j,  t,
            imageFromRules,  iii,  mmm,  inverseVec,  inverse2Vec,
            uuu,  cacheVec,  index,  rule,  variants,  unresolved,
            open,  n,  v,  new,  x1,  ibraid,  inverse,  mark,  inv,  hecke;

    # setup from input data
    simples:= W.generatingReflections;
    H:= Subgroup(W, W.generators{subset});

    T:= Basis(H1, "T");
    q:= X(Rationals);

    # name the constants 0 and 1 in H1
    zero:= Zero(H1);
    one:= One(H1);

    # use orbit algorithm to find coset reps (as words in the generators)
    # and further properties.
    cosets:= [H * ()];  # the actual list of cosets
    words:= [[]];       # coset reps as words in the generators
    rules:= [];         # rules[i] is a list of records with components
                        # s, x1, u stating that i.s = u * x1,
                        # where u = "new" marks an edge in the spanning tree.

    # how to process k.s
    process:= function(k, s)
        local   new,  pos;

        new:= cosets[k] * W.(s);
        pos:= Position(cosets, new);
        if pos = false then
            Add(cosets, new);
            Add(words, Concatenation(words[k], [s]));
            Add(rules[k], rec(s:= s, x1:= Length(cosets), u:= "new"));
        else
            Add(rules[k], rec(s:= s, x1:= pos, u:= "old"));
        fi;
    end;

    # list cosets, or double cosets, according to strategy
    if strategy = 1 then
        i:= 0;
        while i < Length(cosets) do
            i:= i+1;
            rules[i]:= [];
            for s in simples do
                process(i, s);
            od;
        od;
    elif strategy = 2 then
        complement:= Difference(simples, subset);
        i:= 0; j:= 0;
        while j < Length(cosets) do
            j:= j+1;
            for t in complement do
                while i < Length(cosets) do
                    i:= i+1;
                    rules[i]:= [];
                    for s in subset do
                        process(i, s);
                    od;
                od;
                process(j, t);
            od;
        od;
    else
        Error("unknown strategy");
    fi;

    # how to find the index of coset n.s
    imageFromRules:= function(n, s)
        return rules[n][PositionProperty(rules[n], x-> x.s = s)].x1;
    end;

    ##  matrices  (uses globals zero, one, iii, mmm, ...)
    iii:= [1..Length(cosets)]; # matrix indices

    # the matrices, one for each generator
    mmm:= List(simples, x-> []);
    hecke:= HeckeMats(H1, subset, mmm, words);

    # if x_l.s = alpha x_n - (q-1) beta
    # then x_n.s' = alpha^-1 * (x_l + (q-1) beta.s')
    inverseVec:= function(n, s, inverse, l)
        local   vec,  alpha;

        vec:= mmm[s][l];
        alpha:= Get(vec, n);
        if inverse * alpha <> one then
            Error("inverse expected");
        fi;
        vec:= HVec([n], [alpha], zero, hecke) - vec;
        if vec.poss <> [] then
            vec:= HVecUnderWord(vec, [-s], mmm);
            if IsBool(vec) then
                return false;
            fi;
        fi;
#        return inverse * (vec + HVec([l], [q*one], zero, hecke))
        return HVecScaled(inverse, vec + HVec([l], [q*one], zero, hecke))
                      + HVec([n], [(q-1)*one], zero, hecke);
    end;

    # if x_l.s = alpha x_n + (q-1) (x_l + beta)
    # then x_n.s = q alpha' x_l * - (q-1) alpha' beta.s)
    inverse2Vec:= function(n, s, inverse, l)
        local   vec,  alpha;

        vec:= mmm[s][l];
        alpha:= Get(vec, n);
        if inverse * alpha <> one then
            Error("inverse expected");
        fi;
        vec:= HVec([n, l], [alpha, (q-1)*one], zero, hecke) - vec;
        if vec.poss <> [] then
            vec:= HVecUnderMat(vec, mmm[s]);
            if IsBool(vec) then
                return false;
            fi;
        fi;
#        return inverse * (vec + HVec([l], [q*one], zero, hecke));
        return HVecScaled(inverse, vec + HVec([l], [q*one], zero, hecke));
    end;

    ##  Compute the Matrices

    uuu:= [];   # the entries cache

    # how to cache the entries of  vec
    cacheVec:= function(vec)
        local   uu1,  val,  i;

        if IsBound(vec.vals[1].hecke.uuu) then
            uu1:= vec.vals[1].hecke.uuu;
            for val in vec.vals do
                UniteSet(uu1, Set(val.vals));
                for i in [1..Length(val.vals)] do
                    val.vals[i]:= uu1[Position(uu1, val.vals[i])];
                od;
            od;
        fi;

        UniteSet(uuu, Set(vec.vals));
        for i in [1..Length(vec.vals)] do
            vec.vals[i]:= uuu[Position(uuu, vec.vals[i])];
        od;
    end;


    # 1. spanning tree and inverses.
    for index in iii do
        for rule in rules[index] do
            if rule.u = "new" then
                mmm[rule.s][index]:= HVec([rule.x1], [one], zero, hecke);
                cacheVec(mmm[rule.s][index]);
                Print("spanning: ", index, ".", rule.s, " = ", rule.x1, "\n");
                mmm[rule.s][rule.x1]:= inverseVec(rule.x1, rule.s, one, index);
                cacheVec(mmm[rule.s][rule.x1]);
                Print("inversion: ", rule.x1, ".", rule.s, "\n");
            fi;
        od;
    od;

    # 1a. images of coset 1
    for i in [1..Length(subset)] do
        s:= subset[i];
        mmm[s][1]:= HVec([1], [T(i)], zero, hecke);
        cacheVec(mmm[s][1]);
    od;

    # 2. relations.
    variants:= VariantsRelations(relations);
    unresolved:= [];

    repeat
        open:= Sum(simples, s-> Number(iii, n-> not IsBound(mmm[s][n])));
        Print("OPEN: ", open, "\n");

        # loop over empty slots.
        for n in iii do
            for s in simples do

                # loop over all variants
                for v in variants.(s) do
                    if not IsBound(mmm[s][n]) then
                        new:= HVecUnderWord(HVec([n], [one], zero, hecke), v, mmm);
                        if not IsBool(new) then
                            Print("\n", n, ".", s, " = ", n, ".", v, "\n");
#                            mmm[s][n]:= q^(-Number(v, x-> x < 0)) * new;
                            mmm[s][n]:= HVecScaled(q^(-Number(v, x-> x < 0)), new);
                            cacheVec(mmm[s][n]);

                            # check for inverse edge
                            x1:= imageFromRules(n, s);
                            if x1 <> n then
                                if IsBound(mmm[s][x1]) then
                                    Print(n, ".", s, " has a known inverse ", x1, "\n");
                                else
                                    ibraid:= FindInverse(Get(mmm[s][n], x1));
                                    if IsBool(ibraid) then
                                        Print("cannot invert ", x1, ".", s, ": ", Get(mmm[s][n], x1), "\n");
                                    else
                                        if Get(mmm[s][n], n) = (q-1)*one then
                                            Print("***\c");
                                            inverse:= inverse2Vec(x1, s, ibraid, n);
                                            mark:= 2;
                                        else
                                            inverse:= inverseVec(x1, s, ibraid, n);
                                            mark:= 1;
                                        fi;
                                        if IsBool(inverse) then
                                            Print("inversion ", x1, ".", s, " failed\n");
                                            Add(unresolved, [x1, s, ibraid, n, mark]);
                                        else
                                            mmm[s][x1]:= inverse;
                                            cacheVec(mmm[s][x1]);
                                            Print("inversion: ", x1, ".", s, "\n");
                                        fi;
                                    fi;
                                fi;
                            fi;
                        else
                            Print(".\c");
                        fi;
                    fi;
                od;
            od;
        od;

        #  retry unresolved inverses
        Print("\nRETRY:\n");
        for inv in unresolved do
            x1:= inv[1];
            s:= inv[2];
            if not IsBound(mmm[s][x1]) then
                if inv[5] = 2 then
                    inverse:= inverse2Vec(x1, s, inv[3], inv[4]);
                else
                    inverse:= inverseVec(x1, s, inv[3], inv[4]);
                fi;
                if IsBool(inverse) then
                    Print("inversion ", x1, ".", s, " failed again\n");
                else
                    mmm[s][x1]:= inverse;
                    cacheVec(mmm[s][x1]);
                    Print("inversion: ", x1, ".", s, "\n");
                fi;
            fi;
        od;
    until Sum(simples, s-> Number(iii, n-> not IsBound(mmm[s][n]))) = open;

    # show me the entries cache
    hecke.uuu:= uuu;

    return hecke;
end;


# how to test whether the matrix m = mmm[l] satisfies m^2 = (q-1)*m + q*1
testMat:= function(mat)
    local   q,  hecke,  zero,  one,  test,  n,  img,  sum;

    q:= X(Rationals);
    hecke:= mat[1].hecke;
    zero:= Zero(hecke.baseRing);
    one:= One(hecke.baseRing);
    test:= true;
    for n in [1..Length(mat)] do
        Print(".\c");
        if IsBound(mat[n]) then
            img:= HVecUnderMat(mat[n], mat);
            sum:= (q-1) * mat[n] + HVec([n], [q*one], zero, hecke);
            if img <> sum then
                test:= false;
                Print(n);
            fi;
        fi;
    od;
    Print("\n");
    return test;
end;

# how to test whether the matrices mmm satisfy the relation rel
testRel:= function(rel, mmm)
    local   hecke,  zero,  one,  test,  n,  lhs,  rhs;

    hecke:= mmm[1][1].hecke;
    zero:= Zero(hecke.baseRing);
    one:= One(hecke.baseRing);
    test:= true;
    for n in [1..Length(mmm[1])] do
        Print(".\c");
        lhs:= HVecUnderWord(HVec([n], [one], zero, hecke), rel[1], mmm);
        rhs:= HVecUnderWord(HVec([n], [one], zero, hecke), rel[2], mmm);
        if lhs <> rhs then
            test:= false;
            Print(n);
        fi;
    od;
    Print("\n");
    return test;
end;

# how to test whether the matrices mmm define a representation.
testAll:= function(mmm, relations)
    local   test,  l,  rel;

    test:= true;
    for l in [1..Length(mmm)] do
        Print("testing matrix ", l, "...\n");
        test:= test and testMat(mmm[l]);
    od;
    for rel in relations do
        Print("testing ", rel, "...\n");
        test:= test and testRel(rel, mmm);
    od;
    return test;
end;
