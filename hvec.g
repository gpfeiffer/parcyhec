#############################################################################
##
##  hvec.g                                        goetz.pfeiffer@nuigalway.ie
##
##  Hecke algebra elements as sparse vectors over a parabolic subalgebra.
##

#############################################################################
HVecOps:= OperationsRecord("HVecOps");

c:= ComplexReflectionGroup;

# type check
IsHVec:= function(obj)
    return IsRec(obj) and IsBound(obj.isHVec) and obj.isHVec = true;
end;


# constructor
HVec:= function(poss, vals, zero, hecke)
    return rec(poss:= poss, vals:= vals, zero:= zero, hecke:= hecke,
               isHVec:= true, operations:= HVecOps);
end;


# how to print a HVec
HVecOps.Print:= function(self)
    Print("HVec(", self.poss, ", ", self.vals, ", ", self.zero, ", ", self.hecke, ")");
end;


# get the n-th coefficient
Get:= function(obj, n)
    return obj.operations.Get(obj, n);
end;


# how to get the n-th coefficient
HVecOps.Get:= function(self, n)
    local   pos;
    pos:= Position(self.poss, n);
    if pos = false then
        return self.zero;
    else
        return self.vals[pos];
    fi;
end;


# sum
HVecOps.\+:= function(l, r)
    local   vec,  vals,  poss,  i,  pos,  val;

    vec:= HVec(Concatenation(l.poss, r.poss), Concatenation(l.vals, r.vals),
              l.zero, l.hecke);

    # normalize vec in place (sort poss, ignore zero vals)
    l:= Length(vec.poss);
    SortParallel(vec.poss, vec.vals);

    poss:= [];
    vals:= [];

    i:= 1;
    while i <= l do
        pos:= vec.poss[i];
        val:= vec.zero;
        while i <= l and vec.poss[i] = pos do
            val:= val + vec.vals[i];
            i:= i+1;
        od;
        if val <> vec.zero then
            Add(poss, pos);
            Add(vals, val);
        fi;
    od;

    IsSet(poss);
    vec.poss:= poss;
    vec.vals:= vals;
    return vec;
end;

# scalar multiplication
HVecScaled:= function(alpha, vec)
    if alpha = 0 * alpha then
        return Zero(vec.hecke);
    fi;
    return HVec(Copy(vec.poss), List(vec.vals, x-> alpha * x),
                vec.zero, vec.hecke);
end;

# how to subtract
HVecOps.\-:= function(l, r)
#    return l + (-1)*r;
    return l + HVecScaled(-1, r);
end;


# how to compute the image of a vec under a matrix
HVecUnderMat:= function(vec, mat)
    local   img,  k,  pos,  val;
    img:= Zero(vec.hecke);
    for k in [1..Length(vec.poss)] do
        pos:= vec.poss[k];
        val:= vec.vals[k];
        if not IsBound(mat[pos]) then
#            Print("missing: ", l, ".", i, "\n");
            return false;
        else
#            img:= img + val * mat[pos];
            img:= img + HVecScaled(val, mat[pos]);
        fi;
    od;
    return img;
end;


# how to turn a word (in generators 1..4 and their q-inverses -1..-4,
# ie. -s stands for q s') into a vec
HVecUnderWord:= function(vec, word, mmm)
    local   q,  a,  img;
    q:= X(Rationals);
    for a in word do
        img:= HVecUnderMat(vec, mmm[AbsInt(a)]);
        if IsBool(img) then
            return false;
        fi;
        if a < 0 then
#            img:= img + (1-q) * vec;
            img:= img + HVecScaled(1-q, vec);
        fi;
        vec:= img;
    od;
    return vec;
end;


##  Multiplication.

ProdHVecH:= function(l, r)  ## called if IsHVec(l)
    local   sum,  W,  i,  w,  prod;

    sum:= Zero(l.hecke);
    if IsHVec(r) then
        for i in [1..Length(r.poss)] do
            prod:= ProdHVecH(l, r.vals[i]);  # recurse
            w:= r.hecke.words[r.poss[i]];
            sum:= sum + HVecUnderWord(prod, w, l.hecke.mats);
        od;
    elif IsHeckeElt(r) then
        W:= r.hecke.reflectionGroup;
        for i in [1..Length(r.elm)] do
            prod:= r.coeff[i] * l;
            w:= CoxeterWord(W, r.elm[i]);
            sum:= sum + HVecUnderWord(prod, l.hecke.subset{w}, l.hecke.mats);
        od;
    else # scalar: commute!
        return r * l;
    fi;
    return sum;
end;


HVecOps.\*:= function(l, r)

    # presumably we're called because r is a HVec
    if IsHVec(r) then
        if IsHVec(l) and l.hecke = r.hecke then
            return ProdHVecH(l, r);
        elif l = 0*l then
            return r.zero;
        else # scalar
            return HVec(r.poss, List(r.vals, x-> l * x), r.zero, r.hecke);
        fi;
    else
        Error("don't know how to <l> * <r>");
    fi;
end;


##  Inverses.
FindInverseHeckeElt:= function(h)
    local   W,  T,  weight,  list,  j,  word,  qlen, a,  len,  inv;

    # trivial case first
    if Length(h.coeff) = 0 then
        return false;
    fi;

    # easy case next
    if Length(h.coeff) = 1 then
        if h.coeff[1]^-1 <> false then
            return h^-1;
        else
            return false;
        fi;
    fi;

    # otherwise recurse
    W:= Group(h.hecke);
    T:= Basis(h.hecke, h.basis);

    # how to assign a weight to h
    weight:= function(h)
        return Maximum(List(h.elm, x-> CoxeterLength(W, x)));
    end;

    # find a coeff that does not specialize to 0.
    list:= List(h.coeff, x-> Value(x, 1));
    if Number(list, x-> x <> 0) <> 1 then
        return false;
    fi;

    j:= PositionProperty(list, x -> x <> 0);

    word:= CoxeterWord(W, h.elm[j]);
    if word = [] then
        return false;
    fi;
    a:= T(word[1]);

    len:= weight(h);
    if weight(a * h) < len then
        inv:= FindInverseHeckeElt(a * h);
        if inv = false then
            return false;
        else
            return inv * a;
        fi;
    fi;
    if weight(a^-1 * h) < len then
        inv:= FindInverseHeckeElt(a^-1 * h);
        if inv = false then
            return false;
        else
            return inv / a;
        fi;
    fi;

    # if we get here we're out of luck.
    return false;
end;

FindInverse:= function(obj)
    if IsHeckeElt(obj) then
        return FindInverseHeckeElt(obj);
    else
        return obj.operations.FindInverse(obj);
    fi;
end;

HVecOps.FindInverse:= function(self)
    local   inv,  word,  qlen;
    if Length(self.poss) = 1 then
        inv:= FindInverse(self.vals[1]);
        if inv <> false then
            word:= self.hecke.words[self.poss[1]];
            qlen:= X(Rationals)^(-Length(word))*One(self.hecke.baseRing);
            return qlen * HVecUnderWord(One(self.hecke), -Reversed(word), self.hecke.mats)
                   * (inv * One(self.hecke));
        else
            return false;
        fi;
    else
        return false;
    fi;
end;


#############################################################################
##
##  HeckeMats: a Hecke algebra represented as matrices over a parabolic
##

# operations
HeckeMatsOps:= OperationsRecord("HeckeMatsOps");


# type check
IsHeckeMats:= function(obj)
    return IsRec(obj) and IsBound(obj.isHeckeMats) and obj.isHeckeMats = true;
end;


# constructor
HeckeMats:= function(baseRing, subset, mats, words)
    local   new;
    new:= rec(mats:= mats, words:= words, subset:= subset, baseRing:= baseRing);
    new.isHeckeMats:= true;
    new.operations:= HeckeMatsOps;
    return new;
end;


# print
HeckeMatsOps.Print:= function(self)
    if IsBound(self.name) then
        Print(self.name);
    else
        Print("<hecke>");
    fi;
end;


# Equality
HeckeMatsOps.\= := function(l,r)
    return IsIdentical(l, r);
end;


# Zero and One
HeckeAlgebraOps.One:= function(self)
    return Basis(self, "T")();
end;

HeckeAlgebraOps.Zero:= function(self)
    local   one;
    one:= One(self);
    return one - one;
end;

HeckeMatsOps.One:= function(self)
    return HVec([1], [One(self.baseRing)], Zero(self.baseRing), self);
end;

HeckeMatsOps.Zero:= function(self)
    return HVec([], [], Zero(self.baseRing), self);
end;


# T
HeckeMatsOps.Basis:= function(self, name)
    if name = "T" then
        return function(arg)
            return HVecUnderWord(One(self), arg, self.mats);
        end;
    else
        Error("don't know how to Basis( <hecke>, <name> )");
    fi;
end;
