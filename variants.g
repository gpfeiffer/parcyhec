#############################################################################
##
#A  variants.g                                    goetz.pfeiffer@nuigalway.ie
##
##  how to compute variants of relations
##

# invert a word
InverseWord:= function(word)
    return Reversed(-word);
end;

# turn a relation lft = rgt into a relator lft/rgt (=1)
RelatorRelation:= function(lft, rgt)
    return Concatenation(lft, InverseWord(rgt));
end;

# turn a relator into a list of all its variants s = word.
VariantsRelator:= function(relator)
    local   variants,  i,  s,  u,  v;

    variants:= rec();
    for i in [1..Length(relator)] do
        s:= relator[i];                       #  rel = u s v = 1,
        u:= relator{[1..i-1]};                #  s = (v u)^{-1},
        v:= relator{[i+1..Length(relator)]};  #  s^-1 = v u
        if s < 0 then
            if not IsBound(variants.(-s)) then
                variants.(-s):= [];
            fi;
            Add(variants.(-s), Concatenation(v, u));
        else
            if not IsBound(variants.(s)) then
                variants.(s):= [];
            fi;
            Add(variants.(s), InverseWord(Concatenation(v, u)));
        fi;
    od;
    return variants;
end;

# how to merge a list of variants records
# i.e., a list of records all of whose values are lists.
# the resulting record has keys pointing to list which are
# the unions of all the lists found under that key in the various records.
MergeVariants:= function(list)
    local   merge,  record,  key;

    merge:= rec();
    for record in list do
        for key in RecFields(record) do
            if IsBound(merge.(key)) then
                UniteSet(merge.(key), record.(key));
            else
                merge.(key):= Set(ShallowCopy(record.(key)));
            fi;
        od;
    od;
    for key in RecFields(merge) do
        Sort(merge.(key), function(a, b) return Length(a) < Length(b); end);
    od;
    return merge;
end;


VariantsRelations:= function(relations)
    local   relators;
    relators:= List(relations, x-> RelatorRelation(x[1], x[2]));
    return MergeVariants(List(relators, VariantsRelator));
end;
