GAP=../../gap-64/gap3-jm5/bin/gap.sh64
#GAP=../../gap-64/gap3-jm5/bin/gap.sh
MEM=-m 4000m
LIB=hvec.g variants.g parcyhec.g

LOGS = \
g12log.txt\
g22log.txt\
g24log.txt\
g27log.txt\
g29log.txt\
g33log.txt\
g33d4log.txt\
g31g423log.txt\
#g31log.txt\

all: ${LOGS}

g12log.txt: g12input.g
	${GAP} ${MEM} ${LIB} < $< | tee $@

g22log.txt: g22input.g
	${GAP} ${MEM} ${LIB} < $< | tee $@

g24log.txt: g24input.g
	${GAP} ${MEM} ${LIB} < $< | tee $@

g27log.txt: g27input.g
	${GAP} ${MEM} ${LIB} < $< | tee $@

g29log.txt: g29input.g
	${GAP} ${MEM} ${LIB} < $< | tee $@

g31log.txt: g31input.g
	${GAP} ${MEM} ${LIB} < $< | tee $@

g31g423log.txt: g31g423input.g
	${GAP} ${MEM} ${LIB} < $< | tee $@

g33log.txt: g33input.g
	${GAP} ${MEM} ${LIB} < $< | tee $@

g33d4log.txt: g33d4input.g
	${GAP} ${MEM} ${LIB} < $< | tee $@

g34g33log.txt: g34g33input.g
	${GAP} ${MEM} ${LIB} < $< | tee $@


clean:
	rm -fv *~

distclean: clean
	rm -fv *.txt
