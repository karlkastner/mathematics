# export all variables to recursive make calls
export

FC = gfortran
#FC = f77
#FC = ifort
FFLAGS = -u -O3
SHELL=/bin/sh

LIBPATH = -L../jdlib
LIBS = -ljdqr -llapack -lblas
SUBDIRS = jdlib
TEST    = jdtest


.PHONY: all
all: 	libs jdqr doc

.PHONY: jdqr
jdqr: 	
	@for i in $(TEST); do (cd $$i; $(MAKE)); done

.PHONY: libs
libs:
	@for i in $(SUBDIRS); do (cd $$i; $(MAKE)); done

.PHONY: doc
doc : manual.pdf

%.pdf : %.tex
	pdflatex $*.tex
	-rm manual.aux
	-rm manual.log

.PHONY: clean
clean:
	-rm -f jdqz
	-rm -f `find . -name '*.[ao]' -print`
	-rm -f manual.pdf
	@for i in $(TEST); do make -C $(TEST) clean; done

