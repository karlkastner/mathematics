
LIBPATH = -L../jdlib
LIBS = -ljdqz -llapack -lblas
SUBDIRS = ../jdlib
OBJ = example.o

.PHONY: all
all: 	$(OBJ)
	$(FC) $(LDFLAGS) -o example $(OBJ) $(LIBPATH) $(LIBS)

.PHONY: clean
clean:
	-rm -f example
	-rm -f `find . -name '*.[ao]' -print`

