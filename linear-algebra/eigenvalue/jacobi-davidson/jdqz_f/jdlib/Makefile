#
#     Coded by Diederik R. Fokkema
#
#     $Id$
# 
#

LIB = libjdqz.a
OBJS = $(patsubst %.f,%.o,$(wildcard *.f))

$(LIB): $(LIB)($(OBJS)) 
	ranlib $@

$(LIB)(jdqz.o): jdqz.f
	$(FC) $(FFLAGS) -c jdqz.f
	$(AR) rv $(LIB) jdqz.o
	-rm jdqz.o

