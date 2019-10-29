# Fri Aug  3 12:27:58 MSK 2012
# Karl Kästner, Berlin

DIR=class
JC=javac
JFLAGS=-d $(DIR) -Xlint:unchecked

SRC=$(wildcard *.java)
CLASS=$(addprefix $(DIR)/, $(addsuffix .class, $(basename $(SRC))))


.PHONY: all
all : $(CLASS)

$(DIR)/%.class : %.java
	$(JC) $(JFLAGS) $*.java

.PHONY: clean
clean :
	-rm -f $(DIR)/*.class

