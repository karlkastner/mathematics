# Fri Aug  3 12:27:58 MSK 2012
# Karl Kästner, Berlin

#export CLASSPATH='/usr/share/java/*':./class
CLASSPATH='/usr/share/java/*':./class:
#/home/pia/phd/src/lib/class
SOURCEPATH=./

DIR=class
JC=javac
#JFLAGS=-d $(DIR) -Xlint:all -classpath $(CLASSPATH) -sourcepath $(SOURCEPATH) -target 1.6 -source 1.6
JFLAGS=-g -d $(DIR) -Xlint:unchecked -classpath $(CLASSPATH) -sourcepath $(SOURCEPATH) -target 1.6 -source 1.6

SRC=$(wildcard *.java)
CLASS=$(addprefix $(DIR)/, $(addsuffix .class, $(basename $(SRC))))


.PHONY: all
all : $(CLASS)

$(DIR)/%.class : %.java
	$(JC) $(JFLAGS) $*.java

.PHONY: clean
clean :
	-rm -f $(DIR)/*.class

