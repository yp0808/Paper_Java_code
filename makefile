CPLEXDIR  = /apps/cplex121

JC       = javac
J        = java -Djava.library.path
P    	 = $(CPLEXDIR)/bin/x86-64_debian4.0_4.1 -Xmx4000m -classpath $(CPLEXDIR)/lib/cplex.jar:.:.. SALPMPE/ApproxAlgorithm
.SUFFIXES: .java .class

# Source files to compile
CLASSES =  Parameters.java\
	   state.java\
	   strategy.java\
	   sol_sep.java\
	   Fdist.java\
	   Profit99.java\
	   Profit10.java\
	   Profit.java\
	   ApproxLP.java\
	   ApproxAlgorithm.java

.java.class:
	$(JC) -classpath $(CPLEXDIR)/lib/cplex.jar:.:.. $*.java


# Compile source
default: $(CLASSES:.java=.class)


# Run algorithm with wrapper
wrap:
	sge_run --grid_submit=batch --grid_mem=10G --grid_priority=normal "$(J)=$(P)"

# Run algorithm
run:
	$(J)=$(P);

# Delete compiled files
clean: 
	$(RM) *.class