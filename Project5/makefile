CC          = c++ 

#-----------------------------------------
#Optimization ----------------------------
OPT   =  -Wno-deprecated -g

#GL_LIB = -framework OpenGL -framework GLUT -framework foundation
GL_LIB = -lGL -lglut -lGLU

#-----------------------------------------

TARGETS = springmass_skel

OBJECTS = ../common/jsoncpp.o ../common/slVector.o

#-----------------------------------------

LIBS = $(GL_LIB) 

CCOPTS = $(OPT) $(INCS) 
LDOPTS = $(OPT) $(INCS) $(LIBS) 

#-----------------------------------------
#-----------------------------------------

default: $(TARGETS)


clean:
	/bin/rm -f *.o $(TARGETS)
	rm -f *~ *part *#
#-----------------------------------------
#-----------------------------------------

springmass_skel: $(OBJECTS) advect.o 
	$(CC) $(OBJECTS) $(LDOPTS) advect.o -o advect

#-----------------------------------------
#-----------------------------------------

%.o: %.cpp
	$(CC) -c -DNDEBUG $(OPT) ${INCLUDEFLAGS} -o $@ $< 

#.cpp.o: 
#	$(CC) $(CCOPTS) -c $< 

#-----------------------------------------
#-----------------------------------------















