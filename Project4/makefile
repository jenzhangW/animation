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
	rm -f *.obj
#-----------------------------------------
#-----------------------------------------

springmass_skel: $(OBJECTS) springmass_skel.o 
	$(CC) $(OBJECTS) $(LDOPTS) springmass_skel.o -o springmass_skel

#-----------------------------------------
#-----------------------------------------

%.o: %.cpp
	$(CC) -c -DNDEBUG $(OPT) ${INCLUDEFLAGS} -o $@ $< 

#.cpp.o: 
#	$(CC) $(CCOPTS) -c $< 

#-----------------------------------------
#-----------------------------------------















