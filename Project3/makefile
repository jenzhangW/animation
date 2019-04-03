CC          = c++ 

#-----------------------------------------
#Optimization ----------------------------
OPT   =  -Wno-deprecated -g

# on mac use this
#GL_LIB = -framework OpenGL -framework GLUT -framework foundation
# on linux use this
GL_LIB = -lGL -lglut -lGLU

#-----------------------------------------

TARGETS = skinning

OBJECTS = slVector.o slMatrix.o slRotations.o Joint.o

#-----------------------------------------

LIBS = $(GL_LIB) 

CCOPTS = $(OPT) $(INCS) 
LDOPTS = $(OPT) $(INCS) $(LIBS) 

#-----------------------------------------
#-----------------------------------------

default: $(TARGETS)


clean:
	/bin/rm -f *.o $(TARGETS)
	rm -f *~
	rm -f *obj
#-----------------------------------------
#-----------------------------------------

forwardKinematics: $(OBJECTS) Joint.o skinning.o
	$(CC) $(OBJECTS) $(LDOPTS) skinning.o -o skinning

#-----------------------------------------
#-----------------------------------------

.cpp.o: 
	$(CC) $(CCOPTS) -c $< 

#-----------------------------------------
#-----------------------------------------












