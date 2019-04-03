CC          = c++ 

#-----------------------------------------
#Optimization ----------------------------
OPT   =  -Wno-deprecated -g

# on mac use this
#GL_LIB = -framework OpenGL -framework GLUT -framework foundation
# on linux use this
GL_LIB = -lGL -lglut -lGLU

#-----------------------------------------

TARGETS = inverseKinematics

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
#-----------------------------------------
#-----------------------------------------

inverseKinematics: $(OBJECTS) Joint.o inverseKinematics.o
	$(CC) $(OBJECTS) $(LDOPTS) inverseKinematics.o -o inverseKinematics

#-----------------------------------------
#-----------------------------------------

.cpp.o: 
	$(CC) $(CCOPTS) -c $< 

#-----------------------------------------
#-----------------------------------------












