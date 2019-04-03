CC          = c++ 

#-----------------------------------------
#Optimization ----------------------------
OPT   =  -Wno-deprecated -g

# on mac use this
#GL_LIB = -framework OpenGL -framework GLUT -framework foundation
# on linux use this
GL_LIB = -lGL -lglut -lGLU

#-----------------------------------------

TARGETS = forwardKinematics viewer

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

forwardKinematics: $(OBJECTS) Joint.o forwardKinematics.o
	$(CC) $(OBJECTS) $(LDOPTS) forwardKinematics.o -o forwardKinematics

viewer: $(OBJECTS) viewer.o
	$(CC) $(OBJECTS) $(LDOPTS) viewer.o -o viewer

#-----------------------------------------
#-----------------------------------------

.cpp.o: 
	$(CC) $(CCOPTS) -c $< 

#-----------------------------------------
#-----------------------------------------












