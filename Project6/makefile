CC          = c++ 

#-----------------------------------------
#Optimization ----------------------------
OPT   =  -Wno-deprecated -g -O3

#GL_LIB = -framework OpenGL -framework GLUT -framework foundation
GL_LIB = -lGL -lglut -lGLU

#-----------------------------------------

TARGETS = advect fluid

OBJECTS = 

#-----------------------------------------

INCS = -I../common
LIBS = $(GL_LIB) -L../common -lcommon 

CCOPTS = $(OPT) $(INCS) 
LDOPTS = $(OPT) $(INCS) $(LIBS) 

#-----------------------------------------
#-----------------------------------------

default: $(TARGETS)


clean:
	/bin/rm -f *.o $(TARGETS)

#-----------------------------------------
#-----------------------------------------

advect: $(OBJECTS) advect.o
	$(CC) $(OBJECTS) $(LDOPTS) advect.o -o advect 

fluid: $(OBJECTS) fluid.o
	$(CC) $(OBJECTS) $(LDOPTS) fluid.o -o fluid 

sample: $(OBJECTS) sample.o
	$(CC) $(OBJECTS) $(LDOPTS) sample.o -o sample

viewer: $(OBJECTS) viewer.o
	$(CC) $(OBJECTS) $(LDOPTS) viewer.o -o viewer

#-----------------------------------------
#-----------------------------------------

.cpp.o: 
	$(CC) $(CCOPTS) -c $< 

#-----------------------------------------
#-----------------------------------------















