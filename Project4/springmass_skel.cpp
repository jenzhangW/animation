/*
Name: Jenny Zhang
Date: 11/6/17
Email: jennyz1@umbc.edu

Description:  Creates an .obj of a cube like objects falling. Computer the movement through the length of springs and the force acting on the cube.

*/

#include <iostream>
#include <fstream>
#include <vector>
#include "../common/slVector.H"
#include "../common/json/json.h"
using namespace std;
static const double FEPS = 1e-6;

struct Mass {
	SlVector3 pos, vel, frc;
	double vol;
};

struct Spring {
  double l;
  int i, j;
};

class SlTri {
public:
  int indices[3];
  inline int &operator[](const unsigned int &i) { return indices[i];};
  inline int operator[](const unsigned int &i) const { return indices[i];};
  inline void set(int x, int y, int z) {indices[0] = x; indices[1] = y; indices[2] = z;};
  inline SlTri(int x, int y, int z) {indices[0] = x; indices[1] = y; indices[2] = z;};
  inline SlTri() {};
  inline SlTri &operator=(const SlTri &that);
};

inline SlTri &SlTri::operator=(const SlTri &that) {
  this->indices[0] = that.indices[0];
  this->indices[1] = that.indices[1];
  this->indices[2] = that.indices[2];
  return (*this);
};

struct SimulationParameters {
  double dt, total_time, k, d, density;
};

struct Object {
  std::vector<Mass> masses;
  std::vector<Spring> springs;
  std::vector<SlTri> triangles;
};

bool readObject(const char *fname, Object &object) {
  char ch;
  Mass m;
  Spring s;
  SlTri t;
  
  std::ifstream in(fname, std::ios::in);
  while (in>>ch) {
    if (ch == 'm') {
	  in>>m.pos[0]>>m.pos[1]>>m.pos[2]>>m.vol;
	  object.masses.push_back(m);
      continue;
    }
    if (ch == 's') {
	  in>>s.i>>s.j>>s.l;
	  object.springs.push_back(s);
      continue;
    }
    if (ch == 't') {
	  in>>t[0]>>t[1]>>t[2];
	  object.triangles.push_back(t);
      continue;
    }
  }
  in.close();
  std::cout<<"inputfile " << fname <<" read"<<std::endl;
  return true;
}

bool readInputFile(const char *fname, 
	SimulationParameters &params,
	std::vector<Object> &objects) {

  std::ifstream in(fname, std::ios::in);

  Json::Value root;
  Json::Reader jReader;

  if(!jReader.parse(in, root)){
	std::cout << "couldn't read input file: " << fname << '\n'
			  << jReader.getFormattedErrorMessages() << std::endl;
	exit(1);
  }

  params.dt = root.get("dt", 1.0/300.0).asDouble();
  params.total_time = root.get("total_time", 1.0).asDouble();
  params.k = root.get("stiffness", 1.0/300.0).asDouble();
  params.d = root.get("damping", 1.0/300.0).asDouble();
  params.density = root.get("density", 1.0/300.0).asDouble();

  Json::Value objectsIn = root["objects"];
  objects.resize(objectsIn.size());
  for (unsigned int i=0; i<objectsIn.size(); i++) {
	readObject((objectsIn[i]["filename"]).asString().c_str(), objects[i]);
  }
  return true;
}

void writeObj(char *fname, const std::vector<Mass> &meshPts, const std::vector<SlTri> &triangles) {
  std::cout<<"writing "<<fname<<std::endl;
  std::ofstream out;
  std::vector<Mass>::const_iterator p;
  std::vector<SlTri>::const_iterator t;

  out.open(fname);
  
  for (p=meshPts.begin(); p!=meshPts.end(); p++) 
	out<<"v "<<p->pos[0]<<" "<<p->pos[1]<<" "<<p->pos[2]<<std::endl;
  
  for (t=triangles.begin(); t!=triangles.end(); t++) 
	out<<"f "<<(*t)[0]+1<<" "<<(*t)[1]+1<<" "<<(*t)[2]+1<<std::endl;
  
  out.close();
}

int main(int argc, char *argv[]) {
  char fname[80];
  SimulationParameters params;
  std::vector<Object> objects;
  readInputFile(argv[1], params, objects);

  double time = 0;
  int frame = 0;
  int inbetween = 0;
  vector<double> mass;
  SlVector3 empty(0,0,0);
  //force accumulator
  for (unsigned int i=0; i<objects.size(); i++) {
    for (int f = 0; f < objects[i].masses.size(); f++) {
      double m = objects[i].masses[f].vol * params.density;
      mass.push_back(m);	  
      objects[i].masses[f].vel = empty;
    }
  }
  while (time < params.total_time) {
    for (unsigned int i=0; i<objects.size(); i++) {
      sprintf(fname, argv[2], i, frame);
      //calculalte the masses, volume * density
      for (int a = 0; a < objects[i].masses.size(); a++) {
	objects[i].masses[a].frc = empty;
	//set force of gravity
	objects[i].masses[a].frc[2] = -9.81 * mass[a];
      }
      for (int b = 0; b < objects[i].springs.size(); b++) {
	int m_i = objects[i].springs[b].i;
	int m_j = objects[i].springs[b].j;

	//position of node 1 - postition of node 2
	SlVector3 pj_pi = objects[i].masses[m_j].pos - objects[i].masses[m_i].pos; 
	SlVector3 pi_pj = objects[i].masses[m_i].pos - objects[i].masses[m_j].pos; 
	//force on i from j
	SlVector3 force_ij = params.k 
	  * ((mag(pj_pi)/objects[i].springs[b].l) - 1) 
	  * (pj_pi/mag(pj_pi)) 
	  + (params.d * (objects[i].masses[m_j].vel - objects[i].masses[m_i].vel));
	
	//force on j from i
	SlVector3 force_ji = params.k 
	  * ((mag(pi_pj)/objects[i].springs[b].l) - 1) 
	  * (pi_pj/mag(pi_pj)) 
	  + (params.d * (objects[i].masses[m_i].vel - objects[i].masses[m_j].vel));

	objects[i].masses[m_i].frc += force_ij;
	objects[i].masses[m_j].frc += force_ji;
      }
      for (int c = 0; c < objects[i].masses.size(); c++) {
	//calculate accleration
	SlVector3 acceleration = objects[i].masses[c].frc / mass[c];
	//calculate new velocity
	objects[i].masses[c].vel = objects[i].masses[c].vel + (acceleration * params.dt);
	//calculate new position
	objects[i].masses[c].pos = objects[i].masses[c].pos + ( objects[i].masses[c].vel * params.dt);
	if (objects[i].masses[c].pos[2] < 0) {
	  objects[i].masses[c].pos[2] = 0;
	}
      }	  
      //prints the obj when a certain amount of time has passed
      if (inbetween % 10 == 0) {
	writeObj(fname, objects[i].masses, objects[i].triangles);
	frame++;
      }
    }
    // takes a timestep
    time += params.dt;
    inbetween++;
  }
}
