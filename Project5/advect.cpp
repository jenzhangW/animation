/*
Name: Jenny Zhang
Date: 11/21/17
Email: jennyz1@umbc.edu

Description: Simulates fluid motion through the use of grid based fluid simulation. Fluid is consisting on particles that move with a certain velocity. 
The particles are put through a MAC grid where new velocities are computed and then positions are updated.
 */
#include <fstream>
#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <stdio.h>
#include <string.h>
#include "../common/slVector.H"
#include "../common/slArray.H"
#include "../common/json/json.h"
using namespace std;

struct Particle {
  SlVector3 pos;
  SlVector3 vel;
};

struct SimulationParameters {
  double dt, total_time, density, h, flip;
  SlVector3 res, lc;
};
bool readParticle(const char* fname, vector<Particle> &particles, unsigned int &num) {
  char ch;
  Particle p;
  std::ifstream in(fname, std::ios::in);
  in>>num;
  
  for (int i = 0; i < num; i++) {
  in>>p.pos[0]>>p.pos[1]>>p.pos[2]>>p.vel[0]>>p.vel[1]>>p.vel[2];
    particles.push_back(p);
    continue;
  }
  in.close();
  return true;
}

bool writeParticles(const char *fname, unsigned int &nFlipParticles, std::vector<SlVector3> &flipPos, std::vector<SlVector3> &flipVel) {
  std::ofstream out(fname, std::ios::out);
  out<<nFlipParticles<<std::endl;
  for (unsigned int i=0; i<nFlipParticles; i++) {
    out<<flipPos[i][0]<<" "<<flipPos[i][1]<<" "<<flipPos[i][2]<<" "<<flipVel[i][0]<<" "<<flipVel[i][1]<<" "<<flipVel[i][2]<<std::endl;
  }
  out.close();
  std::cout<<"outputfile " << fname <<" written"<<std::endl;
  return true;
}

bool readInputFile(const char *fname,
		   SimulationParameters &params,
		   vector<Particle> &particles,
		   unsigned int &num) {

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
  params.density = root.get("density", 1.0/300.0).asDouble();
  SlVector3 r = (root.get("res", 1.0)[0].asInt(),root.get("res", 1.0)[1].asInt(),root.get("res", 1.0)[2].asInt());
  params.res = r;
  params.h = root.get("h", 1.0/300.0).asDouble();
  SlVector3 l = (root.get("lc", 1.0)[0].asDouble(),root.get("lc", 1.0)[1].asDouble(),root.get("lc", 1.0)[2].asDouble());
  params.lc = l;
  params.flip = root.get("flipRatio", 1.0).asDouble();
 
  //for opening particle file
  string name = root["particles"].asString();
  char file[15];
  strcpy(file, name.c_str());
  readParticle(file, particles,num);
  return true;
}

/*  //front bottom left
    SlVector3 fbl(x0,y0,z0);
    //front bottom right;
    SlVector3 fbr(x1,y0,z0);
    //front top left;
    SlVector3 ftl(x0,y1,z0);
    //front top right;
    SlVector3 ftr(x1,y1,z0);
    //back bottom left
    SlVector3 bbl(x0,y0,z1);
    //back bottom right
    SlVector3 bbr(x1,y0,z1);
    //back top left
    SlVector3 btl(x0,y1,z1);
    //back top right
    SlVector3 btr(x1,y1,z1);
*/

double interpolateX(vector<Particle> &particles, 
		  SimulationParameters &params,
		  SlArray3D<double> &u_v, 
		  int current) {
  double ucx = params.lc[0] + (params.h * params.res[0]);
  double ucy = params.lc[1] + (params.h * params.res[1]);
  double ucz = params.lc[2] + (params.h * params.res[1]);
  SlVector3 upperCorner(ucx,ucy,ucz);

  //find which cell it's in
  SlVector3 cell = (particles[current].pos - params.lc) / params.h;

  //find the weights
  double w_x, w_y, w_z;
  w_x = fmod(cell.data[0],1);
  w_y = fmod(cell.data[1],1);
  w_z = fmod(cell.data[2],1);

  //get the whole number which cell it's in
  int temp = 0;
  for (int a = 0; a < 3; a++) {
    temp = cell.data[a];
    cell.data[a] = temp;
  }

  //finding offset block
  int x0 = cell.data[0];
  if (x0 < 0) {
    x0 = 0;
  }
  int x1 = x0 + 1;
  if (x1 > 25) {
    x1 = 25;
  }

  int y0 = cell.data[1] - (1/2);
  if (y0 < 0) {
    y0 = 0;
  }
  int y1 = y0 + 1;
  int z0 = cell.data[2] - (1/2);
  int z1 = z0 + 1;
  //actual interpolation
  double interX = (1 - w_x)*(1 - w_y)*(1 - w_z)*(u_v(x0,y0,z0)) +
    (1 - w_x)*(1 - w_y)*(w_z)* (u_v(x0,y0,z1)) +
    (1 - w_x)*(w_y)*(1 - w_z)* (u_v(x0,y1,z0))+
    (1 - w_x)*(w_y)*(w_z)* (u_v(x0,y1,z1)) +
    (w_x)*(1 - w_y)*(1 - w_z)* (u_v(x1,y0,z0)) +
    (w_x)*(1 - w_y)*(w_z)* (u_v(x1,y0,z1)) +
    (w_x)*(w_y)*(1 - w_z)* (u_v(x1,y1,z0)) +
    (w_x)*(w_y)*(w_z)* (u_v(x1,y1,z1));
  return interX;
}

double interpolateY(vector<Particle> &particles, 
		  SimulationParameters &params,
		  SlArray3D<double> &v_v, 
		  int current) {
  SlVector3 cell = (particles[current].pos - params.lc) / params.h;
  double w_x, w_y, w_z;
  w_x = fmod(cell.data[0],1);
  w_y = fmod(cell.data[1],1);
  w_z = fmod(cell.data[2],1);

  int temp = 0;
  for (int a = 0; a < 3; a++) {
    temp = cell.data[a];
    cell.data[a] = temp;
  }

  //finding the offset cell
  int x0 = cell.data[0] - (1/2);
  if (x0 < 0) {
    x0 = 0;
  }
  int x1 = x0 + 1;
  if (x1 > 25) {
    x1 = 25;
  }
  int y0 = cell.data[1];
  if (y0 < 0) {
    y0 = 0;
  }

  int y1 = y0 + 1;
  int z0 = cell.data[2] - (1/2);
  int z1 = z0 + 1;

  //interpolating the y velocities
  double interY = (1 - w_x)*(1 - w_y)*(1 - w_z)*(v_v(x0,y0,z0)) +
    (1 - w_x)*(1 - w_y)*(w_z)* (v_v(x0,y0,z1)) +
    (1 - w_x)*(w_y)*(1 - w_z)* (v_v(x0,y1,z0))+
    (1 - w_x)*(w_y)*(w_z)* (v_v(x0,y1,z1)) +
    (w_x)*(1 - w_y)*(1 - w_z)* (v_v(x1,y0,z0)) +
    (w_x)*(1 - w_y)*(w_z)* (v_v(x1,y0,z1)) +
    (w_x)*(w_y)*(1 - w_z)* (v_v(x1,y1,z0)) +
    (w_x)*(w_y)*(w_z)* (v_v(x1,y1,z1));
  return interY;
}

double interpolateZ(vector<Particle> &particles, 
		  SimulationParameters &params,
		  SlArray3D<double> &w_v, 
		  int current) {
  SlVector3 cell = (particles[current].pos - params.lc) / params.h;
  double w_x, w_y, w_z;
  w_x = fmod(cell.data[0],1);
  w_y = fmod(cell.data[1],1);
  w_z = fmod(cell.data[2],1);

  int temp = 0;
  for (int a = 0; a < 3; a++) {
    temp = cell.data[a];
    cell.data[a] = temp;
  }

  int x0 = cell.data[0] - (1/2);
  if (x0 < 0) {
    x0 = 0;
  }
  int x1 = x0 + 1;
  if (x1 > 25) {
    x1 = 25;
  }
  int y0 = cell.data[1] - (1/2);
  if (y0 < 0) {
    y0 = 0;
  }
  int y1 = y0 + 1;
  int z0 = cell.data[2];
  int z1 = z0 + 1;

  double interZ = (1 - w_x)*(1 - w_y)*(1 - w_z)*(w_v(x0,y0,z0)) +
    (1 - w_x)*(1 - w_y)*(w_z)* (w_v(x0,y0,z1)) +
    (1 - w_x)*(w_y)*(1 - w_z)* (w_v(x0,y1,z0))+
    (1 - w_x)*(w_y)*(w_z)* (w_v(x0,y1,z1)) +
    (w_x)*(1 - w_y)*(1 - w_z)* (w_v(x1,y0,z0)) +
    (w_x)*(1 - w_y)*(w_z)* (w_v(x1,y0,z1)) +
    (w_x)*(w_y)*(1 - w_z)* (w_v(x1,y1,z0)) +
    (w_x)*(w_y)*(w_z)* (w_v(x1,y1,z1));
  return interZ;
}

//call the each interpolation for each component of the velocity
SlVector3 interpolateXYZ(vector<Particle> &particles, 
			 SimulationParameters &params,
			 SlArray3D<double> &u_v,
			 SlArray3D<double> &v_v,
			 SlArray3D<double> &w_v,
			 int current) {
  double u = interpolateX(particles, params, u_v, current);
  double v = interpolateY(particles, params, v_v, current);
  double w = interpolateZ(particles, params, w_v, current);
  SlVector3 interXYZ(u,v,w);
  return interXYZ;
}

//Take the velocities and set the faces of the grid to certain velocities
void particleToGrid(vector<Particle> &particles,
		    SimulationParameters &params,
		    SlArray3D<double> &u_v,
		    SlArray3D<double> &v_v,
		    SlArray3D<double> &w_v) {
  double ucx = params.lc[0] + (params.h * params.res[0]);
  double ucy = params.lc[1] + (params.h * params.res[1]);
  double ucz = params.lc[2] + (params.h * params.res[1]);
  SlVector3 upperCorner(ucx,ucy,ucz);

  for (int a = 0; a < particles.size(); a++) { 
    SlVector3 cell = (particles[a].pos - params.lc) / params.h;
    //find the weights
    double w_x, w_y, w_z;
    w_x = fmod(cell.data[0],1);
    w_y = fmod(cell.data[1],1);
    w_z = fmod(cell.data[2],1);

    //find the cell it is in
    int temp = 0;
    for (int b = 0; b < 3; b++) {
      temp = cell.data[b];
      cell.data[b] = temp;
    }
    //for u velocity 
    int x0 = cell.data[0];
    if (x0 < 0) {
      x0 = 0;
    }
    int x1 = x0 + 1;
    int y0 = cell.data[1] - (1/2);
    if (y0 < 0 ) {
      y0 = 0;
    }
    int y1 = y0 + 1;
    int z0 = cell.data[2] - (1/2);
    if (z0 < 0 ) {
      z0 = 0;
    }

    int z1 = z0 + 1;

    u_v(x0,y0,z0) = (1 - w_x)*(1 - w_y)*(1 - w_z)*particles[a].vel[0];
    u_v(x0,y0,z1) = (1 - w_x)*(1 - w_y)*(w_z)*particles[a].vel[0];
    u_v(x0,y1,z0) = (1 - w_x)*(w_y)*(1 - w_z)*particles[a].vel[0];
    u_v(x0,y1,z1) = (1 - w_x)*(w_y)*(w_z)*particles[a].vel[0];
    u_v(x1,y0,z0) = (w_x)*(1 - w_y)*(1 - w_z)*particles[a].vel[0];
    u_v(x1,y0,z1) = (w_x)*(1 - w_y)*(w_z)*particles[a].vel[0];
    u_v(x1,y1,z0) = (w_x)*(w_y)*(1 - w_z)*particles[a].vel[0];
    u_v(x1,y1,z1) = (w_x)*(w_y)*(w_z)*particles[a].vel[0];

    //for v velovity
    x0 = cell.data[0] - (1/2);
    if (x0 < 0) {
      x0 = 0;
    }
    x1 = x0 + 1;
    y0 = cell.data[1];
    if (y0 < 0) {
      y0 = 0;
    }
    y1 = y0 + 1;
    z0 = cell.data[2] - (1/2);
    if (z0 < 0) {
      z0 = 0;
    } 
    z1 = z0 + 1;

    v_v(x0,y0,z0) = (1 - w_x)*(1 - w_y)*(1 - w_z)*particles[a].vel[1];
    v_v(x0,y0,z1) = (1 - w_x)*(1 - w_y)*(w_z)*particles[a].vel[1];
    v_v(x0,y1,z0) = (1 - w_x)*(w_y)*(1 - w_z)*particles[a].vel[1];
    v_v(x0,y1,z1) = (1 - w_x)*(w_y)*(w_z)*particles[a].vel[1];
    v_v(x1,y0,z0) = (w_x)*(1 - w_y)*(1 - w_z)*particles[a].vel[1];
    v_v(x1,y0,z1) = (w_x)*(1 - w_y)*(w_z)*particles[a].vel[1];
    v_v(x1,y1,z0) = (w_x)*(w_y)*(1 - w_z)*particles[a].vel[1];
    v_v(x1,y1,z1) = (w_x)*(w_y)*(w_z)*particles[a].vel[1];

    //for w velovity
    x0 = cell.data[0] - (1/2);
    if (x0 < 0) {
      x0 = 0;
    }
    x1 = x0 + 1;
    y0 = cell.data[1] - (1/2);
    if (y0 < 0) {
      y0 = 0;
    }
    y1 = y0 + 1;
    if (y1 > 25) {
      y1 = 25;
    }
    
    z0 = cell.data[2];
    if (z0 < 0) {
      z0 = 0;
    }
    z1 = z0 + 1;
    if (z1 > 25) {
      z1 = 25;
    }

    w_v(x0,y0,z0) = (1 - w_x)*(1 - w_y)*(1 - w_z)*particles[a].vel[2];
    w_v(x0,y0,z1) = (1 - w_x)*(1 - w_y)*(w_z)*particles[a].vel[2];
    w_v(x0,y1,z0) = (1 - w_x)*(w_y)*(1 - w_z)*particles[a].vel[2];
    w_v(x0,y1,z1) = (1 - w_x)*(w_y)*(w_z)*particles[a].vel[2];
    w_v(x1,y0,z0) = (w_x)*(1 - w_y)*(1 - w_z)*particles[a].vel[2];
    w_v(x1,y0,z1) = (w_x)*(1 - w_y)*(w_z)*particles[a].vel[2];
    w_v(x1,y1,z0) = (w_x)*(w_y)*(1 - w_z)*particles[a].vel[2];
    w_v(x1,y1,z1) = (w_x)*(w_y)*(w_z)*particles[a].vel[2];
  }
}

//taking the velocities of the grid and finding the velocity for the particle. Then use that velovity to find the final position
/*psuedo-code
  o = interpolate(p.pos, u, v, w);
  n = iterpolate(p.pos, fu, fv, fw);
  p.vel = p.vel + 
  (flipRatio * (p.vel + (n - o))) + 
  ((1 - flipRation) * n); 
*/
void gridToParticle(vector<Particle> &particles,
		    SimulationParameters &params,
 		    SlArray3D<double> &u_v,
		    SlArray3D<double> &v_v,
		    SlArray3D<double> &w_v,
 		    SlArray3D<double> &fu_v,
		    SlArray3D<double> &fv_v,
		    SlArray3D<double> &fw_v){

  for (int i = 0; i < particles.size(); i++) {
    SlVector3 o = interpolateXYZ(particles, params, u_v, v_v, w_v, i);
    SlVector3 n = interpolateXYZ(particles, params, fu_v, fv_v, fw_v, i);
    particles[i].vel = particles[i].vel +
      (params.flip * (particles[i].vel + (n - o)))+
      ((1 - params.flip) * n);
  }
}

//just trilinear interpolation. Where the particle is moved through the grid to see how it interacts with the other particles
void advect(vector<Particle> &particles,
	    SimulationParameters &params,
	    SlArray3D<double> &u_v,
	    SlArray3D<double> &v_v,
	    SlArray3D<double> &w_v) {
  for( int i = 0; i < particles.size(); i++) {
    interpolateXYZ(particles, params, u_v, v_v, w_v, i);
    particles[i].pos = particles[i].pos + (particles[i].vel * params.dt);
  }
}

int main(int argc, char *argv[]) {
  char fname[80];
  SimulationParameters params;
  vector<Particle> particles;
  unsigned int numPart;
  readInputFile(argv[1], params, particles, numPart);

  //grid set up
  double ucx = params.lc[0] + (params.h * params.res[0]);
  double ucy = params.lc[1] + (params.h * params.res[1]);
  double ucz = params.lc[2] + (params.h * params.res[1]);
  SlVector3 upperCorner(ucx,ucy,ucz);

  double time = 0;
  int frame = 0;
  int inbetween = 0;
  SlArray3D<double> u(params.res[0]+1,params.res[1],params.res[2]);
  SlArray3D<double> v(params.res[0],params.res[1]+1,params.res[2]);
  SlArray3D<double> w(params.res[0],params.res[1],params.res[2]+1);
  for (int q = 0; q < u.nx(); q++) {
    for (int r = 0; r < u.ny(); r++) {
      for (int e = 0; e < u.nz(); e++) {
	u(q,r,e) = 0;
      }
    }
  }
  for (int q = 0; q < v.nx(); q++) {
    for (int r = 0; r < v.ny(); r++) {
      for (int e = 0; e < v.nz(); e++) {
	v(q,r,e) = 0;
      }
    }
  }
  for (int q = 0; q < w.nx(); q++) {
    for (int r = 0; r < w.ny(); r++) {
      for (int e = 0; e < w.nz(); e++) {
	w(q,r,e) = 0;
      }
    }
  }
  //with forces
  SlArray3D<double> fu(params.res[0]+1,params.res[1],params.res[2]);
  SlArray3D<double> fv(params.res[0],params.res[1]+1,params.res[2]);
  SlArray3D<double> fw(params.res[0],params.res[1],params.res[2]+1);
  particleToGrid(particles, params, u, v, w);
 
  //run until the total time is reached
  while (time < params.total_time) {
    sprintf(fname, argv[2], frame);
    vector<SlVector3> positions;
    vector<SlVector3> velocities;
    //computation of velocities
    advect(particles, params, u, v, w);
    particleToGrid(particles, params, u, v, w);
    gridToParticle(particles, params, u, v, w, fu, fv, fw);
    //set up for the write and calculate new postion
    for (int i = 0; i < particles.size(); i++) {
      particles[i].pos = particles[i].pos + (particles[i].vel * params.dt); 
      for (int h = 0; h < 3; h++) {
	if (particles[i].pos[h] > upperCorner.data[h]) {
	  particles[i].pos[h] = upperCorner.data[h];
	}
      }
      positions.push_back(particles[i].pos);
      velocities.push_back(particles[i].vel);       
    }
    writeParticles(fname, numPart, positions, velocities);
    time+= params.dt;
    frame++;
  }
  return 0;
}
