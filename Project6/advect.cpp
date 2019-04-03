#include <iostream>
#include <fstream>
#include <vector>
#include <slVector.H>
#include <json/json.h>
#include <slArray.H>

static const double FEPS = 1e-6;

struct SimulationParameters {
  double dt, total_time, density;
};

class StaggeredGrid {
public:
  SlArray3D<double> u,v,w;
  SlArray3D<unsigned char> cellLabels;
  unsigned int nx, ny, nz, nynz;
  double h,halfh;
  SlVector3 lc, uc;
  
  bool allocate(unsigned int nx, unsigned int ny, unsigned int nz, double h, SlVector3 lc, SlVector3 uc);
  
  SlVector3 interpVel(const SlVector3 &x) const;
  double interpxface(const SlVector3 &x, const SlArray3D<double> &f) const;
  double interpyface(const SlVector3 &x, const SlArray3D<double> &f) const;
  double interpzface(const SlVector3 &x, const SlArray3D<double> &f) const;
  
  bool advectParticles(double dt);
  bool particlesToGrid();
  bool gridToParticles();
  
  // FLIP stuff
  std::vector<SlVector3> flipPos, flipVel;
  unsigned int nFlipParticles;
  double flipRatio;
  
  static const unsigned char AIR = 1;
  static const unsigned char LIQUID = 2;
  static const unsigned char SURFACE = 4;
  static const unsigned char OBSTACLE = 8;
  static const unsigned char FLUID = 6;
  
private:
  bool particleToGrid(SlArray3D<double> &u, SlArray3D<double> &fu,
	  double w0, double w1, double w2,
	  int i, int j, int k, double v);
  bool setBoundaryVelocity();
  inline bool clipToGrid(SlVector3 &x) const;
  inline bool boundary(unsigned int i, unsigned int j, unsigned int k) const;
  SlArray3D<double> nu,nv,nw,fu,fv,fw;
};


inline bool StaggeredGrid::clipToGrid(SlVector3 &x) const {
	//static const double FEPS = 1e-6;
	static const double FEPS = 1e-4;
	if (x[0] <= lc[0]+h+FEPS) x[0] = lc[0]+h+FEPS;
	if (x[0] >= uc[0]-h-FEPS) x[0] = uc[0]-h-FEPS;
	if (x[1] <= lc[1]+h+FEPS) x[1] = lc[1]+h+FEPS;
	if (x[1] >= uc[1]-h-FEPS) x[1] = uc[1]-h-FEPS;
	if (x[2] <= lc[2]+h+FEPS) x[2] = lc[2]+h+FEPS;
	if (x[2] >= uc[2]-h-FEPS) x[2] = uc[2]-h-FEPS;
	return true;
};

inline bool StaggeredGrid::boundary(unsigned int i, unsigned int j, unsigned int k) const {
	if (i == 0 || i == nx-1 || j == 0 || j == ny-1 || k == 0 || k == nz-1)
		return true;
	return false;
};

inline double trilinearInterp(const SlArray3D<double> &x, int i, int j, int k,
	double w0, double w1, double w2) {
	return (1.0-w0)*(1.0-w1)*(1.0-w2)*x(i,j,k)+
		(1.0-w0)*(1.0-w1)*(w2)*x(i,j,k+1)+
		(1.0-w0)*(w1)*(1.0-w2)*x(i,j+1,k)+
		(1.0-w0)*(w1)*(w2)*x(i,j+1,k+1)+
		(w0)*(1.0-w1)*(1.0-w2)*x(i+1,j,k)+
		(w0)*(1.0-w1)*(w2)*x(i+1,j,k+1)+
		(w0)*(w1)*(1.0-w2)*x(i+1,j+1,k)+
		(w0)*(w1)*(w2)*x(i+1,j+1,k+1);
}

bool StaggeredGrid::setBoundaryVelocity() {
	unsigned int i,j,k;

	// set boundary conditions
	for (j=0; j<ny; j++) {
		for (k=0; k<nz; k++) {
			u(0,j,k) = u(1,j,k) = 0.0;
			v(0,j,k) = v(1,j,k);
			w(0,j,k) = w(1,j,k);
			u(nx,j,k) = u(nx-1,j,k) = 0.0;
			v(nx-1,j,k) = v(nx-2,j,k);
			w(nx-1,j,k) = w(nx-2,j,k);
		}
	}
	for (i=0; i<nx; i++) {
		for (k=0; k<nz; k++) {
			u(i,0,k) = u(i,1,k);
			v(i,0,k) = v(i,1,k) = 0.0;
			w(i,0,k) = w(i,1,k);
			u(i,ny-1,k) = u(i,ny-2,k);
			v(i,ny,k) = v(i,ny-1,k) = 0.0;
			w(i,ny-1,k) = w(i,ny-2,k);
		}
	}
	for (i=0; i<nx; i++) {
		for (j=0; j<ny; j++) {
			u(i,j,0) = u(i,j,1);
			v(i,j,0) = v(i,j,1);
			w(i,j,0) = w(i,j,1) = 0.0;
			u(i,j,nz-1) = u(i,j,nz-2);
			v(i,j,nz-1) = v(i,j,nz-2);
			w(i,j,nz) = w(i,j,nz-1) = 0.0;
		}
	}
	return true;
}

bool StaggeredGrid::allocate(unsigned int nx, unsigned int ny, unsigned int nz, double h, SlVector3 lc, SlVector3 uc){
	this->nx = nx;
	this->ny = ny;
	this->nz = nz;
	this->nynz = ny*nz;
	this->h = h;
	this->halfh = h/2.0;
	this->lc = lc;
	this->uc = uc;
	u.allocate(nx+1, ny, nz);
	v.allocate(nx, ny+1, nz);
	w.allocate(nx, ny, nz+1);
	//p.allocate(nx, ny, nz);
	nu.allocate(nx+1, ny, nz);
	nv.allocate(nx, ny+1, nz);
	nw.allocate(nx, ny, nz+1);
	fu.allocate(nx+1, ny, nz);
	fv.allocate(nx, ny+1, nz);
	fw.allocate(nx, ny, nz+1);
	cellLabels.allocate(nx, ny, nz);
	//laplacian.allocate(nx, ny, nz);
	//precon.allocate(nx, ny, nz);
	//q.allocate(nx, ny, nz);
	//r.allocate(nx, ny, nz);
	//z.allocate(nx, ny, nz);
	//s.allocate(nx, ny, nz);

	return true;
}

double StaggeredGrid::interpxface(const SlVector3 &x, const SlArray3D<double> &f) const {
	SlVector3 y(x-lc);
	int i,j,k;
	y[1] -= halfh;
	y[2] -= halfh;
	i = (int)floor(y[0] / h);
	j = (int)floor(y[1] / h);
	k = (int)floor(y[2] / h);
	double w0, w1, w2;
	w0 = (y[0] - i*h)/h;
	w1 = (y[1] - j*h)/h;
	w2 = (y[2] - k*h)/h;

	return trilinearInterp(f, i, j, k, w0, w1, w2);
}

double StaggeredGrid::interpyface(const SlVector3 &x, const SlArray3D<double> &f) const {
	SlVector3 y(x-lc);
	int i,j,k;
	y[0] -= halfh;
	y[2] -= halfh;
	i = (int)floor(y[0] / h);
	j = (int)floor(y[1] / h);
	k = (int)floor(y[2] / h);
	double w0, w1, w2;
	w0 = (y[0] - i*h)/h;
	w1 = (y[1] - j*h)/h;
	w2 = (y[2] - k*h)/h;

	return trilinearInterp(f, i, j, k, w0, w1, w2);
}

double StaggeredGrid::interpzface(const SlVector3 &x, const SlArray3D<double> &f) const {
	SlVector3 y(x-lc);
	int i,j,k;
	y[0] -= halfh;
	y[1] -= halfh;
	i = (int)floor(y[0] / h);
	j = (int)floor(y[1] / h);
	k = (int)floor(y[2] / h);
	double w0, w1, w2;
	w0 = (y[0] - i*h)/h;
	w1 = (y[1] - j*h)/h;
	w2 = (y[2] - k*h)/h;

	return trilinearInterp(f, i, j, k, w0, w1, w2);
}

SlVector3 StaggeredGrid::interpVel(const SlVector3 &x) const {
	return SlVector3(interpxface(x,u), interpyface(x,v), interpzface(x,w));
}

bool StaggeredGrid::advectParticles(double dt) {
  int i=0;
  for (std::vector<SlVector3>::iterator o = flipPos.begin(); o != flipPos.end(); o++, i++) {
	SlVector3 n = (*o) + dt*interpVel((*o));
	clipToGrid(n);
    (*o) = n;
  }
  return true;
}

bool StaggeredGrid::gridToParticles() {
  std::vector<SlVector3>::iterator pos = flipPos.begin();
  std::vector<SlVector3>::iterator vel = flipVel.begin();
  for (unsigned int i=0; i<nFlipParticles; i++, pos++, vel++) {
	SlVector3 oVel(interpxface((*pos),fu), interpyface((*pos),fv), interpzface((*pos),fw));
	SlVector3 nVel(interpxface((*pos),u), interpyface((*pos),v), interpzface((*pos),w));
	//std::cout<<(*vel)<<" "<<oVel<<" "<<nVel<<std::endl;
	(*vel) = flipRatio*((*vel)-oVel) + nVel;
  }
  return true;
}

bool StaggeredGrid::particleToGrid(SlArray3D<double> &u, SlArray3D<double> &fu,
	double w0, double w1, double w2,
	int i, int j, int k, double v) {

  double w = (1.0-w0)*(1.0-w1)*(1.0-w2);
  u(i,j,k) += w*v;
  fu(i,j,k) += w;
  
  w = (1.0-w0)*(1.0-w1)*(w2);
  u(i,j,k+1) += w*v;
  fu(i,j,k+1) += w;
  
  w = (1.0-w0)*w1*(1.0-w2);
  u(i,j+1,k) += w*v;
  fu(i,j+1,k) += w;
  
  w = (1.0-w0)*w1*w2;
  u(i,j+1,k+1) += w*v;
  fu(i,j+1,k+1) += w;
  
  w = w0*(1.0-w1)*(1.0-w2);
  u(i+1,j,k) += w*v;
  fu(i+1,j,k) += w;

  w = w0*(1.0-w1)*w2;
  u(i+1,j,k+1) += w*v;
  fu(i+1,j,k+1) += w;

  w = w0*w1*(1.0-w2);
  u(i+1,j+1,k) += w*v;
  fu(i+1,j+1,k) += w;
  
  w = w0*w1*w2;
  u(i+1,j+1,k+1) += w*v;
  fu(i+1,j+1,k+1) += w;
  
  return true;
}

bool StaggeredGrid::particlesToGrid() {
  for (unsigned int i=0; i<nx+1; i++) {
	for (unsigned int j=0; j<ny; j++) {
	  for (unsigned int k=0; k<nz; k++) {
		u(i,j,k) = 0.0;
		fu(i,j,k) = 0.0;
	  }
	}
  }
  for (unsigned int i=0; i<nx; i++) {
	for (unsigned int j=0; j<ny+1; j++) {
	  for (unsigned int k=0; k<nz; k++) {
		v(i,j,k) = 0.0;
		fv(i,j,k) = 0.0;
	  }
	}
  }
  for (unsigned int i=0; i<nx; i++) {
	for (unsigned int j=0; j<ny; j++) {
	  for (unsigned int k=0; k<nz+1; k++) {
		w(i,j,k) = 0.0;
		fw(i,j,k) = 0.0;
	  }
	}
  }
  
  for (unsigned int i=0; i<nx; i++) {
	for (unsigned int j=0; j<ny; j++) {
	  for (unsigned int k=0; k<nz; k++) {
		if (boundary(i,j,k)) {
		  cellLabels(i,j,k) = OBSTACLE;
		  continue;
		}
		cellLabels(i,j,k) = AIR;
	  }
	}
  }
  
  std::vector<SlVector3>::iterator pos = flipPos.begin();
  std::vector<SlVector3>::iterator vel = flipVel.begin();
  for (unsigned int p=0; p<nFlipParticles; p++, pos++, vel++) {
	SlVector3 y((*pos)-lc);
	int i,j,k;
	double w0, w1, w2;
	
	i = (int)floor(y[0] / h);
	j = (int)floor(y[1] / h);
	k = (int)floor(y[2] / h);
	cellLabels(i,j,k) = LIQUID;
		
	y[1] -= halfh;
	y[2] -= halfh;
	i = (int)floor(y[0] / h);
	j = (int)floor(y[1] / h);
	k = (int)floor(y[2] / h);
	w0 = (y[0] - i*h)/h;
	w1 = (y[1] - j*h)/h;
	w2 = (y[2] - k*h)/h;
	particleToGrid(u, fu, w0, w1, w2, i, j, k, (*vel)[0]);
	
	y[0] -= halfh;
	y[1] += halfh;
	i = (int)floor(y[0] / h);
	j = (int)floor(y[1] / h);
	k = (int)floor(y[2] / h);
	w0 = (y[0] - i*h)/h;
	w1 = (y[1] - j*h)/h;
	w2 = (y[2] - k*h)/h;
	particleToGrid(v, fv, w0, w1, w2, i, j, k, (*vel)[1]);
	
	y[1] -= halfh;
	y[2] += halfh;
	i = (int)floor(y[0] / h);
	j = (int)floor(y[1] / h);
	k = (int)floor(y[2] / h);
	w0 = (y[0] - i*h)/h;
	w1 = (y[1] - j*h)/h;
	w2 = (y[2] - k*h)/h;
	particleToGrid(w, fw, w0, w1, w2, i, j, k, (*vel)[2]);
  }
  
  for (unsigned int i=0; i<nx+1; i++) {
	for (unsigned int j=0; j<ny; j++) {
	  for (unsigned int k=0; k<nz; k++) {
		if (i == 0 || i == 1 || i == nx || i == nx-1) {
		  u(i,j,k) = 0.0;
		} else {
		  if (fu(i,j,k) > FEPS) u(i,j,k) /= fu(i,j,k);
		}
		fu(i,j,k) = u(i,j,k);
	  }
	}
  }
  for (unsigned int i=0; i<nx; i++) {
	for (unsigned int j=0; j<ny+1; j++) {
	  for (unsigned int k=0; k<nz; k++) {
		if (j == 0 || j == 1 || j == ny || j == ny-1) {
		  v(i,j,k) = 0.0;
		} else {
		  if (fv(i,j,k) > FEPS) v(i,j,k) /= fv(i,j,k);
		}
		fv(i,j,k) = v(i,j,k);
	  }
	}
  }
  for (unsigned int i=0; i<nx; i++) {
	for (unsigned int j=0; j<ny; j++) {
	  for (unsigned int k=0; k<nz+1; k++) {
		if (k == 0 || k == 1 || k == nz || k == nz-1) {
		  w(i,j,k) = 0.0;
		} else {
		  if (fw(i,j,k) > FEPS) w(i,j,k) /= fw(i,j,k);
		}
		fw(i,j,k) = w(i,j,k);
	  }
	}
  }
  
  for (unsigned int i=0; i<nx; i++) {
	for (unsigned int j=0; j<ny; j++) {
	  for (unsigned int k=0; k<nz; k++) {
		if (cellLabels(i,j,k) == FLUID) {
		  if (cellLabels(i-1,j,k) == AIR ||
			  cellLabels(i+1,j,k) == AIR ||
			  cellLabels(i,j-1,k) == AIR ||
			  cellLabels(i,j+1,k) == AIR ||
			  cellLabels(i,j,k-1) == AIR ||
			  cellLabels(i,j,k+1) == AIR) {
			cellLabels(i,j,k) = SURFACE;
		  }
		}
	  }
	}
  }

  setBoundaryVelocity();
  return true;
}

bool readParticles(const char *fname, unsigned int &nFlipParticles, std::vector<SlVector3> &flipPos, std::vector<SlVector3> &flipVel) {
  std::ifstream in(fname, std::ios::in);
  in>>nFlipParticles;
  flipPos.resize(nFlipParticles);
  flipVel.resize(nFlipParticles);
  for (unsigned int i=0; i<nFlipParticles; i++) {
	in>>flipPos[i][0]>>flipPos[i][1]>>flipPos[i][2]>>flipVel[i][0]>>flipVel[i][1]>>flipVel[i][2];
  }
  in.close();
  std::cout<<"inputfile " << fname <<" read"<<std::endl;
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
	StaggeredGrid &grid) {
  
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
  grid.flipRatio = root.get("flipRatio", 0.95).asDouble();

  int nx = root["res"][0].asInt();
  int ny = root["res"][1].asInt();
  int nz = root["res"][2].asInt();
  double h = root["h"].asDouble();
  SlVector3 lc, uc;
  lc.set(root["lc"][0].asDouble(), root["lc"][1].asDouble(), root["lc"][2].asDouble());
  uc.set(root["uc"][0].asDouble(), root["uc"][1].asDouble(), root["uc"][2].asDouble());
  uc = lc + SlVector3(h*nx, h*ny, h*nz);
  grid.allocate(nx, ny, nz, h, lc, uc);

  readParticles(root["particles"].asString().c_str(), grid.nFlipParticles, grid.flipPos, grid.flipVel);
  return true;
}

int main(int argc, char *argv[]) {
  char fname[80];
  SimulationParameters params;
  StaggeredGrid grid;
  double time = 0;
  int frame = 0;
  double frameTime = -1.0;

  readInputFile(argv[1], params, grid);
  grid.particlesToGrid();

  std::cout<<"Domain is "<<grid.lc<<"x"<<grid.uc<<std::endl;
  
  while (time < params.total_time) {
	if (frameTime < 0) {
	  sprintf(fname, argv[2], frame);
	  writeParticles(fname, grid.nFlipParticles, grid.flipPos, grid.flipVel);
	  frameTime = 1.0/30.0-0.0001;
	  frame++;
	}

	grid.advectParticles(params.dt);
	grid.particlesToGrid();

	// other fluid simulation steps

	grid.gridToParticles();

	time += params.dt;
	frameTime -= params.dt;
  }
}
