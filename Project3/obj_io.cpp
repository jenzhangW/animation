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

bool loadObj(char *fname, std::vector<SlVector3> &pts, std::vector<SlTri> &triangles) {
  char c[500];
  
  int numVertices=0, numFaces=0;
  bool normals = false, texture = false;
  int tint;
  char ch;
  int p, q, r;
  double x, y, z;
  std::vector<SlVector3>::iterator v;
  std::vector<SlTri>::iterator t;
  std::ifstream in1(fname, std::ios::in);
  if (!in1.is_open()) {
	return false;
  }
  in1.flags(in1.flags() & ~std::ios::skipws);
  
  while (in1>>ch) {
    if (ch == 'v') {
      in1>>ch;
      if (ch == ' ') numVertices++;
      else if (ch == 'n') normals = true;
      else if (ch == 't') texture = true;
      else std::cerr<<"error \'"<<ch<<"\'"<<std::endl;
    } else if (ch == '#') {
	  while (in1 >> ch && ch != '\n') ; // Read to the end of the line.
    } else if (ch == 'f') numFaces++;
  }
  in1.close();
  
  pts.resize(numVertices);
  triangles.resize(numFaces);
  v = pts.begin();
  t = triangles.begin();
  
  std::ifstream in(fname, std::ios::in);
  if (!in.is_open()) {
	return false;
  }
  
  while (in>>ch) {
    if (ch == '#') {
      in.getline(c,500);
      continue;
    }
    if (ch == 'g') {
      in.getline(c,500);
      continue;
    }
    if (ch == 's') {
      in.getline(c,500);
      continue;
    }
    if (ch == 'm') {
      in.getline(c,500);
      continue;
    }
    if (ch == 'u') {
      in.getline(c,500);
      continue;
    }
    if (ch == 'v') {
      ch = in.peek();
      if (ch != 't' && ch != 'n') {
		in>>x>>y>>z;
		(*v).set(x,y,z);
		v++;
      } else {
		in.getline(c, 500);
      }
      continue;
    }
    if (ch == 'f') {
      if (normals && texture) {
		in>>p>>ch>>tint>>ch>>tint>>q>>ch>>tint>>ch>>tint>>r>>ch>>tint>>ch>>tint;
      } else if (normals) {
		in>>p>>ch>>ch>>tint>>q>>ch>>ch>>tint>>r>>ch>>ch>>tint;
      } else if (texture) {
		in>>p>>ch>>tint>>q>>ch>>tint>>r>>ch>>tint;
      } else {
		in>>p>>q>>r;
      }
      (*t)[0] = p-1;
      (*t)[1] = q-1;
      (*t)[2] = r-1;
      t++;
      continue;
    }
  }
  in.close();
  return true;
}

void writeObj(char *fname, const std::vector<SlVector3> &meshPts, const std::vector<SlTri> &triangles) {
  std::ofstream out;
  std::vector<SlVector3>::const_iterator p;
  std::vector<SlTri>::const_iterator t;

  out.open(fname);
  
  for (p=meshPts.begin(); p!=meshPts.end(); p++) 
	out<<"v "<<(*p)[0]<<" "<<(*p)[1]<<" "<<(*p)[2]<<std::endl;
  
  for (t=triangles.begin(); t!=triangles.end(); t++) 
	out<<"f "<<(*t)[0]+1<<" "<<(*t)[1]+1<<" "<<(*t)[2]+1<<std::endl;
  
  out.close();
}
