#include <slVector.H>
#include <vector>
int main(int argc, char *argv[]) {
  std::vector<SlVector3> pos, vel;
  for (double x = -0.1; x <= 0.1; x+=0.005) {
	for (double y = -0.1; y <= 0.1; y+=0.005) {
	  for (double z = -0.1; z <= 0.1; z+=0.005) {
		SlVector3 p(x,y,z);
		//std::cout<<sqrMag(p)<<std::endl;
		if (sqrMag(p) > 0.01) continue;
		pos.push_back(p);
		vel.push_back(cross(p,SlVector3(0.0,0.0,1.0)));
	  }
	}
  }
  std::cout<<pos.size()<<std::endl;
  for (signed int i=0; i<pos.size(); i++) {
	std::cout<<pos[i][0]<<" "<<pos[i][1]<<" "<<pos[i][2]<<" "<<vel[i][0]<<" "<<vel[i][1]<<" "<<vel[i][2]<<std::endl;
  }
}
