#ifndef JOINT_H
#define JOINT_H

#include "slVector.H"
#include "slMatrix.H"
#include <vector>

class Joint {
 public:
  Joint(int weight, int parent, double x, double y, double z);
  Joint(const Joint &j);
  int getWeight();
  int getParent();
  double getX();
  double getY();
  double getZ();
  //calculate the outboard/inboard;
  SlVector3 calculateInBoard(std::vector<Joint> &joints, int currentJoint);
  SlVector3 calculateOutBoard(std::vector<Joint> &joints, int currentJoint);
  //calculate the rotation matrix, multiply the joint's rotation matrix and it's parent's rotation matrix
  SlMatrix3x3 calculateRotationMatrix(std::vector< std::vector<double> > &poses, std::vector<Joint> &joints, int currentJoint, int currentPose);

  //calculate the cordinates in world space.
  SlVector3 calculateOffset(std::vector< std::vector<double> > &poses, 
			    std::vector<Joint> &joints,
			    int currentJoint,
			    int currentPose);
 private:
  int m_weight;
  int m_parent;
  double m_x;
  double m_y;
  double m_z;
};
#endif
