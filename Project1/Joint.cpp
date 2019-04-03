#include <iostream>
#include <vector>
using namespace std;

#include "Joint.h"

Joint::Joint(int weight, int parent, double x, double y, double z){
  m_weight = weight;
  m_parent = parent;
  m_x = x;
  m_y = y;
  m_z = z;
}

Joint::Joint(const Joint &j) {
  m_weight = j.m_weight;
  m_parent = j.m_parent;
  m_x = j.m_x;
  m_y = j.m_y;
  m_z = j.m_z;
}
int Joint::getWeight() {
  return m_weight;
}

int Joint::getParent() {
  return m_parent;
}
 
double Joint::getX() {
  return m_x;
}

double Joint::getY() {
  return m_y;
}

double Joint::getZ() {
  return m_z;
}

SlMatrix3x3 Joint::calculateRotationMatrix(std::vector< std::vector<double> > &poses, std::vector<Joint> &joints, int currentJoint, int currentPose) {
  double pi = 3.1415926535897;
  
  //findd the joint's respective Euler angles.
  Joint joint = joints[currentJoint];
  SlVector3 eulerVector(poses[currentPose][currentJoint*3],
			poses[currentPose][currentJoint*3 + 1],
			poses[currentPose][currentJoint*3 + 2]);
  //converts to radians
  for (int i = 0; i < 3; i++) {
    eulerVector.data[i] = (eulerVector.data[i]/360) * 2 * pi;
  } 
  
  //computes the joints rotation matrix based on the past 
  SlMatrix3x3 rotationMatrix;
  SlEulerAngToMatrixXYZ(eulerVector, rotationMatrix);
  int parent = joint.getParent();
  if (currentJoint > 0) {
    SlMatrix3x3 parentRotationMatrix = calculateRotationMatrix(poses, joints, parent, currentPose); 
    return rotationMatrix * parentRotationMatrix;  
  }
  else {
    return rotationMatrix;  
  }
}

SlVector3 Joint::calculateOffset(std::vector< std::vector<double> > &poses, std::vector<Joint> &joints, int currentJoint, int currentPose) {
  Joint joint = joints[currentJoint];
  SlVector3 jointVector(joint.getX(), joint.getY(), joint.getZ());
  
  //returns the rotation matrix, then multiply the current joint by it
  SlMatrix3x3 rotationMatrix = calculateRotationMatrix(poses, joints, currentJoint, currentPose);
  int parent = joint.getParent();  
  SlVector3 newJointVector = rotationMatrix * jointVector;

  //adds up the new joint with it's parents offsets
  if (currentJoint > 0) {
  SlVector3 parentJointVector = calculateOffset(poses, joints, parent, currentPose); 

  double x = newJointVector.data[0] + parentJointVector.data[0];
  double y = newJointVector.data[1] + parentJointVector.data[1];
  double z = newJointVector.data[2] + parentJointVector.data[2];
  
  SlVector3 finalJointVector(x, y, z);
  return finalJointVector;
  }

  else {
    return newJointVector;

  }
}
