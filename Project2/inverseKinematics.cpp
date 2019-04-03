/*Description
  -Takes in a .bf and a .con file
  -Changes the angles of joints to match the position in world space of the joints to the constraints.
  -Calculates a Jacobian, then the transpose and then calculates theta
  -outputs a file with world coordinates and the error.
 */
#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>
#include <sstream>
#include <vector>
#include <cmath>  

#include "Joint.h"
#include "slVector.H"
#include "slRotations.H"
#include "slMatrix.H"
#include "slIO.H"
using namespace std;

void calculateJacobian(std::vector< std::vector<double> > &poses, std::vector< std::vector <double> > &jacobian, std::vector<Joint> &endeffector, std::vector<Joint> &joints) {
  //finds the endefector joint
  //calculate the vector between the inboard of a joint and the endeffector
  for (int i = 0; i < endeffector.size(); i++) {
    vector<int> parents;
    int index = endeffector[i].getWeight();
    //find the joints that are affect by the endeffector
    while (index > 0) {
      parents.push_back(joints[index].getParent());
      index = joints[index].getParent();
    }
    int end = endeffector[i].getWeight();
    SlVector3 endJoint = joints[end].calculateOffset(poses, joints, end, 0); 
    //looping through the joints and seeing which is affected by the endeffector
    for (int j = 1; j < joints.size(); j++) {
      for (int k = 0; k < parents.size(); k++) {
	if ( j == parents[k]) {
	  //inboard rotation and position
	  int currentParent = joints[j].getParent();
	  SlVector3 currentJoint = joints[currentParent].calculateOffset(poses,joints, currentParent, 0);
	  SlMatrix3x3 rotation = joints[currentParent].calculateRotationMatrix(poses,joints,currentParent,0);
	  //calculate jacobian entries
	  SlVector3 difference = endJoint - currentJoint;
	  for (int l = 0; l < 3; l++) {
	    switch(l){
	    case 0: {
	      SlVector3 xAxis(rotation.data[0][0],rotation.data[0][1],rotation.data[0][2]);
	      SlVector3 crossX = cross(xAxis,difference);
	      //adding into jacobian
	      jacobian[i * 3][j * 3 + l] = crossX[0];
	      jacobian[i * 3 + 1][j * 3 + l] = crossX[1];
	      jacobian[i * 3 + 2][j * 3 + l] = crossX[2];
	    }
	      break;
	    case 1: {
	      //computing x rotation
	      double xTheta = poses[0][currentParent * 3];
	      SlMatrix3x3 yRotation = rotation * SlMatrix3x3(1,0,0,0,cos(xTheta),-sin(xTheta),0,sin(xTheta),cos(xTheta)); 
	      SlVector3 yAxis(yRotation.data[1][0],yRotation.data[1][1],yRotation.data[1][2]);
	      SlVector3 crossY = cross(yAxis, difference);
	      //adding into jacobian
	      jacobian[i * 3][j * 3 + l] = crossY[0];
	      jacobian[i * 3 + 1][j * 3 + l] = crossY[1];
	      jacobian[i * 3 + 2][j * 3 + l] = crossY[2];
	    }
	      break;
	    case 2: {
	      //computing x rotation
	      double xTheta2 = poses[0][currentParent * 3];
	      SlMatrix3x3 yRotation2 = rotation * SlMatrix3x3(1,0,0,0,cos(xTheta2),-sin(xTheta2),0,sin(xTheta2),cos(xTheta2)); 
	      //computing y rotation
	      double yTheta = poses[0][currentParent * 3 + 1];
	      SlMatrix3x3 zRotation = yRotation2 * SlMatrix3x3(cos(yTheta),0,sin(yTheta), 0, 1, 0, -sin(yTheta), 0, cos(yTheta)); 
	      SlVector3 zAxis(zRotation.data[2][0],zRotation.data[2][1],zRotation.data[2][2]);
	      SlVector3 crossZ = cross(zAxis, difference);
	      //adding into jacobian
	      jacobian[i * 3][j * 3 + l] = crossZ[0];
	      jacobian[i * 3 + 1][j * 3 + l] = crossZ[1];
	      jacobian[i * 3 + 2][j * 3 + l] = crossZ[2];
	      }
	      break;
	    }
	  }
	}
      }
    }
  }
}

//just inverses indices to create the transpose
void calculateTranspose(std::vector< std::vector<double> > &jacobian, std::vector< std::vector<double> > &tpose)  {
  for (int i = 0; i < jacobian.size(); i++) {
    for (int j = 0; j < jacobian[0].size(); j++) {
      tpose[j][i] = jacobian[i][j];
    }
  }  
  /*  cout << "[";
  for (int x = 0; x < tpose.size(); x++) {
    for (int y = 0; y < tpose[x].size(); y++) {
      cout << tpose[x][y] << ",";
    }
    cout << "]" << endl ; 
  }
  cout << endl;*/
}

//subtracting the position goal from the current endeffetors to get error
double calculateError(std::vector< std::vector<double> > &poses, std::vector<Joint> &joints, int currentJoint, int currentPose, std::vector<int> &endeffector, std::vector<Joint> &endeffectorJoints) {
  double errorSum;
  for (int i = 0; i < endeffector.size(); i++) { 
    int end = endeffector[i] + 1;
    Joint joint = joints[end];
    SlVector3 offsetVector = joint.calculateOffset(poses, joints, end,0); 
    double x = endeffectorJoints[i].getX();
    double y = endeffectorJoints[i].getY();
    double z = endeffectorJoints[i].getZ();
    errorSum += (x - offsetVector[0]) + (y - offsetVector[1]) + (z - offsetVector[2]);
  }
}

//find the change in angle and then change the vector that holds the euler  angles.
vector<double> calculateDeltaTheta(std::vector< std::vector<double> > &poses, std::vector< std::vector<double> > &transpose,  std::vector<Joint> &joints,std::vector<int> &endeffector, std::vector<Joint> &endeffectorJoints ) {
  vector<double> difference;
  //finds the difference between the goal and the endeffectors current position.
  for (int i = 0; i < endeffector.size(); i++) {
    int end = endeffector[i] + 1;
    Joint joint = joints[end];
    SlVector3 offsetVector = joint.calculateOffset(poses, joints, end, 0);
    double x = endeffectorJoints[i].getX();
    double y = endeffectorJoints[i].getY();
    double z = endeffectorJoints[i].getZ();
    difference.push_back(joint.getX() - x);
    difference.push_back(joint.getY() - y);
    difference.push_back(joint.getY() - z);
  }
  //calculates change in theta
  vector<double> deltaTheta;
  for (int j = 0; j < transpose.size(); j++) {
    double sum  = 0;
    for (int k = 0; k < transpose[j].size(); k++) {
      sum = sum + transpose[j][k] * difference[k];
    }
    deltaTheta.push_back(sum);
  }
  return deltaTheta;
}

//changes the euler angles
void calculateTheta(std::vector< std::vector<double> > &poses,std::vector< std::vector<double> > &transpose, std::vector<Joint> &joints, std::vector<int> &endeffector, std::vector<Joint> &endeffectorJoints, double step) {
  vector<double> deltaTheta = calculateDeltaTheta(poses,transpose, joints, endeffector, endeffectorJoints);
  for (int i = 0; i < deltaTheta.size(); i++) {
    deltaTheta[i] = deltaTheta[i]*step;
    poses[0][i] += deltaTheta[i];
  }
}

int main(int argc, char**argv) {
  vector<Joint> jointList;
  ifstream skeletonFile;
  ifstream constraintFile;
  string line;
  string buf;
  // open the skeleteon file to get the joints at rest and parse the file
  skeletonFile.open("ogre-skeleton.bf");
  if (skeletonFile.fail()) {
    cout << "The skeleton file was not successfully opened." << endl;
    exit(1);
  }
  while (getline(skeletonFile,line)) {
    stringstream ss(line);
    vector<string> data;
    while (ss >> buf) {
      data.push_back(buf);
    }
    int weight = atoi(data[0].c_str());
    int parent = atoi(data[1].c_str());
    double x = atof(data[2].c_str());
    double y = atof(data[3].c_str());
    double z = atof(data[4].c_str());
    Joint newJoint(weight, parent, x, y, z);
    jointList.push_back(newJoint);
  }
  skeletonFile.close();
  // finished formating the skeleton hierachy.

  // formating constraints
  vector<Joint> constraintJoints;
  constraintFile.open(argv[1]);
  if (constraintFile.fail()) {
    cout << "The constraint file was not successfully opened." << endl;
    exit(1);
  }    
  vector<string> con;
  while (getline(constraintFile,line)) {
    stringstream ss(line);
    while (ss >> buf) {
      con.push_back(buf);
    }
  }
  //holds the weight of the endeffector
  vector<int> end;
  int numConstraints = atoi(con[0].c_str());
  con.erase(con.begin());
  for (int a = 0; a < con.size(); a = a+4) {
    int weight = atoi(con[a].c_str());
    double x = atof(con[a + 1].c_str());
    double y = atof(con[a + 2].c_str());
    double z = atof(con[a + 3].c_str());
    Joint newJoint(weight, 0, x, y, z);
    constraintJoints.push_back(newJoint);
    end.push_back(weight);
  }
  // finished formating constraints
  
  //finding the jacobian and the jacobian transpose
  vector< vector<double> > jacobian;
  for (int c = 0; c < (numConstraints * 3); c++) {
    vector< double> jacobianIn;
    for (int d = 0; d < 69; d++) {
      jacobianIn.push_back(0);
    }
    jacobian.push_back(jacobianIn);
  }

  //initializing the tranpose matrix and the matrix of euelr angles
  vector< vector<double> > transpose;
  for (int f = 0; f < 69; f++) {
    vector<double> transposeIn;
    for (int g = 0; g < (numConstraints * 3); g++) {
      transposeIn.push_back(0);    
    }
    transpose.push_back(transposeIn);
  }
  //the vector that holds the angles of rotation
  vector<double> rotations;
  vector< vector<double> > poses;
  for (int e = 0; e < 69; e++) {
    rotations.push_back(0);
  }
  poses.push_back(rotations);
  
  //calculating the difference and seeing it it's close enough to the endeffector
  double stepsize = 1;
  double error = (calculateError(poses, jointList, 0, 0, end, constraintJoints));

  //continues to compute until the error is less than a threshold.
  while (error > .0001) {
    cout << error << endl;
    calculateJacobian(poses, jacobian, constraintJoints, jointList);
    calculateTranspose(jacobian, transpose);
    calculateTheta(poses, transpose, jointList, end, constraintJoints,stepsize);
    error = abs(calculateError(poses, jointList, 0, 0, end, constraintJoints));
  }
  
  //printing world coordinates to file using forward kinematics for checking
  ofstream outFile;
  char buffer[20];
  sprintf(buffer, "output-00000.pose");
  outFile.open(buffer);
  //calculating offsets
  int currentEulerAngle = 0;
  for (int j = 0; j < jointList.size(); j++) {
    Joint joint = jointList[j];
    int weight = joint.getWeight();
    SlVector3 offsetVector = joint.calculateOffset(poses, jointList, j, 0);
    outFile << weight << " " << offsetVector.data[0] << " " << offsetVector.data[1] << " " << offsetVector.data[2]  << endl;
  }
  outFile << endl;
  
  return 0;
}

