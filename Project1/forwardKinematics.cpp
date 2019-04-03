#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>
#include <sstream>
#include <vector>

#include "Joint.h"
#include "slVector.H"
#include "slRotations.H"
#include "slMatrix.H"
#include "slIO.H"
using namespace std;


int main() {
  vector<Joint> jointList;

  ifstream skeletonFile;
  ifstream poseFile;
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

  vector<double> rotations;
  vector< vector <double> > poses;

  // open and parse pose file
  poseFile.open("pose.dmat");
  if (poseFile.fail()) {
    cout << "The pose file was not successfully opened." << endl;
    exit(1);
  }
  while (getline(poseFile,line)) {
    stringstream ss(line);
    if (line == "" ) {
      poses.push_back(rotations);
      rotations.clear();
    }
    while (ss >> buf) {
      double angle = atof(buf.c_str());
      rotations.push_back(angle);
    }
  }
  int numPoses = poses[0][0];
  int numAngles = poses[0][1];
  poses[0].erase(poses[0].begin());
  poses[0].erase(poses[0].begin());
  //finished formating the rotations.
  
  //applying rotations and finding offset in world space.
  for (int p = 0; p < numPoses; p++) {
    //open file for writing
    ofstream outFile;
    char buffer[20];
    sprintf(buffer, "output-%04d.pose",p); 
    outFile.open(buffer);
    
    //calculating offsets
    int currentEulerAngle = 0;
    for (int j = 0; j < jointList.size(); j++) {
      Joint joint = jointList[j];
      int weight = joint.getWeight();
      SlVector3 offsetVector = joint.calculateOffset(poses, jointList, j, p);
      outFile << weight << " " << offsetVector.data[0] << " " << offsetVector.data[1] << " " << offsetVector.data[2]  << endl;
    }    
    outFile << endl;
  } 
  return 0;
}
