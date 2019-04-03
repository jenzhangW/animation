//-------------------------------------------------------------------
//-------------------------------------------------------------------
//
// Simple Spring Mass System
// -- rotations lib
//
// Primary Author: James F. O'Brien (obrienj@cc.gatech.edu)
//
//-------------------------------------------------------------------
//-------------------------------------------------------------------
//
// Copyright (c) 2003-2005, Regents of the University of California.  All
// rights reserved.
//
// This software is part of the Berkeley Fluid Animation & Simulation
// Toolkit.  The name "Berkeley Fluid Animation & Simulation Toolkit" is
// a trademark of the Regents of the University of California.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
//   Redistributions of source code must retain the above copyright
//  notice, this list of conditions and the following disclaimer.
//
//  Redistributions in binary form must reproduce the above copyright
//  notice, this list of conditions and the following disclaimer in the
//  documentation and/or other materials provided with the distribution.
//
//  Redistributions in binary form as an executable program, or part of
//  an executable program must, when run for the first time by a given
//  user, prominently display the above copyright notice, this list of
//  conditions and the following disclaimer.
//
//  Neither the name of the University of California, Berkeley nor the
//  names of its contributors may be used to endorse or promote products
//  derived from this software without specific prior written
//  permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
//-------------------------------------------------------------------
//-------------------------------------------------------------------



#include "slVector.H"
#include "slMatrix.H"
#include "slRotations.H"

#ifndef M_PI
#define M_PI		3.14159265358979323846
#endif

//-------------------------------------------------------------------

void SlMomentToPosMatrix(const SlVector3 &v, SlMatrix3x3 &m) {
  // from Graphics Gems Vol I, page 466
  // With a transpose!

  register double p = mag(v);

  if (p!=0.0) {

    register double d = 1.0/p;
    register double x = v[0] * d;
    register double y = v[1] * d;
    register double z = v[2] * d;

    register double s = sin(p);
    register double c = cos(p);
    register double t = 1.0-c;

    m(0,0) = t*x*x+c;
    m(1,0) = t*x*y+s*z;
    m(2,0) = t*x*z-s*y;
    
    m(0,1) = t*x*y-s*z;
    m(1,1) = t*y*y+c;
    m(2,1) = t*y*z+s*x;

    m(0,2) = t*x*z+s*y;
    m(1,2) = t*y*z-s*x;
    m(2,2) = t*z*z+c;

  }else{
    m.setIdentity();
  }
}
//-------------------------------------------------------------------

void SlMomentToVelMatrix(const SlVector3 &v, SlMatrix3x3 &m) {
  // adapted from Graphics Gems Vol I, page 466
  // Set m so that m*b == vxb
  m(0,0) =   0.0;
  m(0,1) = -v[2];
  m(0,2) =  v[1];

  m(1,0) =  v[2]; 
  m(1,1) =  0.0;
  m(1,2) = -v[0];

  m(2,0) = -v[1];
  m(2,1) =  v[0];
  m(2,2) =   0.0;

}


//-------------------------------------------------------------------

void SlPosMatrixToMoment(const SlMatrix3x3 &m, SlVector3 &v) {
  // from Graphics Gems Vol I, page 466
  // With a transpose!

  // ** The "gem" has a bug in that it does not deal with
  // the singularity at +/-Pi in a reasonable way.
  // The singularity at 0 (aka +/-2Pi) occurs when theta
  // vanishes and v->[0,0,0].  At +/-Pi the limit does
  // not exist because approach from + gives a vector v1
  // approach from - gives v2==-v1.  This is where the
  // "wraps" overlap in an infinitessimal area.
  // So I have added a fix to detect this and fix it.
  // Unfortunatly, I am not sure WHY the fix I came up
  // works... I basicly used some intuition and made it
  // up.  So while it works, I can not prove it WILL
  // always work....

  register double cosT = (m(0,0)+m(1,1)+m(2,2)-1)*0.5;
  if (cosT >  1) cosT =  1;
  if (cosT < -1) cosT = -1;
  register double t = acos(cosT);
  if (t!=t) t = 0.0;
  register double sinT2 = 2*sin(t);

  if ((sinT2<1e-6) && (sinT2>-1e-6)) {
    if ((t<1e-6) && (t>-1e-6)) {
      v = 0.0;  // Singularity at 0,2Pi,4Pi....
    }else{
      SlMatrix3x3 s = m;
      s(0,0) += 1.0;
      s(1,1) += 1.0;
      s(2,2) += 1.0;
      v(0) = s(0,0)>=0.0 ? s(0,0) : -s(0,0);
      v(1) = s(1,0)>=0.0 ? s(1,1) : -s(1,1);
      v(2) = s(2,0)>=0.0 ? s(2,2) : -s(2,2);
      v *= t / mag(v);
    }
  }else{
    register double adj = t/sinT2;
    v[0] = (m(2,1)-m(1,2))*adj; 
    v[1] = (m(0,2)-m(2,0))*adj; 
    v[2] = (m(1,0)-m(0,1))*adj;
  }
}

//-------------------------------------------------------------------

void SlVelMatrixToMoment(const SlMatrix3x3 &m, SlVector3 &v) {
  // Set m so that vxb == m*b 
  v[0] = (m(2,1)-m(1,2))/2.0;
  v[1] = (m(0,2)-m(2,0))/2.0;  
  v[2] = (m(1,0)-m(0,1))/2.0;  
}


//-------------------------------------------------------------------
//-------------------------------------------------------------------

void SlEulerAngToMatrixXYZ(const SlVector3 &ea, SlMatrix3x3 &m) {
  
  // Adapted from SDFast generated routine sdang2dc

  double cos1 = cos(ea[0]);
  double cos2 = cos(ea[1]);
  double cos3 = cos(ea[2]);
  double sin1 = sin(ea[0]);
  double sin2 = sin(ea[1]);
  double sin3 = sin(ea[2]);
  m(0,0) = (cos2*cos3);
  m(0,1) = -(cos2*sin3);
  m(0,2) = sin2;
  m(1,0) = ((cos1*sin3)+(sin1*(cos3*sin2)));
  m(1,1) = ((cos1*cos3)-(sin1*(sin2*sin3)));
  m(1,2) = -(cos2*sin1);
  m(2,0) = ((sin1*sin3)-(cos1*(cos3*sin2)));
  m(2,1) = ((cos1*(sin2*sin3))+(cos3*sin1));
  m(2,2) = (cos1*cos2);
}
  
//-------------------------------------------------------------------

void SlMatrixToEulerAngXYZ(const SlMatrix3x3 &m, SlVector3 &ea) {

  // Adapted from SDFast generated routine sddc2eng

  double quot,angle,th1,th2,th3,costh2,temp[10];
  
  if ((fabs((fabs(m(0,2))-1.)) <= 1e-10)  ) {
    if ((m(0,2) > 0.)  ) {
      temp[0] = 1.5707963267949;
    } else {
      temp[0] = -1.5707963267949;
    }
    th2 = temp[0];
    if ((m(2,1) > 1.)  ) {
      temp[0] = 1.;
    } else {
      if ((m(2,1) < -1.)  ) {
	temp[1] = -1.;
      } else {
	temp[1] = m(2,1);
      }
      temp[0] = temp[1];
    }
    angle = asin(temp[0]);
    if ((m(1,1) >= 0.)  ) {
      temp[0] = angle;
    } else {
      temp[0] = (3.14159265358979-angle);
    }
    th1 = temp[0];
    th3 = 0.;
  } else {
    if ((m(0,2) > 1.)  ) {
      temp[0] = 1.;
    } else {
      if ((m(0,2) < -1.)  ) {
	temp[1] = -1.;
      } else {
	temp[1] = m(0,2);
      }
      temp[0] = temp[1];
    }
    th2 = asin(temp[0]);
    costh2 = cos(th2);
    quot = ((-m(1,2))/costh2);
    if ((quot > 1.)  ) {
      temp[0] = 1.;
    } else {
      if ((quot < -1.)  ) {
	temp[1] = -1.;
      } else {
	temp[1] = quot;
      }
      temp[0] = temp[1];
    }
    angle = asin(temp[0]);
    if ((m(2,2) >= 0.)  ) {
      temp[0] = angle;
    } else {
      temp[0] = (3.14159265358979-angle);
    }
    th1 = temp[0];
    quot = ((-m(0,1))/costh2);
    if ((quot > 1.)  ) {
      temp[0] = 1.;
    } else {
      if ((quot < -1.)  ) {
	temp[1] = -1.;
      } else {
	temp[1] = quot;
      }
      temp[0] = temp[1];
    }
    angle = asin(temp[0]);
    if ((m(0,0) >= 0.)  ) {
      temp[0] = angle;
    } else {
      temp[0] = (3.14159265358979-angle);
    }
    th3 = temp[0];
  }
  if ((th1 > 3.14159265358979)  ) {
    temp[0] = (th1-6.28318530717959);
  } else {
    temp[0] = th1;
  }
  ea[0] = temp[0];
  ea[1] = th2;
  if ((th3 > 3.14159265358979)  ) {
    temp[0] = (th3-6.28318530717959);
  } else {
    temp[0] = th3;
  }
  ea[2] = temp[0];
}
  
//-------------------------------------------------------------------
//-------------------------------------------------------------------

void SlQuaternionToMatrix(const SlVector3 &xyz, double w, SlMatrix3x3 &m) {

  // Adapted from SDFast generated routine sdquat2dc

  double e1,e2,e3,e4,e11,e22,e33,e44,norm;
  
  e11 = xyz[0]*xyz[0];
  e22 = xyz[1]*xyz[1];
  e33 = xyz[2]*xyz[2];
  e44 = w   *w   ;
  norm = sqrt(e11+e22+e33+e44);
  if (norm == 0.) {
    e4 = 1.;
    norm = 1.;
  } else {
    e4 = w;
  }
  norm = 1./norm;
  e1 = xyz[0]*norm;
  e2 = xyz[1]*norm;
  e3 = xyz[2]*norm;
  e4 = w     *norm;
  e11 = e1*e1;
  e22 = e2*e2;
  e33 = e3*e3;
  m(0,0) = 1.-(2.*(e22+e33));
  m(0,1) = 2.*(e1*e2-e3*e4);
  m(0,2) = 2.*(e1*e3+e2*e4);
  m(1,0) = 2.*(e1*e2+e3*e4);
  m(1,1) = 1.-(2.*(e11+e33));
  m(1,2) = 2.*(e2*e3-e1*e4);
  m(2,0) = 2.*(e1*e3-e2*e4);
  m(2,1) = 2.*(e2*e3+e1*e4);
  m(2,2) = 1.-(2.*(e11+e22));
}
  
//-------------------------------------------------------------------

void SlMatrixToQuaternion(const SlMatrix3x3 &m, SlVector3 &xyz, double &w) {

  // Adapted from SDFast generated routine sddc2quat  

  double tmp,tmp1,tmp2,tmp3,tmp4,temp[10];
  
  tmp = (.25*(1.-(m(0,0)+(m(1,1)+m(2,2)))));
  tmp4 = (.5-tmp);
  if ((tmp4 <= 0.)  ) {
    temp[0] = 0.;
  } else {
    temp[0] = sqrt(tmp4);
  }
  tmp4 = temp[0];
  tmp1 = (tmp+(.5*m(0,0)));
  if ((tmp1 <= 0.)  ) {
    temp[0] = 0.;
  } else {
    temp[0] = sqrt(tmp1);
  }
  tmp1 = temp[0];
  tmp2 = (tmp+(.5*m(1,1)));
  if ((tmp2 <= 0.)  ) {
    temp[0] = 0.;
  } else {
    temp[0] = sqrt(tmp2);
  }
  tmp2 = temp[0];
  tmp3 = (tmp+(.5*m(2,2)));
  if ((tmp3 <= 0.)  ) {
    temp[0] = 0.;
  } else {
    temp[0] = sqrt(tmp3);
  }
  tmp3 = temp[0];
  if (((tmp1 >= tmp2) && (tmp1 >= tmp3))  ) {
    if ((m(2,1) < m(1,2))  ) {
      temp[0] = -1.;
    } else {
      temp[0] = 1.;
    }
    tmp1 = (tmp1*temp[0]);
    if (((tmp1*(m(0,1)+m(1,0))) < 0.)  ) {
      temp[0] = -1.;
    } else {
      temp[0] = 1.;
    }
    tmp2 = (tmp2*temp[0]);
    if (((tmp1*(m(0,2)+m(2,0))) < 0.)  ) {
      temp[0] = -1.;
    } else {
      temp[0] = 1.;
    }
    tmp3 = (tmp3*temp[0]);
  } else {
    if (((tmp2 >= tmp1) && (tmp2 >= tmp3))  ) {
      if ((m(0,2) < m(2,0))  ) {
				temp[0] = -1.;
      } else {
				temp[0] = 1.;
      }
      tmp2 = (tmp2*temp[0]);
      if (((tmp2*(m(1,2)+m(2,1))) < 0.)  ) {
				temp[0] = -1.;
      } else {
				temp[0] = 1.;
      }
      tmp3 = (tmp3*temp[0]);
      if (((tmp2*(m(0,1)+m(1,0))) < 0.)  ) {
				temp[0] = -1.;
      } else {
				temp[0] = 1.;
      }
      tmp1 = (tmp1*temp[0]);
    } else {
      if ((m(1,0) < m(0,1))  ) {
				temp[0] = -1.;
      } else {
				temp[0] = 1.;
      }
      tmp3 = (tmp3*temp[0]);
      if (((tmp3*(m(0,2)+m(2,0))) < 0.)  ) {
				temp[0] = -1.;
      } else {
				temp[0] = 1.;
      }
      tmp1 = (tmp1*temp[0]);
      if (((tmp3*(m(1,2)+m(2,1))) < 0.)  ) {
				temp[0] = -1.;
      } else {
				temp[0] = 1.;
      }
      tmp2 = (tmp2*temp[0]);
    }
  }
  tmp = (1./sqrt(((tmp1*tmp1)+((tmp2*tmp2)+((tmp3*tmp3)+(tmp4*tmp4))))));
  xyz[0] = (tmp*tmp1);
  xyz[1] = (tmp*tmp2);
  xyz[2] = (tmp*tmp3);
  w      = (tmp*tmp4);
}
  
//-------------------------------------------------------------------

void SlQuaternionToMoment(const SlVector3 &xyz, double w, SlVector3 &m) {
  double l = sqrt(dot(xyz,xyz) + w * w);
  double g = sqrt(dot(xyz,xyz)        );
  if ((l==0.0) || (g==0.0)) {
    m = 0.0;
  }else{
    double theta = (2.0 * acos(w/l));
    if (theta > M_PI) theta =  theta - 2.0 * M_PI;
    m = xyz * theta / g;
  }
}
    
void SlMomentToQuaternion(const SlVector3 &m, SlVector3 &xyz, double &w) {
  double theta = mag(m);
  if (theta == 0.0) {
    xyz = 0.0;
    w   = 1.0;
  }else{
    xyz = m * ( sin(theta*0.5) / theta );
    w   = cos(theta*0.5);
  }
}

//-------------------------------------------------------------------
//-------------------------------------------------------------------
//-------------------------------------------------------------------


void SlSlerpQuaternion(SlVector3  xyz0, double  w0,
		       SlVector3  xyz1, double  w1, 
		       SlVector3 &xyzT, double &wT, double t) {

  //void slerp(SlVector3  xyz0, double  w0,
  //	   SlVector3  xyz1, double  w1, 
  //	   SlVector3 &xyzT, double &wT, double t) {
  
  double omega,cosom,sinom,sclp,sclq, pminusq, pplusq;

  static const double EPS = ( 1.0e-10 );
  static const double HPI = ( M_PI / 2.0);

  pminusq = sqrMag(xyz0-xyz1)+(w0-w1)*(w0-w1);
  pplusq  = sqrMag(xyz0+xyz1)+(w0+w1)*(w0+w1); 

  if (pplusq < pminusq) {
	xyz1 = -xyz1;
	w1   = -w1;
  }
    
  cosom = dot(xyz0,xyz1) + w0*w1;
  
  if ((1.0+cosom) > EPS) { 
	if ( (1.0-cosom) > EPS) { 
	  omega = acos(cosom); 
	  sinom =  sin(omega); 
	  sclp  =  sin( ( 1.0 - t ) * omega ) /sinom; 
	  sclq  =  sin(         t   * omega ) /sinom; 
	}else{ 
	  sclp  = 1.0 - t; 
	  sclq  =       t; 
	} 
  
	xyzT = sclp*xyz0 + sclq*xyz1; 
	wT = sclp*w0 + sclq*w1; 
  }else{ 
    
	xyzT[0] = -xyz0[1];   
	xyzT[1] =  xyz0[0]; 
	
	xyzT[2] = -w0   ;   
	wT    =  xyz0[2]; 
	
	sclp = sin((1.0-t)*HPI); 
	sclq = sin(t*HPI); 
	
	xyzT = sclp*xyz0 + sclq*xyzT; 
	wT   = sclp*w0   + sclq*wT; 
  } 
  double m = sqrt(dot(xyzT,xyzT) + wT*wT);
  xyzT /= m;
  wT /= m;
} 


//-------------------------------------------------------------------
//-------------------------------------------------------------------
//-------------------------------------------------------------------
//-------------------------------------------------------------------






