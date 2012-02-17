#ifndef mathutils_h
#define mathutils_h 
#include <cmath>
#include <cstdlib> 
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

template <class T>
T deltaPhi (T phi1, T phi2) {
  T result = phi1 - phi2;
  while (result > M_PI) result -= 2*M_PI;
  while (result <= -M_PI) result += 2*M_PI;
  return result;
}

template <class T>
T deltaR2 (T eta1, T phi1, T eta2, T phi2) {
  T deta = eta1 - eta2;
  T dphi = deltaPhi (phi1, phi2);
  return deta*deta + dphi*dphi;
}

template <class T>
T deltaR (T eta1, T phi1, T eta2, T phi2) {
  return sqrt (deltaR2 (eta1, phi1, eta2, phi2));
}

template<typename T1, typename T2>
double deltaR2( const T1 & t1, const T2 & t2 ) {
  return deltaR2( t1.eta(), t1.phi(), t2.eta(), t2.phi() );
}

template<typename T1, typename T2>
double deltaR( const T1 & t1, const T2 & t2 ) {
  return deltaR( t1.eta(), t1.phi(), t2.eta(), t2.phi() );
}

template <class T>
double invm2 ( T e1,T e2,T pt1,T pt2,T phi1,T phi2,T pz1, T pz2 ) {
  return (2*((e1*e2)-(pt1*pt2*cos(phi1-phi2))-(pz1*pz2)) ) ; 
}

template <class T>
double invm ( T e1,T e2,T pt1,T pt2,T phi1,T phi2,T pz1, T pz2 ) {
  return sqrt(invm2(e1,e2,pt1,pt2,phi1,phi2,pz1,pz2)) ; 
}

template <class T>
double mT2 ( T pt1, T pt2, T phi1, T phi2 ) {
  return ( 2*pt1*pt2*(1 - cos(deltaPhi(phi1,phi2))) ) ; 
}

template <class T>
double mT ( T pt1, T pt2, T phi1, T phi2 ) {
  return sqrt(mT2(pt1,pt2,phi1,phi2)) ; 
}

template <class T>
inline T pow2(T& x) {return x*x;} 
  
#include "rand.h"

#endif
