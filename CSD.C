#include "TMatrixD.h"
#include "TVectorD.h"
#include "TMath.h"
#include "TDecompSVD.h"
#include "TDecompQRH.h"
#include <iostream>
using namespace std;

void ReverseColumns(TMatrixD& A);
void ReverseRows(TMatrixD& A);
void ReverseVector(TVectorD& v);
void SwapColumns(TMatrixD &A, int col1, int col2);
void SwapElements(TVectorD& v, int j1, int j2);

struct House
{
  int n;
  double beta;
  TVectorD v;
};

struct QRDecomp
{
  TMatrixD Q;
  TMatrixD R;
};
QRDecomp QR(TMatrixD& A);

void CSD()
{
  // Decompose Q = [Q1; Q2] into Q1 = UCZ' and Q2 = VSZ'
  // Q (m+p x l, m+p >= l) is partitioned as Q1(mxl), Q2(pxl).
  // Dimensions: U(mxm), V(pxp), Z(lxl), C(mxl), S(pxl).

  int m,p,l,q1,q2,r,n;

  // Example Q matrix
  TMatrixD Q(7,4);
  for (int i=0; i<7; i++) Q(i,0) = 1./TMath::Sqrt(7.);
  for (int i=0; i<7; i++) Q(i,1) = 1./TMath::Sqrt(10.);
  for (int i=0; i<7; i++) Q(i,2) = 1./4;
  for (int i=0; i<7; i++) Q(i,3) = 1./TMath::Sqrt(3.);
  Q(0,1) *=  0; Q(0,2) *=  0;
  Q(1,1) *= -2; Q(1,2) *= -2; Q(1,3) *= -1./2;
  Q(2,1) *= -1; Q(2,2) *=  1; Q(2,3) *=  3./4;
  Q(3,1) *=  0; Q(3,2) *=  3; Q(3,3) *= -3./4;
  Q(4,1) *=  0; Q(4,2) *=  0; Q(4,3) *=  0.;
  Q(5,1) *=  1; Q(5,2) *= -1; Q(5,3) *= -3./4;
  Q(6,1) *=  2; Q(6,2) *= -1; Q(6,3) *=  1./4;
  Q.Print();

  // ***************************************
  // m <= p case
  // ***************************************
  m=3; p=4; l=4;
  q1 = TMath::Min(m,l); // 3
  q2 = TMath::Min(p,l); // 4
  TMatrixD C(m,l);
  TMatrixD S(p,l);

  // Partition Q
  TMatrixD Q1 = Q.GetSub(0, m-1,   0, l-1); // m x l
  TMatrixD Q2 = Q.GetSub(m, m+p-1, 0, l-1); // p x l
  cout << "Q1: ";  Q1.Print();
  cout << "Q2: ";  Q2.Print();

  // Conventional SVD of Q2: VSZ'
  TDecompSVD svdQ2(Q2);
  TMatrixD V    = svdQ2.GetU();    // p x p
  TVectorD beta = svdQ2.GetSig();  // l
  TMatrixD Z    = svdQ2.GetV();    // l x l

  // cout << "V: ";        V.Print();
  // cout << "beta: ";  beta.Print();
  // cout << "Z: ";        Z.Print();

  // Re-order V, S, Z cols
  ReverseColumns(V);
  ReverseVector(beta);
  ReverseColumns(Z);

  for (int i=0; i<q2; i++) {
    S(i, l-q2+i) = beta(l-q2+i);
  }

  cout << "V: ";  V.Print();
  cout << "S: ";  S.Print();
  cout << "Z: ";  Z.Print();

  // Find r: beta(r) <= 1/sqrt(2) < beta(r+1) 
  r=0;
  for (int i=0; i<m-1; i++) {
    if (beta(i)  <= 1./TMath::Sqrt(2.) && 
	beta(i+1) > 1./TMath::Sqrt(2.)) {
      r = i;
      break;
    }
  }
  cout << "r: " << r << endl;

// Compute T = Q1(mxl) * Z(lxl) --> mxl
TMatrixD T = Q1*Z;

// Pad with extra rows to enable QR decomp
 if (m < l)
   T.ResizeTo(l,l);

cout << "T: ";  T.Print();

 // Check:
 TMatrixD TTT(T, TMatrixD::kTransposeMult, T);
 TMatrixD STS(S, TMatrixD::kTransposeMult, S);
 TMatrixD TS(TTT);
 TS+=STS;
 cout << "T\'T + S\'S: ";  TS.Print();

 // QR decomp of T: T = UR
 TMatrixD Tcopy(T);
 QRDecomp qrT = QR(Tcopy);
 TMatrixD U = qrT.Q;
 TMatrixD R = qrT.R;

 cout << "U and R: ";
 U.Print();
 R.Print();
 TMatrixD UTU(U, TMatrixD::kTransposeMult, U);
 UTU.Print();

 // TDecompQRH qrt(T);
 // TMatrixD U = qrt.GetQ();
 // TMatrixD R = qrt.GetR();
 // TVectorD hBeta = qrt.GetW();
 // TVectorD hV = qrt.GetUp();

// cout << "U: ";  U.Print();
// cout << "R: ";  R.Print();
// cout << "Householder beta: ";  hBeta.Print();

//  TMatrixD QQ(U); QQ.UnitMatrix();

//  for (int j=0; j<l; j++) {
//    TVectorD v(l);// = TMatrixDColumn(U, j);
//    for (int i=0; i<l; i++) {
//      if (i<j)  v(i) = 0;
//      if (i==j) v(i) = 1.;
//      if (i>j)  v(i) = U(i-1, j);
//    }
//    TMatrixD H = OuterProduct(v, v);
//    H *= hBeta(j);
//    TMatrixD I(H); I.UnitMatrix();
//    H = I - H;
//    QQ *= H;

//    Printf("v %d", j); v.Print();
//  }
//  QQ.Print();
//  TMatrixD QTQ(QQ, TMatrixD::kTransposeMult, QQ);
//  QTQ.Print();



// Get the sub-matrix from R: [R33 R34]
 int r1,r2,c1,c2;
 r1=l-q2;
 r2=r-l+q2;
 c1=r1;
 c2=r2;
 TMatrixD R22 = R.GetSub(r1,r2,c1,c2);
 cout << "R22: ";  R22.Print();

 r1=r2+1;
 r2=r1+q1-r-2;
 c1=c2+1;
 c2=r2;
 TMatrixD R33 = R.GetSub(r1,r2,c1,c2);
 cout << "R33: ";  R33.Print();

 c1=c2+1;
 c2=R.GetNcols()-1;
 TMatrixD R34 = R.GetSub(r1,r2,c1,c2);
 cout << "R34: ";  R34.Print();

 int r33 = R33.GetNrows(), c33 = R33.GetNcols();
 int c34 = R34.GetNcols();

 TMatrixD R3334(r33, c33+c34);
 if (r33 < c33+c34)
   R3334.ResizeTo(c33+c34, c33+c34);
 
 TMatrixDSub(R3334, 0, r33-1, 0,   c33-1) += R33;
 TMatrixDSub(R3334, 0, r33-1, c33, c33+c34-1) += R34;
 cout << "[R33 R34]: ";  R3334.Print();

 // Compute SVD of [R33 R34]
 TDecompSVD svd2(R3334);
 TMatrixD Ur = svd2.GetU();
 TMatrixD Zr = svd2.GetV();
 TVectorD alpha = svd2.GetSig();

 cout << "Ur: ";  Ur.Print();
 cout << "Zr: ";  Zr.Print();
 cout << "alpha: ";  alpha.Print();

 // Assign C
 // int t = q1 + q2 - l; // 3
 for (int i=0; i<q1; i++) {
   if (i<l-q2) C(i,i) = 1.;
   else if (i<alpha.GetNrows()) C(i, i) = alpha(i);
 }

 cout << "C: ";  C.Print();
 cout << "S: ";  S.Print();
 //#if(0)#endif

 return;
}

void ReverseColumns(TMatrixD& A)
{
  TMatrixD B(A);
  int n = B.GetNcols();

  for (int j=0; j<n; j++) {
    TMatrixDColumn(B, j) = TMatrixDColumn(A, n-j-1);
  }
  A = B;
}

void ReverseRows(TMatrixD& A)
{
  TMatrixD B(A);
  int m = B.GetNrows();

  for (int i=0; i<m; i++) {
    TMatrixDColumn(B, i) = TMatrixDColumn(A, m-i-1);
  }
  A = B;
}

void ReverseVector(TVectorD& v)
{
  TVectorD vtmp(v);
  int m = v.GetNrows();
  for (int i=0; i<m; i++) {
    vtmp(i) = v(m-i-1);
  }
  v = vtmp;
}

TMatrixD OuterProduct(TVectorD a, TVectorD b)
{
  // Return a*b'
  int na = a.GetNrows();
  int nb = b.GetNrows();
  TMatrixD M(na,nb);
  for (int i=0; i<na; i++) {
    for (int j=0; j<nb; j++) {
      M(i,j) = a(i)*b(j);
    }
  }
  return M;
}


House Householder(TVectorD& x, int m)
{
  // Golub & Van Loan, 3rd ed, Alg. 5.1.1
  House h;
  int n = x.GetNrows();
  TVectorD x1 = x.GetSub(1,n-1);
  double sigma = x1*x1; 
  TVectorD v(x);
  v(0) = 1.;
  double beta = 0;
  if (sigma != 0) {
    double mu = TMath::Sqrt(x(0)*x(0) + sigma);
    if (x(0) <= 0)
      v(0) = x(0) - mu;
    else
      v(0) = -sigma/(x(0) + mu);

    //    beta = 2./TMath::Sqrt(v*v);
    beta = 2*v(0)*v(0) / (sigma + v(0)*v(0));
    v *= 1./v(0);
  }

  h.n = n;
  h.beta = beta;
  h.v.ResizeTo(m);
  for (int i=0; i<n; i++)
    h.v(i+m-n) = v(i);
  return h;
}


QRDecomp QR(TMatrixD& A)
{
  
  QRDecomp qr;
  int m = A.GetNrows();
  int n = A.GetNcols();
  TMatrixD Q(m,m); Q.UnitMatrix();
  TMatrixD R(A);
  TMatrixD Qj(m,m);

  int nIter = TMath::Min(m-1, n);
  for (int j=0; j<nIter; j++) {
    TVectorD col = TMatrixDColumn(R,j);
    TVectorD x = col.GetSub(j,m-1);
    int sign = (x(0)<0.)? 1. : -1.;
    double alpha = sign*TMath::Sqrt(x*x);
    TVectorD u(x);
    u(0) += alpha;
    
    // Compute Householder vector v and matrix H
    double unorm = TMath::Sqrt(u*u);
    TVectorD v(u); v *= (unorm==0)? 0. : 1./unorm;  
    TMatrixD H = OuterProduct(v,v);
    H *= 2;
    TMatrixD I(H); I.UnitMatrix();
    H = I - H;
    
    Qj.UnitMatrix();
    Qj.SetSub(j,j,H);
    
    Q = Q*Qj;    
    R = Qj*R;
    
  }








  /*
  // G & VL,  ex. 5.4.1, p250
  QRDecomp qr;
  const int m = A.GetNrows();
  const int n = A.GetNcols();
  TMatrixD Q(m,m);
  Q.UnitMatrix();
  TMatrixD R(m,m);

  int r=0; // Column rank of A
  int k=0; // Permutation index
  double tau=0;
  int piv[n];
  TVectorD c(n); // vector of squared column norms

  for (int j=0; j<n; j++) {
    TVectorD aj= TMatrixDColumn(A, j);
    c(j) = aj*aj;
  }
  tau = c.Max();

  // Find smallest column index k where c(k) = tau
  for (int j=0; j<n; j++) {
    if (c(k) < tau)
      k++;
    else
      break;
  }

  while (tau > 0.) {
    piv[r] = k;
    SwapColumns(A, r, k);
    SwapElements(c, r, k);
    TVectorD xfull = TMatrixDColumn(A,r); 
    TVectorD x = xfull.GetSub(r, xfull.GetNrows()-1);
    House h = Householder(x, m);
    TMatrixD H = OuterProduct(h.v, h.v); H *= h.beta;
    TMatrixD I(H); I.UnitMatrix();
    H = I - H;
    //    TMatrixDSub(A, r, m-1, r, n-1) = H*A.GetSub(r, m-1, r, n-1);
    TMatrixDSub(A, r, m-1, r, n-1) = 
      H.GetSub(r, m-1, r, m-1)*A.GetSub(r, m-1, r, n-1);

    Q *= H;

    TMatrixDSub(R, r, m-1, r, n-1) = A.GetSub(r, m-1, r, n-1);

    // Printf("H_%d (%d x %d)", r, H.GetNrows(), H.GetNcols());
    // H.Print();

    for (int i=0; i<m-r-1; i++) 
      A(i+r+1,r) = h.v(i+1);

    for (int i=r+1; i<n; i++) 
      c(i) -= A(r,i)*A(r,i);

    if (r < n-1) {
    // Set tau as max { c(r+1)..c(n-1) }
    tau = 0.;
    for (int j=r+1; j<n; j++) {
      if (c(j) > tau)       
	tau = c(j);
    }
    // Find smallest column index k where c(k) = tau
    for (int j=r+1; j<n; j++) {
      if (c(k) < tau)
	k++;
      else
	break;
    }
    }
    else
      tau = 0;
    r++;

    //    Printf("A"); A.Print();
  }
    // Q.Print();
    // R.Print();

    // TMatrixD QR = Q*R;
    // QR.Print();
    // for (int i=0; i<3; i++)
    //   Printf("piv[%d]= %d", i, piv[i]);

    // TMatrixD QTQ(Q, TMatrixD::kTransposeMult, Q);
    // QTQ.Print();
    */

  qr.Q.ResizeTo(Q);
  qr.R.ResizeTo(R);
  qr.Q = Q;
  qr.R = R;
  
  return qr;
}


void SwapColumns(TMatrixD &A, int col1, int col2)
{
  int nc = A.GetNcols();
  if (col1 >= nc || col2 >= nc)
    Error("SwapColumns", "col1 or col2 index out of bounds");

  TMatrixD B(A);
  
  TMatrixDColumn(B, col1) = TMatrixDColumn(A, col2);
  TMatrixDColumn(B, col2) = TMatrixDColumn(A, col1);

  A = B;
}

void SwapElements(TVectorD& v, int j1, int j2)
{
  int nr = v.GetNrows();
  if (j1 >= nr || j2 >= nr)
    Error("SwapElements", "an index is out of bounds");

  TVectorD v2(v);

  v2(j1) = v(j2);
  v2(j2) = v(j1);

  v = v2;
}
