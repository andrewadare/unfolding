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

struct QRDecompResult
{
  TMatrixD Q;
  TMatrixD R;
};

struct CSDecompResult
{
  TMatrixD C;
  TMatrixD S;
  TMatrixD U;
  TMatrixD V;
  TMatrixD Z;
};

QRDecompResult QRDecomp(TMatrixD& A);
QRDecompResult QL(TMatrixD& A);
CSDecompResult CSDecompQ1Taller(TMatrixD& Q1, TMatrixD& Q2);
CSDecompResult CSDecomp(TMatrixD& Q1, TMatrixD& Q2);

void CSD()
{
  // Decompose Q = [Q1; Q2] into Q1 = UCZ' and Q2 = VSZ'
  // Q (m+p x l, m+p >= l) is partitioned as Q1(mxl), Q2(pxl).
  // Dimensions: U(mxm), V(pxp), Z(lxl), C(mxl), S(pxl).

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

  // Partition Q
  int l=4;
    int m=3; int p=4; 
    //int m=5; int p=2;
  TMatrixD Q1 = Q.GetSub(0, m-1,   0, l-1); // m x l
  TMatrixD Q2 = Q.GetSub(m, m+p-1, 0, l-1); // p x l
  cout << "Q1: ";  Q1.Print();
  cout << "Q2: ";  Q2.Print();


  if (m > p) {
    CSDecompResult cs = CSDecompQ1Taller(Q1, Q2);
    cout << "+++++++++++++++ m > p result ++++++++++++++++: "<<endl;
    cout << "C: ";  cs.C.Print();
    cout << "S: ";  cs.S.Print();
    cout << "U: ";  cs.U.Print();
    cout << "V: ";  cs.V.Print();
    cout << "Z: ";  cs.Z.Print();
    TMatrixD Q1test(cs.C, TMatrixD::kMultTranspose, cs.Z);
    TMatrixD Q2test(cs.S, TMatrixD::kMultTranspose, cs.Z);
    Q1test = Q1 - cs.U*Q1test;
    Q2test = Q2 - cs.V*Q2test;
    cout << "Q1 - UCZ': ";  Q1test.Print();
    cout << "Q2 - VSZ': ";  Q2test.Print();
  }
  else {
    CSDecompResult cs = CSDecomp(Q1, Q2);
    cout << "+++++++++++++++ m <= p result +++++++++++++++: "<<endl;
    cout << "C: ";  cs.C.Print();
    cout << "S: ";  cs.S.Print();
    cout << "U: ";  cs.U.Print();
    cout << "V: ";  cs.V.Print();
    cout << "Z: ";  cs.Z.Print();
    TMatrixD Q1test(cs.C, TMatrixD::kMultTranspose, cs.Z);
    TMatrixD Q2test(cs.S, TMatrixD::kMultTranspose, cs.Z);
    Q1test = Q1 - cs.U*Q1test;
    Q2test = Q2 - cs.V*Q2test;
    cout << "Q1 - UCZ': ";  Q1test.Print();
    cout << "Q2 - VSZ': ";  Q2test.Print();
  }
 
  return;
}

CSDecompResult 
CSDecompQ1Taller(TMatrixD& Q1, TMatrixD& Q2)
{
  // ***************************************
  // m > p case
  // ***************************************
  int m,p,l,q1,q2,r;
  m  = Q1.GetNrows();
  p  = Q2.GetNrows();
  l  = Q1.GetNcols();
  q1 = TMath::Min(m,l);
  q2 = TMath::Min(p,l);
  r  = 0;
  TMatrixD C(m,l);
  TMatrixD S(p,l);
  CSDecompResult cs;

  // SVD of Q1: UCZ'
  TDecompSVD svdQ1(Q1);
  TMatrixD U     = svdQ1.GetU();    // m x m
  TVectorD alpha = svdQ1.GetSig();  // l
  TMatrixD Z     = svdQ1.GetV();    // l x l

  for (int j=0; j<l; j++) C(j,j) = alpha(j);

  cout << "U: ";  U.Print();
  cout << "alpha: ";  alpha.Print();
  cout << "Z: ";  Z.Print();

  // Find r where alpha(r) >= 1/sqrt(2) > alpha(r+1) 
  double thr = 1./TMath::Sqrt(2.);
  for (int i=0; i<l-1; i++) {
    if (alpha(i) >= thr && alpha(i+1) < thr) {
      r = i;
      break;
    }
  }
  cout << "r: " << r << endl;

  TMatrixD T = Q2*Z;
  cout << "T: ";  T.Print();

  // QL decomp of T: VL
  QRDecompResult vl = QL(T);

  TMatrixD V = vl.Q;   // p x p
  TMatrixD L = vl.R;

  cout << "V: ";  V.Print();
  cout << "L: ";  L.Print();

  // Create permutation matrix Pi; L = Pi*L.
  TMatrixD Pi(p,p);
  TMatrixD Iq2(q2,q2); Iq2.UnitMatrix();
  Pi.SetSub(0,p-q2,Iq2);
  if (p>q2) {
    int d = p-q2;
    TMatrixD Id(d,d); Id.UnitMatrix();
    Pi.SetSub(p-1,0,Id);
  }
  // And its inverse (for later)
  TMatrixD PiInv(TMatrixD::kInverted, Pi);
  cout << "Pi: ";  Pi.Print();
  cout << "PiInv: ";  PiInv.Print();

  L = Pi*L;
  cout << "Pi*L: ";  L.Print();

  TMatrixD L1 = L.GetSub(0,r-1,0,l-q2+r-1);

  if (l > q2) // So TDecompSVD works
    L1.ResizeTo(l-q2+r, l-q2+r);

  cout << "[L_11 L_12]: ";  L1.Print();

  TDecompSVD svdl(L1);
  TMatrixD Vlbig = svdl.GetU();
  TVectorD bl = svdl.GetSig();
  TMatrixD Zl = svdl.GetV(); 
  TMatrixD Vl = Vlbig.GetSub(0,r-1,0,r-1);

  ReverseColumns(Vl);
  ReverseColumns(Zl);
  ReverseVector(bl);

  for (int i=0; i<p; i++) 
    S(i,l-p+i) = bl(l-p+i);

  cout << "Vl: ";  Vl.Print();
  cout << "Vlbig: ";  Vlbig.Print();
  cout << "Zl: ";  Zl.Print();
  cout << "bl: ";  bl.Print();

  // If the dimensions don't work out, these Xl's may need to be
  // padded with unit block matrices on the diagonals.
  V = V*Vl*PiInv;
  Z = Z*Zl;

  TMatrixD W(r+l-q2, r+l-q2);
  for (int i=0; i<r+l-q2; i++) 
    W(i,i) = alpha(i);
  W = W*Zl;
  cout << "W: ";  W.Print();

  QRDecompResult qrw = QRDecomp(W);
  TMatrixD Qw = qrw.Q;


  int ndiff =  U.GetNcols() - Qw.GetNrows(); 

  if (ndiff > 0) {
    int nr = Qw.GetNrows();
    Qw.ResizeTo(nr+ndiff, Qw.GetNcols()+ndiff);
    for (int j=0; j<ndiff; j++)
      Qw(nr+j, nr+j) = 1.;
  }
  cout << "Qw: ";  Qw.Print();

  U = U*Qw;

  cs.C.ResizeTo(C);
  cs.S.ResizeTo(S);
  cs.U.ResizeTo(U);
  cs.V.ResizeTo(V);
  cs.Z.ResizeTo(Z);

  cs.C = C;
  cs.S = S;
  cs.U = U;
  cs.V = V;
  cs.Z = Z;

  return cs;
}

CSDecompResult 
CSDecomp(TMatrixD& Q1, TMatrixD& Q2)
{
  // ***************************************
  // m <= p case
  // ***************************************
  int m,p,l,q1,q2,r,n;
  m  = Q1.GetNrows();
  p  = Q2.GetNrows();
  l  = Q1.GetNcols();
  r  = 0;

  TMatrixD C(m,l);
  TMatrixD S(p,l);
  CSDecompResult csd;
  TVectorD alpha(l);
  TVectorD beta(l);

  // 1.
  q1 = TMath::Min(m,l); // 3
  q2 = TMath::Min(p,l); // 4

  // 2. SVD of Q2: Q2 = VSZ'
  TDecompSVD svdQ2(Q2);
  TMatrixD V    = svdQ2.GetU();    // p x p
  TMatrixD Z    = svdQ2.GetV();    // l x l
  beta = svdQ2.GetSig();           // l

  // 3-5. Re-order V, S, Z cols
  ReverseColumns(V);
  ReverseColumns(Z);
  ReverseVector(beta);

  // 6.
  for (int i=0; i<l-q2; i++) {
    alpha(i) = 1.0;
    beta(i)  = 0.0; 
  }

  // // Assign S.
  // // Q2 has l-q2 zero s.v.'s. Make them zero.
  // for (int i=0; i<q2; i++) {
  //   if (i<l-q2) beta(i) = 0.;
  //   S(i,i) = beta(i);
  // }
  
  cout << "V: ";  V.Print();
  cout << "S: ";  S.Print();
  cout << "Z: ";  Z.Print();
  cout << "beta:"; beta.Print(); // non-decreasing
  
  // 7.
  // Find r where beta(r) <= 1/sqrt(2) < beta(r+1).  C++ problem:
  // index != dimension! Need two variables, r and rdim, to resolve
  // the ambiguity when r=0. rdim is the number of betas below 0.707,
  // r is the index.
  double thr = 1./TMath::Sqrt(2.);
  int rdim = 0;
  r = -1;
  for (int i=0; i<m-1; i++) {
    if (beta(i) <= thr && beta(i+1) > thr) {
      r = i;
      break;
    }
  }
  if (r == -1) {
    r = 0;
  }
  else rdim = r+1;
  Printf("r = %d, rdim = %d.",r,rdim);
  
  // 8.
  TMatrixD T = Q1*Z; // (m x l)
  cout << "T: ";  T.Print();

  // 9.
  // QR decomp of T: T = UR
  QRDecompResult qrT = QRDecomp(T);
  TMatrixD U = qrT.Q;
  TMatrixD R = qrT.R;
  
  cout << "U and R: ";
  U.Print();
  R.Print();

  // Get R2 and R3 from R
  TMatrixD R2 = R.GetSub(l-q2,r,l-q2,r);
  TMatrixD R3 = R.GetSub(rdim,q1-1,rdim,l-1);
  int r3r = R3.GetNrows();
  int r3c = R3.GetNcols();
  if (r3r < r3c)
    R3.ResizeTo(r3c, r3c);

  cout << "R2: ";  R2.Print();
  cout << "R3: ";  R3.Print();

  // 10.
  // Compute SVD of R3: R3 = Ur*Cr*Zr'
  TDecompSVD svd2(R3);
  TMatrixD Ur = svd2.GetU();
  TMatrixD Zr = svd2.GetV();
  TVectorD a3 = svd2.GetSig();

  for (int i=0; i<a3.GetNrows(); i++)
    alpha(rdim+i) = a3(i);

  // 11.
  for (int i=q1; i<l; i++) {
    alpha(i) = 0.0;
    beta(i)  = 1.0;
  }

  // 12.
  for (int i=l-q2; i<rdim; i++) {
    alpha(i) = R2(i,i);
  }
  cout << "alpha: ";  alpha.Print();

  // 13.
  // Form final U matrix
  // First resize U to un-do TDecompSVD-required modification
  if (r3r < r3c) {
    Ur.ResizeTo(r3r,r3r);
  }
  cout << "Ur: ";  Ur.Print();
  TMatrixD Ur1(U);
  Ur1.UnitMatrix();
  Ur1.SetSub(rdim,rdim,Ur);
  cout << "Ur1: ";  Ur1.Print();

  // 14.
  // Form final Z matrix
  TMatrixD Zrt(Zr);
  Zrt.ResizeTo(q2-rdim,q2-rdim);

  TMatrixD Zr1(Z);
  Zr1.UnitMatrix();  
  Zr1.SetSub(rdim,rdim,Zr);
  cout << "Zr1: ";  Zr1.Print();

  /*
  // Fix sign?
  for (int i=0; i<Ur1.GetNcols(); i++)
    if (alpha(i)<0) {
      alpha(i) *= -1.0;
      Ur1(i,i) *= -1.0;
      Zr1(i,i) *= -1.0;
    }
  cout << "Ur1 after sign adjustments: ";  Ur1.Print();
  cout << "Zr1 after sign adjustments: ";  Zr1.Print();
  */

  U = U*Ur1;
  Z = Z*Zr1;

  // 15.
  TMatrixD St(Zrt);
  St.Zero();
  for (int i=0; i<St.GetNcols(); i++)
    St(i,i) = beta(i+rdim);

  TMatrixD W = St*Zrt;
  
  // 16.
  QRDecompResult qrw = QRDecomp(W);
  TMatrixD Qw = qrw.Q;
  n = TMath::Min(rdim,l-q2);

  // 17.
  TMatrixD Vpost(V);
  Vpost.UnitMatrix();
  Vpost.SetSub(rdim-n, rdim-n, Qw);
  V = V*Vpost;
  
  cout << "St: ";  St.Print();
  cout << "Zrt: ";  Zrt.Print();
  cout << "W: ";  W.Print();
  cout << "Qw: ";  Qw.Print();
  cout << "Vpost: ";  Vpost.Print();

  // Construct C from alpha
  for (int i=0; i<q1; i++)
    C(i,i) = alpha(i);

  // And S from beta
  for (int i=0; i<q2; i++)
    S(i,i) = beta(i);
  
  csd.C.ResizeTo(C);
  csd.S.ResizeTo(S);
  csd.U.ResizeTo(U);
  csd.V.ResizeTo(V);
  csd.Z.ResizeTo(Z);

  csd.C = C;
  csd.S = S;
  csd.U = U;
  csd.V = V;
  csd.Z = Z;
  Printf("m=%d, p=%d, l=%d, q1=%d, q2=%d, r=%d, n=%d",  m,p,l,q1,q2,r,n);
  return csd;
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
  // Not currently used.
  // Probably has a bug.
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

QRDecompResult QL(TMatrixD& A)
{
  // Compute QL decomposition of A using Householder transformations
  QRDecompResult ql;
  int m = A.GetNrows();
  int n = A.GetNcols();
  TMatrixD Q(m,m); Q.UnitMatrix();
  TMatrixD L(A);
  TMatrixD Qj(m,m);

  int nIter = TMath::Min(m-1, n);
  for (int j=0; j<nIter; j++) {
    TVectorD col = TMatrixDColumn(L,n-j-1);
    TVectorD x = col.GetSub(0,nIter-j);

    int sign = (x(nIter-j)<0.)? 1. : -1.;
    double alpha = sign*TMath::Sqrt(x*x);
    TVectorD u(x);
    u(nIter-j) += alpha;
    
    // Compute Householder vector v and matrix H
    double unorm = TMath::Sqrt(u*u);
    TVectorD v(u); v *= (unorm==0)? 0. : 1./unorm;  
    TMatrixD H = OuterProduct(v,v);
    H *= 2;
    TMatrixD I(H); I.UnitMatrix();
    H = I - H;

    // Full-dimension Householder matrix
    Qj.UnitMatrix();
    Qj.SetSub(0,0,H);
    
    // Update Q and L
    Q = Q*Qj;    
    L = Qj*L;
  }

  // Store results
  ql.Q.ResizeTo(Q);
  ql.R.ResizeTo(L);
  ql.Q = Q;
  ql.R = L;
  
  return ql;
}

QRDecompResult 
QRDecomp(TMatrixD& A)
{
  // Compute QR decomposition of A using Householder transformations
  QRDecompResult qr;
  int m = A.GetNrows();
  int n = A.GetNcols();
  TMatrixD Q(m,m); Q.UnitMatrix();
  TMatrixD R(A);
  TMatrixD Qj(m,m);

  int nIter = TMath::Min(m-1, n);
  for (int j=0; j<nIter; j++) {
    TVectorD col = TMatrixDColumn(R,j);
    TVectorD x = col.GetSub(j,m-1);
 
    int sign = (col(j)<0.)? -1. : 1.;
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
    
    // Full-dimension Householder matrix
    Qj.UnitMatrix();
    Qj.SetSub(j,j,H);

    // Update Q and R    
    Q = Q*Qj;    
    R = Qj*R;
    
  }

  // *** Extra feature ***
  TMatrixD U1(Q);
  U1.UnitMatrix();
  TMatrixD R1(R.GetNcols(),R.GetNcols());
  R1.UnitMatrix();

  for (int i=0; i<R.GetNrows(); i++) {
    if(R(i,i)<0) {
      U1(i,i) = -1;
      R1(i,i) = -1;
    }
  }

  Q = Q*U1;
  R = R*R1;

  // Store results
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
