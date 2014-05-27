void TestGSVD()
{
  // Quotient SVD of A(m,n), B(p,n) where m >= n >= p
  
  gROOT->LoadMacro("UnfoldingUtils.C+");
  UnfoldingUtils uu;
  const int m=5;
  const int n=4;
  const int p=3;
  int r;
  TMatrixD A(m,n);
  TMatrixD B(p,n);
  double a[m][n] = { {1,2,1,0},
		     {2,3,1,1},
		     {3,4,1,2},
		     {4,5,1,3},
		     {5,6,1,4}};
  double b[p][n] = { {6,7,1,5},
		     {7,1,-6,13},
		     {-4,8,9,-2}};
  for (int i=0; i<m; i++)
    for (int j=0; j<n; j++)
      A(i,j) = a[i][j];
  for (int i=0; i<p; i++)
    for (int j=0; j<n; j++)
      B(i,j) = b[i][j];
  
  GSVDecompResult g = uu.GSVD(A,B);

  TMatrixD upper = g.U*g.C*g.XT;
  TMatrixD lower = g.V*g.S*g.XT;

  cout << "U1*C*XT: ";  upper.Print();
  cout << "U2*S*XT: ";  lower.Print();

  cout << "C: "; g.C.Print();
  cout << "S: "; g.S.Print();

}
