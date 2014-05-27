void StressTestGSVD()
{
  // Quotient SVD of A(m,n), B(p,n) where m >= n >= p
  
  gROOT->LoadMacro("UnfoldingUtils.C+");
  UnfoldingUtils uu;
  TRandom3 ran;

  ran.SetSeed(0);
  
  int t=1, ntests=5;
  while (t<=ntests) {

    int p = 2 + ran.Integer(5);
    int n = p + ran.Integer(5);
    int m = n + ran.Integer(5);

    m=    n = p;

    Printf("\n\nGSVD Test %d: m=%d, n=%d, p=%d.", t,m,n,p);

    TMatrixD A(m,n);
    TMatrixD B(p,n);
    for (int i=0; i<m; i++)
      for (int j=0; j<n; j++)
	A(i,j) = ran.Integer(9);
    for (int i=0; i<p; i++)
      for (int j=0; j<n; j++)
	B(i,j) = ran.Integer(9);
  
    GSVDecompResult g = uu.GSVD(A,B);

    Printf("U (%d x %d), V (%d x %d), C (%d x %d), S (%d x %d), X' (%d x %d)",
	   g.U.GetNrows(),  g.U.GetNcols(),
	   g.V.GetNrows(),  g.V.GetNcols(),
	   g.C.GetNrows(),  g.C.GetNcols(),
	   g.S.GetNrows(),  g.S.GetNcols(),
	   g.XT.GetNrows(), g.XT.GetNcols());
	   
    TMatrixD upper = g.U*g.C*g.XT;
    TMatrixD lower = g.V*g.S*g.XT;

    // // Printing -----------------------------------
    // A.Print();
    // B.Print();

    g.C.Print();
    g.S.Print();

    // cout << "U*C*XT: ";  upper.Print();
    // cout << "V*S*XT: ";  lower.Print();
    // // Printing -----------------------------------

    TMatrixD M1(m+p,n);
    M1.SetSub(0,0,A);
    M1.SetSub(m,0,B);

    TMatrixD M2(M1);
    M2.SetSub(0,0,upper);
    M2.SetSub(upper.GetNrows(),0,lower);
    
    M1 -= M2;

    double sum = M1.Sum()/M1.GetNoElements();
    if (sum < 1e-13) 
      sum = 0.;
    // cout << "A: "; A.Print();
    // cout << "B: "; B.Print();

    // cout << "U*C*XT: ";  upper.Print();
    // cout << "V*S*XT: ";  lower.Print();

    // cout << "M1: "; M1.Print();
    // cout << "M2: "; M2.Print();

    Printf("M1-M2 = %g", sum);
    if (sum > 0) 
      M1.Print();
    
    t++;
  }

}
