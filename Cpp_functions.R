library(Rcpp)

cppFunction('NumericVector simAdp(double conc, NumericVector receptorPars, NumericVector LIFPars, 
NumericVector fitPars, NumericVector times, NumericVector initial, NumericVector valveStates) {
  double ki = receptorPars[0];
  double k1 = receptorPars[1];
  double km1 = receptorPars[2];
  double k2 = receptorPars[3];
  double km2 = receptorPars[4];
  double k3 = receptorPars[5];
  double km3 = receptorPars[6];
  double k4 = receptorPars[7];
  double Rtot = receptorPars[8];
  double Ntot = receptorPars[9];

  double Cm = LIFPars[0];
  double gL = LIFPars[1]; 
  double EL = LIFPars[2];
  double ER = LIFPars[3]; 
  double Vreset = LIFPars[4]; 
  double theta0 = LIFPars[5];

  double tau = fitPars[0];
  double Delta = fitPars[1];
  double gamma = fitPars[2];
  double nhill = fitPars[3];

  int n = times.size();
  NumericVector L(n);
  NumericVector R(n);
  NumericVector Rstar(n);
  NumericVector N(n);
  L[0] = 0;
  R[0] = Rtot-initial[0]-initial[1];
  Rstar[0] = initial[1];
  N[0] = Ntot;
  double deltat = times[1] - times[0];
  NumericVector voltage(n);
  voltage[0] = EL;
  NumericVector wth(n);
  wth[0] = 0;
  NumericVector spikestimes(n);
  
  for (int i = 1; i < n; ++i){

    // RECEPTOR DYNAMICS
    L[i] = L[i-1] + (ki*conc*valveStates[i-1] - k1*nhill*R[i-1]*pow(L[i-1], nhill) + 
      km1*nhill*(Rtot - R[i-1] - Rstar[i-1]) - k3*N[i-1]*L[i-1] + km3*(Ntot - N[i-1]))*deltat;
    if (L[i] < 0) {
      L[i] = 0;
    }
    R[i] = R[i-1] + (- k1*R[i-1]*pow(L[i-1], nhill) + km1*(Rtot - Rstar[i-1] - R[i-1]))*deltat;
    Rstar[i] = Rstar[i-1] + (k2*(Rtot - R[i-1] - Rstar[i-1]) - km2*Rstar[i-1])*deltat;
    N[i] = N[i-1] + (-k3*L[i-1]*N[i-1] + (km3 + k4)*(Ntot - N[i-1]))*deltat;
    
    // VOLTAGE AND THRESHOLD DYNAMICS
    voltage[i] = voltage[i-1] + (- gL/Cm*(voltage[i-1] - EL) - gamma*Rstar[i-1]/Cm*(voltage[i-1] - ER))*deltat; 
    wth[i] = wth[i-1]*exp(-deltat/tau);
    if (voltage[i] > theta0 + wth[i]) {
      voltage[i] = Vreset;
      spikestimes[i] = 1;
      wth[i] = wth[i] + Delta/tau;
    } else {
      spikestimes[i] = 0;
    } 
  } 

  return spikestimes;
}')