#include <iostream>
#include <iomanip>
#include <cmath>
#include <random>
#include <vector>
#include <chrono>

// Benutzte "Black Box" Funktion
template<typename T>
T f(const T& x){
    if(x<=4) return 12;
    else if (x<=8) return -4*cos(3/4. * x * M_PI) + 8;
    else if (x<=12) return 4*sin((x-8.0) * M_PI) + 4;
    else return pow(x-13., 2) + 3.;
}

template<typename T>
T integral_mc(const int n, const T& xlower, const T& xupper, const T& ylower, const T& yupper) {  
  // TODO 2: Monte Carlo Integration mit std::random_engine_generator und std::uniform_real_distribution 
  std::default_random_engine generator;
  std::uniform_real_distribution<double> distribution_x(xlower,xupper);
  std::uniform_real_distribution<double> distribution_y(ylower,yupper);
  
  double m = 0;
  double untereTeil = 0;
  
  for(int i=0;i<n;i++){
  double random_x = distribution_x(generator);
  double random_y = distribution_y(generator);
  
    if(f(random_x)>random_y) {m++ ;}

  }
  untereTeil = m/n;
  
  double integralWerte = untereTeil * ((yupper-ylower)*(xupper-xlower));

  return integralWerte;

}


template<typename T>
T integral_riemann(const int n, const T& xlower, const T& xupper) {
  // TODO 3: Recheckabschätzung
  
  double breite = (xupper-xlower)/n;
  double aktuelle_x = xlower + (breite/2);
  
  double rechteck_sum =0;
  for(int i=0;i<n-1;i++){
    aktuelle_x += breite;

    rechteck_sum += f(aktuelle_x)*breite;
  }
return rechteck_sum;
}

// Exakter Integralwert für X=[0,16]
double integral_exact() {
    // TODO 1: Bereche den exakten Wert
    double sum = 0;

    sum += 4*12;
    sum += (8*8 - 16*sin((3*M_PI*8)/4)/3*M_PI) - (8*4 - 16*sin((3*M_PI*4)/4)/3*M_PI);
    sum += (4*12 - (4*cos(M_PI*12)/M_PI)) - (4*8 - (4*cos(M_PI*8)/M_PI));
    sum += (3*16 + pow((16-13),3)/3) - (3*12 + pow((12-13),3)/3);
    
    return sum;
}

template<typename T>
T error(const T& approx, const T& exact) {
  // TODO 4: berechne den relativen Fehler
  double fehler = fabs((approx-exact)/exact);

  return fehler;
}

int main() {

  // TODO 5: Abfrage von Benutzereingaben und Ausgabe der gewünschten Ergebnisse
  bool method1,method2;
  std::cout<< "Integral mit Monte Carlo Methode berechnen? " << "(0: N, 1: Y): ";
  std::cin >> method1 ; std::cout<<std::endl;
  std::cout<< "Integral mit Recheckverfahren berechnen? " << "(0: N, 1: Y): ";
  std::cin >> method2 ; std::cout<<std::endl;

  std::vector<std::vector<double>> tabelle(5,std::vector<double>(6));
  
  int n_approx[5] = {10,100,1000,10000,100000};

  
  for(int i=0;i<5;i++){
      tabelle[i][0] = n_approx[i];
      for(int j=1;j<6;j++){
        if(j==1)  tabelle[i][j] = integral_exact() ; 
        if(method1&&method2){
        if(j==2) tabelle[i][j] = integral_mc(n_approx[i],0.0,16.0,0.0,16.0);
        else if(j==3) tabelle[i][j] = integral_riemann(n_approx[i],0.0,16.0);
        else if(j==4) tabelle[i][j] = error(integral_mc(n_approx[i],0.0,16.0,0.0,16.0),integral_exact());
        else if(j==5) tabelle[i][j] = error(integral_riemann(n_approx[i],0.0,16.0),integral_exact());
        }
        if(method1&&!method2)
        {
          if(j==2) tabelle[i][j] = integral_mc(n_approx[i],0.0,16.0,0.0,16.0);
          else if(j==4) tabelle[i][j] = error(integral_mc(n_approx[i],0.0,16.0,0.0,16.0),integral_exact());
        }
        if(!method1&&method2){
         if(j==3) tabelle[i][j] = integral_riemann(n_approx[i],0.0,16.0);
        else if(j==5) tabelle[i][j] = error(integral_riemann(n_approx[i],0.0,16.0),integral_exact());

        }
      }
    }
  std::cout<<"n"<<"\t"<<"SYMB"<<"\t"<<"MC"<<"\t"<<"RR"<<"\t"<<"MC(err)"<<"\t"<<"\t"<<"RR(err)"<<std::endl;
  for (int i=0;i<5;i++){
    for(int j=0;j<6;j++){
      std::cout<< tabelle[i][j] <<"\t";
    }
    std::cout<<std::endl;
  }
  return 0;
}
