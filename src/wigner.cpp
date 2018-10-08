#include "WIG_369j.h"
#include <iostream>
#include <sstream>
/*
Quick routine that outputs numeric values for Wigner 3,6,9J symbols.
*/

int caseCG(void){
  double j1, m1, j2, m2, J, M;
  std::cout<<"Clebsh-Gordon coef. Enter in order: <j1 m1, j2 m2 | J M>:\n";
  std::cin>>j1>>m1>>j2>>m2>>J>>M;
  double cgc = WIG::cg(j1,m1,j2,m2,J,M);
  printf("<%.1f %.1f, %.1f %.1f| %.1f %.1f> = %f\n",j1,m1,j2,m2,J,M,cgc);
  return 0;
}

int case3j(void){
  double j1, j2, j3, m1, m2, m3;
  std::cout<<"3j symbol. Enter in order: j1 j2 j3 m1 m2 m3:\n";
  std::cin>>j1>>j2>>j3>>m1>>m2>>m3;
  double cgc = WIG::threej(j1,j2,j3,m1,m2,m3);
  printf("(%5.1f %5.1f %5.1f) = %f\n",j1,j2,j3,cgc);
  printf("(%5.1f %5.1f %5.1f)\n",m1,m2,m3);
  return 0;
}

int case6j(void){
  double j1, j2, j3, j4, j5, j6;
  std::cout<<"6j symbol. Enter in order: j1 j2 j3 j4 j5 j6:\n";
  std::cin>>j1>>j2>>j3>>j4>>j5>>j6;
  double cgc = WIG::sixj(j1,j2,j3,j4,j5,j6);
  printf("{%5.1f %5.1f %5.1f} = %f\n",j1,j2,j3,cgc);
  printf("{%5.1f %5.1f %5.1f}\n",j4,j5,j6);
  return 0;
}

int case9j(void){
  double j1, j2, j3, j4, j5, j6, j7, j8, j9;
  std::cout<<"9j symbol. Enter in order: j1 j2 j3 j4 j5 j6 j7 j8 j9:\n";
  std::cin>>j1>>j2>>j3>>j4>>j5>>j6>>j7>>j8>>j9;
  double cgc = WIG::ninej(j1,j2,j3,j4,j5,j6,j7,j8,j9);
  printf("{%5.1f %5.1f %5.1f}\n",j1,j2,j3);
  printf("{%5.1f %5.1f %5.1f} = %f\n",j4,j5,j6,cgc);
  printf("{%5.1f %5.1f %5.1f}\n",j7,j8,j9);
  return 0;
}

int main(int num_in, char* argv[]){

  BEG:

  int i_type=-1;
  if(num_in>1){
    i_type = std::stoi(argv[1]);
  }else{
    std::cout<<"Which symbol do you want? 3,6,9? (0 for CG)\n";
    std::cin>>i_type;
  }

  if(i_type==0){
    return caseCG();
  }else if(i_type==3){
    return case3j();
  }else if(i_type==6){
    return case6j();
  }else if(i_type==9){
    return case9j();
  }else{
    std::cout<<i_type<<" is not valid. Enter 0, 3, 6, or 9.\n";
    goto BEG;
    // return 0;
  }

  return 0;
}
