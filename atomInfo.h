#ifndef _ATINFO_H
#define _ATINFO_H
#include <string>
#include <iostream>

// List of "most-common" A's, for each Z.
//ans list of atom names
// Used Mathematica to generate list!

//NOTE: First entry is zero, so the array index matches Z!
//For now, goes up to E120 [121 etries]

//XXX change to ATINF ??

//LATER: can add magnetic moments, "shells" etc..
//Also: rn(?), skin thickness, half-density radius?

//Default values for A for each atom.
//Note: array index matches Z, so first entry is blank.
//Goes up to E120 (Z=120)
const int atinfo_a[121]={0,
    1,  4,  7,  9, 11, 12, 14, 16, 19, 20,
   23, 24, 27, 28, 31, 32, 35, 40, 39, 40,
   45, 48, 51, 52, 55, 56, 59, 59, 64, 65,
   70, 73, 75, 79, 80, 84, 85, 88, 89, 91,
   93, 96, 97,101,103,106,108,112,115,119,
  122,128,127,131,133,137,139,140,141,144,
  145,150,152,157,159,162,165,167,169,173,
  175,178,181,184,186,190,192,195,197,201,
  204,207,209,209,210,222,223,226,227,232,
  231,238,237,244,243,247,247,251,252,257,
  258,259,262,267,270,269,270,270,278,281,
  281,285,286,289,289,293,293,294,315,320};

const std::string atinfo_sym[121]={"0",
   "H","He","Li","Be", "B", "C", "N", "O", "F","Ne",
  "Na","Mg","Al","Si", "P", "S","Cl","Ar", "K","Ca",
  "Sc","Ti", "V","Cr","Mn","Fe","Co","Ni","Cu","Zn",
  "Ga","Ge","As","Se","Br","Kr","Rb","Sr", "Y","Zr",
  "Nb","Mo","Tc","Ru","Rh","Pd","Ag","Cd","In","Sn",
  "Sb","Te", "I","Xe","Cs","Ba","La","Ce","Pr","Nd",
  "Pm","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm","Yb",
  "Lu","Hf","Ta", "W","Re","Os","Ir","Pt","Au","Hg",
  "Tl","Pb","Bi","Po","At","Rn","Fr","Ra","Ac","Th",
  "Pa", "U","Np","Pu","Am","Cm","Bk","Cf","Es","Fm",
  "Md","No","Lr","Rf","Db","Sg","Bh","Hs","Mt","Ds",
  "Rg","Cn","Nh","Fl","Mc","Lv","Ts","Og","E119","E120"};

//Short function that returns orbital term given l
const std::string atinfo_l_array[5]={"s","p","d","f","g"};
inline std::string atinfo_l(int l){
  if(l<5) return atinfo_l_array[l];
  else return std::to_string(l);
}

//Given an atomic symbol (H, He, etc.), will return Z
//Note: Symbol must be exact, including capitalisation
inline int atinfo_get_z(std::string at){
  for(int z=0; z<121; z++){
    if(at==atinfo_sym[z]) return z;
  }
  std::cout<<"\n FAILURE 47 in atomInfo: "<<at<<" not found.\n";
  return 0;
}

#endif
