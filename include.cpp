//[this should be placed in a file on its own.. to be included in compiling! (not here!)
#include <iostream>
#include <cmath>
#include <stdio.h>
#include <time.h>

// order is important!  
#include "./data.cpp"
#include "./params.cpp"		//for parameters
#include "./funs.cpp"

int dodebug=0;
void debug(int lineno){
	if (dodebug==1){
		printf("DEBUG: Line no. %i\n",lineno);
	}
}
