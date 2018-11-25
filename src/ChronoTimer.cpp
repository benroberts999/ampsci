#include "ChronoTimer.h"
#include <chrono>
#include <string>
#include <iomanip> // setprecision
#include <sstream> // stringstream

//******************************************************************************
ChronoTimer::ChronoTimer(bool auto_start)
{
  running = false;
  total_time_ms = 0;
  if(auto_start) start(); //false by default
}

//******************************************************************************
void ChronoTimer::start()
{
  //note: will over-ride any existing reading!
  if(running) stop();
  running = true;
  tstart = std::chrono::high_resolution_clock::now();
}

//******************************************************************************
void ChronoTimer::stop()
{
  if(!running) return;

  double current_time = lap_reading_ms();
  running = false; //"turn off" stopwatch

  //update total time
  total_time_ms += current_time;
}

//******************************************************************************
void ChronoTimer::reset(){
  running = false;
  total_time_ms = 0;
}

//******************************************************************************
double ChronoTimer::lap_reading_ms()
/*
Returns value for current riming run (lap)
Returns double (milliseconds)
*/
{

  if(!running) return 0;

  std::chrono::high_resolution_clock::time_point tcurrent
    = std::chrono::high_resolution_clock::now();

  auto duration =
    std::chrono::duration_cast<std::chrono::microseconds>
      (tcurrent - tstart).count();

  return ((double) duration)*1.e-3;
}

//******************************************************************************
double ChronoTimer::reading_ms()
/*
Returns total value for timnig run
Returns double (milliseconds)
*/
{
  return lap_reading_ms() + total_time_ms;
}
//******************************************************************************

std::string ChronoTimer::reading_str(){
  return convertHR(reading_ms());
}
//******************************************************************************
std::string ChronoTimer::lap_reading_str(){
  return convertHR(lap_reading_ms());
}

//******************************************************************************
std::string ChronoTimer::convertHR(double t)
/*
Convers double (in ms) into formmated 2 d.p. string in units of either
ms, s, mins, or hours, depending on size.
*/
{
  double ot;

  std::string un;
  if(t<1000){
    ot = t;
    un = "ms";
  }else if(t<60000){
    ot = t/1000.;
    un = "s";
  }else if(t<3600000){
    ot = t/60000.;
    un = "mins";
  }else{
    ot = t/3600000.;
    un = "hours";
  }

  std::stringstream ss;
  ss << std::fixed << std::setprecision(2) << ot;
  return ss.str()+" "+un;

}
