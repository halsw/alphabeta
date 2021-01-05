/*
 * This file is part of the alphabeta library
 * Usage: Provide an example use of the library
 * 
 * Version 1.1.0
 * Developed by Evan https://github.com/halsw
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, version 3 of the License
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */
#include "AlphaBeta.h"
#include <limits.h> 

#define NOISE_PROCESS 0.5
#define NOISE_MEASURE 0.25
#define PERIOD_MS 50

GHKFilter <float> filter(sqrt(NOISE_PROCESS),sqrt(NOISE_MEASURE),PERIOD_MS);

unsigned int wait() {
  static unsigned int load = 0;
  static unsigned long loopMS=PERIOD_MS;
  unsigned long now=millis();
  if (now>loopMS) {
    load = 100<<8;
    loopMS = now + PERIOD_MS;
  } else {
    load += (25600 - (loopMS-now)*(25600/PERIOD_MS) - load)>>3;
    while (loopMS>millis());
    loopMS += PERIOD_MS;
  }
  return(load);
}

void setup() {
  Serial.begin(115200);
  randomSeed(analogRead(0)); //assuming A0 is not connected
}

void loop() {
  static double t=0;
  static double w=TWO_PI/PERIOD_MS/10;
  static double a=NOISE_PROCESS*sqrt(5.0);
  static double mean = 0;
  static double variance = 0;
  static unsigned long sample = 1;
  double actual, measure, recover, error, accumulate;
  int load;
  t += PERIOD_MS/1000.0; //update simulation time
  actual = a*sin(w*t); //Generate fast changing sine wave (x10 the period of sampling)
  load = wait()>>8; //Sample at exact intervals
  measure = actual + (NOISE_MEASURE*((LONG_MAX>>1)-random(LONG_MAX)))/LONG_MAX; //Add noise to measurement
  recover = filter.update(measure); // apply filter
  error = recover - actual - mean; // calculate mean and variance of error
  accumulate = error / sample;
  if (sample>1) variance -= variance / ( sample - 1);
  variance += error * accumulate;
  mean += accumulate;
  if (sample++ % 100 == 0) {
    Serial.print("Recovered signal error stdev:");
    Serial.print(sqrt(variance),4);
    Serial.print(" load:");
    Serial.print(load);
    Serial.print(" ");
    Serial.println("%");
  }
}
