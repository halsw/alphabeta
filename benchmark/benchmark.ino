/*
 * This file is part of the alphabeta library
 * Usage: Benchmatk the iteration time of the filter
 * 
 * Version 1.0
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

#define NOISE_PROCESS 1
#define NOISE_MEASURE 1
#define PERIOD 0.001
#define BENCHMARK 1000

int counter;
unsigned long bmtime;
GHKFilter <float> filter(NOISE_PROCESS,NOISE_MEASURE,PERIOD);

void setup()  { 
  Serial.begin(115200);
  while (! Serial);
  Serial.println("AlphaBetaGamma Filter benchmark");
  counter = 0;
  bmtime = millis();
}

void loop() {
  if (counter++ < BENCHMARK)
    filter.update(1.0);
  else {
    float t=(millis()-bmtime) / (BENCHMARK*1.0);
    Serial.print("Filter update in ");
    Serial.print(t);
    Serial.println(" ms");
    counter = 0;
    bmtime = millis();
  }
}
