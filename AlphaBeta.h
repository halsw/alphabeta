/*
 * This file is part of the alphabeta library
 * Usage: A template library for the implementation
 *        of alpha beta gamma filters for arduino/teensy
 * Dependencies: MathFixed lbrary (https://github.com/halsw/MathFixed)       
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
 * 
 * Classes:
 *  GFilter the alpha filter
 *  GHFilter the alpha beta filter
 *  GHKFilter the alpha beta gamma filter
 *  
 *  Types:
 *    GFState a designator to the state variables of the filters
 */
#include "MathFixed.h"

#ifndef ALPHABETA_H
#define ALPHABETA_H

enum GFState
{
    GFDisplacement = 0,
    GFVelocity = 1,
    GFAcceleration = 2
};

template <class GF>
class GFilter {
  protected:  
    unsigned long tp;
    GF err;
    GF g;
    GF nP;
    GF nM;
    GF dt;
    GF x[3];
    virtual void coefs() {
      GF v;
      v =  this->nP / this->nM * this->dt * this->dt;
      v *= v;
      this->g = (fxsqrt(v*(v+16.0))-v)/8.0;
    }
  public:
    GFilter (GF processDev, GF noiseDev, GF period, GF x0 = 0.0) {
      this->tp = millis();
      this->nP = 0.0;  
      this->nM = 0.0;  
      this->dt = 0.0;
      setPeriod(period);
      setProcessDEV(processDev);
      setMeasurementDEV(noiseDev);
      this->x[GFDisplacement] = x0;
      this->x[GFVelocity] = 0;
      this->x[GFAcceleration] = 0;
    }
    virtual GF update( GF measure) {
      this->tp = millis();
      if (this->g == 0.0) coefs();
      this->err = measure - this->x[GFDisplacement];
      return(this->x[GFDisplacement] += this->g*err);
    }    
    inline void setState( GF x, GFState s) {
      this->x[s] = x;
    }    
    inline GF getState(GFState s = GFDisplacement) {
      return(this->x[s]);
    }    
    inline GF getError() {
      return(this->err);
    }    
    virtual void setProcessDEV( GF n ) {
      this->nP = n;
    }    
    inline void setProcessVAR( GF n ) {
      setProcessDEV(fxsqrt(n));
    }    
    inline GF getProcessDEV() {
      return(this->nP);
    }    
    inline GF getProcessVAR( GF n ) {
      return(this->nP*this->nP);
    }    
    virtual void setMeasurementDEV( GF n ) {
      this->nM = n;
    }    
    inline void setMeasurementVAR( GF n ) {
      setMeasurementDEV(fxsqrt(n));
    }    
    inline GF getMeasurementDEV() {
      return(this->nM);
    }    
    inline GF getMeasurementVAR( GF n ) {
      return(this->nM*this->nM);
    }    
    virtual void setPeriod( GF t ) {
      if (t<1)
        this->dt = t;
      else  
        this->dt = t/1000.0;
      if ( this->dt != 0 )
        coefs();
    }    
    inline GF getPeriod() {
      return(this->dt);
    }    
    virtual GF getCoef(GFState s = GFDisplacement) {
      return s == GFDisplacement?this->g:0;
    }    
  };


template <class GF>
class GHFilter: public GFilter<GF> {
  protected:  
    GF h;
    virtual void coefs() {
      GF v;  
      v =  this->nP / this->nM * this->dt * this->dt;
      v = 1.0 + (v - fxsqrt(v*(8.0+v)))/4.0;
      this->g = 1.0 - v*v;
      this->h = 2.0*(2.0 - this->g)-4.0*fxsqrt(1.0-this->g);
    }
  public:
    GHFilter (GF processDev, GF noiseDev, GF period, GF x0 = 0.0, GF v0 = 0.0) : GFilter<GF>(processDev, noiseDev, 0, x0) {
      this->x[GFVelocity] = v0;
      setPeriod(period);
    }
    virtual GF update( GF measure) {
      this->tp = millis();
      if (this->g == 0.0) coefs();
      this->err = this->x[GFVelocity] * this->dt;  
      this->x[GFDisplacement] += this->err;
      this->err = measure - this->x[GFDisplacement];
      this->x[GFDisplacement] += this->g * this->err;
      this->x[GFVelocity] += this->h * this->err / this->dt;         
      return( this->x[GFDisplacement] );
    }
    virtual GF getCoef(GFState s = GFDisplacement) {
      return s == GFVelocity?this->h:GFilter<GF>::getCoef(s);
    }        
};

template <class GF>
class GHKFilter: public GHFilter<GF> {
  protected:  
    GF k;
    virtual void coefs() {
      GF v;  
      v =  this->nP / this->nM * this->dt * this->dt; //l
      this->g = v / 2.0 - 3.0; //b = l/2.0 - 3.0
      this->k = this->g + 6.0; //c = l/2.0 + 3.0
      this->h = this->k - this->g*this->g/3.0; //p = c - sq(b)/3.0
      this->k = this->g*(this->g*this->g*(2.0/27.0) - this->k/3.0) - 1.0; //q = b*sq(b)*(2.0/27.0) - b*c/3.0 - 1
      v = fxsqrt(this->k*this->k + this->h*this->h*this->h*(4.0/27.0)); //v = sqrt(sq(q) + p*sq(p)*(4.0/27.0))
      v = -fxcbrt(this->k + v/2.0); //z = -pow(q + v/2.0, 0.33333333333333333333333)
      v -= (this->h/v + this->g)/3.0; //s = z - p/(3.0*z) - b/3.0;
      this->g = 1.0 - v*v; //g = 1.0 - sq(s);
      v = 1.0 - v;
      this->h = 2.0*v*v; //h = 2*sq(1.0 - s);
      this->k = this->h*this->h/(2.0*this->g); //k=sq(h)/(2.0*g);
    }
  public:
    GHKFilter (GF processDev, GF noiseDev, GF period, GF x0 = 0.0, GF v0 = 0.0, GF a0 = 0.0) : GHFilter<GF>(processDev, noiseDev, 0, x0, v0) {
      this->x[GFAcceleration] = a0;
      setPeriod(period);
    }
    virtual GF update( GF measure) {
      GF e = this->x[GFAcceleration]*this->dt;  
      this->tp = millis();
      if (this->g == 0.0) coefs();
      this->x[GFDisplacement] += (this->x[GFVelocity] + 0.5*e) * this->dt;
      this->x[GFVelocity] += e;
      this->err = measure - this->x[GFDisplacement];
      this->x[GFDisplacement] += this->g * this->err;
      e = this->err / this->dt;
      this->x[GFVelocity] += this->h * e;
      this->x[GFAcceleration] += this->k * (e / this->dt);         
      return( this->x[GFDisplacement] );
    }    
    virtual GF getCoef(GFState s = GFDisplacement) {
      return s == GFAcceleration?this->k:GHFilter<GF>::getCoef(s);
    }        
};
#endif
