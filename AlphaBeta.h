#ifndef ALPHABETA_H
#define ALPHABETA_H
//#define GF_ACCURACY 0.00001
//if GF_ACCURACY is defined then square & cubic root are calculated within this accuracy for datatypes not supported 
//by Arduino builtin functions defined below, which can also be overriden
#define GFSQRT(x) sqrt(x)
#define GFCBRT(x) pow((x),0.3333333333333333)
enum GFState
{
    GFDisplacement = 0,
    GFVelocity = 1,
    GFAcceleration = 2
};

template <class GF>
#ifdef GF_ACCURACY
  GF gfsqrt (GF x) {
    GF r, e, t;  
    r = 0.5*x;
    do {
      t = x / r;
      e = t - x;
      if (e<0) e=-e;
      r = (r + t) / 2.0; 
    } while (e>GF_ACCURACY);
    return(r);
};
#else
  inline GF gfsqrt (GF x) { return(GFSQRT(x)); };
#endif

template <class GF>
#ifdef GF_ACCURACY
  GF gfcbrt (GF x) {
    GF r, e, t;  
    r = x / 3.0;
    do {
      t = (x / r) / r;
      e = t - x;
      if (e<0) e=-e;
      r = (2.0*r + t) / 3.0; 
    } while (e>GF_ACCURACY);
    return(r);
};
#else
  inline GF gfcbrt (GF x) { return(GFCBRT(x)); };
#endif


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
      this->g = (gfsqrt(v*(v+16.0))-v)/8.0;
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
      setProcessDEV(gfsqrt(n));
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
      setMeasurementDEV(gfsqrt(n));
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
      v = 1.0 + (v - gfsqrt(v*(8.0+v)))/4.0;
      this->g = 1.0 - v*v;
      this->h = 2.0*(2.0 - this->g)-4.0*gfsqrt(1.0-this->g);
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
      v = gfsqrt(this->k*this->k + this->h*this->h*this->h*(4.0/27.0)); //v = sqrt(sq(q) + p*sq(p)*(4.0/27.0))
      v = -gfcbrt(this->k + v/2.0); //z = -pow(q + v/2.0, 0.33333333333333333333333)
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
