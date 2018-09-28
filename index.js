/*
  ======================================================================
   solver03.c --- simple fluid solver
  ----------------------------------------------------------------------
   Author : Jos Stam (jstam@aw.sgi.com)
   Creation Date : Jan 9 2003

   Description:

  This code is a simple prototype that demonstrates how to use the
  code provided in my GDC2003 paper entitles "Real-Time Fluid Dynamics
  for Games". This code uses OpenGL and GLUT for graphics and interface

  =======================================================================
*/

function IX(N,i,j) { return ((i)+(N+2)*(j)) }

function add_source ( N, x, s, dt ) {
  const size=(N+2)*(N+2);
  for ( let i=0 ; i<size ; i++ ) x[i] += dt*s[i];
}

function set_bnd ( N, b, x ) {
  for ( let i=1 ; i<=N ; i++ ) {
    x[IX(N, 0  ,i)] = b==1 ? -x[IX(N, 1,i)] : x[IX(N, 1,i)];
    x[IX(N, N+1,i)] = b==1 ? -x[IX(N, N,i)] : x[IX(N, N,i)];
    x[IX(N, i,0  )] = b==2 ? -x[IX(N, i,1)] : x[IX(N, i,1)];
    x[IX(N, i,N+1)] = b==2 ? -x[IX(N, i,N)] : x[IX(N, i,N)];
  }
  x[IX(N, 0  ,0  )] = 0.5*(x[IX(N, 1,0  )]+x[IX(N, 0  ,1)]);
  x[IX(N, 0  ,N+1)] = 0.5*(x[IX(N, 1,N+1)]+x[IX(N, 0  ,N)]);
  x[IX(N, N+1,0  )] = 0.5*(x[IX(N, N,0  )]+x[IX(N, N+1,1)]);
  x[IX(N, N+1,N+1)] = 0.5*(x[IX(N, N,N+1)]+x[IX(N, N+1,N)]);
}

function lin_solve ( N, b, x, x0, a, c ) {
  for ( let k=0 ; k<20 ; k++ ) {
    for ( let i=1 ; i<=N ; i++ ) { for ( let j=1 ; j<=N ; j++ ) {
      x[IX(N, i,j)] = (x0[IX(N, i,j)] + a*(x[IX(N, i-1,j)]+x[IX(N, i+1,j)]+x[IX(N, i,j-1)]+x[IX(N, i,j+1)]))/c;
    } }
    set_bnd ( N, b, x );
  }
}

function diffuse ( N, b, x, x0, diff, dt ) {
  const a=dt*diff*N*N;
  lin_solve ( N, b, x, x0, a, 1+4*a );
}

function advect ( N, b, d, d0, u, v, dt ) {
  let i0, j0, i1, j1;
  let x, y, s0, t0, s1, t1;

  const dt0 = dt*N;
  for ( let i=1 ; i<=N ; i++ ) { for ( let j=1 ; j<=N ; j++ ) {
    x = i-dt0*u[IX(N, i,j)]; y = j-dt0*v[IX(N, i,j)];
    if (x<0.5) x=0.5; if (x>N+0.5) x=N+0.5; i0=x|0; i1=i0+1;
    if (y<0.5) y=0.5; if (y>N+0.5) y=N+0.5; j0=y|0; j1=j0+1;
    s1 = x-i0; s0 = 1-s1; t1 = y-j0; t0 = 1-t1;
    d[IX(N, i,j)] = s0*(t0*d0[IX(N, i0,j0)]+t1*d0[IX(N, i0,j1)])+
           s1*(t0*d0[IX(N, i1,j0)]+t1*d0[IX(N, i1,j1)]);
  } }
  set_bnd ( N, b, d );
}

function project ( N, u, v, p, div ) {
  for ( let i=1 ; i<=N ; i++ ) { for ( let j=1 ; j<=N ; j++ ) {
    div[IX(N, i,j)] = -0.5*(u[IX(N, i+1,j)]-u[IX(N, i-1,j)]+v[IX(N, i,j+1)]-v[IX(N, i,j-1)])/N;
    p[IX(N, i,j)] = 0;
  } }
  set_bnd ( N, 0, div );
  set_bnd ( N, 0, p );

  lin_solve ( N, 0, p, div, 1, 4 );

  for ( let i=1 ; i<=N ; i++ ) { for ( let j=1 ; j<=N ; j++ ) {
    u[IX(N, i,j)] -= 0.5*N*(p[IX(N, i+1,j)]-p[IX(N, i-1,j)]);
    v[IX(N, i,j)] -= 0.5*N*(p[IX(N, i,j+1)]-p[IX(N, i,j-1)]);
  } }
  set_bnd ( N, 1, u );
  set_bnd ( N, 2, v );
}

function dens_step ( N, x, x0, u, v, diff, dt ) {
  add_source ( N, x, x0, dt );
  diffuse ( N, 0, x0, x, diff, dt );
  advect ( N, 0, x, x0, u, v, dt );
}

function vel_step ( N, u, v, u0, v0, visc, dt ) {
  add_source ( N, u, u0, dt );
  add_source ( N, v, v0, dt );
  diffuse ( N, 1, u0, u, visc, dt );
  diffuse ( N, 2, v0, v, visc, dt );
  project ( N, u0, v0, u, v );
  advect ( N, 1, u, u0, u0, v0, dt );
  advect ( N, 2, v, v0, u0, v0, dt );
  project ( N, u, v, u0, v0 );
}

class Fluid {
  constructor(N, visc, diff) {
    this.N = N;
    this.visc = visc;
    this.diff = diff;
    this.u = new Float32Array((N+2)*(N+2));
    this.u0 = new Float32Array((N+2)*(N+2));
    this.v = new Float32Array((N+2)*(N+2));
    this.v0 = new Float32Array((N+2)*(N+2));
    this.x = new Float32Array((N+2)*(N+2));
    this.x0 = new Float32Array((N+2)*(N+2));
  }

  step(dt) {
    vel_step(this.N, this.u, this.v, this.u0, this.v0, this.visc, dt);
    dens_step(this.N, this.x, this.x0, this.u, this.v, this.diff, dt);
    this.u0.fill(0)
    this.v0.fill(0)
    this.x0.fill(0)
  }

  addDensity(x, y, dd) {
    this.x0[IX(this.N, x, y)] += dd
  }
  addForce(x, y, fx, fy) {
    this.u0[IX(this.N, x, y)] += fx
    this.v0[IX(this.N, x, y)] += fy
  }

  densityAt(x, y) {
    return this.x[IX(this.N, x, y)]
  }

  velocityAt(x, y) {
    const ix = IX(this.N, x, y)
    return [this.u[ix], this.v[ix]]
  }
}

exports.Fluid = Fluid
