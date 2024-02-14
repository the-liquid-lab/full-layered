/**
# Vertical diffusion

We consider the vertical diffusion of a tracer $s$ with a diffusion
coefficient $D$ for the multilayer solver.

For stability, we discretise the vertical diffusion equation implicitly as
$$
\frac{(hs_l)^{n + 1} - (hs_l)^{\star}}{\Delta t} =
D \left( \frac{s_{l + 1} - s_l}{h_{l + 1 / 2}} -
\frac{s_l - s_{l - 1}}{h_{l - 1 / 2}} \right)^{n + 1}
$$
which can be expressed as the linear system
$$
\mathbf{Ms}^{n + 1} = \mathrm{rhs}
$$
where $\mathbf{M}$ is a 
[tridiagonal matrix](https://en.wikipedia.org/wiki/Tridiagonal_matrix). 
The lower, principal and upper diagonals are *a*, *b* and *c* respectively. 

Boundary conditions on the top and bottom layers need to be added to close the
system. We chose to impose a Neumann condition on the free-surface i.e.
$$
\partial_z s |_t = \dot{s}_t
$$
and a Navier slip condition on the bottom i.e.
$$
s|_b = s_b + \lambda_b \partial_z s|_b
$$ */

bool h_diffusion = false;

void vertical_diffusion_NeumanNeuman (Point point, scalar h, scalar s, double dt, double D,
				        double dst, double dsb)
{
  double a[nl], b[nl], c[nl], rhs[nl];

  /**
  The *rhs* of the tridiagonal system is $h_l s_l$. */
      
  foreach_layer()
    rhs[_layer] = s[]*h[];

  /**
  The lower, principal and upper diagonals $a$, $b$ and $c$ are given by
  $$
  a_{l > 0} = - \left( \frac{D \Delta t}{h_{l - 1 / 2}} \right)^{n + 1}
  $$
  $$
  c_{l < \mathrm{nl} - 1} = - \left( \frac{D \Delta t}{h_{l + 1 / 2}}
  \right)^{n + 1}
  $$
  $$
  b_{0 < l < \mathrm{nl} - 1} = h_l^{n + 1} - a_l - c_l
  $$
  */

  for (int l = 1; l < nl - 1; l++) {
    a[l] = - 2.*D*dt/(h[0,0,l-1] + h[0,0,l]);
    c[l] = - 2.*D*dt/(h[0,0,l] + h[0,0,l+1]);
    b[l] = h[0,0,l] - a[l] - c[l];
    //fprintf(stderr,"%g ",dt);
  }
  //fprintf(stderr,"\n");
    
  /**
  For the top layer the boundary conditions give the (ghost)
  boundary value
  $$
  s_{\mathrm{nl}} = s_{\mathrm{nl} - 1} + \dot{s}_t h_{\mathrm{nl} - 1},
  $$
  which gives the diagonal coefficient and right-hand-side
  $$
  b_{\mathrm{nl} - 1} = h_{\mathrm{nl} - 1}^{n + 1}
  - a_{\mathrm{nl} - 1}
  $$
  $$
  \mathrm{rhs}_{\mathrm{nl} - 1} = 
  (hs)_{\mathrm{nl} - 1}^{\star} + D \Delta t \dot{s}_t
  $$
  */

  a[nl-1] = - 2.*D*dt/(h[0,0,nl-2] + h[0,0,nl-1]);
  b[nl-1] = h[0,0,nl-1] - a[nl-1];
  rhs[nl-1] += D*dt*dst;

  /**
  For the bottom layer a third-order discretisation of the Navier slip
  condition gives
  $$
  \begin{aligned}
  b_0 & = h_0 + 2 \Delta t D \left( \frac{1}{h_0 + h_1} + \frac{h^2_1 + 3
  h_0 h_1 + 3 h^2_0}{\det} \right),\\
  c_0 & = - 2 \Delta t D \left( \frac{1}{h_0 + h_1} + \frac{h^2_0}{\det}
  \right),\\
  \text{rhs}_0 & = (hs_0)^{\star} + 2 \Delta t D s_b  \frac{h^2_1 + 3 h_0
  h_1 + 2 h^2_0}{\det},\\
  \det & = h_0 (h_0 + h_1)^2  + 2\lambda (3\,h_0 h_1 + 2\,h_0^2 + h_1^2),
  \end{aligned}
  $$
  */
  /**

  double den = h[]*sq(h[] + h[0,0,1]) 
    + 2.*lambda_b*(3.*h[]*h[0,0,1] + 2.*sq(h[]) + sq(h[0,0,1]));
  b[0] = h[] + 2.*dt*D*(1./(h[] + h[0,0,1]) +
			  (sq(h[0,0,1]) + 3.*h[]*h[0,0,1] + 3.*sq(h[]))/den);
  c[0] = - 2.*dt*D*(1./(h[] + h[0,0,1]) + sq(h[])/den);
  rhs[0] += 2.*dt*D*s_b*(sq(h[0,0,1]) + 3.*h[]*h[0,0,1] + 2.*sq(h[0]))/den;
  */

  c[0]=- 2.*D*dt/(h[0,0,0] + h[0,0,1]);
  b[0]= h[] - c[0];
  rhs[0] -= D*dt*dsb;

  //TO DO change Bc for 1 Layer
  if (nl == 1) {
    b[0] += c[0];
    rhs[0] += (- c[0]*h[] - D*dt) * dst;
  }
    
  /**
  We can now solve the tridiagonal system using the [Thomas
  algorithm](https://en.wikipedia.org/wiki/Tridiagonal_matrix_algorithm). */
  
  for (int l = 1; l < nl; l++) {
    b[l] -= a[l]*c[l-1]/b[l-1];
    //fprintf(stderr,"%g ",b[l]);
    rhs[l] -= a[l]*rhs[l-1]/b[l-1];
  }
  //fprintf(stderr,"\n");
  a[nl-1] = rhs[nl-1]/b[nl-1];
  s[0,0,nl-1] = a[nl-1];
  for (int l = nl - 2; l >= 0; l--)
    s[0,0,l] = a[l] = (rhs[l] - c[l]*a[l+1])/b[l];
  
}


void vertical_diffusion_NeumanNavier (Point point, scalar h, scalar s, double dt, double D,
				        double dst, double s_b, double lambda_b)
{
  double a[nl], b[nl], c[nl], rhs[nl];

  /**
  The *rhs* of the tridiagonal system is $h_l s_l$. */
      
  foreach_layer()
    rhs[_layer] = s[]*h[];

  /**
  The lower, principal and upper diagonals $a$, $b$ and $c$ are given by
  $$
  a_{l > 0} = - \left( \frac{D \Delta t}{h_{l - 1 / 2}} \right)^{n + 1}
  $$
  $$
  c_{l < \mathrm{nl} - 1} = - \left( \frac{D \Delta t}{h_{l + 1 / 2}}
  \right)^{n + 1}
  $$
  $$
  b_{0 < l < \mathrm{nl} - 1} = h_l^{n + 1} - a_l - c_l
  $$
  */

  for (int l = 1; l < nl - 1; l++) {
    a[l] = - 2.*D*dt/(h[0,0,l-1] + h[0,0,l]);
    c[l] = - 2.*D*dt/(h[0,0,l] + h[0,0,l+1]);
    b[l] = h[0,0,l] - a[l] - c[l];
  }
    
  /**
  For the top layer the boundary conditions give the (ghost)
  boundary value
  $$
  s_{\mathrm{nl}} = s_{\mathrm{nl} - 1} + \dot{s}_t h_{\mathrm{nl} - 1},
  $$
  which gives the diagonal coefficient and right-hand-side
  $$
  b_{\mathrm{nl} - 1} = h_{\mathrm{nl} - 1}^{n + 1}
  - a_{\mathrm{nl} - 1}
  $$
  $$
  \mathrm{rhs}_{\mathrm{nl} - 1} = 
  (hs)_{\mathrm{nl} - 1}^{\star} + D \Delta t \dot{s}_t
  $$
  */

  a[nl-1] = - 2.*D*dt/(h[0,0,nl-1]);
  b[nl-1] = h[0,0,nl-1] - a[nl-1];
  rhs[nl-1] += D*dt*dst;
  //fprintf(stderr,"%g\n",lambda_b);

  /**
  For the bottom layer a third-order discretisation of the Navier slip
  condition gives
  $$
  \begin{aligned}
  b_0 & = h_0 + 2 \Delta t D \left( \frac{1}{h_0 + h_1} + \frac{h^2_1 + 3
  h_0 h_1 + 3 h^2_0}{\det} \right),\\
  c_0 & = - 2 \Delta t D \left( \frac{1}{h_0 + h_1} + \frac{h^2_0}{\det}
  \right),\\
  \text{rhs}_0 & = (hs_0)^{\star} + 2 \Delta t D s_b  \frac{h^2_1 + 3 h_0
  h_1 + 2 h^2_0}{\det},\\
  \det & = h_0 (h_0 + h_1)^2  + 2\lambda (3\,h_0 h_1 + 2\,h_0^2 + h_1^2),
  \end{aligned}
  $$
  */

  double den = h[]*sq(h[] + h[0,0,1]) 
    + 2.*lambda_b*(3.*h[]*h[0,0,1] + 2.*sq(h[]) + sq(h[0,0,1]));
  b[0] = h[] + 2.*dt*D*(1./(h[] + h[0,0,1]) +
			  (sq(h[0,0,1]) + 3.*h[]*h[0,0,1] + 3.*sq(h[]))/den);
  c[0] = - 2.*dt*D*(1./(h[] + h[0,0,1]) + sq(h[])/den);
  rhs[0] += 2.*dt*D*s_b*(sq(h[0,0,1]) + 3.*h[]*h[0,0,1] + 2.*sq(h[0]))/den;
  

  if (nl == 1) {
    b[0] += c[0];
    rhs[0] += (- c[0]*h[] - D*dt) * dst;
  }
    
  /**
  We can now solve the tridiagonal system using the [Thomas
  algorithm](https://en.wikipedia.org/wiki/Tridiagonal_matrix_algorithm). */
  
  for (int l = 1; l < nl; l++) {
    b[l] -= a[l]*c[l-1]/b[l-1];
    rhs[l] -= a[l]*rhs[l-1]/b[l-1];
  }
  a[nl-1] = rhs[nl-1]/b[nl-1];
  s[0,0,nl-1] = a[nl-1];
  for (int l = nl - 2; l >= 0; l--)
    s[0,0,l] = a[l] = (rhs[l] - c[l]*a[l+1])/b[l];
}


/**
## Horizontal diffusion

This approximates
$$
h \partial_t s = D \nabla \cdot (h \nabla s)
$$
with $D$ the diffusion coefficient. Note that metric terms linked to
the slope of the layers are not taken into account. Note also that the
time discretisation is explicit so that the timestep must be limited
(manually) by $\min(\Delta^2/D)$. */


void horizontal_diffusion (scalar s, double D, double dt, scalar dst)
{
  if (D > 0.) {
    scalar d2s[];
    foreach_layer() {
      foreach() {
	double a = 0.;
	foreach_dimension(){
	  a += (s[-1] -2.*s[] + s[1]);
	}
	d2s[] = a/(sq(Delta));
      }
    }
  
    
    scalar zl[];
    foreach()
      zl[]=zb[];
    scalar d2sz[];
    for (int l = 0; l < nl; l++) {
      foreach(){
	double b = 0;
	foreach_dimension(){
	  if (l<nl-1){
	    b += (s[1,0,l]-s[-1,0,l]-s[1,0,l+1]+s[-1,0,l+1])*(h[1,0,l]-h[-1,0,l])/4.;
	    b += (s[0,0,l]-s[0,0,l+1])*(h[1,0,l]-2.*h[0,0,l]+h[-1,0,l])/2.;
	    if(l>0){
	      b -= (s[1,0,l+1]-s[-1,0,l+1]-s[1,0,l-1]+s[-1,0,l-1])*(zl[1,0]-zl[-1,0])/4.;
	      b -= (s[0,0,l+1]-s[0,0,l-1])*(zl[1,0]-2.*zl[0,0]+zl[-1,0])/2.;
	    }
	  }
	  else{
	    b += (-dst[1,0]*h[1,0,l]+dst[-1,0]*h[-1,0,l])*(h[1,0,l]-h[-1,0,l])/4.;
	    b += (-dst[]*h[0,0,l])*(h[1,0,l]-2.*h[0,0,l]+h[-1,0,l])/2.;
	    if(l>0){
	      b -= (s[1,0,l]+dst[1,0]*h[1,0,l]-s[-1,0,l]-dst[-1,0]*h[-1,0,l]-s[1,0,l-1]+s[-1,0,l-1])*(zl[1,0]-zl[-1,0])/4.;
	      b -= (s[0,0,l]+dst[1,0]*h[1,0,l]-s[0,0,l-1])*(zl[1,0]-2.*zl[0,0]+zl[-1,0])/2.;
	    }
	  }
	}
	d2sz[0,0,l]= b/sq(Delta);
      }
      foreach()
	zl[]+=h[0,0,l];
    } 
    foreach_layer() 
      foreach()
	if (h[] > dry) 
	    s[] += dt*D*d2s[];
    for (int l = 0; l < nl; l++) 
      foreach()
	if (h[0,0,l] > dry)
	  s[0,0,l] += dt*D*d2sz[0,0,l]/h[0,0,l];
  }
}


/**
void horizontal_diffusion (scalar * list, double D, double dt)
{
  if (D > 0.) {
    scalar * d2sl = list_clone (list);
    foreach_layer() {
      foreach() {
	scalar s, d2s;
	for (s,d2s in list,d2sl) {
	  double a = 0.;
	  foreach_dimension()
	    a += (hf.x[]*fm.x[]/(cm[-1] + cm[])*(s[-1] - s[]) +
		  hf.x[1]*fm.x[1]/(cm[1] + cm[])*(s[1] - s[]));
	  d2s[] = 2.*a/(cm[]*sq(Delta));
        }
      }
      foreach()
	if (h[] > dry) {
	  scalar s, d2s;
	  for (s,d2s in list,d2sl)
	    s[] += dt*D*d2s[]/h[];
	}
    }
    delete (d2sl);
    free (d2sl);
  }
}
*/

/**
event stability (i++)
{
  double rhom = 1. [-3];

  foreach_face (reduction (min:dtmax)) {
    double Hf = 0.;
    foreach_layer()
      Hf += hf.x[];
    if (Hf > dry) {
      Hf /= fm.x[];
      if (sigma[] > 0.) {
	double cp = sqrt((pi*sigma[]/(rhom*Delta))*tanh(Hf/Delta));            
	foreach_layer() {
          double c = fabs(hu.x[]/hf.x[])/CFL + cp/CFL_H;
          if (c > 0){
            double dt = min(cm[], cm[-1])*Delta/(c*fm.x[]);
	    if (dt < dtmax)
	      dtmax = dt;
          }
      	}
      }
    }
  }
}
*/

/**
# Viscous friction between layers

By default the viscosity is zero and we impose free-slip on the
free-surface and no-slip on the bottom boundary
i.e. $\dot{\mathbf{u}}_t = 0$, $\mathbf{\lambda}_b = 0$, $\mathbf{u}_b
= 0$. */

double nu = 0.;
/*
fixme: should just be:
(const) vector lambda_b[] = {0,0,0}, dut[] = {0,0,0}, u_b[] = {0,0,0};
*/
const vector lambda0[] = {0,0,0}, dut0[] = {0,0,0}, u_b0[] = {0,0,0}, dub0[] = {0,0,0};
(const) vector lambda_b = lambda0, dut = dut0, u_b = u_b0, dub = dub0 ;


/**
In the [layered solver](hydro.h), vertical viscosity is applied to the
velocity field just after advection, but before the pressure
gradient/acceleration term is applied. To take the pressure gradient
into account, we first apply the acceleration of the previous
timestep, apply vertical viscosity and then substract the previous
acceleration. */

event viscous_term (i++,last)
{
  if (nu > 0.) {    

    foreach() {
      foreach_layer()
	foreach_dimension()
	  u.x[] += dt*(ha.x[] + ha.x[1])/(hf.x[] + hf.x[1] + dry);
      foreach_dimension()
	//vertical_diffusion_NeumanNeuman (point, h, u.x, dt, nu,
      	//	    dut.x[], dub.x[]);
	vertical_diffusion_NeumanNavier (point, h, u.x, dt, nu,
        				 dut.x[], u_b.x[], dub.x[]);
    }
    if (h_diffusion){
      vector dup[];
      foreach()
	foreach_dimension()
	dup.x[]=dut.x[];
      foreach_dimension()
    	horizontal_diffusion (u.x, nu, dt, dup.x);
    }
    foreach() {
      foreach_layer()
	foreach_dimension()
	  u.x[] -= dt*(ha.x[] + ha.x[1])/(hf.x[] + hf.x[1] + dry);
    }
  }
}

/**
## References

~~~bib
@hal{popinet2020, hal-02365730}

@hal{devita2019, hal-02295398}
~~~
*/
