#include "grid/octree.h"

#include "navier-stokes/centered.h"
#include "navier-stokes/perfs.h"

#define mu(f)  (1./(clamp(f,0,1)*(1./mu1 - 1./mu2) + 1./mu2))
#include "two-phase.h"//自带vof.h与embed.h公用变量

#include "tension.h"
#if dimension == 3
# include "lambda2.h"
#endif
#include "view.h"

//#include "embed.h"不要添加如下这两个头文件，其中定义的一些东西会与two-phase中的东西相重合，造成不可逆转的bug
//#include "vof.h"
/**
We can control the maximum runtime. */

#include "maxruntime.h"

/**
The density ratio is 1000 and the dynamic viscosity ratio 100. */

#define RHOR 10.
#define MUR 10.

# define Ga 29.9
# define Bo 2.

//# define Ga 31.51737
//# define Bo 2.2222
# define MAXTIME 200
/**
We choose as length unit the diameter of the bubble. The domain is
$120^3$. *Z0* is the initial position of the bubble relative to the
bottom wall. The acceleration of gravity is set to unity, which gives
a characteristic rise velocity also of order unity, which gives a
maximum time for the simulation comparable to the domain size. */

#define WIDTH 1.61290
#define R0 0.5

int LEVEL = 4;

int main (int argc, char * argv[]) {
  maxruntime (&argc, argv);
  if (argc > 1)
    LEVEL = atoi (argv[1]);
  
  /**
  We set the domain geometry and initial refinement. */
  
  size (WIDTH);
  origin (-L0/2.,-L0/2.,-L0/2.);
  init_grid (16);
  periodic(right);
  periodic(front);
  periodic(top);
 

  /**
  We set the physical parameters: densities, viscosities and surface
  tension. */
  
  rho1 = 1.;//第一相中f为1，第二相为0，是故为1的相必须为水
  rho2 = 1./RHOR;
  mu1 = 1./Ga;
  mu2 = 1./(MUR*Ga);
  f.sigma = 1./Bo;

  /**
  We reduce the tolerance on the divergence of the flow. This is
  important to minimise mass conservation errors for these simulations
  which are very long. */
  
  TOLERANCE = 1e-4;
  run();
}


event init (t = 0) {
  if (!restore (file = "restart")) {
  fraction (f, sq(x) + sq(y) + sq(z) - sq(R0));//气泡初始位置
  }
}

/**
We add the acceleration of gravity (unity) in the downward (-y)
direction. */

event acceleration (i++) {
face vector av = a;
foreach_face(y){
double ff1 = (f[0,0,0]+f[0,-1,0])/2.;

double rho_m = fm.y[]/rho(ff1);

av.y[] -= (1. - (0.875212*rho1+0.124788*rho2)*rho_m);
}
}


//用于检测运行时间
event logfile (i = 0 ; t <= MAXTIME; i += 10) {
  double sb = 0.;
  double vby = 0., v_y = 0.;
  double Re = 0.;
  foreach(
	  reduction(+:vby) reduction(+:v_y)
	  reduction(+:sb)) {
    double dv = (1. - f[])*dv();
    vby += u.y[]*dv;
    v_y += u.y[]*dv();
    sb += dv;
  }
  
  Re = (vby/sb - v_y/(pow(WIDTH,3)))*Ga;
  
  fprintf (stderr,
	   "%.8f %.8f %.8f %.8f %.8f\n", 
	   t, sb,
	   vby/sb, v_y/(pow(WIDTH,3)), Re);
  fflush (stderr);
}

