#include <acado_toolkit.hpp>                    // Include the ACADO toolkit
#include <acado/utils/matlab_acado_utils.hpp>   // Include specific Matlab utils

USING_NAMESPACE_ACADO

#define PI 3.1415926536
#define deg2rad(d) (d/180.0*PI)

const int controlHorizon = 50;

//using namespace std;

void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] )  // Start the MEX function. Do NOT change the header of this function.
{
  USING_NAMESPACE_ACADO


  DifferentialEquation f;

  DifferentialState x, y, z;
  DifferentialState vx, vy, vz;
  DifferentialState phi, theta, psi;
  DifferentialState xl, yl;
  DifferentialState vxl, vyl;

  Control p, q, r;
  Control fT;
  
  OnlineData xr, yr, zr;
  OnlineData g;
  OnlineData mF;
  OnlineData CDL, rho, AL, 
  OnlineData l;
          
  IntermediateState fHx, fHy, fHz;
  
  auto dot_x = cos(psi) * cos(theta) * vx + (cos(psi) * sin(phi) * sin(theta) - cos(phi) * sin(psi)) * vy + (sin(phi)*sin(psi) + cos(phi)*cos(psi)*sin(theta)) * vz;
  auto dot_y = cos(theta) * sin(psi) * vx + (cos(phi) * cos(psi) + sin(phi) * sin(psi) * sin(theta)) * vy + (cos(phi) * sin(psi) * sin(theta) - cos(psi) * sin(phi)) * vz;
  auto dot_z = -sin(theta) * vx + cos(theta) * sin(phi) * vy + cos(phi) * cos(theta) * vz;
  auto delta_x = x - xl;
  auto delta_y = y - yl;
  auto zl = z + sqrt(l * l - delta_x * delta_x - delta_y * delta_y);
  auto delta_z = z - zl;
  auto dot_delta_x = dot_x - vxl;
  auto dot_delta_y = dot_y - vyl;
  auto H = dot_delta_x * delta_x + dot_delta_y * delta_y;
  auto dot_delta_z = H / delta_z;
  auto dot_zl = dot_delta_z + dot_z;
  auto vl_square = vxl * vxl + vyl * vyl + dot_zl * dot_zl;
  auto fDL_m = -0.5 * cDL * rho * AL * vl_square;
  auto fDLx = 
  
  // Equations of motion
  f << dot(x) == dot_x; //cos(psi) * cos(theta) * vx + (cos(psi) * sin(phi) * sin(theta) - cos(phi) * sin(psi)) * vy + (sin(phi)*sin(psi) + cos(phi)*cos(psi)*sin(theta)) * vz;
  f << dot(y) == dot_y; //cos(theta) * sin(psi) * vx + (cos(phi) * cos(psi) + sin(phi) * sin(psi) * sin(theta)) * vy + (cos(phi) * sin(psi) * sin(theta) - cos(psi) * sin(phi)) * vz;
  f << dot(z) == dot_z; //-sin(theta) * vx + cos(theta) * sin(phi) * vy + cos(phi) * cos(theta) * vz;
  f << dot(vx) == r * vy - q * vz - g * sin(theta) + (fHz * cos(psi) * cos(theta) - fHz * cos(theta) * sin(psi)) / mF;
  f << dot(vy) == p * vz - r * vx + g * cos(theta) * sin(phi) +  (fHy * (cos(phi) * cos(psi) + sin(phi) * sin(psi) * sin(theta)) - fHz * (cos(phi) * sin(psi) - cos(psi) * sin(phi) * sin(theta)) + fHz * cos(theta) * sin(phi)) / mF;
  f << dot(vz) == q * vx - p * vy + g * cos(phi) * cos(theta) -(fT + fHy * (cos(psi) * sin(phi) - cos(phi) * sin(psi) * sin(theta)) - fHz * (sin(phi) * sin(psi) + cos(phi) * cos(psi) * sin(theta)) - fHz * cos(phi) * cos(theta)) / mF;
  f << dot(phi) == p + (r * cos(phi) * sin(theta)) / cos(theta) + (q * sin(phi) * sin(theta)) / cos(theta);
  f << dot(theta) == q * cos(phi) - r * sin(phi);
  f << dot(psi) == (r * cos(phi)) / cos(theta) + (q * sin(phi)) / cos(theta);
  f << dot(xl) == vxl;
  f << dot(yl) == vyl;
  f << dot(vxl) == ;
  f << dot(vyl) == ;

  // Running cost
  Function h;

  // Distance errors
  h << x - xr;
  h << y - yr;
  h << z - zr;
  h << xl - xr;
  h << yl - yr;
  h << p;
  h << q;
  h << r;
  h << fT;

  DMatrix Q(9,9); 
//Q.setAll(false);
  Q(0,0) = 0.1;
  Q(1,1) = 0.1;
  Q(2,2) = 1.0;
  Q(3,3) = 1.0;
  Q(4,4) = 1.0;
  Q(5,5) = 0.01;
  Q(6,6) = 0.01;
  Q(7,7) = 0.01;
  Q(8,8) = 0.01;


  // Terminal cost
  Function hN;

  // Distance errors
  h << x - xr;
  h << y - yr;
  h << z - zr;
  h << xl - xr;
  h << yl - yr;

  DMatrix QN(5,5); 
//QN.setAll(false);
  QN(0,0) = 0.1;
  QN(1,1) = 0.1;
  QN(2,2) = 1.0;
  QN(3,3) = 1.0;
  QN(4,4) = 1.0;

  // Non uniform time grid
  // First 5 timesteps are 0.05, after that it's 0.15
  DMatrix numSteps(20, 1);
  for (int i = 0; i < 5; i++){
    numSteps(i) = 1;
  }
  for (int i = 5; i < 20; i++){
    numSteps(i) = 3;
  }

  // Setup Optimal Control Problem
  const double tStart = 0.0;
  const double tEnd   = 2.5;

  OCP ocp( tStart, tEnd, numSteps);
  ocp.subjectTo(f);

  ocp.minimizeLSQ(Q, h);
  ocp.minimizeLSQEndTerm(QN, hN);

  // car can't go backward to avoid "circles"
  ocp.subjectTo( deg2rad(-90) <= psi <= deg2rad(90));
  // more than absolute max steer angle
  ocp.subjectTo( deg2rad(-50) <= delta <= deg2rad(50));
  ocp.setNOD(18);

  OCPexport mpc(ocp);
  mpc.set( HESSIAN_APPROXIMATION, GAUSS_NEWTON );
  mpc.set( DISCRETIZATION_TYPE, MULTIPLE_SHOOTING );
  mpc.set( INTEGRATOR_TYPE, INT_RK4 );
  mpc.set( NUM_INTEGRATOR_STEPS, 1 * controlHorizon);
  mpc.set( MAX_NUM_QP_ITERATIONS, 500);
  // mpc.set( CG_USE_VARIABLE_WEIGHTING_MATRIX, YES);

  mpc.set( SPARSE_QP_SOLUTION, CONDENSING );
  mpc.set( QP_SOLVER, QP_QPOASES );
  mpc.set( HOTSTART_QP, YES );
  mpc.set( GENERATE_TEST_FILE, YES);
  mpc.set( GENERATE_MAKE_FILE, YES );
  mpc.set( GENERATE_MATLAB_INTERFACE, YES );
  mpc.set( GENERATE_SIMULINK_INTERFACE, YES );

  if (mpc.exportCode( "mpc_export" ) != SUCCESSFUL_RETURN)
    exit( EXIT_FAILURE );

  mpc.printDimensionsQP( );

 // return EXIT_SUCCESS;
}
