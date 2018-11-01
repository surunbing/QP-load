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

  DifferentialState xl, yl, zl;
  DifferentialState dxl, dyl,dzl;
  DifferentialState qx, qy, qz;
  DifferentialState wx, wy, wz;
  DifferentialState r11, r12, r13, r21, r22, r23, r31, r21, r33;
  DifferentialState wwx, wwy, wwz;

  Control mx, my, mz;
  Control f;
  
  OnlineData xr, yr, zr;
  OnlineData g;
  OnlineData mf, ml;
  OnlineData L;
  OnlineData j11, j12, j13, j21, j22, j23, j31, j32, j33; 

  auto dqx = qz*wy - qy*wz;
  auto dqy = qx*wz - qz*wx;
  auto dqz = qy*wx - qx*wy;
  auto qfre3mqlqq = f * r13 * qx + f * r23 * qy + f * r33 * qz - L * mq * (dqx * dqx + dqy * dqy + dqz * dqz);
            
  // Equations of motion
  f << dot(xl) == dxl; //cos(psi) * cos(theta) * vx + (cos(psi) * sin(phi) * sin(theta) - cos(phi) * sin(psi)) * vy + (sin(phi)*sin(psi) + cos(phi)*cos(psi)*sin(theta)) * vz;
  f << dot(yl) == dyl; //cos(theta) * sin(psi) * vx + (cos(phi) * cos(psi) + sin(phi) * sin(psi) * sin(theta)) * vy + (cos(phi) * sin(psi) * sin(theta) - cos(psi) * sin(phi)) * vz;
  f << dot(zl) == dzl; //-sin(theta) * vx + cos(theta) * sin(phi) * vy + cos(phi) * cos(theta) * vz;
  f << dot(dxl) == qfre3mqlqq * qx / (mq + ml);
  f << dot(dyl) == qfre3mqlqq * qy / (mq + ml);
  f << dot(dzl) == qfre3mqlqq * qz / (mq + ml) - g;
  f << dot(qx) == dqx;
  f << dot(qy) == dqy;
  f << dot(qz) == dqz;
  f << dot(wx) ==  -(f * qy * r33 - f * qz * r23) / (L * mq);
  f << dot(wy) == (f * qx * r33 - f * qz * r13) / (L * mq);
  f << dot(wz) == -(f * qx * r23 - f * qy * r13)/(L * mq);
  f << dot(r11) == r12 * wwz - r13 * wwy;
  f << dot(r12) == r13 * wwx - r11 * wwz;
  f << dot(r13) == r11 * wwy - r12 * wwx;
  f << dot(r21) == r22 * wwz - r23 * wwy;
  f << dot(r22) == r23 * wwx - r21 * wwz;
  f << dot(r23) == r21 * wwy - r22 * wwx;
  f << dot(r31) == r32 * wwz - r33 * wwy;
  f << dot(r32) == r33 * wwx - r31 * wwz;
  f << dot(r33) == r31 * wwy - r32 * wwx;
  f << 0 == dot(wwx)*j11 - mx + dot(wwy)*j12 + dot(wwz)*j13 + wy*(j31*wx + j32*wy + j33*wz) - wz*(j21*wx + j22*wy + j23*wz);
  f << 0 == dot(wwx)*j21 - my + dot(wwy)*j22 + dot(wwz)*j23 - wx*(j31*wx + j32*wy + j33*wz) + wz*(j11*wx + j12*wy + j13*wz);
  f << 0 == dot(wwx)*j21 - my + dot(wwy)*j22 + dot(wwz)*j23 - wx*(j31*wx + j32*wy + j33*wz) + wz*(j11*wx + j12*wy + j13*wz);
  // Running cost
  
  Function h;

  // Distance errors
  h << xl - xr;
  h << yl - yr;
  h << zl - zr;
  h << mx;
  h << my;
  h << mz;
  h << f;

  DMatrix Q(7,7); 
//Q.setAll(false);
  Q(0,0) = 0.1;
  Q(1,1) = 0.1;
  Q(2,2) = 1.0;
  Q(3,3) = 1.0;
  Q(4,4) = 0.01;
  Q(5,5) = 0.01;
  Q(6,6) = 0.01;


  // Terminal cost
  Function hN;

  // Distance errors
  hN << xl - xr;
  hN << yl - yr;
  hN << zl - zr;
  
  DMatrix QN(3,3); 
//QN.setAll(false);
  QN(0,0) = 0.1;
  QN(1,1) = 0.1;
  QN(2,2) = 1.0;

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
  const double tEnd   = 5.0;

  OCP ocp( tStart, tEnd, 50);
  ocp.subjectTo(f);

  ocp.minimizeLSQ(Q, h);
  ocp.minimizeLSQEndTerm(QN, hN);

  // car can't go backward to avoid "circles"
  ocp.subjectTo( 0 <= f <= 15.0);
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
