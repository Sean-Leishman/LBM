// This is an example 2D lattice Boltzmann method code.
// It uses the D2Q9 lattice including forcing term.
//
// This code may be changed and distributed freely.
//
// The lattice velocities are defined according to the following scheme:
// index:   0  1  2  3  4  5  6  7  8
// ----------------------------------
// x:       0 +1 -1  0  0 +1 -1 +1 -1
// y:       0  0  0 +1 -1 +1 -1 -1 +1
//
// 8 3 5  ^y
//  \|/   |   x
// 2-0-1   --->
//  /|\
// 6 4 7

/// *********************
/// PREPROCESSOR COMMANDS
/// *********************

#include <vector> // vector containers
#include <cmath> // mathematical library
#include <iostream> // for the use of 'cout'
#include <fstream> // file streams
#include <sstream> // string streams
#include <cstdlib> // standard library
#define SQ(x) ((x) * (x)) // square function; replaces SQ(x) by ((x) * (x)) in the code

using namespace std; // permanently use the standard namespace

/// *********************
/// SIMULATION PARAMETERS
/// *********************

// These are the relevant simulation parameters.
// They are set in the parameter file.

// Geometry
int Nx = 0; // number of lattice nodes along the x-axis (periodic)
int Ny = 0; // number of lattice nodes along the y-axis (including two wall nodes)
int fingerNum = 0; // number of fingers
int fingerDepth = 0; // depth of finger
int fingerWidth = 0; // width of finger
int fingerOpening = 0; // width of finger opening
// Fluid and flow properties
double tau = 0.; // relaxation time for Navier-Stokes equation (viscosity = (tau - 0.5) / 3)
double tauG = 0.; // relaxation time for advection-diffusion equation (diffusion = (tauG - 0.5) / 3)
double gravity = 0.; // force density due to gravity (in positive x-direction)
double velWall = 0.; // velocity of "fluid wall" in main channel
double reactionRate = 0.; // rate of solute consumption
// Simulation control
int tNum = 0; // number of time steps (running from 1 to tNum)
int tDiffOn = 0; // number of time steps after which diffusion is switched on
int tInfo = 0; // info time step (screen message will be printed every tInfo step)
int tVTK = 0; // disk write time step for VTK files
int tData = 0; // disk write time step for data output

/// *****************
/// DECLARE VARIABLES
/// *****************

double omega = 0.; // relaxation frequency for Navier-Stokes equation
double omegaG = 0.; // relaxation frequency for advection-diffusion equation
double ***pop, ***pop_old; // LBM populations (old and new) for Navier-Stokes equation
double ***popG, ***popG_old; // LBM populations (old and new) for advection-diffusion equation
int **obstacle; // obstacle sites
int **catalyst; // catalyst sites
double **density; // fluid density
double **conc; // concentration
double **velocity_x; // fluid velocity (x-component)
double **velocity_y; // fluid velocity (y-component)
double **force_x; // fluid force (x-component)
double **force_y; // fluid force (y-component)
double *reaction_profile; // reaction profile along finger axis
double pop_eq[9]; // equilibrium populations for Navier-Stokes equation
double popG_eq[9]; // equilibrium populations for advection-diffusion equation
const double weight[9] = {4./9., 1./9., 1./9., 1./9., 1./9., 1./36., 1./36., 1./36., 1./36.}; // lattice weights
double reacted = 0.; // amount of substance reacted each time step
double convergence_mass_old = 0.;
double convergence_mass_new = 0.;
double convergence_vel_old = 0.;
double convergence_vel_new = 0.;

/// *****************
/// DECLARE FUNCTIONS
/// *****************

void read_parameters(); // read parameters from file
void initialise(); // allocate memory and initialise variables
void inletOutletBCs(); // impose inlet and outlet boundary conditions
void LBM(int); // perform LBM operations
void momenta(); // compute fluid density and velocity from the populations
void reaction(); // reaction model
void equilibrium(double, double, double); // compute the equilibrium populations for Navier-Stokes equation
void equilibriumG(double, double, double); // compute the equilibrium populations for advection-diffusion equation
void write_fluid_vtk(int); // write the fluid state to the disk as VTK file
void write_fluid_data(int); // write the fluid data to the disk
void check_convergence(int); // checks for convergence criteria
void write_reaction_profile(); // writes reaction profile to disk

/// *************
/// MAIN FUNCTION
/// *************

// This is the main function, containing the simulation initialization and the simulation loop.

int main() {

  /// ************
  /// PREPARATIONS
  /// ************

  read_parameters(); // read parameter file
  initialise(); // allocate memory and initialise variables

  /// ***************
  /// SIMULATION LOOP
  /// ***************

  cout << "Starting simulation..." << endl;
  for(int t = 1; t <= tNum; ++t) { // run over all times between 1 and tNum
    inletOutletBCs(); // impose inlet and outlet boundary conditions
    LBM(t); // perform collision, propagation, and bounce-back
    reaction(); // perform reaction model
    /// Write fluid to VTK files
    if(t % tVTK == 0) {
      write_fluid_vtk(t);
    }
    // /// Write data files
    if(t % tData == 0) {
      write_fluid_data(t);
    }
    /// Report end of time step
    if(t % tInfo == 0) {
      cout << "  Completed time step " << t << " in [1, " << tNum << "]" << endl;
    }
    /// Check convergence criteria
    if(t % tData == 0 && t > 5000) {
      check_convergence(t);
    }
  }
  /// Report successful end of simulation
  cout << "Simulation complete without convergence reached!" << endl;

  return 0;
}

/// *************************
/// READ PARAMETERS FROM FILE
/// *************************

// The parameter file is read, and all parameters and set.

void read_parameters() {
  cout << "Reading parameter file..." << endl;
  ifstream data_file("parameters.dat");
  if(!data_file) {
    cout << "  parameters.dat could not be found!" << endl;
    exit(1);
  }
  char dummy_char[1000]; // dummy string for reading lines
  // Read geometry
  data_file.getline(dummy_char, 999);
  // Read Nx
  data_file >> Nx;
  data_file.getline(dummy_char, 999);
  cout << "  Nx = " << Nx << endl;
  // Read Ny
  data_file >> Ny;
  data_file.getline(dummy_char, 999);
  cout << "  Ny = " << Ny << endl;
  // Read fingerNum
  data_file >> fingerNum;
  data_file.getline(dummy_char, 999);
  cout << "  fingerNum = " << fingerNum << endl;
  // Read fingerDepth
  data_file >> fingerDepth;
  data_file.getline(dummy_char, 999);
  cout << "  fingerDepth = " << fingerDepth << endl;
  // Read fingerWidth
  data_file >> fingerWidth;
  data_file.getline(dummy_char, 999);
  cout << "  fingerWidth = " << fingerWidth << endl;
  // Read fingerOpening
  data_file >> fingerOpening;
  data_file.getline(dummy_char, 999);
  cout << "  fingerOpening = " << fingerOpening << endl;
  // Read fluid and flow properties
  data_file.getline(dummy_char, 999);
  // Read tau
  data_file >> tau ;
  data_file.getline(dummy_char, 999);
  cout << "  tau = " << tau << endl;
  omega = 1. / tau;
  // Read tauG
  data_file >> tauG ;
  data_file.getline(dummy_char, 999);
  cout << "  tauG = " << tauG << endl;
  omegaG = 1. / tauG;
  // Read gravity
  data_file >> gravity ;
  data_file.getline(dummy_char, 999);
  cout << "  gravity = " << gravity << endl;
  // Read velWall
  data_file >> velWall ;
  data_file.getline(dummy_char, 999);
  cout << "  velWall = " << velWall << endl;
  // Read reactionRate
  data_file >> reactionRate ;
  data_file.getline(dummy_char, 999);
  cout << "  reactionRate = " << reactionRate << endl;
  // Read simulation control
  data_file.getline(dummy_char, 999);
  // Read tNum
  data_file >> tNum ;
  data_file.getline(dummy_char, 999);
  cout << "  tNum = " << tNum << endl;
  // Read tDiffOn
  data_file >> tDiffOn ;
  data_file.getline(dummy_char, 999);
  cout << "  tDiffOn = " << tDiffOn << endl;
  // Read tInfo
  data_file >> tInfo ;
  data_file.getline(dummy_char, 999);
  cout << "  tInfo = " << tInfo << endl;
  // Read tVTK
  data_file >> tVTK ;
  data_file.getline(dummy_char, 999);
  cout << "  tVTK = " << tVTK << endl;
  // Read tData
  data_file >> tData ;
  data_file.getline(dummy_char, 999);
  cout << "  tData = " << tData << endl;

  return;
}

/// ****************************************
/// ALLOCATE MEMORY AND INITIALIZE VARIABLES
/// ****************************************

// The memory for lattice variables (populations, density, velocity, forces) is allocated.
// The variables are initialised.

void initialise() {
  cout << "Initialising simulation..." << endl;
  /// Create folders, delete data file
  // Make sure that the VTK folders exist.
  int ignore; // ignore return value of system calls
  ignore = system("mkdir -p vtk_fluid"); // create folder if not existing
  /// Allocate memory.
  obstacle = new int*[Nx];
  catalyst = new int*[Nx];
  density = new double*[Nx];
  conc = new double*[Nx];
  velocity_x = new double*[Nx];
  velocity_y = new double*[Nx];
  force_x = new double*[Nx];
  force_y = new double*[Nx];
  reaction_profile = new double[Ny];
  for(int X = 0; X < Nx; ++X) {
    obstacle[X] = new int[Ny];
    catalyst[X] = new int[Ny];
    density[X] = new double[Ny];
    conc[X] = new double[Ny];
    velocity_x[X] = new double[Ny];
    velocity_y[X] = new double[Ny];
    force_x[X] = new double[Ny];
    force_y[X] = new double[Ny];
  }
  /// Initialise the fluid density and velocity.
  // Start with unit density and zero velocity.
  for(int X = 0; X < Nx; ++X) {
    for(int Y = 0; Y < Ny; ++Y) {
      obstacle[X][Y] = 0;
      catalyst[X][Y] = 0;
      density[X][Y] = 1.;
      conc[X][Y] = 0.;
      velocity_x[X][Y] = 0.;
      velocity_y[X][Y] = 0.;
      force_x[X][Y] = 0.;
      force_y[X][Y] = 0.;
    }
  }
  /// Initialise geometry.
  // Set side walls.
  for(int X = 0; X < Nx; ++X) {
    obstacle[X][0] = 1;
    obstacle[X][Ny - 1] = 1;
  }
  // Set finger obstacles and catalyst.
  const int unitLength = Nx / fingerNum; // length of a single repeating finger unit
  for(int i = 0; i < fingerNum; ++i) {
    // Find positions.
    const int startUnit = i * unitLength; // first node along x-axis in unit
    const int startFinger = startUnit + (unitLength - fingerWidth) / 2; // first node along x-axis in finger
    const int startOpening = startUnit + (unitLength - fingerOpening) / 2; // first node along x-axis in finger opening
    const int endOpening = startOpening + fingerOpening - 1; // last node along x-axis in finger opening
    const int endFinger = startFinger + fingerWidth - 1; // last node along x-axis in finger
    const int endUnit = (i + 1) * unitLength - 1; // last node along x-axis in unit
    // Add upstream finger obstacle.
    for(int X = startUnit; X < startFinger; ++X) {
      for(int Y = Ny - fingerDepth - 1; Y < Ny - 1; ++Y) {
        obstacle[X][Y] = 1;
      }
    }
    // Add downstream finger obstacle.
    for(int X = endFinger + 1; X <= endUnit; ++X) {
      for(int Y = Ny - fingerDepth - 1; Y < Ny - 1; ++Y) {
        obstacle[X][Y] = 1;
      }
    }
    // Add upstream finger opening obstacle.
    for(int X = startUnit; X < startOpening; ++X) {
      for(int Y = Ny - fingerDepth - 1; Y <= Ny - fingerDepth; ++Y) {
        obstacle[X][Y] = 1;
      }
    }
    // Add downstream finger opening obstacle.
    for(int X = endOpening + 1; X <= endUnit; ++X) {
      for(int Y = Ny - fingerDepth - 1; Y <= Ny - fingerDepth; ++Y) {
        obstacle[X][Y] = 1;
      }
    }
    // Set upstream/downstream catalyst sites.
    for(int Y = Ny - fingerDepth; Y < Ny - 1; ++Y) {
      catalyst[startFinger][Y] = 1;
      catalyst[endFinger][Y] = 1;
    }
    // Set end catalyst sites.
    for(int X = startFinger; X <= endFinger; ++X) {
      catalyst[X][Ny - 2] = 1;
    }
  }
  /// Initialise concentration.
  // for(int X = 0; X < Nx; ++X) {
  //   for(int Y = 0; Y < Ny; ++Y) {
  //     if(Y > 0 && Y < Ny - fingerDepth - 1) {
  //       conc[X][Y] = 1.;
  //     }
  //     else {
  //       conc[X][Y] = 0.;
  //     }
  //   }
  // }
  /// Allocate memory for populations.
  pop = new double**[9];
  pop_old = new double**[9];
  popG = new double**[9];
  popG_old = new double**[9];
  for(int c_i = 0; c_i < 9; ++c_i) {
    pop[c_i] = new double*[Nx];
    pop_old[c_i] = new double*[Nx];
    popG[c_i] = new double*[Nx];
    popG_old[c_i] = new double*[Nx];
    for(int X = 0; X < Nx; ++X) {
      pop[c_i][X] = new double[Ny];
      pop_old[c_i][X] = new double[Ny];
      popG[c_i][X] = new double[Ny];
      popG_old[c_i][X] = new double[Ny];
    }
  }
  /// Initialise populations.
  // Navier-Stokes
  for(int X = 0; X < Nx; ++X) {
    for(int Y = 0; Y < Ny; ++Y) {
      equilibrium(density[X][Y], velocity_x[X][Y], velocity_y[X][Y]);
      for(int c_i = 0; c_i < 9; ++c_i) {
        pop_old[c_i][X][Y] = pop_eq[c_i];
        pop[c_i][X][Y] = pop_eq[c_i];
      }
    }
  }
  // Advection diffusion
  for(int X = 0; X < Nx; ++X) {
    for(int Y = 0; Y < Ny; ++Y) {
      equilibriumG(conc[X][Y], velocity_x[X][Y], velocity_y[X][Y]);
      for(int c_i = 0; c_i < 9; ++c_i) {
        popG_old[c_i][X][Y] = popG_eq[c_i];
        popG[c_i][X][Y] = popG_eq[c_i];
      }
    }
  }

  return;
}

/// *******************
/// COMPUTE EQUILIBRIUM
/// *******************

// This function computes the equilibrium populations for the Navier-Stokes equation.

void equilibrium(double den, double vel_x, double vel_y) {
  pop_eq[0] = weight[0] * den * (1                                                     - 1.5 * (SQ(vel_x) + SQ(vel_y)));
  pop_eq[1] = weight[1] * den * (1 + 3 * (  vel_x        ) + 4.5 * SQ(  vel_x        ) - 1.5 * (SQ(vel_x) + SQ(vel_y)));
  pop_eq[2] = weight[2] * den * (1 + 3 * (- vel_x        ) + 4.5 * SQ(- vel_x        ) - 1.5 * (SQ(vel_x) + SQ(vel_y)));
  pop_eq[3] = weight[3] * den * (1 + 3 * (          vel_y) + 4.5 * SQ(          vel_y) - 1.5 * (SQ(vel_x) + SQ(vel_y)));
  pop_eq[4] = weight[4] * den * (1 + 3 * (        - vel_y) + 4.5 * SQ(        - vel_y) - 1.5 * (SQ(vel_x) + SQ(vel_y)));
  pop_eq[5] = weight[5] * den * (1 + 3 * (  vel_x + vel_y) + 4.5 * SQ(  vel_x + vel_y) - 1.5 * (SQ(vel_x) + SQ(vel_y)));
  pop_eq[6] = weight[6] * den * (1 + 3 * (- vel_x - vel_y) + 4.5 * SQ(- vel_x - vel_y) - 1.5 * (SQ(vel_x) + SQ(vel_y)));
  pop_eq[7] = weight[7] * den * (1 + 3 * (  vel_x - vel_y) + 4.5 * SQ(  vel_x - vel_y) - 1.5 * (SQ(vel_x) + SQ(vel_y)));
  pop_eq[8] = weight[8] * den * (1 + 3 * (- vel_x + vel_y) + 4.5 * SQ(- vel_x + vel_y) - 1.5 * (SQ(vel_x) + SQ(vel_y)));
  return;
}

// This function computes the equilibrium populations for the advection-diffusion equation.

void equilibriumG(double conc, double vel_x, double vel_y) {
  popG_eq[0] = weight[0] * conc * (1                                                     - 1.5 * (SQ(vel_x) + SQ(vel_y)));
  popG_eq[1] = weight[1] * conc * (1 + 3 * (  vel_x        ) + 4.5 * SQ(  vel_x        ) - 1.5 * (SQ(vel_x) + SQ(vel_y)));
  popG_eq[2] = weight[2] * conc * (1 + 3 * (- vel_x        ) + 4.5 * SQ(- vel_x        ) - 1.5 * (SQ(vel_x) + SQ(vel_y)));
  popG_eq[3] = weight[3] * conc * (1 + 3 * (          vel_y) + 4.5 * SQ(          vel_y) - 1.5 * (SQ(vel_x) + SQ(vel_y)));
  popG_eq[4] = weight[4] * conc * (1 + 3 * (        - vel_y) + 4.5 * SQ(        - vel_y) - 1.5 * (SQ(vel_x) + SQ(vel_y)));
  popG_eq[5] = weight[5] * conc * (1 + 3 * (  vel_x + vel_y) + 4.5 * SQ(  vel_x + vel_y) - 1.5 * (SQ(vel_x) + SQ(vel_y)));
  popG_eq[6] = weight[6] * conc * (1 + 3 * (- vel_x - vel_y) + 4.5 * SQ(- vel_x - vel_y) - 1.5 * (SQ(vel_x) + SQ(vel_y)));
  popG_eq[7] = weight[7] * conc * (1 + 3 * (  vel_x - vel_y) + 4.5 * SQ(  vel_x - vel_y) - 1.5 * (SQ(vel_x) + SQ(vel_y)));
  popG_eq[8] = weight[8] * conc * (1 + 3 * (- vel_x + vel_y) + 4.5 * SQ(- vel_x + vel_y) - 1.5 * (SQ(vel_x) + SQ(vel_y)));
  return;
}

/// ********************************
/// INLET/OUTLET BOUNDARY CONDITIONS
/// ********************************

void inletOutletBCs() {
  for(int Y = 0; Y < Ny - fingerDepth - 1; Y++) {
    // Treat concentration on inlet boundary
    conc[0][Y] = 1.;
    equilibriumG(conc[0][Y], velocity_x[0][Y], velocity_y[0][Y]);
    for(int i = 0; i < 9; ++i) {
      popG[i][0][Y] = popG_eq[i];
    }
    // Treat concentration on outlet boundary
    conc[Nx - 1][Y] = conc[Nx - 2][Y];
    equilibriumG(conc[Nx - 1][Y], velocity_x[Nx - 1][Y], velocity_y[Nx - 1][Y]);
    for(int i = 0; i < 9; ++i) {
      popG[i][Nx - 1][Y] = popG_eq[i];
    }
  }
  return;
}

/// **********************
/// PERFORM LBM OPERATIONS
/// **********************

void LBM(int time) {
  /// Swap populations
  // The present code used old and new populations which are swapped at the beginning of each time step.
  // This is sometimes called 'double-buffered' or 'ping-pong' algorithm.
  // This way, the old populations are not overwritten during propagation.
  // The resulting code is easier to write and to debug.
  // The memory requirement for the populations is twice as large.
  double ***swap_temp = pop_old;
  pop_old = pop;
  pop = swap_temp;
  double ***swapG_temp = popG_old;
  popG_old = popG;
  popG = swapG_temp;
  /// Lattice Boltzmann equation
  // The lattice Boltzmann equation is solved in the following.
  // The algorithm includes
  // - computation of the lattice force
  // - combined collision and propagation (faster than first collision and then propagation)
  for(int X = 0; X < Nx; ++X) {
    for(int Y = 1; Y < Ny - 1; ++Y) {
      /// Compute new populations
      // This is the lattice Boltzmann equation for the Navier-Stokes equation
      equilibrium(density[X][Y], (velocity_x[X][Y] + (force_x[X][Y] + gravity) * tau / density[X][Y]), (velocity_y[X][Y] + (force_y[X][Y]) * tau / density[X][Y]));
      pop[0][X]                [Y]     = pop_old[0][X][Y] * (1 - omega) + pop_eq[0] * omega;
      pop[1][(X + 1) % Nx]     [Y]     = pop_old[1][X][Y] * (1 - omega) + pop_eq[1] * omega;
      pop[2][(X - 1 + Nx) % Nx][Y]     = pop_old[2][X][Y] * (1 - omega) + pop_eq[2] * omega;
      pop[3][X]                [Y + 1] = pop_old[3][X][Y] * (1 - omega) + pop_eq[3] * omega;
      pop[4][X]                [Y - 1] = pop_old[4][X][Y] * (1 - omega) + pop_eq[4] * omega;
      pop[5][(X + 1) % Nx]     [Y + 1] = pop_old[5][X][Y] * (1 - omega) + pop_eq[5] * omega;
      pop[6][(X - 1 + Nx) % Nx][Y - 1] = pop_old[6][X][Y] * (1 - omega) + pop_eq[6] * omega;
      pop[7][(X + 1) % Nx]     [Y - 1] = pop_old[7][X][Y] * (1 - omega) + pop_eq[7] * omega;
      pop[8][(X - 1 + Nx) % Nx][Y + 1] = pop_old[8][X][Y] * (1 - omega) + pop_eq[8] * omega;
      // This is the lattice Boltzmann equation for the Navier-Stokes equation
      if(time > tDiffOn) {
        equilibriumG(conc[X][Y], velocity_x[X][Y], velocity_y[X][Y]);
        popG[0][X]                [Y]     = popG_old[0][X][Y] * (1 - omegaG) + popG_eq[0] * omegaG;
        popG[1][(X + 1) % Nx]     [Y]     = popG_old[1][X][Y] * (1 - omegaG) + popG_eq[1] * omegaG;
        popG[2][(X - 1 + Nx) % Nx][Y]     = popG_old[2][X][Y] * (1 - omegaG) + popG_eq[2] * omegaG;
        popG[3][X]                [Y + 1] = popG_old[3][X][Y] * (1 - omegaG) + popG_eq[3] * omegaG;
        popG[4][X]                [Y - 1] = popG_old[4][X][Y] * (1 - omegaG) + popG_eq[4] * omegaG;
        popG[5][(X + 1) % Nx]     [Y + 1] = popG_old[5][X][Y] * (1 - omegaG) + popG_eq[5] * omegaG;
        popG[6][(X - 1 + Nx) % Nx][Y - 1] = popG_old[6][X][Y] * (1 - omegaG) + popG_eq[6] * omegaG;
        popG[7][(X + 1) % Nx]     [Y - 1] = popG_old[7][X][Y] * (1 - omegaG) + popG_eq[7] * omegaG;
        popG[8][(X - 1 + Nx) % Nx][Y + 1] = popG_old[8][X][Y] * (1 - omegaG) + popG_eq[8] * omegaG;
      }
    }
  }
  /// Bounce-back
  // The populations have to be bounced back at rigid walls.
  // Periodicity of the lattice in x-direction is taken into account via the %-operator.
  for(int X = 0; X < Nx; ++X) {
    // Standard bounce-back everywhere, except at the bottom and top walls
    for(int Y = 1; Y < Ny - 1; ++Y) {
      if(obstacle[X][Y]) {
        // Bounce-back for Navier-Stokes equation
        pop[1][(X + 1) % Nx]     [Y]     = pop[2][X][Y];
        pop[2][(X - 1 + Nx) % Nx][Y]     = pop[1][X][Y];
        pop[3][X]                [Y + 1] = pop[4][X][Y];
        pop[4][X]                [Y - 1] = pop[3][X][Y];
        pop[5][(X + 1) % Nx]     [Y + 1] = pop[6][X][Y];
        pop[6][(X - 1 + Nx) % Nx][Y - 1] = pop[5][X][Y];
        pop[7][(X + 1) % Nx]     [Y - 1] = pop[8][X][Y];
        pop[8][(X - 1 + Nx) % Nx][Y + 1] = pop[7][X][Y];
        // Bounce-back for advection-diffusion equation
        if(time > tDiffOn) {
          popG[1][(X + 1) % Nx]     [Y]     = popG[2][X][Y];
          popG[2][(X - 1 + Nx) % Nx][Y]     = popG[1][X][Y];
          popG[3][X]                [Y + 1] = popG[4][X][Y];
          popG[4][X]                [Y - 1] = popG[3][X][Y];
          popG[5][(X + 1) % Nx]     [Y + 1] = popG[6][X][Y];
          popG[6][(X - 1 + Nx) % Nx][Y - 1] = popG[5][X][Y];
          popG[7][(X + 1) % Nx]     [Y - 1] = popG[8][X][Y];
          popG[8][(X - 1 + Nx) % Nx][Y + 1] = popG[7][X][Y];
        }
      }
    }
    // Standard bounce-back at the top wall
      // Bounce-back for Navier-Stokes equation
      pop[4][X]                [Ny - 2] = pop[3][X][Ny - 1];
      pop[6][(X - 1 + Nx) % Nx][Ny - 2] = pop[5][X][Ny - 1];
      pop[7][(X + 1) % Nx]     [Ny - 2] = pop[8][X][Ny - 1];
      // Bounce-back for advection-diffusion equation
      if(time > tDiffOn) {
        popG[4][X]                [Ny - 2] = popG[3][X][Ny - 1];
        popG[6][(X - 1 + Nx) % Nx][Ny - 2] = popG[5][X][Ny - 1];
        popG[7][(X + 1) % Nx]     [Ny - 2] = popG[8][X][Ny - 1];
      }
    // Velocity bounce-back at the bottom "wall"
      // Bounce-back for Navier-Stokes equation
      pop[3][X]                [1] = pop[4][X][0];
      pop[5][(X + 1) % Nx]     [1] = pop[6][X][0] + 6 * weight[6] * density[X][1] * velWall;
      pop[8][(X - 1 + Nx) % Nx][1] = pop[7][X][0] - 6 * weight[7] * density[X][1] * velWall;
      // Bounce-back for advection-diffusion equation
      if(time > tDiffOn) {
        popG[3][X]                [1] = popG[4][X][0];
        popG[5][(X + 1) % Nx]     [1] = popG[6][X][0] + 6 * weight[6] * conc[X][1] * velWall;
        popG[8][(X - 1 + Nx) % Nx][1] = popG[7][X][0] - 6 * weight[7] * conc[X][1] * velWall;
      }
  }
  /// Compute fluid density and velocity
  // The fluid density and velocity are obtained from the populations.
  momenta();

  return;
}

/// **********************************
/// COMPUTE FLUID DENSITY AND VELOCITY
/// **********************************

// This function computes the fluid density and velocity from the populations.
// The velocity correction due to body force is not included here.
// It must be taken into account whenever the physical velocity is required.

void momenta() {
  // Density and velocity
  for(int X = 0; X < Nx; ++X) {
    for(int Y = 1; Y < Ny - 1; ++Y) {
      if(!obstacle[X][Y]) {
        density[X][Y] = pop[0][X][Y] + pop[1][X][Y] + pop[2][X][Y] + pop[3][X][Y] + pop[4][X][Y] + pop[5][X][Y] + pop[6][X][Y] + pop[7][X][Y] + pop[8][X][Y];
        velocity_x[X][Y] = (pop[1][X][Y] - pop[2][X][Y] + pop[5][X][Y] - pop[6][X][Y] + pop[7][X][Y] - pop[8][X][Y]) / density[X][Y];
        velocity_y[X][Y] = (pop[3][X][Y] - pop[4][X][Y] + pop[5][X][Y] - pop[6][X][Y] - pop[7][X][Y] + pop[8][X][Y]) / density[X][Y];
      }
    }
  }
  // Concentration
  for(int X = 1; X < Nx - 1; ++X) {
    for(int Y = 1; Y < Ny - 1; ++Y) {
      if(!obstacle[X][Y]) {
        conc[X][Y] = popG[0][X][Y] + popG[1][X][Y] + popG[2][X][Y] + popG[3][X][Y] + popG[4][X][Y] + popG[5][X][Y] + popG[6][X][Y] + popG[7][X][Y] + popG[8][X][Y];
      }
    }
  }

  return;
}

/// **************
/// REACTION MODEL
/// **************

void reaction() {
  /// Reset reaction profile.
  for(int Y = 0; Y < Ny; ++Y) {
    reaction_profile[Y] = 0.;
  }
  /// Calculate reactions.
  reacted = 0.;
  for(int X = 0; X < Nx; ++X) {
    for(int Y = 0; Y < Ny; ++Y) {
      if(catalyst[X][Y]) {
        reacted += conc[X][Y] * reactionRate;
        reaction_profile[Y] += conc[X][Y] * reactionRate;
        conc[X][Y] *= (1 - reactionRate);
        for(int i = 0; i < 9; ++i) {
          popG[i][X][Y] *= (1 - reactionRate);
        }
      }
    }
  }
  
  return;
}

/// *****************************
/// WRITE FLUID STATE TO VTK FILE
/// *****************************

// The fluid state is writen to a VTK file at each tVTK step.
// The following data is written:
// - Obstacles
// - Catalyst
// - Density difference (density - 1)
// - Concentration
// - Velocity

void write_fluid_vtk(int time) {
  /// Create filename.
  stringstream output_filename;
  output_filename << "vtk_fluid/fluid_t" << time << ".vtk";
  ofstream output_file;
  /// Open file.
  output_file.open(output_filename.str().c_str());
  /// Write VTK header.
  output_file << "# vtk DataFile Version 3.0\n";
  output_file << "fluid_state\n";
  output_file << "ASCII\n";
  output_file << "DATASET RECTILINEAR_GRID\n";
  output_file << "DIMENSIONS " << Nx << " " << Ny - 1 << " 1" << "\n";
  output_file << "X_COORDINATES " << Nx << " float\n";
  for(int X = 0; X < Nx; ++X) {
    output_file << X + 0.5 << " ";
  }
  output_file << "\n";
  output_file << "Y_COORDINATES " << Ny - 1 << " float\n";
  for(int Y = 1; Y < Ny; ++Y) {
    output_file << Y + 0.5 << " ";
  }
  output_file << "\n";
  output_file << "Z_COORDINATES " << 1 << " float\n";
  output_file << 0 << "\n";
  output_file << "POINT_DATA " << Nx * (Ny - 1) << "\n";
  /// Write obstacle sites.
  output_file << "SCALARS obstacles float 1\n";
  output_file << "LOOKUP_TABLE default\n";
  for(int Y = 1; Y < Ny; ++Y) {
    for(int X = 0; X < Nx; ++X) {
      output_file << obstacle[X][Y] << "\n";
    }
  }
  /// Write catalyst sites.
  output_file << "SCALARS catalyst float 1\n";
  output_file << "LOOKUP_TABLE default\n";
  for(int Y = 1; Y < Ny; ++Y) {
    for(int X = 0; X < Nx; ++X) {
      output_file << catalyst[X][Y] << "\n";
    }
  }
  /// Write density difference
  output_file << "SCALARS density_difference float 1\n";
  output_file << "LOOKUP_TABLE default\n";
  for(int Y = 1; Y < Ny; ++Y) {
    for(int X = 0; X < Nx; ++X) {
      output_file << (density[X][Y] - 1.) * (1 - obstacle[X][Y]) << "\n";
    }
  }
  /// Write concentration
  output_file << "SCALARS concentration float 1\n";
  output_file << "LOOKUP_TABLE default\n";
  for(int Y = 1; Y < Ny; ++Y) {
    for(int X = 0; X < Nx; ++X) {
      output_file << conc[X][Y] * (1 - obstacle[X][Y])<< "\n";
    }
  }
  /// Write velocity
  output_file << "VECTORS velocity_vector float\n";
  for(int Y = 1; Y < Ny; ++Y) {
    for(int X = 0; X < Nx; ++X) {
      output_file << (velocity_x[X][Y] + 0.5 * (force_x[X][Y] + gravity) / density[X][Y]) * (1 - obstacle[X][Y]) << " " << (velocity_y[X][Y] + 0.5 * (force_y[X][Y]) / density[X][Y]) * (1 - obstacle[X][Y]) << " 0\n";
    }
  }
  /// Close file
  output_file.close();
  
  return;
}

/// ****************
/// WRITE FLUID DATA
/// ****************

// The fluid data is writen to the disk every tData step.
// The following data is written:
// TODO

void write_fluid_data(int time) {
  /// Open file.
  ofstream output_file;
  output_file.open("Data.dat", ios::app);
  /// Write header first time function is called.
  if(time == tData) {
    output_file << "time vel_channel vel_finger vel_lid reacted flux_finger_advection flux_finger_diffusion flux_inlet_advection flux_inlet_diffusion flux_outlet_advection mass_total mass_finger" << endl;
  }
  /// Calculate average channel velocity.
  double vel_channel = 0.;
  int counter = 0;
  for(int X = 0; X < Nx; ++X) {
    for(int Y = 0; Y < Ny - fingerDepth - 1; ++Y) {
      vel_channel += (velocity_x[X][Y] + 0.5 * (force_x[X][Y] + gravity) / density[X][Y]);
      counter++;
    }
  }
  vel_channel /= counter;
  /// Calculate average finger velocity.
  double vel_finger = 0.;
  counter = 0;
  for(int X = 0; X < Nx; ++X) {
    for(int Y = Ny - fingerDepth - 1; Y < Ny; ++Y) {
      if(!obstacle[X][Y]) {
        const double veloX = (velocity_x[X][Y] + 0.5 * (force_x[X][Y] + gravity) / density[X][Y]);
        const double veloY = (velocity_y[X][Y] + 0.5 * (force_y[X][Y]) / density[X][Y]);
        vel_finger += sqrt(veloX * veloX + veloY * veloY);
        counter++;
      }
    }
  }
  vel_finger /= counter;
  /// Calculate average lid velocity (velocity of the fluid layer defining the finger opening).
  const int length_without_finger = (Nx - fingerWidth) / 2;
  double vel_lid = 0.;
  for(int X = length_without_finger; X < Nx - length_without_finger; ++X) {
    vel_lid += (velocity_x[X][Ny - fingerDepth - 2] + 0.5 * (force_x[X][Ny - fingerDepth - 2] + gravity) / density[X][Ny - fingerDepth - 2]);
  }
  vel_lid /= fingerWidth;  
  /// Calculate advective concentration flux into finger.
  double flux_finger_advection = 0.;
  const int Y_finger = Ny - fingerDepth - 1;
  for(int X = 0; X < Nx; ++X) {
    flux_finger_advection += (conc[X][Y_finger] * (velocity_y[X][Y_finger] + 0.5 * (force_y[X][Y_finger]) / density[X][Y_finger]));
  }
  /// Calculate diffusive concentration flux into finger.
  double flux_finger_diffusion = 0.;
  const double diffusion_constant = (tauG - 0.5) / 3.;
  for(int X = 0; X < Nx; ++X) {
    flux_finger_diffusion += -0.5 * diffusion_constant * (conc[X][Y_finger + 2] - conc[X][Y_finger]);
  }
  /// Calculate diffusive inlet concentration flux.
  double flux_inlet_diffusion = 0.;
  for(int Y = 0; Y < Ny; ++Y) {
    flux_inlet_diffusion += -diffusion_constant * (conc[2][Y] - conc[1][Y]);
  }
  /// Calculate advective inlet concentration flux.
  double flux_inlet_advection = 0.;
  const int X_inlet = 1;
  for(int Y = 0; Y < Ny; ++Y) {
    flux_inlet_advection += (conc[X_inlet][Y] * (velocity_x[X_inlet][Y] + 0.5 * (force_x[X_inlet][Y] + gravity) / density[X_inlet][Y]));
  }
  /// Calculate advective outlet concentration flux.
  double flux_outlet_advection = 0.;
  const int X_outlet = Nx - 2;
  for(int Y = 0; Y < Ny; ++Y) {
    flux_outlet_advection += (conc[X_outlet][Y] * (velocity_x[X_outlet][Y] + 0.5 * (force_x[X_outlet][Y] + gravity) / density[X_outlet][Y]));
  }
  /// Calculate total mass in system.
  double mass_total = 0.;
  for(int X = 0; X < Nx; ++X) {
    for(int Y = 0; Y < Ny; ++Y) {
      if(!obstacle[X][Y]) {
        mass_total += conc[X][Y];
      }
    }
  }
  /// Calculate mass in finger.
  double mass_finger = 0.;
  for(int X = 0; X < Nx; ++X) {
    for(int Y = Ny - fingerDepth - 1; Y < Ny; ++Y) {
      if(!obstacle[X][Y]) {
        mass_finger += conc[X][Y];
      }
    }
  }
  /// Update observables for convergence.
  convergence_mass_old = convergence_mass_new;
  convergence_mass_new = mass_finger;
  convergence_vel_old = convergence_vel_new;
  convergence_vel_new = vel_finger;
  /// Write data.
  output_file
    << time << " "
    << vel_channel << " "
    << vel_finger << " "
    << vel_lid << " "
    << reacted << " "
    << flux_finger_advection << " "
    << flux_finger_diffusion << " "
    << flux_inlet_advection << " "
    << flux_inlet_diffusion << " "
    << flux_outlet_advection << " "
    << mass_total << " "
    << mass_finger << endl;
  /// Close file.
  output_file.close();
  
  return;
}

/// **************************
/// CHECK CONVERGENCE CRITERIA
/// **************************

// It is checked whether the simulation has converged.
// The criteria are that between t-tData and t both the
// average velocity and average mass in the finger do
// not change more than 1e-9 * tData (relative).

void check_convergence(int time) {
  bool converged = false;
  /// Check for convergence.
  const double convergence_mass_change = (convergence_mass_new - convergence_mass_old) / convergence_mass_old;
  const double convergence_vel_change = (convergence_vel_new - convergence_vel_old) / convergence_vel_old;
  if(convergence_mass_change < 1e-9 * tData && convergence_vel_change < 1e-9 * tData) {
    converged = true;
  }
  /// Instruct simulation outcome.
  if(converged) {
    write_reaction_profile();
    cout << "Simulation complete with convergence reached!" << endl;
    exit(0);
  }
  else {
    return;
  }
}

/// **********************
/// WRITE REACTION PROFILE
/// **********************

// The reaction profile is written to the disk.

void write_reaction_profile() {
  /// Open file.
  ofstream output_file;
  output_file.open("ReactionProfile.dat", ios::app);
  /// Write header.
  output_file << "Y reacted" << endl;
  /// Write reaction profile.
  for(int Y = 0; Y < Ny; ++Y) {
    output_file << Y << " " << reaction_profile[Y] << endl;
  }
  /// Close file.
  output_file.close();
  return;
}
