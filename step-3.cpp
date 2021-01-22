// Translate this file with
//
// g++ -O3 assignment-code.cpp -o assignment-code
//
// Run it with
//
// ./demo-code
//
// There should be a result.pvd file that you can open with Paraview.
// Sometimes, Paraview requires to select the representation "Point Gaussian"
// to see something meaningful.
//
// (C) 2018-2020 Tobias Weinzierl

#include <fstream>
#include <sstream>
#include <iostream>
#include <string>
#include <math.h>
#include <limits>
#include <iomanip>



#include <cmath>


double t          = 0;
double tFinal     = 0;
double tPlot      = 0;
double tPlotDelta = 0;

int NumberOfBodies = 0;

// Collision constant
double c;

// array to keep track of merged particles
bool *alive;

/**
 * Pointer to pointers. Each pointer in turn points to three coordinates, i.e.
 * each pointer represents one molecule/particle/body.
 */
double** x;
// backups of positions for rk(2)
double** x_back;

/**
 * Equivalent to x storing the velocities.
 */
double** v;

// backup of velocities for rk(2)
double** v_back;

/**
 * One mass entry per molecule/particle.
 */
double*  mass;

/**
 * Global time step size used.
 */
double   timeStepSize = 0.0;

/**
 * Maximum velocity of all particles.
 */
double   maxV;

/**
 * Minimum distance between two elements.
 */
double   minDx;


/**
 * Set up scenario from the command line.
 *
 * If you need additional helper data structures, you can
 * initialise them here. Alternatively, you can introduce a
 * totally new function to initialise additional data fields and
 * call this new function from main after setUp(). Either way is
 * fine.
 *
 * This operation's semantics is not to be changed in the assignment.
 */
void setUp(int argc, char** argv) {
  NumberOfBodies = (argc-4) / 7;

  x    = new double*[NumberOfBodies];
  v    = new double*[NumberOfBodies];
  mass = new double [NumberOfBodies];

  // init c
  c = 0.01/NumberOfBodies;

  // init alive
  alive = new bool[NumberOfBodies];
  x_back = new double*[NumberOfBodies];
  v_back = new double*[NumberOfBodies];

  int readArgument = 1;

  tPlotDelta   = std::stof(argv[readArgument]); readArgument++;
  tFinal       = std::stof(argv[readArgument]); readArgument++;
  timeStepSize = std::stof(argv[readArgument]); readArgument++;

  for (int i=0; i<NumberOfBodies; i++) {

    // init alive to True for all particles at start
    alive[i] = true;

    x[i] = new double[3];
    v[i] = new double[3];
    x_back[i] = new double[3];
    v_back[i] = new double[3];

    x[i][0] = std::stof(argv[readArgument]); readArgument++;
    x[i][1] = std::stof(argv[readArgument]); readArgument++;
    x[i][2] = std::stof(argv[readArgument]); readArgument++;

    v[i][0] = std::stof(argv[readArgument]); readArgument++;
    v[i][1] = std::stof(argv[readArgument]); readArgument++;
    v[i][2] = std::stof(argv[readArgument]); readArgument++;

    mass[i] = std::stof(argv[readArgument]); readArgument++;

    if (mass[i]<=0.0 ) {
      std::cerr << "invalid mass for body " << i << std::endl;
      exit(-2);
    }
  }

  std::cout << "created setup with " << NumberOfBodies << " bodies" << std::endl;
  
  if (tPlotDelta<=0.0) {
    std::cout << "plotting switched off" << std::endl;
    tPlot = tFinal + 1.0;
  }
  else {
    std::cout << "plot initial setup plus every " << tPlotDelta << " time units" << std::endl;
    tPlot = 0.0;
  }
}


std::ofstream videoFile;


/**
 * This operation is not to be changed in the assignment.
 */
void openParaviewVideoFile() {
  videoFile.open( "result.pvd" );
  videoFile << "<?xml version=\"1.0\"?>" << std::endl
            << "<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\" compressor=\"vtkZLibDataCompressor\">" << std::endl
            << "<Collection>";
}





/**
 * This operation is not to be changed in the assignment.
 */
void closeParaviewVideoFile() {
  videoFile << "</Collection>"
            << "</VTKFile>" << std::endl;
}


/**
 * The file format is documented at http://www.vtk.org/wp-content/uploads/2015/04/file-formats.pdf
 *
 * This operation is not to be changed in the assignment.
 */
void printParaviewSnapshot() {
  static int counter = -1;
  counter++;
  std::stringstream filename;
  filename << "result-" << counter <<  ".vtp";
  std::ofstream out( filename.str().c_str() );
  out << "<VTKFile type=\"PolyData\" >" << std::endl
      << "<PolyData>" << std::endl
      << " <Piece NumberOfPoints=\"" << NumberOfBodies << "\">" << std::endl
      << "  <Points>" << std::endl
      << "   <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">";
//      << "   <DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">";

  for (int i=0; i<NumberOfBodies; i++) {
    out << x[i][0]
        << " "
        << x[i][1]
        << " "
        << x[i][2]
        << " ";
  }

  out << "   </DataArray>" << std::endl
      << "  </Points>" << std::endl
      << " </Piece>" << std::endl
      << "</PolyData>" << std::endl
      << "</VTKFile>"  << std::endl;

  videoFile << "<DataSet timestep=\"" << counter << "\" group=\"\" part=\"0\" file=\"" << filename.str() << "\"/>" << std::endl;
}



// This function carries out the rk(2) prediction shot using explicit euler
void future_shot(){
  double* force0 = new double[NumberOfBodies];
  double* force1 = new double[NumberOfBodies];
  double* force2 = new double[NumberOfBodies];
  double* distances = new double[NumberOfBodies];
  for (int j = 0;j<NumberOfBodies;j++){
      force0[j] = 0;
      force1[j] = 0;
      force2[j] = 0;
  }

  for (int j = 0; j < NumberOfBodies; j++){
      for (int i = j+1; i < NumberOfBodies; i++) {
        distances[i] = sqrt((x[j][0]-x[i][0]) * (x[j][0]-x[i][0]) +(x[j][1]-x[i][1]) * (x[j][1]-x[i][1]) + (x[j][2]-x[i][2]) * (x[j][2]-x[i][2]));
        // x,y,z forces acting on particle j from i
        const double f0 = (x[i][0]-x[j][0]) * mass[i]*mass[j] / distances[i] / distances[i] / distances[i] ;
        const double f1 = (x[i][1]-x[j][1]) * mass[i]*mass[j] / distances[i] / distances[i] / distances[i] ;
        const double f2 = (x[i][2]-x[j][2]) * mass[i]*mass[j] / distances[i] / distances[i] / distances[i] ;
        force0[j] += f0;
        force1[j] += f1;
        force2[j] += f2;
        // x,y,z forces acting on i from j
        force0[i] += -f0;
        force1[i] += -f1;
        force2[i] += -f2;
      }
  }
      #pragma omp simd 
      for (int j = 0;j<NumberOfBodies;j++){
          x[j][0] = x[j][0] + timeStepSize/2 * v[j][0];
          x[j][1] = x[j][1] + timeStepSize/2 * v[j][1];
          x[j][2] = x[j][2] + timeStepSize/2 * v[j][2];
      }
      #pragma omp simd
      for(int j = 0;j<NumberOfBodies;j++){
          v[j][0] = v[j][0] + timeStepSize/2 * force0[j] / mass[j];
          v[j][1] = v[j][1] + timeStepSize/2 * force1[j] / mass[j];
          v[j][2] = v[j][2] + timeStepSize/2 * force2[j] / mass[j];
      }
      delete[] force0;
      delete[] force1;
      delete[] force2;
      delete[] distances; 
}




/**
 * This is the main operation you should change in the assignment. You might
 * want to add a few more variables or helper functions, but this is where the
 * magic happens.
 */
void updateBody() {
  // backup current y values for rk(2)
  for (int j = 0;j<NumberOfBodies;j++){
      x_back[j][0] = x[j][0];
      x_back[j][1] = x[j][1];
      x_back[j][2] = x[j][2];
      v_back[j][0] = v[j][0];
      v_back[j][1] = v[j][1];
      v_back[j][2] = v[j][2];
  }
  // perform future shot store y_hat in x and v
  future_shot();
  maxV   = 0.0;
  minDx  = std::numeric_limits<double>::max();

  // force0 = force along x direction
  // force1 = force along y direction
  // force2 = force along z direction
  double* force0 = new double[NumberOfBodies];
  double* force1 = new double[NumberOfBodies];
  double* force2 = new double[NumberOfBodies];
  double* distances = new double[NumberOfBodies];
  double* velocities = new double[NumberOfBodies];
  for (int j = 0;j<NumberOfBodies;j++){
      force0[j] = 0;
      force1[j] = 0;
      force2[j] = 0;
  }

  for (int j = 0; j < NumberOfBodies; j++){
      for (int i = j+1; i < NumberOfBodies; i++) {
        // calculate forces using y_hat as this is the derivative for velocities update
        distances[i] = sqrt((x[j][0]-x[i][0]) * (x[j][0]-x[i][0]) +(x[j][1]-x[i][1]) * (x[j][1]-x[i][1]) + (x[j][2]-x[i][2]) * (x[j][2]-x[i][2]));
        // x,y,z forces acting on particle j from i
        const double f0 = (x[i][0]-x[j][0]) * mass[i]*mass[j] / distances[i] / distances[i] / distances[i] ;
        const double f1 = (x[i][1]-x[j][1]) * mass[i]*mass[j] / distances[i] / distances[i] / distances[i] ;
        const double f2 = (x[i][2]-x[j][2]) * mass[i]*mass[j] / distances[i] / distances[i] / distances[i] ;
        force0[j] += f0;
        force1[j] += f1;
        force2[j] += f2;
        // x,y,z forces acting on i from j
        force0[i] += -f0;
        force1[i] += -f1;
        force2[i] += -f2;
      }
  }
      #pragma omp simd 
      for (int j = 0;j<NumberOfBodies;j++){
          x[j][0] = x_back[j][0] + timeStepSize * v[j][0];
          x[j][1] = x_back[j][1] + timeStepSize * v[j][1];
          x[j][2] = x_back[j][2] + timeStepSize * v[j][2];
      }
      #pragma omp simd
      for(int j = 0;j<NumberOfBodies;j++){
          v[j][0] = v_back[j][0] + timeStepSize * force0[j] / mass[j];
          v[j][1] = v_back[j][1] + timeStepSize * force1[j] / mass[j];
          v[j][2] = v_back[j][2] + timeStepSize * force2[j] / mass[j];
      }
      for(int j = 0; j<NumberOfBodies;j++){
          for (int i = j+1;i<NumberOfBodies;i++){
            distances[i] = sqrt((x[j][0]-x[i][0]) * (x[j][0]-x[i][0]) +(x[j][1]-x[i][1]) * (x[j][1]-x[i][1]) + (x[j][2]-x[i][2]) * (x[j][2]-x[i][2]));
          }
          for (int i = j+1;i<NumberOfBodies;i++){
            minDx = std::min(minDx,distances[i]);
          }
      }
      for (int j = 0;j<NumberOfBodies;j++){
        velocities[j] = std::sqrt( v[j][0]*v[j][0] + v[j][1]*v[j][1] + v[j][2]*v[j][2]);
      }

  for (int j = 0;j<NumberOfBodies;j++){
      maxV = std::max(maxV, velocities[j]);
  }
  t += timeStepSize;
   
  int num_cols = 0;
  int elements_rem = 0;
  do {
      num_cols = 0;
      for(int i = 0;i<NumberOfBodies;i++){
          for(int j = i+1;j<NumberOfBodies;j++){
              distances[j] = sqrt((x[i][0]-x[j][0]) * (x[i][0]-x[j][0]) +
                                   (x[i][1]-x[j][1]) * (x[i][1]-x[j][1]) +
                                   (x[i][2]-x[j][2]) * (x[i][2]-x[j][2]));
          }
          for(int j = i+1;j<NumberOfBodies;j++){
              if(alive[i] && alive[j]){ 
                    if(distances[j] <= c*(mass[i] + mass[j])) {
                        for (int k = 0;k<3;k++){
                        // update velocity of merged particle
                        v[i][k] = ((mass[i]*v[i][k])/(mass[i]+mass[j]))+((mass[j]*v[j][k])/(mass[i]+mass[j]));
                        // update position of merged particle
                        x[i][k] = ((mass[i]*x[i][k])+(mass[j]*x[j][k]))/(mass[i] + mass[j]);
                        }
                        // update mass of merged particle
                        mass[i] = mass[i] + mass[j];
                        alive[j] = false; // J has been merged so is no longer valid
                        num_cols += 1;
                        elements_rem += 1;
                    }
              }
          }
      }
  }while(num_cols > 0);
  NumberOfBodies -= elements_rem;
  int tmp_new_index = 0;
  for (int i = 0; i<(NumberOfBodies+elements_rem);i++){
      if(alive[i]){
        alive[i] = false;
        mass[tmp_new_index] = mass[i];
        for (int j = 0;j<3;j++){
            v[tmp_new_index][j] = v[i][j];
            x[tmp_new_index][j] = x[i][j];
        }
        alive[tmp_new_index] = true;
        tmp_new_index += 1;
      }
  }




  delete[] force0;
  delete[] force1;
  delete[] force2;
  delete[] distances; 
  delete[] velocities; 
}


/**
 * Main routine.
 *
 * No major changes in assignment. You can add a few initialisation
 * or stuff if you feel the need to do so. But keep in mind that you
 * may not alter what the program plots to the terminal.
 */
int main(int argc, char** argv) {
  if (argc==1) {
    std::cerr << "usage: " + std::string(argv[0]) + " snapshot final-time dt objects" << std::endl
              << "  snapshot        interval after how many time units to plot. Use 0 to switch off plotting" << std::endl
              << "  final-time      simulated time (greater 0)" << std::endl
              << "  dt              time step size (greater 0)" << std::endl
              << std::endl
              << "Examples:" << std::endl
              << "0.01  100.0  0.001    0.0 0.0 0.0  1.0 0.0 0.0  1.0 \t One body moving form the coordinate system's centre along x axis with speed 1" << std::endl
              << "0.01  100.0  0.001    0.0 0.0 0.0  1.0 0.0 0.0  1.0     0.0 1.0 0.0  1.0 0.0 0.0  1.0  \t One spiralling around the other one" << std::endl
              << "0.01  100.0  0.001    3.0 0.0 0.0  0.0 1.0 0.0  0.4     0.0 0.0 0.0  0.0 0.0 0.0  0.2     2.0 0.0 0.0  0.0 0.0 0.0  1.0 \t Three body setup from first lecture" << std::endl
              << "0.01  100.0  0.001    3.0 0.0 0.0  0.0 1.0 0.0  0.4     0.0 0.0 0.0  0.0 0.0 0.0  0.2     2.0 0.0 0.0  0.0 0.0 0.0  1.0     2.0 1.0 0.0  0.0 0.0 0.0  1.0     2.0 0.0 1.0  0.0 0.0 0.0  1.0 \t Five body setup" << std::endl
              << std::endl
              << "In this naive code, only the first body moves" << std::endl;

    return -1;
  }
  else if ( (argc-4)%7!=0 ) {
    std::cerr << "error in arguments: each planet is given by seven entries (position, velocity, mass)" << std::endl;
    std::cerr << "got " << argc << " arguments (three of them are reserved)" << std::endl;
    std::cerr << "run without arguments for usage instruction" << std::endl;
    return -2;
  }

  std::cout << std::setprecision(15);

  setUp(argc,argv);

  openParaviewVideoFile();

  int snapshotCounter = 0;
  if (t > tPlot) {
    printParaviewSnapshot();
    std::cout << "plotted initial setup" << std::endl;
    tPlot = tPlotDelta;
  }

  int timeStepCounter = 0;
  while (t<=tFinal) {
    updateBody();
    timeStepCounter++;
    if (t >= tPlot) {
      printParaviewSnapshot();
      std::cout << "plot next snapshot"
    		    << ",\t time step=" << timeStepCounter
    		    << ",\t t="         << t
				<< ",\t dt="        << timeStepSize
				<< ",\t v_max="     << maxV
				<< ",\t dx_min="    << minDx
				<< std::endl;

      tPlot += tPlotDelta;
    }
  }

  std::cout << "Number of remaining objects: " << NumberOfBodies << std::endl;
  std::cout << "Position of first remaining object: " << x[0][0] << ", " << x[0][1] << ", " << x[0][2] << std::endl;

  closeParaviewVideoFile();

  return 0;
}
