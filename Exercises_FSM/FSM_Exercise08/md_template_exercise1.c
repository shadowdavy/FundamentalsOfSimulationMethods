/*
This is a template version for Problem Set 8.
All parts that need to be completed are marked with "TODO"
*/


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <unistd.h>
#include <math.h>
#include <time.h>
#define MAXPART 10000

// Modified data structure for neighbour lists
typedef struct {
    double pos[3];
    double vel[3];
    double acc[3];
    double acc_prev[3];
    double pot;
    int n_neighbors;        // length of neighbor list
    int neighbors[MAXPART]; // list of neighbors (with index larger than current particle)
} particle;


// auxiliary function to create a Gaussian random deviate
double gaussian_rnd(void) {
    double x, y;

    do {
      x = drand48();
      y = drand48();
    } while(x == 0);

    return sqrt(-2.0 * log(x)) * cos(2 * M_PI * y);
}


// This function initializes our particle set.
void initialize(particle *p, int nperdim, double boxsize, double temp, double epsilon_in_Kelvin) {
    int i, j, k, n = 0;
    double vel_rms = sqrt(temp / epsilon_in_Kelvin);

    for(i = 0; i < nperdim; i++) {
        for(j = 0; j < nperdim; j++) {
            for(k = 0; k < nperdim; k++) {
              p[n].pos[0] = (i + 0.5) * boxsize / nperdim;
              p[n].pos[1] = (j + 0.5) * boxsize / nperdim;
              p[n].pos[2] = (k + 0.5) * boxsize / nperdim;

              p[n].vel[0] = vel_rms * gaussian_rnd();
              p[n].vel[1] = vel_rms * gaussian_rnd();
              p[n].vel[2] = vel_rms * gaussian_rnd();

              n++;
            }
        }
    }
}


// This function updates the velocities by applying the accelerations for the given time interval.
void kick(particle * p, int ntot, double dt) {
    int i, k;

    for(i = 0; i < ntot; i++) {
        for(k = 0; k < 3; k++) {
            p[i].vel[k] += p[i].acc[k] * dt;
        }
    }
}


// This function drifts the particles with their velocities for the given time interval.
// Afterwards, the particles are mapped periodically back to the box if needed.
void drift(particle * p, int ntot, double boxsize, double dt) {
    int i, k;

    for(i = 0; i < ntot; i++) {
        for(k = 0; k < 3; k++) {
            p[i].pos[k] += p[i].vel[k] * dt;

            if (p[i].pos[k] >= boxsize) {
                p[i].pos[k] -= boxsize;
            }

            if (p[i].pos[k] < 0) {
                p[i].pos[k] += boxsize;
            }
        }
    }
}


// This function calculates the potentials and forces for all particles. For simplicity,
// we do this by going through all particle pairs i-j, and then adding the contributions both to i and j.
void calc_forces(particle * p, int ntot, double boxsize, double rcut)
{
    int i, j, k, n;
    double rcut2 = rcut * rcut;
    double r2, r, r6, r12, dr[3], acc[3], pot;

    // first, set all the accelerations and potentials to zero
    for (i = 0; i < ntot; i++) {
       p[i].pot = 0;

       for (k = 0; k < 3; k++) {
           p[i].acc[k] = 0;
       }
    }

    // sum over all distinct pairs
    for (i = 0; i < ntot; i++) {
        
        for (int jp = 0; jp < p[i].n_neighbors; jp++) {        // TODO: Neighbour list modification has to be added here
            j = p[i].neighbors[jp];
            for (k = 0, r2 = 0; k < 3; k++) {
                dr[k] = p[i].pos[k] - p[j].pos[k];

                if (dr[k] > 0.5 * boxsize) {
                    dr[k] -= boxsize;
                }

                if (dr[k] < -0.5 * boxsize) {
                    dr[k] += boxsize;
                }
                r2 += dr[k] * dr[k];
             }

            if (r2 < rcut2) {
                r = sqrt(r2);
                r6 = r2 * r2 * r2;
                r12 = r6 * r6;

                // now calculate the Lennard-Jones potential for the pair
                pot = 4.0 * (1.0 / r12 - 1.0 / r6);

                p[i].pot += pot;
                p[j].pot += pot;

                // now calculate the Lennard-Jones force between the particles
                for (k = 0; k < 3; k++) {
                    acc[k] = 4.0 * (12.0 / r12 - 6.0 / r6) * dr[k] / r2;

                    p[i].acc[k] += acc[k];
                    p[j].acc[k] -= acc[k];
                }
            }
        }
    }
}


// This function calculates the total kinetic and total potential energy, averaged per particle.
// It also returns the instantanous kinetic temperature.
void calc_energies(particle *p, int ntot, double epsilon_in_Kelvin, double *ekin, double *epot, double *temp) {
    int i, k;
    double sum_pot = 0, sum_kin = 0;

    for (i = 0; i < ntot; i++) {
      sum_pot += p[i].pot;

      for(k = 0; k < 3; k++) {
          sum_kin += p[i].vel[k] * p[i].vel[k];
      }
    }

    *ekin = 0.5 * sum_kin / ntot;
    *epot = 0.5 * sum_pot / ntot;

    *temp = 2.0 / 3 * (*ekin) * epsilon_in_Kelvin;
}



// This function rescales the velocities by a prescribed factor 'fac'
void rescale_velocities(particle * p, int ntot, double fac) {
    int i, k;

    for(i = 0; i < ntot; i++) {
      for(k = 0; k < 3; k++) {
          p[i].vel[k] *= fac;
      }
    }
}

// Updates the radial distribution function
void update_rdf(particle *p, int ntot, double boxsize, int *rdf_grid, int n_rdf_grid) {
    int i, j, k, n;
    double r, r2;
    double dr[3];

    // sum over all distinct pairs
    for (i = 0; i < ntot; i++) {
        for (j = i + 1; j < ntot; j++) {

            /* TODO:
            Find the corresponding bin in rdf_grid
            and increase the counter.
            */

        }
    }
}

// Updates the neighbour lists
void update_neighbors(particle *p, int ntot, double boxsize, double rcut) {

    /* TODO:
    Find particles within the cutoff radius
    and add them to the prepared list in the particle struct.
    Consequently, the loop in calc_forces() has to be modified.
    */
    double rl, dr[3];
    for (int i = 0; i < ntot; i++) {
        p[i].n_neighbors = 0;

        for (int j = i + 1; j < ntot; j++) {
            //if (i == j) continue;
            /*dr[0] = p[i].pos[0] - p[j].pos[0];  // for three the calculation is faster than a loop
            dr[1] = p[i].pos[1] - p[j].pos[1];
            dr[2] = p[i].pos[2] - p[j].pos[2];
            rl = (dr[0] * dr[0]) + (dr[1] * dr[1]) + (dr[2] * dr[2]);*/
            for (int k = 0, rl = 0; k < 3; k++) {
                dr[k] = p[i].pos[k] - p[j].pos[k];

                if (dr[k] > 0.5 * boxsize) {
                    dr[k] -= boxsize;
                }

                if (dr[k] < -0.5 * boxsize) {
                    dr[k] += boxsize;
                }
                rl += dr[k] * dr[k];
            }

            if (sqrt(rl) < rcut) { 
                p[i].neighbors[p[i].n_neighbors] = j;
                p[i].n_neighbors += 1;
            }
            
        }
    }

}

/*
 * main driver routine
 */
int main(int argc, char **argv) {

    // Temperature and density can be passed as parameters from the command line 
    // add command line parameter N_1d 3 -> 4
    if (argc != 3) {
        printf("Error! Please call as: %s (temperature in kelvin) (density)\n", argv[0]);
        return 1;
    }

    double epsilon_in_Kelvin = 120.0;               // energy scale in Lennard Jones potential
    double rho = atof(argv[2]);                     // initial density
    int N1d = 8;                                    // particles per dimension
    //int N1d = atof(argv[3]);
    int N = N1d * N1d * N1d;                        // total particle number
    double boxsize = pow(N / rho, 1.0/3.0);         // dimensionless boxsize
    double mean_spacing = boxsize / N1d;            // mean spacing in dimensionless internal units
    double target_temp_in_Kelvin = atof(argv[1]);   // target temperature
    int nsteps = 10000;                             // number of steps to take
    int start_rdf = 5000;                           // steps to start RDF calculation
    double dt = 0.01;                               // timestep size
    double cutoff = 5.0;                            // cutoff
    int n_rdf_grid = 1000;                          // grid size for radial distr. function
    double cutoff_neighbors = 10.0;
    int update_neighbors_freq = 100;

    double ekin, epot, temp;
    int i, step;
    double r, rbin;

    // allocate storage for our particles
    particle *p = malloc(N * sizeof(particle));
    int *rdf_grid = malloc(n_rdf_grid * sizeof(int));          // radial distribution function counts

    // let's initialize the particles
    initialize(p, N1d, boxsize, target_temp_in_Kelvin, epsilon_in_Kelvin);
    update_neighbors(p, N, boxsize, cutoff_neighbors);

    // calculate the forces at t=0
    calc_forces(p, N, boxsize, cutoff);

    // create an output file
    char fname[100];
    sprintf(fname, "output_%d.txt", (int)target_temp_in_Kelvin);
    FILE *fd = fopen(fname, "w");

    sprintf(fname, "output_rdf_%d.txt", (int)target_temp_in_Kelvin);
    FILE *frdf = fopen(fname, "w");

    // measure energies at beginning, and output this to the file and screen
    calc_energies(p, N, epsilon_in_Kelvin, &ekin, &epot, &temp);
    fprintf(fd, "%6d   %10g   %10g   %10g       %10g\n", 0, ekin, epot, ekin + epot, temp);
    printf("%6d   %10g   %10g   %10g       %10g\n", 0, ekin, epot, ekin + epot, temp);

    for (i = 0; i < n_rdf_grid; i++) {
        rdf_grid[i] = 0;
    }

    // now we carry out time integration steps using the leapfrog algorithm
    for (step = 0; step < nsteps; step++) {
        kick(p, N, 0.5 * dt);
        drift(p, N, boxsize, dt);

        if (step>=start_rdf) {
          update_rdf(p, N, boxsize, rdf_grid, n_rdf_grid);     //rdf
        }

        calc_forces(p, N, boxsize, cutoff);
        kick(p, N, 0.5 * dt);

        /* exercise 1 b) turn of velocity rescaling
        // every 100-th step, we print out the coordinates and rescale velocities
        if (step % 100 == 0) {
             double fac = sqrt(target_temp_in_Kelvin / temp);
             rescale_velocities(p, N, fac);
        }*/ 

        if (step % update_neighbors_freq == 0) {
           update_neighbors(p, N, boxsize, cutoff_neighbors);
        }

        // measure energies and output this to the file and screen
        calc_energies(p, N, epsilon_in_Kelvin, &ekin, &epot, &temp);
        fprintf(fd, "%6d   %10g   %10g   %10g       %10g\n", step + 1, ekin, epot, ekin + epot, temp);
        printf("%6d   %10g   %10g   %10g       %10g\n", step + 1, ekin, epot, ekin + epot, temp);
    }


    // TODO: Print radial distribution function data to the output file
    for (i = 1; i < n_rdf_grid; i++) {
        /* TODO:
        So far the array rdf_grid only contains particle counts.
        Compute the radial distribution function and print it to the output file frdf.
        Alternatively, one can print the particle counts and do the postprocessing
        with python.
        */
    }


    fclose(fd);
    fclose(frdf);
    return 0;
}
