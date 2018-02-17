#ifndef FLUIDSIM_H
#define FLUIDSIM_H

#include "array2.h"
#include "vec.h"
#include "pcgsolver/sparse_matrix.h"
#include "pcgsolver/pcg_solver.h"
#include "rigidbody.h"

#include <vector>

void extrapolate(Array2f& grid, Array2c& valid);

class FluidSim {

public:
   void initialize(float width, int ni_, int nj_);
   void set_boundary(float (*phi)(const Vec2f&));
   void advance(float dt);

   //Grid dimensions
   int ni,nj;
   float dx;

   //Fluid velocity
   Array2f u, v;
   Array2f temp_u, temp_v;

   //Static geometry representation
   Array2f nodal_solid_phi;

   RigidBody* rbd;
   RigidGeometry* rigidgeom;
   Array2f nodal_rigid_phi;
   Array2f rigid_u_weights, rigid_v_weights;
   float rigid_u_mass, rigid_v_mass;

   Array2f solid_u, solid_v;

   //Data for pressure solve and extrapolation
   Array2c u_valid, v_valid;
   Array2f liquid_phi; //extracted from particles
   Array2f u_weights, v_weights; //solid v.s. fluid face weights

   //Data for viscosity solve
   Array2f u_vol, v_vol, c_vol, n_vol;
   Array2f viscosity;

   std::vector<Vec2f> particles; //For marker particle simulation
   std::vector<Vec2f> particles_velocity;
   float particle_radius;
   Array2f sum;

   //Data arrays for extrapolation
   Array2c valid, old_valid;

   //Viscosity solver data
   PCGSolver<double> solver;
   SparseMatrixd matrix;
   std::vector<double> rhs;
   std::vector<double> pressure;
   std::vector<double> velocities;

   Vec2f get_velocity(const Vec2f& position);
   Vec2f get_solid_velocity(const Vec2f& position);
   void add_particle(const Vec2f& position);

   Vec2f trace_rk3(const Vec2f& position, float dt);

   void advect_particles(float dt);

   void accumulate(Array2f &accum, float q, int i, int j, float fx, float fy);
   void transfer_to_grid();
   void update_from_grid();

   void compute_phi();
   void extrapolate_phi();

   float cfl();

   //fluid velocity operations
   void advect(float dt);
   void add_force(float dt);

   void apply_projection(float dt);
   void compute_pressure_weights();

   void solve_pressure(float dt);

   void constrain_velocity();
   void process_collisions();
   void update_rigid_body_grids();
   void recompute_solid_velocity();
};

#endif
