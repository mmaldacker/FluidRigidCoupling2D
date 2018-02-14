#include "fluidsim.h"

#include "array2_utils.h"
#include "box2dgeometry.h"
#include "levelset.h"

#include "pcgsolver/sparse_matrix.h"
#include "pcgsolver/pcg_solver.h"
#include <numeric>
#include <functional>
#include <fstream>

float fraction_inside(float phi_left, float phi_right);
void extrapolate(Array2f& grid, Array2c& valid);

float circle_phi(const Vec2f& pos) {
    Vec2f centre(0.5f, 0.75f);
    float rad = 0.1f;
    Vec2f centre1(0.4f, 0.3f);
    float rad1 = 0.15f;
    float phi0 = dist(centre, pos) - rad;
    float phi1 = dist(centre1, pos) - rad1;
    return min(phi0, phi1);
}

void FluidSim::initialize(float width, int ni_, int nj_) {
    ni = ni_;
    nj = nj_;
    dx = width / (float)ni;
    u.resize(ni + 1, nj); temp_u.resize(ni + 1, nj); u_weights.resize(ni + 1, nj); u_valid.resize(ni + 1, nj); u_vol.resize(ni + 1, nj);
    v.resize(ni, nj + 1); temp_v.resize(ni, nj + 1); v_weights.resize(ni, nj + 1); v_valid.resize(ni, nj + 1); v_vol.resize(ni, nj + 1);
    solid_u.resize(ni + 1, nj); solid_v.resize(ni, nj + 1);
    solid_u.set_zero();
    solid_v.set_zero();
    rigid_u_weights.resize(ni + 1, nj); rigid_v_weights.resize(ni, nj + 1);
    rigid_u_weights.set_zero();
    rigid_v_weights.set_zero();
    c_vol.resize(ni, nj);
    n_vol.resize(ni + 1, nj + 1);
    u.set_zero();
    v.set_zero();
    sum.resize(ni+1,nj+1);
    nodal_solid_phi.resize(ni + 1, nj + 1);
    nodal_rigid_phi.resize(ni + 1, nj + 1);
    nodal_rigid_phi.assign(3*dx);
    valid.resize(ni + 1, nj + 1);
    old_valid.resize(ni + 1, nj + 1);
    liquid_phi.resize(ni, nj);
    particle_radius = dx / sqrt(2.0f);
    viscosity.resize(ni, nj);
    viscosity.assign(0.005f);
    rigid_u_mass = 0.0f;
    rigid_v_mass = 0.0f;

    /*
    rigidgeom = new Box2DGeometry(0.3f, 0.2f);
    rbd = new RigidBody(0.4f, *rigidgeom);
    rbd->setCOM(Vec2f(0.5f, 0.65f));
    rbd->setAngle(0.0);
    rbd->setAngle(-(float)M_PI / 2.0f);
    rbd->setAngularMomentum(0.0f);
    rbd->setLinearVelocity(Vec2f(0, 0));
    */

    rbd = NULL;
}

//Initialize the grid-based signed distance field that dictates the position of the solid boundary
void FluidSim::set_boundary(float(*phi)(const Vec2f&)) {

    for (int j = 0; j < nj + 1; ++j) for (int i = 0; i < ni + 1; ++i) {
        Vec2f pos(i*dx, j*dx);
        nodal_solid_phi(i, j) = phi(pos);
    }

}

float FluidSim::cfl() {
    float maxvel = 0;
    for (unsigned int i = 0; i < u.a.size(); ++i)
        maxvel = std::max(maxvel, (float)fabs(u.a[i]));
    for (unsigned int i = 0; i < v.a.size(); ++i)
        maxvel = std::max(maxvel, (float)fabs(v.a[i]));
    return dx / maxvel;
}

//The main fluid simulation step
void FluidSim::advance(float dt) {
    float t = 0;

    while (t < dt) {
        float substep = cfl();
        if (t + substep > dt)
            substep = dt - t;

        //Passively advect particles
        advect_particles(substep);


        //Time integration of rigid bodies
        if (rbd) rbd->advance(substep);

        //Collisions between rbd and static solid boundary
        process_collisions();

        //Recompute the distance fields and face areas
        update_rigid_body_grids();

        //Estimate the liquid signed distance
        compute_phi();

        //Advance the velocity
        advect(substep);
        add_force(substep);

        //apply_viscosity(substep);

        apply_projection(substep);

        //Pressure projection only produces valid velocities in faces with non-zero associated face area.
        //Because the advection step may interpolate from these invalid faces,
        //we must extrapolate velocities from the fluid domain into these zero-area faces.
        extrapolate(u, u_valid);
        extrapolate(v, v_valid);

        //get grid-based approximate solid velocities for use in
        //constrained extrapolation
        recompute_solid_velocity();

        //For extrapolated velocities, replace the normal component with
        //that of the object.
        constrain_velocity();

        t += substep;
    }
}

void FluidSim::add_force(float dt) {

    for (int j = 0; j < nj + 1; ++j) for (int i = 0; i < ni; ++i) {
        v(i, j) -= 0.1f*dt;
    }
}

//For extrapolated points, replace the normal component
//of velocity with the object velocity (in this case zero).
void FluidSim::constrain_velocity() {
    temp_u = u;
    temp_v = v;

    //(At lower grid resolutions, the normal estimate from the signed
    //distance function is poor, so it doesn't work quite as well.
    //An exact normal would do better.)

    //constrain u
    for (int j = 0; j < u.nj; ++j) for (int i = 0; i < u.ni; ++i) {
        if (u_weights(i, j) == 0) {
            //apply constraint
            Vec2f pos(i*dx, (j + 0.5f)*dx);
            Vec2f vel = get_velocity(pos);
            Vec2f solid_vel = get_solid_velocity(pos);

            Vec2f solid_normal(0, 0);
            double solid_phi_val = interpolate_gradient(solid_normal, pos / dx, nodal_solid_phi);
            normalize(solid_normal);

            Vec2f rigid_normal(0, 0);
            double rigid_phi_val = interpolate_gradient(rigid_normal, pos / dx, nodal_rigid_phi);
            normalize(rigid_normal);

            if (solid_phi_val < rigid_phi_val) {
                vel -= dot(vel, solid_normal) * solid_normal;
                vel += dot(solid_vel, solid_normal) * solid_normal;
                temp_u(i, j) = vel[0];
            }
            else {
                vel -= dot(vel, rigid_normal) * rigid_normal;
                vel += dot(solid_vel, rigid_normal) * rigid_normal;
                temp_u(i, j) = vel[0];
            }
        }
    }

    //constrain v
    for (int j = 0; j < v.nj; ++j) for (int i = 0; i < v.ni; ++i) {
        if (v_weights(i, j) == 0) {
            //apply constraint
            Vec2f pos((i + 0.5f)*dx, j*dx);
            Vec2f vel = get_velocity(pos);
            Vec2f solid_vel = get_solid_velocity(pos);

            Vec2f solid_normal(0, 0);
            float solid_phi_val = interpolate_gradient(solid_normal, pos / dx, nodal_solid_phi);
            normalize(solid_normal);

            Vec2f rigid_normal(0, 0);
            float rigid_phi_val = interpolate_gradient(rigid_normal, pos / dx, nodal_rigid_phi);
            normalize(rigid_normal);

            if (solid_phi_val < rigid_phi_val) {
                vel -= dot(vel, solid_normal) * solid_normal;
                vel += dot(solid_vel, solid_normal) * solid_normal;
                temp_v(i, j) = vel[1];
            }
            else {
                vel -= dot(vel, rigid_normal) * rigid_normal;
                vel += dot(solid_vel, rigid_normal) * rigid_normal;
                temp_v(i, j) = vel[1];
            }
        }
    }

    //update
    u = temp_u;
    v = temp_v;

}


void FluidSim::process_collisions() {

    if (!rbd) return;
    //Handle rigid body/rigid body collisions.
    //This pretty hacky since it's not really the focus here.

    //Should really probably find first collision point and back up to that point,
    //rather than trying to project out of walls.

    //Handle velocities of collision first,
    //then process positions afterwards

    //detect all penetrating, non-separating points
    //This is roughly following the Guendelman '05
    //description of impulse response
    std::vector<Vec2f> vertices;
    rbd->get2DVertices(vertices);
    float min_phi = dx;
    int min_ind = -1;
    Vec2f min_normal;
    std::vector<unsigned int> penetrating_points;
    for (unsigned int i = 0; i < vertices.size(); ++i) {

        float phi;
        Vec2f normal;
        Vec2f vertex(vertices[i][0], vertices[i][1]);
        phi = interpolate_normal(normal, vertex, nodal_solid_phi, Vec2f(0, 0), dx);
        if (phi < 0) {
            Vec2f pt_vel = rbd->getPointVelocity(vertex);
            if (dot(pt_vel, normal) < 1e-2) {
                penetrating_points.push_back(i);
            }
            if (phi < min_phi) {
                min_phi = phi;
                min_ind = i;
                min_normal = normal;
            }
        }
    }
    int iteration = 0;


    while (penetrating_points.size() > 0) {
        iteration++;
        //Process collisions on all interfering, non-separating points

        //Take last one in order, should really grab deepest instead
        int p = penetrating_points[0];
        penetrating_points.pop_back();

        Vec2f point = vertices[p];
        Vec2f com;
        rbd->getCOM(com);
        Vec2f rad = point - com;

        Vec2f velocity = rbd->getPointVelocity(point);
        Vec2f normal;
        interpolate_normal(normal, point, nodal_solid_phi, Vec2f(0, 0), dx);

        //Project out normal velocity (inelastic collision)
        float normal_vel = dot(velocity, normal);

        //Compute impulse (go up to 3D since I'm too lazy to work out
        //what the 2D simplification is)
        Vec3f rad3d(rad[0], rad[1], 0);
        Vec3f normal3d(normal[0], normal[1], 0);
        Mat33f K;
        zero(K);
        K(0, 0) = 1 / rbd->getMass(); //why this term at all?
        K(1, 1) = 1 / rbd->getMass(); //why this term at all?
        K(2, 2) = 1 / rbd->getMass();
        Mat33f R = star_matrix(rad3d);
        Mat33f Iinv;
        float imodinv;
        imodinv = rbd->getInvInertiaModulus();
        Iinv(0, 0) = 0; //do I actually need these 2 dimensions' values?
        Iinv(1, 1) = 0;
        Iinv(2, 2) = imodinv;
        K = K + R.transpose()*Iinv*R;
        float j = -normal_vel / dot(normal3d, K*normal3d);

        //update linear vel
        Vec2f linear_vel;
        rbd->getLinearVelocity(linear_vel);
        linear_vel = linear_vel + j / rbd->getMass() * normal;
        rbd->setLinearVelocity(linear_vel);
        //update angular vel
        Vec3f momentum;
        rbd->getAngularMomentum(momentum[2]);
        momentum = momentum + j*cross(rad3d, normal3d);
        rbd->setAngularMomentum(momentum[2]);

        //recompute list of separating points, remove any separating points
        std::vector<int> new_nonseparating;
        for (unsigned int i = 0; i < penetrating_points.size(); ++i) {
            float phi;
            Vec2f normal;
            int index = penetrating_points[i];
            Vec2f vertex(vertices[index][0], vertices[index][1]);
            phi = interpolate_normal(normal, vertex, nodal_solid_phi, 0.5f*dx*Vec2f(1, 1), dx);

            Vec2f pt_vel = rbd->getPointVelocity(vertex);
            if (dot(pt_vel, normal) < 1e-2) {
                new_nonseparating.push_back(index);
            }
        }
        penetrating_points.clear();
        for (unsigned int i = 0; i < new_nonseparating.size(); ++i) {
            penetrating_points.push_back(new_nonseparating[i]);
        }


    }


    //project the whole object to be outside the wall.
    //again, really not the ideal way to deal with this...
    bool found_penetrating = false;
    do {
        found_penetrating = false;
        vertices.clear();
        rbd->get2DVertices(vertices);
        min_phi = 10 * dx;
        min_ind = -1;
        min_normal = Vec2f(0, 0);
        for (unsigned int i = 0; i < vertices.size(); ++i) {

            float phi;
            Vec2f normal;
            Vec2f vertex(vertices[i][0], vertices[i][1]);
            phi = interpolate_normal(normal, vertex, nodal_solid_phi, Vec2f(0, 0), dx);

            if (phi < min_phi) {
                min_phi = phi;
                min_normal = normal;
            }
        }

        if (min_phi < 1e-3) {
            found_penetrating = true;
            Vec2f com0;
            rbd->getCOM(com0);
            rbd->setCOM(com0 - (min_phi - 1e-2f)*min_normal);
        }
    } while (found_penetrating);

}


void FluidSim::update_rigid_body_grids()
{
    if (!rbd) return;

	//update level set from current position
	for (int j = 0; j < nj + 1; ++j) {
		for (int i = 0; i < ni + 1; ++i) {
            Vec2f location((i+0.5f)*dx, (j+0.5f)*dx);
			nodal_rigid_phi(i, j) = rbd->getSignedDist(location);
		}
	}

    //compute face area fractions from distance field
    rigid_u_weights.set_zero();
    rigid_v_weights.set_zero();
    for (int i = 0; i < rigid_u_weights.ni; ++i) {
        for (int j = 0; j < rigid_u_weights.nj; ++j) {
            float rigid_bottom = nodal_rigid_phi(i, j);
            float rigid_top = nodal_rigid_phi(i, j + 1);
            rigid_u_weights(i, j) = fraction_inside(rigid_bottom, rigid_top);
        }
    }

    for (int i = 0; i < v_vol.ni; ++i) {
        for (int j = 1; j < v_vol.nj - 1; ++j) {
            float rigid_left = nodal_rigid_phi(i, j);
            float rigid_right = nodal_rigid_phi(i + 1, j);
            rigid_v_weights(i, j) = fraction_inside(rigid_left, rigid_right);
        }
    }

	//recompute the grid-based "effective" masses per axis, so that we can *exactly*
	//balance in hydrostatic scenarios.
	double u_sum = std::accumulate(rigid_u_weights.a.begin(), rigid_u_weights.a.end(), 0.0);
	double v_sum = std::accumulate(rigid_v_weights.a.begin(), rigid_v_weights.a.end(), 0.0);
	rigid_u_mass = rbd->getDensity() * (float)u_sum;
	rigid_v_mass = rbd->getDensity() * (float)v_sum;
}

void FluidSim::recompute_solid_velocity()
{
    if (!rbd) return;

    for (int i = 0; i < ni + 1; ++i) {
        for (int j = 0; j < nj - 1; ++j) {
            Vec2f pos(i*dx, (j + 0.5f)*dx);
            if (0.5f*(nodal_solid_phi(i, j) + nodal_solid_phi(i, j + 1)) < 0.5f*(nodal_rigid_phi(i, j) + nodal_rigid_phi(i, j + 1)))
                solid_u(i, j) = 0;
            else {
                solid_u(i, j) = rbd->getPointVelocity(pos)[0];
            }
        }
    }
    for (int i = 0; i < ni - 1; ++i) {
        for (int j = 0; j < nj + 1; ++j) {
            Vec2f pos((i + 0.5f)*dx, j*dx);
            if (0.5f*(nodal_solid_phi(i, j) + nodal_solid_phi(i + 1, j)) < 0.5f*(nodal_rigid_phi(i, j) + nodal_rigid_phi(i + 1, j)))
                solid_v(i, j) = 0;
            else {
                solid_v(i, j) = rbd->getPointVelocity(pos)[1];
            }
        }
    }

}

//Add a tracer particle for visualization
void FluidSim::add_particle(const Vec2f& position) {
    particles.push_back(position);
    particles_velocity.push_back(Vec2f(0.0f, 0.0f));
}

//Basic first order semi-Lagrangian advection of velocities
void FluidSim::advect(float dt) {

    //semi-Lagrangian advection on u-component of velocity
    for (int j = 0; j < nj; ++j) for (int i = 0; i < ni + 1; ++i) {
        Vec2f pos(i*dx, (j + 0.5f)*dx);
        pos = trace_rk2(pos, -dt);
        temp_u(i, j) = get_velocity(pos)[0];
    }

    //semi-Lagrangian advection on v-component of velocity
    for (int j = 0; j < nj + 1; ++j) for (int i = 0; i < ni; ++i) {
        Vec2f pos((i + 0.5f)*dx, j*dx);
        pos = trace_rk2(pos, -dt);
        temp_v(i, j) = get_velocity(pos)[1];
    }

    //move update velocities into u/v vectors
    u = temp_u;
    v = temp_v;
}

//Perform 2nd order Runge Kutta to move the particles in the fluid
void FluidSim::advect_particles(float dt) {

    for (unsigned int p = 0; p < particles.size(); ++p) {
        Vec2f before = particles[p];
        Vec2f start_velocity = get_velocity(before);
        Vec2f midpoint = before + 0.5f*dt*start_velocity;
        Vec2f mid_velocity = get_velocity(midpoint);
        particles[p] += dt*mid_velocity;
        Vec2f after = particles[p];

        //Particles can still occasionally leave the domain due to truncation errors,
        //interpolation error, or large timesteps, so we project them back in for good measure.

        //Try commenting this section out to see the degree of accumulated error.
        float phi_value = interpolate_value(particles[p] / dx, nodal_solid_phi);
        if (phi_value < 0) {
            Vec2f normal;
            interpolate_gradient(normal, particles[p] / dx, nodal_solid_phi);
            normalize(normal);
            particles[p] -= phi_value*normal;
        }

        if (rbd)
        {
            rbd->testCollisionAndProject(particles[p], particles[p]);
        }
    }
}

void FluidSim::accumulate(Array2f& accum, float q, int i, int j, float fx, float fy)
{
    float weight;

    weight=(1-fx)*(1-fy);
    accum(i,j)+=weight*q;
    sum(i,j)+=weight;

    weight=fx*(1-fy);
    accum(i+1,j)+=weight*q;
    sum(i+1,j)+=weight;

    weight=(1-fx)*fy;
    accum(i,j+1)+=weight*q;
    sum(i,j+1)+=weight;

    weight=fx*fy;
    accum(i+1,j+1)+=weight*q;
    sum(i+1,j+1)+=weight;
}

void FluidSim::transfer_to_grid()
{
    u.set_zero();
    sum.set_zero();
    for(int p = 0; p < particles.size(); ++p)
    {
        Vec2f point = particles[p];
        int i,j;
        float fx,fy;

        //determine containing cell;
        get_barycentric((point[0])/dx, i, fx, 0, ni);
        get_barycentric((point[1])/dx-0.5f, j, fy, 0, nj);

        accumulate(u, particles_velocity[p][0], i, j, fx, fy);
    }

    for(int j = 0; j < u.nj; ++j) for(int i = 0; i < u.ni; ++i)
    {
        if(sum(i,j) != 0) u(i,j) /= sum(i,j);
    }

    v.set_zero();
    sum.set_zero();
    for(int p = 0; p < particles.size(); ++p)
    {
        Vec2f point = particles[p];
        int i,j;
        float fx,fy;

        //determine containing cell;
        get_barycentric((point[0])/dx-0.5f, i, fx, 0, ni);
        get_barycentric((point[1])/dx, j, fy, 0, nj);

        accumulate(v, particles_velocity[p][1], i, j, fx, fy);
    }

    for(int j = 0; j < v.nj; ++j) for(int i = 0; i < v.ni; ++i)
    {
        if(sum(i,j) != 0) v(i,j) /= sum(i,j);
    }
}

void FluidSim::update_from_grid()
{
    for(int p = 0; p < particles.size(); ++p)
    {
        Vec2f point = particles[p];

        // FLIP
        //u[p]+=Vec2f(grid.du.bilerp(ui, j, ufx, fy), grid.dv.bilerp(i, vj, fx, vfy));

        // PIC
        particles_velocity[p] = get_velocity(point);
        //u[p]=Vec2f(grid.u.bilerp(ui, j, ufx, fy), grid.v.bilerp(i, vj, fx, vfy));
    }
}

void FluidSim::compute_phi() {

    //Estimate from particles
    liquid_phi.assign(3 * dx);
    for (unsigned int p = 0; p < particles.size(); ++p) {
        Vec2f point = particles[p];
        int i, j;
        float fx, fy;
        //determine containing cell;
        get_barycentric((point[0]) / dx - 0.5f, i, fx, 0, ni);
        get_barycentric((point[1]) / dx - 0.5f, j, fy, 0, nj);

        //compute distance to surrounding few points, keep if it's the minimum
        for (int j_off = j - 2; j_off <= j + 2; ++j_off) for (int i_off = i - 2; i_off <= i + 2; ++i_off) {
            if (i_off < 0 || i_off >= ni || j_off < 0 || j_off >= nj)
                continue;

            Vec2f pos((i_off + 0.5f)*dx, (j_off + 0.5f)*dx);
            float phi_temp = dist(pos, point) - 1.02f*particle_radius;
            liquid_phi(i_off, j_off) = min(liquid_phi(i_off, j_off), phi_temp);
        }
    }
}

void FluidSim::extrapolate_phi()
{
    //"extrapolate" phi into solids if nearby
    for (int j = 0; j < nj; ++j) {
        for (int i = 0; i < ni; ++i) {
            if (liquid_phi(i, j) < 0.5*dx) {
                float solid_phi_val = 0.25f*(nodal_solid_phi(i, j) + nodal_solid_phi(i + 1, j) + nodal_solid_phi(i, j + 1) + nodal_solid_phi(i + 1, j + 1));
                if (solid_phi_val < 0)
                    liquid_phi(i, j) = -0.5f*dx;
            }
        }
    }

}



void FluidSim::apply_projection(float dt) {
    //Compute finite-volume type face area weight for each velocity sample.
    compute_pressure_weights();
    solve_pressure(dt); //original Batty '07 style: SPD, but with dense blocks due to J'MJ terms
}

//Apply RK2 to advect a point in the domain.
Vec2f FluidSim::trace_rk2(const Vec2f& position, float dt) {
    Vec2f input = position;
    Vec2f velocity = get_velocity(input);
    velocity = get_velocity(input + 0.5f*dt*velocity);
    input += dt*velocity;
    return input;
}

//Interpolate velocity from the MAC grid.
Vec2f FluidSim::get_velocity(const Vec2f& position) {

    //Interpolate the velocity from the u and v grids
    float u_value = interpolate_value(position / dx - Vec2f(0, 0.5f), u);
    float v_value = interpolate_value(position / dx - Vec2f(0.5f, 0), v);

    return Vec2f(u_value, v_value);
}

Vec2f FluidSim::get_solid_velocity(const Vec2f& position) {

    //Interpolate the velocity from the u and v grids
    float u_value = interpolate_value(position / dx - Vec2f(0, 0.5f), solid_u);
    float v_value = interpolate_value(position / dx - Vec2f(0.5f, 0), solid_v);

    return Vec2f(u_value, v_value);
}


//Compute finite-volume style face-weights for fluid from nodal signed distances
void FluidSim::compute_pressure_weights() {

    //recompute the face area weights, consider static and dynamic solid
    for (int j = 0; j < u_weights.nj; ++j) for (int i = 0; i < u_weights.ni; ++i) {
        u_weights(i, j) = 1 - fraction_inside(nodal_solid_phi(i, j + 1), nodal_solid_phi(i, j)) - rigid_u_weights(i, j);
        u_weights(i, j) = clamp(u_weights(i, j), 0.0f, 1.0f);
    }
    for (int j = 0; j < v_weights.nj; ++j) for (int i = 0; i < v_weights.ni; ++i) {
        v_weights(i, j) = 1 - fraction_inside(nodal_solid_phi(i + 1, j), nodal_solid_phi(i, j)) - rigid_v_weights(i, j);
        v_weights(i, j) = clamp(v_weights(i, j), 0.0f, 1.0f);
    }

}

//An implementation of the variational pressure projection solve for static geometry
void FluidSim::solve_pressure(float dt) {

    // Assemble the data for the J vectors
    Array2d base_trans_x(ni, nj), base_trans_y(ni, nj), base_rot_z(ni, nj);
    Vec2f centre_of_mass;
    if (rbd)
    {
        rbd->getCOM(centre_of_mass);

        for (int j = 0; j != nj; ++j) {
            for (int i = 0; i != ni; ++i) {
                double u_term = (rigid_u_weights(i + 1, j) - rigid_u_weights(i, j)) / dx;
                double v_term = (rigid_v_weights(i, j + 1) - rigid_v_weights(i, j)) / dx;

                // Translation coupling
                base_trans_x(i, j) = u_term;
                base_trans_y(i, j) = v_term;

                // Rotation coupling
                Vec2f position((i + 0.5f) * dx, (j + 0.5f) * dx);
                Vec2f rad = position - centre_of_mass;
                base_rot_z(i, j) = rad[0] * v_term - rad[1] * u_term;
            }
        }
    }

    int ni = v.ni;
    int nj = u.nj;
    int system_size = ni*nj;
    if (rhs.size() != system_size) {
        rhs.resize(system_size);
        pressure.resize(system_size);
        matrix.resize(system_size);
    }
    matrix.zero();

    //Build the linear system for pressure
    for (int j = 1; j < nj - 1; ++j) {
        for (int i = 1; i < ni - 1; ++i) {
            int index = i + ni*j;
            rhs[index] = 0;
            pressure[index] = 0;
            float centre_phi = liquid_phi(i, j);
            if (centre_phi < 0) {

                //right neighbour
                float term = u_weights(i + 1, j) * dt / sqr(dx);
                if (term > 0) {
                    float right_phi = liquid_phi(i + 1, j);
                    if (right_phi < 0) {
                        matrix.add_to_element(index, index, term);
                        matrix.add_to_element(index, index + 1, -term);
                    }
                    else {
                        float theta = fraction_inside(centre_phi, right_phi);
                        if (theta < 0.01f) theta = 0.01f;
                        matrix.add_to_element(index, index, term / theta);
                    }
                    rhs[index] -= u_weights(i + 1, j)*u(i + 1, j) / dx;
                }

                //left neighbour
                term = u_weights(i, j) * dt / sqr(dx);
                if (term > 0) {
                    float left_phi = liquid_phi(i - 1, j);
                    if (left_phi < 0) {
                        matrix.add_to_element(index, index, term);
                        matrix.add_to_element(index, index - 1, -term);
                    }
                    else {
                        float theta = fraction_inside(centre_phi, left_phi);
                        if (theta < 0.01f) theta = 0.01f;
                        matrix.add_to_element(index, index, term / theta);
                    }
                    rhs[index] += u_weights(i, j)*u(i, j) / dx;
                }

                //top neighbour
                term = v_weights(i, j + 1) * dt / sqr(dx);
                if (term > 0) {
                    float top_phi = liquid_phi(i, j + 1);
                    if (top_phi < 0) {
                        matrix.add_to_element(index, index, term);
                        matrix.add_to_element(index, index + ni, -term);
                    }
                    else {
                        float theta = fraction_inside(centre_phi, top_phi);
                        if (theta < 0.01f) theta = 0.01f;
                        matrix.add_to_element(index, index, term / theta);
                    }
                    rhs[index] -= v_weights(i, j + 1)*v(i, j + 1) / dx;
                }

                //bottom neighbour
                term = v_weights(i, j) * dt / sqr(dx);
                if (term > 0) {
                    float bot_phi = liquid_phi(i, j - 1);
                    if (bot_phi < 0) {
                        matrix.add_to_element(index, index, term);
                        matrix.add_to_element(index, index - ni, -term);
                    }
                    else {
                        float theta = fraction_inside(centre_phi, bot_phi);
                        if (theta < 0.01f) theta = 0.01f;
                        matrix.add_to_element(index, index, term / theta);
                    }
                    rhs[index] += v_weights(i, j)*v(i, j) / dx;
                }
            }
        }
    }

    if (rbd)
    {
        Vec2f solidLinearVelocity;
        float angular_velocity;
        rbd->getLinearVelocity(solidLinearVelocity);
        rbd->getAngularVelocity(angular_velocity);

        const float Jinv = rbd->getInvInertiaModulus();

        for (int j = 0; j < nj; ++j) {
            for (int i = 0; i < ni; ++i) {
                int index = i + ni*j;
                float centre_phi = liquid_phi(i, j);
                if (centre_phi < 0) {

                    //RHS contributions...
                    // Translation
                    rhs[index] -= solidLinearVelocity[0] * base_trans_x(i, j);
                    rhs[index] -= solidLinearVelocity[1] * base_trans_y(i, j);

                    //Rotation
                    rhs[index] -= angular_velocity * base_rot_z(i, j);

                    //LHS matrix contributions
                    for (int k = 0; k < ni; ++k) {
                        for (int m = 0; m < nj; ++m) {
                            double val = 0;
                            float other_phi = liquid_phi(k, m);
                            if (other_phi < 0) {
                                //Translation
                                val += dt * base_trans_x(i, j) * base_trans_x(k, m) / rigid_u_mass;
                                val += dt * base_trans_y(i, j) * base_trans_y(k, m) / rigid_v_mass;

                                //Rotation
                                val += dt * base_rot_z(i, j) * base_rot_z(k, m) * Jinv;

                                if (fabs(val) > 1e-10) {
                                    matrix.add_to_element(i + ni*j, k + ni*m, val);
                                }
                            }
                        }
                    }

                }
            }
        }
    }

    //replace empty rows/cols to make (at least) positive semi-definite
   for (unsigned row = 0; row < matrix.n; row++) {
     if (matrix.index[row].empty()) {
       rhs[row] = 0.0;
     }
   }

    //Solve the system using Robert Bridson's incomplete Cholesky PCG solver
    double tolerance;
    int iterations;
    bool success = solver.solve(matrix, rhs, pressure, tolerance, iterations);
    if (!success) {
        printf("WARNING: Pressure solve failed!************************************************\n");
    }

    //Apply the velocity update
    u_valid.assign(0);
    for (int j = 0; j < u.nj; ++j) for (int i = 1; i < u.ni - 1; ++i) {
        int index = i + j*ni;
        if (u_weights(i, j) > 0 && (liquid_phi(i, j) < 0 || liquid_phi(i - 1, j) < 0)) {
            float theta = 1;
            if (liquid_phi(i, j) >= 0 || liquid_phi(i - 1, j) >= 0)
                theta = fraction_inside(liquid_phi(i - 1, j), liquid_phi(i, j));
            if (theta < 0.01f) theta = 0.01f;
            u(i, j) -= dt  * (float)(pressure[index] - pressure[index - 1]) / dx / theta;
            u_valid(i, j) = 1;
        }
        else
            u(i, j) = 0;
    }
    v_valid.assign(0);
    for (int j = 1; j < v.nj - 1; ++j) for (int i = 0; i < v.ni; ++i) {
        int index = i + j*ni;
        if (v_weights(i, j) > 0 && (liquid_phi(i, j) < 0 || liquid_phi(i, j - 1) < 0)) {
            float theta = 1;
            if (liquid_phi(i, j) >= 0 || liquid_phi(i, j - 1) >= 0)
                theta = fraction_inside(liquid_phi(i, j - 1), liquid_phi(i, j));
            if (theta < 0.01f) theta = 0.01f;
            v(i, j) -= dt  * (float)(pressure[index] - pressure[index - ni]) / dx / theta;
            v_valid(i, j) = 1;
        }
        else
            v(i, j) = 0;
    }

    //Get pressure update to apply to solid
    if (rbd)
    {
        Vec2f updated_rigid_linear_velocity;
        float updated_rigid_angular_momentum;
        rbd->getLinearVelocity(updated_rigid_linear_velocity);
        rbd->getAngularMomentum(updated_rigid_angular_momentum);
        for (int j = 0; j < nj; ++j) {
            for (int i = 0; i < ni; ++i) {
                int index = i + ni*j;
                float centre_phi = liquid_phi(i, j);
                if (centre_phi < 0) {
                    updated_rigid_linear_velocity[0] += (float)(dt*base_trans_x(i, j)*pressure[index] / rigid_u_mass);
                    updated_rigid_linear_velocity[1] += (float)(dt*base_trans_y(i, j)*pressure[index] / rigid_v_mass);

                    updated_rigid_angular_momentum += (float)(dt*base_rot_z(i, j) * pressure[index]);
                }
            }
        }

        rbd->setLinearVelocity(updated_rigid_linear_velocity);
        rbd->setAngularMomentum(updated_rigid_angular_momentum);
    }
}

//Apply several iterations of a very simple "Jacobi"-style propagation of valid velocity data in all directions
void extrapolate(Array2f& grid, Array2c& valid) {

    Array2c old_valid(valid.ni, valid.nj);
    for (int layers = 0; layers < 10; ++layers) {
        old_valid = valid;
        Array2f temp_grid = grid;
        for (int j = 1; j < grid.nj - 1; ++j) for (int i = 1; i < grid.ni - 1; ++i) {
            float sum = 0;
            int count = 0;

            if (!old_valid(i, j)) {

                if (old_valid(i + 1, j)) {
                    sum += grid(i + 1, j); \
                    ++count;
                }
                if (old_valid(i - 1, j)) {
                    sum += grid(i - 1, j); \
                    ++count;
                }
                if (old_valid(i, j + 1)) {
                    sum += grid(i, j + 1); \
                    ++count;
                }
                if (old_valid(i, j - 1)) {
                    sum += grid(i, j - 1); \
                    ++count;
                }

                //If any of neighbour cells were valid,
                //assign the cell their average value and tag it as valid
                if (count > 0) {
                    temp_grid(i, j) = sum / (float)count;
                    valid(i, j) = 1;
                }

            }
        }
        grid = temp_grid;

    }

}
