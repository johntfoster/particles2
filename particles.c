#include <math.h>
#include <stdio.h>

void transfer_momentum(double *x, 
                       double *y,
                       double *vx, 
                       double *vy,
                       int number_of_particles,
                       int *neighbor_length_list,
                       int *neighbor_list,
                       double particle_diameter)
{
    int i,j;
    int neighbor_id;
    int number_of_neighbors;
    int first_contact;
    double normal_x, normal_y, tangent_normal_x, tangent_normal_y;
    double vx_normal, vy_normal, vx_tangent, vy_tangent;
    double distance_x, distance_y, distance;

    for(i=0; i<number_of_particles; i++, neighbor_length_list++){

        //Dereference the neighbors pointer
        number_of_neighbors = *neighbor_length_list;


        // Set the 'first contact' bool to false (0)
        first_contact = 0;

        //Loop over neighbors
        for(j=0; j<number_of_neighbors; j++,neighbor_list++){

            //Only transfer the momentum of the first contact pair
            if(first_contact == 0){
            
                //Set the neighbor id to the dereferenced pointer value
                neighbor_id = *neighbor_list;

                
                //Comptute the x and y distance components
                distance_x = x[neighbor_id] - x[i];
                distance_y = y[neighbor_id] - y[i];

                //Compute the distance between particles
                distance = sqrt(distance_x * distance_x + 
                                distance_y * distance_y);

                //Check for contact
                if(distance < particle_diameter){

                    //Compute the normal vectors
                    normal_x = distance_x / distance;
                    normal_y = distance_y / distance;
                    tangent_normal_x = -normal_y;
                    tangent_normal_y =  normal_x;

                    //Transfer the momentum
                    vx_tangent = vx[i] * tangent_normal_x;
                    vy_tangent = vy[i] * tangent_normal_y;
                    vx_normal = vx[neighbor_id] * normal_x;
                    vy_normal = vy[neighbor_id] * normal_y;
                    vx[i] = vx_normal + vx_tangent;
                    vy[i] = vy_normal + vy_tangent;

                    //Set the contact bool to true (1) so the momentum
                    //transfer is skipped for the rest of the neighbors
                    first_contact = 1;
                }
            }
        }
    }
    return;
}

double sine_solve_function(double x,
                           double a,
                           double b,
                           double c,
                           double d,
                           double e,
                           double f)
{
    return e - b * c * cos(c * x + d) * (a + b * sin(c * x + d)) - f - x;
}

double sine_solve_function_derivative(double x,
                                      double a,
                                      double b,
                                      double c,
                                      double d,
                                      double e,
                                      double f)
{
    return (-1.0 - b * b * c * c * cos(d + c * x) * cos(d + c * x) + b * c * c *  
            sin(d + c * x) * (a - f + b * sin(d + c * x)));
}

double newtons_method_for_sin_wave(double a,
                                   double b,
                                   double c,
                                   double d,
                                   double e,
                                   double f,
                                   double x_guess)
{
    double xnew = 0.0;
    double xold = x_guess;
    double fx, fpx;

    for(int i=0; i<1e10; i++)
    {

        fx = sine_solve_function(xold, a, b, c, d, e, f);
        fpx = sine_solve_function_derivative(xold, a, b, c, d, e, f);

        xnew = xold - fx / fpx;

        if(fabs(xnew - xold)/xnew < 1.0e-6){
            return xnew;
        }
    }
    printf("The Newton iteration failed to converge in 1e5 iterations.\n");
    return 0.0;
}

void wall_contact(double *x, 
                  double *y,
                  double *vx, 
                  double *vy,
                  int number_of_particles,
                  double width,
                  double height,
                  int sine_wave_bool,
                  double particle_radius)
{

    int i;
    double x_on_sine_wave;
    double y_on_sine_wave;
    double normal_x, normal_y, tangent_normal_x, tangent_normal_y;
    double vx_normal, vy_normal, vx_tangent, vy_tangent;
    double distance_x, distance_y, distance;

    for(i=0; i<number_of_particles; i++){

        //Uniform reflection off of flat wall
        if(x[i] < particle_radius | 
           (width - x[i]) < particle_radius | 
           y[i] < particle_radius |
           (height - y[i] < particle_radius))
        {
            vx[i] = -vx[i];
            vy[i] = -vy[i];
        }

        //Do the following only if an actual sine wave wall exists.
        if( sine_wave_bool == 0){
            //The reflection off the sinusoid might be a little expensive, so let's
            //not even try unless there is even a remote chance a particle could be
            //in contact.
            if (y[i] < 2.0+particle_radius){

                //Now following:
                //http://math.stackexchange.com/questions/514820/distance-between-point-and-sine-wave
                //we have the point, x, on the sine wave y = a + b * sin(c*x + d), which is closest to a point (e,f) in space, can be found by solving:
                // e - b*c*cos(c*x + d) * (a + b * sin(c * x + d) - f - x == 0
                // we will do this with Newton's method
                
                // These could become nonconstant
                x_on_sine_wave = newtons_method_for_sin_wave(1.0, 1.0, 1.0, 
                        asin(1.0), x[i], y[i], x[i]);
                y_on_sine_wave = sin(x_on_sine_wave + asin(1.0)) + 1.0;

                distance_x = x_on_sine_wave - x[i];
                distance_y = y_on_sine_wave - y[i];
                distance = sqrt(distance_x * distance_x + distance_y * distance_y);

                if(distance < particle_radius){
                    
                    normal_x = distance_x / distance;
                    normal_y = distance_y / distance;

                    tangent_normal_x = -normal_y;
                    tangent_normal_y =  normal_x;

                    //Transfer the momentum
                    vx_tangent = vx[i] * tangent_normal_x;
                    vy_tangent = vy[i] * tangent_normal_y;
                    vx_normal = -vx[i] * normal_x;
                    vy_normal = -vy[i] * normal_y;
                    vx[i] = vx_normal + vx_tangent;
                    vy[i] = vy_normal + vy_tangent;
                }
            }
        }
    }
}
