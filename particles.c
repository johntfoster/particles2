#include <math.h>

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
    double normal_x, normal_y;
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

                    //Transfer the momentum
                    vx[i] = vx[neighbor_id] * normal_x - vx[i] * normal_y;
                    vy[i] = vy[neighbor_id] * normal_x - vy[i] * normal_y;

                    //Set the contact bool to true (1) so the momentum
                    //transfer is skipped for the rest of the neighbors
                    first_contact = 1;
                }
            }
        }
    }
}
