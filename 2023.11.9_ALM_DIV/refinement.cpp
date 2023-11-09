#include "refinement.hpp"


double **ruppert_refinement(int num_segments, int num_seeds, double range_of_seeds, double rho, double error) {

    printf(" refining meshes...\n");

    srand(time(nullptr));
    double **segments = new double* [num_segments];
    segments[0] = new double [5];
    segments[0][4] = num_segments;
    for (int i = 1; i < num_segments; i++) {
        segments[i] = new double [4];
    }
    for (int i = 0; i < num_segments; i++) {
        segments[i][0] = cos(2*pi/num_segments*i);
        segments[i][1] = sin(2*pi/num_segments*i);
    }
    for (int i = 0; i < num_segments - 1; i++) {
        segments[i][2] = segments[i+1][0];
        segments[i][3] = segments[i+1][1];
    }
    segments[num_segments-1][2] = segments[0][0];
    segments[num_segments-1][3] = segments[0][1];

    int num_vertices = num_segments + num_seeds;
    double **vertices = new double* [num_vertices];
    vertices[0] = new double [3];
    vertices[0][2] = num_vertices;
    for (int i = 1; i < num_vertices; i++) {
        vertices[i] = new double [2];
    }
    for (int i = 0; i < num_segments; i++) {
        vertices[i][0] = cos(2*pi/num_segments*i);
        vertices[i][1] = sin(2*pi/num_segments*i);
    }
    double num1;
    if (num_seeds !=0) {
        for (int i = num_segments; i < num_vertices; i++){
            num1 = rand()/double(RAND_MAX)*2*pi;
            vertices[i][0] = range_of_seeds*rand()/double(RAND_MAX)*cos(num1);
            vertices[i][1] = range_of_seeds*rand()/double(RAND_MAX)*sin(num1);
        }
    }

    int num_triangles; 
    double **triangles = new double*;
    printf(" 2D constrained delaunay triangulating...\n");
    triangles = gift_wrapping_2d(segments, vertices, error);

    int termination;

    for (int iter = 0; ; iter++) {

        if ( iter % 5 == 0) {
            printf("iteration = %d\n",iter + 1);
        }

        termination = 0;
// printf("1...\n");
        for (int i = 0; i < (int) vertices[0][2]; i++) {
            for (int j = 0; j < (int) segments[0][4]; j++) {
                
                if (in_min_containment_circle(segments[j][0], segments[j][1], segments[j][2], segments[j][3], vertices[i][0], vertices[i][1], error) > 0) {
// printf("no...\n");   
                    termination++;
// printf("j = %d\n", j);                    
// printf("num_segments = %d\n", (int) segments[0][4]);
// for (int k = 0; k < (int) segments[0][4]; k++) {
//     printf("%f %f; %f %f\n",segments[k][0],segments[k][1],segments[k][2],segments[k][3]);
// }                    
// printf("1.1...\n");                     
                    segments = split_segment(segments[j], segments);

// printf("num_segments = %d\n", (int) segments[0][4]);
// for (int k = 0; k < (int) segments[0][4]; k++) {
//     printf("%f %f; %f %f\n",segments[k][0],segments[k][1],segments[k][2],segments[k][3]);
// }
// abort();
                    for (int k = 0; k < (int) triangles[0][6]; k++) {
                        delete [] triangles[k];
                    }          
                    triangles = gift_wrapping_2d(segments, vertices, error);                    
                    break;
                }
            }
        }  
// printf("2...\n");
        for (int i = 0; i < (int) triangles[0][6]; i++) {

            if (radius_edge_ratio(triangles[i][0], triangles[i][1], triangles[i][2], triangles[i][3], triangles[i][4], triangles[i][5]) > rho) {
// printf("yes...\n");   
// printf("i = %d\n",i);  
// printf("triangles: %0.9f %0.9f ; %0.9f %0.9f ; %0.9f %0.9f;\n",triangles[i][0],triangles[i][1],triangles[i][2],triangles[i][3],triangles[i][4],triangles[i][5]);            
                termination++;
                vertices = enlarge_2d(vertices, num_vertices, num_vertices + 1, 2);
                vertices[0] = enlarge_1d(vertices[0], 2, 3);
                vertices[0][2] = num_vertices + 1;
                vertices[num_vertices] = find_circumcenter(triangles[i][0], triangles[i][1], triangles[i][2], triangles[i][3], triangles[i][4], triangles[i][5]);
                num_vertices++;
                break;
            }
 
        }
// printf("3...\n");
// printf("num_vertices = %d\n",num_vertices);
// printf("termination = %d\n",termination);
        if (termination == 0) {
            printf("iteration = %d\n",iter + 1);
            break;
        }
        // if (temp[1] == temp[0]) {
        //     printf("poor convergence in meshing");
        //     exit(0);  
        // }
// printf("4...\n");
// printf("(int) triangles[0][6] = %d\n",(int) triangles[0][6]);
        for (int k = 0; k < (int) triangles[0][6]; k++) {
            delete [] triangles[k];
        }
        triangles = gift_wrapping_2d(segments, vertices, error);
// printf("5...\n");     
// printf("(int) triangles[0][6] = %d\n",(int) triangles[0][6]);
    }      
    
    for (int i = 0; i < (int) triangles[0][6]; i++) {
        printf("%0.9f %0.9f ; %0.9f %0.9f ; %0.9f %0.9f;\n",triangles[i][0],triangles[i][1],triangles[i][2],triangles[i][3],triangles[i][4],triangles[i][5]);
    }
    
    // printf("(int) triangles[0][6] = %d\n",(int) triangles[0][6]);
    // for (int k = 0; k < (int) triangles[0][6]; k++) {
    //     printf("k = %d\n",k);
    //     delete [] triangles[k];
    //     printf("kk = %d\n",k);
    // }
    // printf("ok...\n");

    return triangles;
}


double **chew_first_refinement(int num_segments, double rho, double scale, double error) {

    printf(" refining meshes...\n");

    srand(time(nullptr));
    double **segments = new double* [num_segments];
    segments[0] = new double [5];
    segments[0][4] = num_segments;
    for (int i = 1; i < num_segments; i++) {
        segments[i] = new double [4];
    }
    for (int i = 0; i < num_segments; i++) {
        segments[i][0] = cos(2*pi/num_segments*i);
        segments[i][1] = sin(2*pi/num_segments*i);
    }
    for (int i = 0; i < num_segments - 1; i++) {
        segments[i][2] = segments[i+1][0];
        segments[i][3] = segments[i+1][1];
    }
    segments[num_segments-1][2] = segments[0][0];
    segments[num_segments-1][3] = segments[0][1];

    int num_vertices = num_segments;
    double **vertices = new double* [num_vertices];
    vertices[0] = new double [3];
    vertices[0][2] = num_vertices;
    for (int i = 1; i < num_vertices; i++) {
        vertices[i] = new double [2];
    }
    for (int i = 0; i < num_vertices; i++) {
        vertices[i][0] = cos(2*pi/num_segments*i);
        vertices[i][1] = sin(2*pi/num_segments*i);
    }

    double num1;
    // double num3 = 0;
    
    // double triangle[6] = {0};
    double *triangle = new double [6];

    int num_triangles; 
    double **triangles;
    printf(" 2D constrained delaunay triangulating...\n");
    triangles = gift_wrapping_2d(segments, vertices, error);

    int termination;

    for (int iter = 0; ; iter++) {

        if (iter % 50 == 0) {
            printf("iteration = %d\n",iter + 1);
        }

        termination = 0;

        for (int i = 0; i < (int) segments[0][4]; i++) {
            if (sqrt((segments[i][0]-segments[i][2])*(segments[i][0]-segments[i][2]) + (segments[i][1]-segments[i][3])*(segments[i][1]-segments[i][3])) > rho + error) {
                termination++;
                segments = split_segment(segments[i], segments); 
            }
        }

        if (termination != 0) {
            for (int i = 0; i < (int) triangles[0][6]; i++) {
                delete [] triangles[i];
            }      
            triangles = gift_wrapping_2d(segments, vertices, error);  
        }
        for (int i = 0; i < (int) triangles[0][6]; i++) {
            if (length_of_radius(triangles[i][0], triangles[i][1], triangles[i][2], triangles[i][3], triangles[i][4], triangles[i][5]) > (rho + error) / scale) {        
                termination++;
                num1 = 1;            
                vertices = enlarge_2d(vertices, num_vertices, num_vertices + 1, 2);                    
                vertices[0] = enlarge_1d(vertices[0], 2, 3);
               
                vertices[0][2] = num_vertices + 1;
               
                vertices[num_vertices] = find_circumcenter(triangles[i][0], triangles[i][1], triangles[i][2], triangles[i][3], triangles[i][4], triangles[i][5]);
                num_vertices++;

                triangle[0] = triangles[i][0];
                triangle[1] = triangles[i][1];
                triangle[2] = triangles[i][2];
                triangle[3] = triangles[i][3];
                triangle[4] = triangles[i][4];
                triangle[5] = triangles[i][5];

                break;

            }
 
        }

        if (termination == 0) {
            if (iter % 50 != 0) {
                printf("iteration = %d\n",iter + 1);
                break;
            }
        }
      

        if (num1 == 1){ 
            // for (int k = 0; k < (int) triangles[0][6]; k++) {
            //     delete [] triangles[k];
            // }
            // triangles = gift_wrapping_2d(segments, vertices, error);

            triangles = bowyer_watson(vertices[num_vertices-1][0], vertices[num_vertices-1][1], \
                        triangle[0], triangle[1], triangle[2], triangle[3], triangle[4], triangle[5], triangles, vertices, segments, error);
            num1 = 0;
        }
              
    }      

    triangles[0][7] = num_vertices;
    
    // for (int i = 0; i < (int) triangles[0][6]; i++) {
    //     printf("%0.9f %0.9f ; %0.9f %0.9f ; %0.9f %0.9f;\n",triangles[i][0],triangles[i][1],triangles[i][2],triangles[i][3],triangles[i][4],triangles[i][5]);
    // }
    // exit(0);

    for (int i = 1; i < num_vertices; i++) {
        delete[] vertices[i];
    }
    delete[] vertices;

    for (int i = 1; i < num_segments; i++) {
        delete[] segments[i];
    }
    delete[] segments;
    delete[] triangle;


    return triangles;
}


double **split_segment(double *segment, double **segments) {

    int num_segments = (int) segments[0][4];
    double x1, y1, x2, y2, angle1, angle2, angle;

    for (int i = 0; i < num_segments; i++) {
        if (segment[0] == segments[i][0] && segment[1] == segments[i][1] && segment[2] == segments[i][2] && segment[3] == segments[i][3]) {

            x1 = segment[0];
            y1 = segment[1];
            x2 = segment[2];
            y2 = segment[3];
            angle1 = atan2(y1,x1);
            angle2 = atan2(y2,x2);

            if (angle2 - angle1 < -pi) {
                angle2 += 2*pi;
            }
            else if (angle2 - angle1 > pi) {
                angle2 -= 2*pi;
            }

            angle = (angle1 + angle2) / 2;

            segments = insert_2d(segments, num_segments, i, 1, 4);

            segments[i][0] = x1;
            segments[i][1] = y1;
            // segments[i][2] = (x1+x2)/2;
            // segments[i][3] = (y1+y2)/2;
            segments[i][2] = cos(angle);
            segments[i][3] = sin(angle);

            // segments[i+1][0] = (x1+x2)/2;
            // segments[i+1][1] = (y1+y2)/2;
            segments[i+1][0] = cos(angle);
            segments[i+1][1] = sin(angle);
            segments[i+1][2] = x2;
            segments[i+1][3] = y2;

            segments[0] = enlarge_1d(segments[0], 4, 5);
            segments[0][4] = num_segments + 1;

            break;

        }
    }

    return segments;

}


double *find_circumcenter(double p1_x, double p1_y, double p2_x, double p2_y, double p3_x, double p3_y) {

    double dx2 = p2_x - p1_x;
    double dy2 = p2_y - p1_y;
    double dx3 = p3_x - p1_x;
    double dy3 = p3_y - p1_y;
    double *p4 = new double[2];

    p4[0] = p1_x + (dx2*dx2*dy3-dx3*dx3*dy2+dy2*dy2*dy3-dy3*dy3*dy2) / (dx2*dy3-dx3*dy2) / 2;
    p4[1] = p1_y + (dx2*dx2*dx3-dx3*dx3*dx2+dy2*dy2*dx3-dy3*dy3*dx2) / (dx3*dy2-dx2*dy3) / 2;

    return p4;

}


double in_min_containment_circle(double p1_x, double p1_y, double p2_x, double p2_y, double p3_x, double p3_y, double error) {

    double cx = (p1_x+p2_x)/2, cy = (p1_y+p2_y)/2;
    double num1 = 0;

    if ( sqrt((p1_x-p2_x)*(p1_x-p2_x)+(p1_y-p2_y)*(p1_y-p2_y))/2 - sqrt((p3_x-cx)*(p3_x-cx)+(p3_y-cy)*(p3_y-cy)) > error ) {
        num1 = 1;
    }

    return num1;

}


double length_of_radius(double p1_x, double p1_y, double p2_x, double p2_y, double p3_x, double p3_y) {

    double dx2 = p2_x - p1_x;
    double dy2 = p2_y - p1_y;
    double dx3 = p3_x - p1_x;
    double dy3 = p3_y - p1_y;
    double xc_x1 = (dx2*dx2*dy3-dx3*dx3*dy2+dy2*dy2*dy3-dy3*dy3*dy2) / (dx2*dy3-dx3*dy2) / 2;
    double yc_y1 = (dx2*dx2*dx3-dx3*dx3*dx2+dy2*dy2*dx3-dy3*dy3*dx2) / (dx3*dy2-dx2*dy3) / 2;
    double radius = sqrt(xc_x1*xc_x1 + yc_y1*yc_y1);

    return radius;
}


double length_of_min_edge(double p1_x, double p1_y, double p2_x, double p2_y, double p3_x, double p3_y) {

    double num1 = sqrt((p1_x-p2_x)*(p1_x-p2_x)+(p1_y-p2_y)*(p1_y-p2_y));
    double num2 = sqrt((p1_x-p3_x)*(p1_x-p3_x)+(p1_y-p3_y)*(p1_y-p3_y));
    double num3 = sqrt((p2_x-p3_x)*(p2_x-p3_x)+(p2_y-p3_y)*(p2_y-p3_y));

    if (num1 > num2) {
        num1 = num2;
    }

    if (num1 > num3) {
        num1 = num3;
    }

    return num1;
}


double radius_edge_ratio(double p1_x, double p1_y, double p2_x, double p2_y, double p3_x, double p3_y) {

    double num1 = sqrt((p1_x-p2_x)*(p1_x-p2_x)+(p1_y-p2_y)*(p1_y-p2_y));
    double num2 = sqrt((p1_x-p3_x)*(p1_x-p3_x)+(p1_y-p3_y)*(p1_y-p3_y));
    double num3 = sqrt((p2_x-p3_x)*(p2_x-p3_x)+(p2_y-p3_y)*(p2_y-p3_y));

    double dx2 = p2_x - p1_x;
    double dy2 = p2_y - p1_y;
    double dx3 = p3_x - p1_x;
    double dy3 = p3_y - p1_y;
    double xc_x1 = (dx2*dx2*dy3-dx3*dx3*dy2+dy2*dy2*dy3-dy3*dy3*dy2) / (dx2*dy3-dx3*dy2) / 2;
    double yc_y1 = (dx2*dx2*dx3-dx3*dx3*dx2+dy2*dy2*dx3-dy3*dy3*dx2) / (dx3*dy2-dx2*dy3) / 2;
    double radius = sqrt(xc_x1*xc_x1 + yc_y1*yc_y1);

    if (num1 > num2) {
        num1 = num2;
    }

    if (num1 > num3) {
        num1 = num3;
    }

    return radius/num1;
}


double **bowyer_watson(double x1, double y1, double x2, double y2, double x3, double y3, double x4, double y4, double** triangles, double **vertices, double **segments, double error) {

    int i = 0, num_triangles = (int) triangles[0][6];
    int num1 = 0;

    for (i = 0; i < num_triangles; i++) {
        if ((value_comparison(triangles[i][0], x2, error) && value_comparison(triangles[i][1], y2, error)) && \
            (value_comparison(triangles[i][2], x3, error) && value_comparison(triangles[i][3], y3, error)) && \
            (value_comparison(triangles[i][4], x4, error) && value_comparison(triangles[i][5], y4, error))) {
                num1 = 1;
                break;
        }
        if ((value_comparison(triangles[i][0], x3, error) && value_comparison(triangles[i][1], y3, error)) && \
            (value_comparison(triangles[i][2], x4, error) && value_comparison(triangles[i][3], y4, error)) && \
            (value_comparison(triangles[i][4], x2, error) && value_comparison(triangles[i][5], y2, error))) {
                num1 = 1;
                break;
        }
        if ((value_comparison(triangles[i][0], x4, error) && value_comparison(triangles[i][1], y4, error)) && \
            (value_comparison(triangles[i][2], x2, error) && value_comparison(triangles[i][3], y2, error)) && \
            (value_comparison(triangles[i][4], x3, error) && value_comparison(triangles[i][5], y3, error))) {
                num1 = 1;
                break;
        }
    }

    if (num1 == 0) {
        printf("bowyer_watson: error in searching for triangles which should be deleted\n");
        printf("%0.9f %0.9f\n", x2, y2);
        printf("%0.9f %0.9f\n", x3, y3);
        printf("%0.9f %0.9f\n", x4, y4);

        printf("%d\n",(int) triangles[0][6]);

        exit(0);
    }

    triangles = reduce_2d(triangles, num_triangles, i, 1, 6);  
    triangles[0] = enlarge_1d(triangles[0], 6, 8);
    triangles[0][6] = num_triangles - 1;
    num_triangles--;

    triangles = dig_cavity(x1, y1, x2, y2, x3, y3, triangles, vertices, segments, error); 
    triangles = dig_cavity(x1, y1, x3, y3, x4, y4, triangles, vertices, segments, error);
    triangles = dig_cavity(x1, y1, x4, y4, x2, y2, triangles, vertices, segments, error);

    return triangles;

}


double **dig_cavity(double x1, double y1, double x2, double y2, double x3, double y3, double** triangles, double **vertices, double **segments, double error) {

    double *p4 = adjacent(x3, y3, x2, y2, triangles, vertices, segments, error);

    int i = 0, j = 0, num_triangles = (int) triangles[0][6], num_vertices = (int) vertices[0][2];
    int num1 = 0;
    double array1[2];

    if (p4[0] != 9999.0) {
        if (in_circle(x1, y1, x2, y2, x3, y3, p4[0], p4[1]) > error) {  

            for (i = 0; i < num_triangles; i++) {
    
                if ((triangles[i][0] == x3 && triangles[i][1] == y3) && \
                    (triangles[i][2] == x2 && triangles[i][3] == y2) && \
                    (triangles[i][4] == p4[0] && triangles[i][5] == p4[1])) {
                        triangles = reduce_2d(triangles, num_triangles, i, 1, 6);
                        triangles[0] = enlarge_1d(triangles[0], 6, 8);
                        triangles[0][6] = num_triangles - 1;
                        num_triangles--;
                        break;
                }
                if ((triangles[i][0] == x2 && triangles[i][1] == y2) && \
                    (triangles[i][2] == p4[0] && triangles[i][3] == p4[1]) && \
                    (triangles[i][4] == x3 && triangles[i][5] == y3)) {
                        triangles = reduce_2d(triangles, num_triangles, i, 1, 6);
                        triangles[0] = enlarge_1d(triangles[0], 6, 8);
                        triangles[0][6] = num_triangles - 1;
                        num_triangles--;
                        break;
                }
                if ((triangles[i][0] == p4[0] && triangles[i][1] == p4[1]) && \
                    (triangles[i][2] == x3 && triangles[i][3] == y3) && \
                    (triangles[i][4] == x2 && triangles[i][5] == y2)) {
                        triangles = reduce_2d(triangles, num_triangles, i, 1, 6);
                        triangles[0] = enlarge_1d(triangles[0], 6, 8);
                        triangles[0][6] = num_triangles - 1;
                        num_triangles--;
                        break;
                } 
                // four co-circle points exist
                if (i == num_triangles - 1)  {
                    for (j = 0; j < num_vertices; j++) {
                        if (p4[0] == vertices[j][0] && p4[1] == vertices[j][1]) {

                            array1[0] = vertices[j][0];
                            array1[1] = vertices[j][1];

                            vertices = reduce_2d(vertices, num_vertices, j, 1, 2);
                            vertices[0] = enlarge_1d(vertices[0], 2, 3);
                            vertices[0][2] = num_vertices - 1;
                            num_vertices--;

                            p4 = adjacent(x3, y3, x2, y2, triangles, vertices, segments, error);

                            vertices = insert_2d(vertices, num_vertices, j-1, 1, 2);
                            vertices[j][0] = array1[0];
                            vertices[j][1] = array1[1];
                            vertices[0] = enlarge_1d(vertices[0], 2, 3);
                            vertices[0][2] = num_vertices + 1;
                            num_vertices++;

                            i = 0;
                            num1 = 1;
                        
                            break;
                        }
                    }
                    if (num1 == 0) {
                        printf("dig_cavity: error in searching for triangles which should be deleted\n");
                        printf("x1 = %f, y1 = %f, x2 = %f, y2 = %f, x3 = %f, y3 = %f, p4[0] = %f, p4[1] = %f\n",x1, y1, x2, y2, x3, y3, p4[0], p4[1]);
                        exit(0);
                    }
                    else {
                        continue;
                    }
                }      

            }

            triangles = dig_cavity(x1, y1, x2, y2, p4[0], p4[1], triangles, vertices, segments, error);
            triangles = dig_cavity(x1, y1, p4[0], p4[1], x3, y3, triangles, vertices, segments, error);

        }
        else {
            triangles = enlarge_2d(triangles, num_triangles, num_triangles + 1, 6);
            triangles[num_triangles][0] = x1;
            triangles[num_triangles][1] = y1;
            triangles[num_triangles][2] = x2;
            triangles[num_triangles][3] = y2;
            triangles[num_triangles][4] = x3;
            triangles[num_triangles][5] = y3;
            triangles[0] = enlarge_1d(triangles[0], 6, 8);
            triangles[0][6] = num_triangles + 1;
            num_triangles++;

        } 
    }
    
    if (p4[0] == 9999.0) {
        for (int i = 0; i < (int) segments[0][4]; i++) {
            if (x2 == segments[i][0] && y2 == segments[i][1] && x3 == segments[i][2] && y3 == segments[i][3]) {
                triangles = enlarge_2d(triangles, num_triangles, num_triangles + 1, 6);
                triangles[num_triangles][0] = x1;
                triangles[num_triangles][1] = y1;
                triangles[num_triangles][2] = x2;
                triangles[num_triangles][3] = y2;
                triangles[num_triangles][4] = x3;
                triangles[num_triangles][5] = y3;
                triangles[0] = enlarge_1d(triangles[0], 6, 8);
                triangles[0][6] = num_triangles + 1;
                num_triangles++;
                break;
            }
        }
    }

    return triangles;
}


double *adjacent(double x1, double y1, double x2, double y2, double** triangles, double **vertices, double **segments, double error) {

    double *p3 = new double[2];
    p3[0] = 9999.0;
    
    double relative_interior_p[2] = {(x1+x2)/2, (y1+y2)/2}; 
    double num1 = 0;
    for (int i = 0; i < (int) vertices[0][2] - 1; i++) {
        num1 = 0;
        if ((orientation_2d(x1, y1, x2, y2, vertices[i][0], vertices[i][1]) > error) && \
        ((p3[0] == 9999.0) || (in_circle(x1, y1, x2, y2, p3[0], p3[1], vertices[i][0], vertices[i][1]) > error))) {
            for (int j = 0; j < (int) segments[0][2]; j++) {
                if (if_itersection(relative_interior_p[0], relative_interior_p[1], vertices[i][0], vertices[i][1], \
                segments[j][0], segments[j][1], segments[j][2], segments[j][3], error) > 0 && ( vertices[i][0] != x1 && vertices[i][1] != y1)) {
                    num1++;
                    break;
                }
            }
            if (num1 == 0) {
                p3[0] = vertices[i][0];
                p3[1] = vertices[i][1];
            }
        }
    }

    return p3;
}