#include "triangulation.hpp"

double **gift_wrapping_2d(double **segments, double **vertices, double error) {

    int num_edges = (int) segments[0][4];    
    double **edges = new double* [4];
    for (int i = 0; i < 4; i++) {
        edges[i] = new double[num_edges];
  
    }
    for (int i = 0; i < num_edges; i++) {
        edges[0][i] = segments[i][0];
        edges[1][i] = segments[i][1];
        edges[2][i] = segments[i][2];
        edges[3][i] = segments[i][3];
    }

    int num_triangles = 0;
    double **triangles = new double*; 
    triangles[0] = new double [6];

    double edge[4];
    double *triangle = new double [6];

    int num1 = 0, num2 = 0, num3 = 0, num4 = 0;
    for (; num_edges != 0; ) {

        num4 = 0;
     
        for (int i = 0; i < 4; i++) {
            edge[i] = edges[i][0];
            edges[i] = reduce_1d(edges[i], num_edges, 0, 1); 
        }
        num_edges--;
       
        triangle = positively_oriented_triangle_2d(edge[0], edge[1], edge[2], edge[3], vertices, segments, error);

        if (triangle[4] == 9999) {
            continue;
        }

        num1 = 0;  
        num3 = 0;
        for (int i = 0; i < num_edges; i++) {
            if (  ((triangle[2] == edges[0][i]) && (triangle[3] == edges[1][i]) && (triangle[4] == edges[2][i]) && (triangle[5] == edges[3][i]))  ) {
                for (int j = 0; j < 4; j++) {
                    edges[j] = reduce_1d(edges[j], num_edges, i, 1);
                } 
                num_edges--;
                num3++;
                num4++;
                break;
            }
            else {
                num1++;    
            }                    
        }
        if (num1 == num_edges + num3) {
            for (int j = 0; j < 4; j++) {
                    edges[j] = enlarge_1d(edges[j], num_edges, num_edges+1);
            }
            edges[0][num_edges] = triangle[4];
            edges[1][num_edges] = triangle[5];
            edges[2][num_edges] = triangle[2];
            edges[3][num_edges] = triangle[3];
            num_edges++;
        }
        
        num2 = 0;  
        num3 = 0;
        for (int i = 0; i < num_edges; i++) {
            if (  ((triangle[4] == edges[0][i]) && (triangle[5] == edges[1][i]) && (triangle[0] == edges[2][i]) && (triangle[1] == edges[3][i]))  ) {
                for (int j = 0; j < 4; j++) {
                    edges[j] = reduce_1d(edges[j], num_edges, i, 1);
                }  
                num_edges--;
                num3++;
                num4++;
                break;
            }
            else {    
                num2++;
            }                    
        }
        if (num2 == num_edges + num3) {
            for (int j = 0; j < 4; j++) {
                    edges[j] = enlarge_1d(edges[j], num_edges, num_edges+1);
            }
            edges[0][num_edges] = triangle[0];
            edges[1][num_edges] = triangle[1];
            edges[2][num_edges] = triangle[4];
            edges[3][num_edges] = triangle[5];
            num_edges++;
        }

        if (num_triangles > (int) segments[0][4] && num4 == 2) {
            continue;
        }

        for (int i = 0; i < 6; i ++) {
            triangles[num_triangles][i] = triangle[i];
        }

        num_triangles++;
        triangles = enlarge_2d(triangles, num_triangles, num_triangles+1, 6);
        
    }

    // triangles = enlarge_2d(triangles, num_triangles, num_triangles, 6);     
    triangles[0] = enlarge_1d(triangles[0], 6, 8);   
    triangles[0][6] = num_triangles; 
    triangles[0][7] = num_edges; 

    delete [] edges[0];
    delete [] edges[1];
    delete [] edges[2];
    delete [] edges[3];
    delete [] edges;
    delete [] triangle;

    return triangles;

}


double *positively_oriented_triangle_2d(double p1_x, double p1_y, double p2_x, double p2_y, double **vertices, double **segments, double error) {

    double p3[2] = {9999};

    double relative_interior_p[2] = {(p1_x+p2_x)/2, (p1_y+p2_y)/2}; 
    
    double num1 = 0;
    for (int i = 0; i < (int) vertices[0][2]; i++) {
        num1 = 0;
        if ((orientation_2d(p1_x, p1_y, p2_x, p2_y, vertices[i][0], vertices[i][1]) > error) && \
        ((p3[0] == 9999) || (in_circle(p1_x, p1_y, p2_x, p2_y, p3[0], p3[1], vertices[i][0], vertices[i][1]) > error))) {
            for (int j = 0; j < (int) segments[0][2]; j++) {
                if (if_itersection(relative_interior_p[0], relative_interior_p[1], vertices[i][0], vertices[i][1], \
                segments[j][0], segments[j][1], segments[j][2], segments[j][3], error) > error) {
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

    double *triangle = new double [6];
    triangle[0] = p1_x;
    triangle[1] = p1_y;
    triangle[2] = p2_x;
    triangle[3] = p2_y;
    triangle[4] = p3[0];
    triangle[5] = p3[1];

    // if (p3[0] != 9999) {
    //     printf("%0.9f %0.9f ; %0.9f %0.9f ; %0.9f %0.9f;\n",triangle[0],triangle[1],triangle[2],triangle[3],triangle[4],triangle[5]);
    // }

    return triangle;

}


double orientation_2d(double a_x, double a_y, double b_x, double b_y, double c_x, double c_y) {

    return triple_product(\
    a_x, a_y, 1, \
    b_x, b_y, 1, \
    c_x, c_y, 1);
}


double in_circle(double a_x, double a_y, double b_x, double b_y, double c_x, double c_y, double d_x, double d_y) {

    return triple_product(\
    a_x-d_x, a_y-d_y, (a_x-d_x)*(a_x-d_x)+(a_y-d_y)*(a_y-d_y), \
    b_x-d_x, b_y-d_y, (b_x-d_x)*(b_x-d_x)+(b_y-d_y)*(b_y-d_y), \
    c_x-d_x, c_y-d_y, (c_x-d_x)*(c_x-d_x)+(c_y-d_y)*(c_y-d_y));

}


double triple_product(double a_x, double a_y, double a_z, double b_x, double b_y, double b_z, double c_x, double c_y, double c_z) {

    return a_x*b_y*c_z + a_y*b_z*c_x + a_z*b_x*c_y - a_z*b_y*c_x - a_x*b_z*c_y- a_y*b_x*c_z;

}


double if_itersection(double p1_x, double p1_y, double p2_x, double p2_y, double p3_x, double p3_y, double p4_x, double p4_y, double error) {

    double num1 = 0.0;

    if (std::max(p1_x, p2_x) < (std::min(p3_x, p4_x) + error) || \
    std::max(p1_y, p2_y) < (std::min(p3_y, p4_y) + error) || \
    (std::min(p1_x, p2_x) + error) > std::max(p3_x, p4_x) || \
    (std::min(p1_y, p2_y) + error) > std::max(p3_y, p4_y)) {
        printf("four-point colinearity");
        abort();
    }

    if (((p3_x-p1_x)*(p3_y-p4_y) - (p3_y-p1_y)*(p3_x-p4_x)) * ((p3_x-p2_x)*(p3_y-p4_y) - (p3_y-p2_y)*(p3_x-p4_x)) <= 0 & \
    ((p1_x-p3_x)*(p1_y-p2_y) - (p1_y-p3_y)*(p1_x-p2_x)) * ((p1_x-p4_x)*(p1_y-p2_y) - (p1_y-p4_y)*(p1_x-p2_x)) <= 0) {
        num1 = 1.0;
    }

    return num1;

}
