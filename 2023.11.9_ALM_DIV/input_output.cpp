#include "input_output.hpp"

double **read_boundary(const char *file_name) {
    printf("reading the fourier coefficients of boundary...\n");

    char char1;
    int position_initial, row = 0, num1 = 0;
    string file_str, str1;

    fstream file1;
    file1.open(file_name,std::ios::in);
    if (!file1.is_open()) {
        printf("read_boundary: file does not exist or has improper access permission\n");
        abort();
    }
    
    // import data to a string variable
    for (;;) {

        char1 = file1.get();
        
        if (char1 == '&') {
            break;
        }
        else if (char1 == '/') {
            char1 = file1.get();
            if (char1 == '/') {
                while (char1 != '\n') {
                    char1 = file1.get();
                }
                char1 = file1.get();
            }
            else {
                file1.seekg(-1,std::ios::cur);
            }
        }

        if (num1 > 99999) {
            printf("read_boundary: wrong format\n");
            exit(0);
        }

        file_str = file_str + char1;
        num1++;
    }

    file1.close();

    row = std::count(file_str.begin(),file_str.end(),'\n');
    
    // construct fourier coefficients
    double **coefficients = new double* [row];
    coefficients[0] = new double [7];
    for (int j = 0; j < 7; j++) {
        coefficients[0][j] = 0;
    }
    for (int i = 1; i < row; i++) {
        coefficients[i] = new double [6];
        for (int j = 0; j < 6; j++) {
            coefficients[i][j] = 0;
        }
    }
    coefficients[0][6] = row;

    for (int i = 0, k = 0; k < file_str.size() - 1; k++) {

        for (int j = 0; j < 6; j++) {
    
            while (file_str[k] == ' ' || file_str[k] == '\t') {
                k++;
            }

            while (file_str[k] != ' ' && file_str[k] != '\t' && file_str[k] != '\n'){
                str1 += file_str[k];
                k++;
            } 
           
            coefficients[i][j] = stod(str1);
            str1.clear(); 

            while (file_str[k] == ' ' || file_str[k] == '	') {
                k++;
            }               
                
            if (file_str[k] == '\n') {
                i++;
                break;
            }
        }
    }
    
    // for (int i = 0; i < row; i++) {
    //     printf("rbc(%d,%d) = %e, zbs(%d,%d) = %e, rbs(%d,%d) = %e, zbc(%d,%d) = %e\n",\
    //     (int) coefficients[i][0], (int) coefficients[i][1], coefficients[i][2], \
    //     (int) coefficients[i][0], (int) coefficients[i][1], coefficients[i][3], \
    //     (int) coefficients[i][0], (int) coefficients[i][1], coefficients[i][4], \
    //     (int) coefficients[i][0], (int) coefficients[i][1], coefficients[i][5]);
    // }
    // exit(0);
    
    return coefficients;
}

double **fourier_coefficients_dzeta(double **fourier_coefficients, int nfp) {
    printf("...\n");
    int num_coefficients = (int) fourier_coefficients[0][6];
    double **fc = new double *[num_coefficients];
    for (int i = 0; i < num_coefficients; i++) {
        if (i == 0) {
            fc[i] = new double [7]();
            fc[i][6] = num_coefficients;
        }
        else {
            fc[i] = new double [6]();
        }
        // m * theta + n * zeta
        fc[i][0] =  fourier_coefficients[i][0];
        fc[i][1] =  fourier_coefficients[i][1];
        fc[i][2] =  fourier_coefficients[i][0] * fourier_coefficients[i][4] * nfp;
        fc[i][3] = -fourier_coefficients[i][0] * fourier_coefficients[i][5] * nfp;
        fc[i][4] = -fourier_coefficients[i][0] * fourier_coefficients[i][2] * nfp;
        fc[i][5] =  fourier_coefficients[i][0] * fourier_coefficients[i][3] * nfp;
    }

    // for (int i = 0; i < num_coefficients; i++) {
    //     printf("rbc(%d,%d) = %f, zbs(%d,%d) = %f, rbs(%d,%d) = %f, zbc(%d,%d) = %f\n",\
    //     (int) fc[i][0], (int) fc[i][1], fc[i][2], \
    //     (int) fc[i][0], (int) fc[i][1], fc[i][3], \
    //     (int) fc[i][0], (int) fc[i][1], fc[i][4], \
    //     (int) fc[i][0], (int) fc[i][1], fc[i][5]);
    // }
    // exit(0);

    return fc;
}

double **fourier_coefficients_dtheta(double **fourier_coefficients, int nfp) {
    printf("...\n");
    int num_coefficients = (int) fourier_coefficients[0][6];
    double **fc = new double *[num_coefficients];
    for (int i = 0; i < num_coefficients; i++) {
        if (i == 0) {
            fc[i] = new double [7]();
            fc[i][6] = num_coefficients;
        }
        else {
            fc[i] = new double [6]();
        }

        fc[i][0] =  fourier_coefficients[i][0];
        fc[i][1] =  fourier_coefficients[i][1];
        fc[i][2] =  fourier_coefficients[i][1] * fourier_coefficients[i][4];
        fc[i][3] = -fourier_coefficients[i][1] * fourier_coefficients[i][5];
        fc[i][4] = -fourier_coefficients[i][1] * fourier_coefficients[i][2];
        fc[i][5] =  fourier_coefficients[i][1] * fourier_coefficients[i][3];
    }

    return fc;
}

double ***read_coils_p(const char *file_name, int inflation, double tao) {
    printf("importing the points of coils...\n");

    std::fstream file1;
    file1.open(file_name, std::ios::in);
    if (!file1.is_open()) {
        printf("read_coils_p: file does not exist or has improper access permission\n");
        abort();
    }

    string file_str;
    char char1;
    int num1 = 0, row = 0;
    for (;;) {
        char1 = file1.get();
        
        if (char1 == '&') {
            break;
        }
        else if (char1 == '/') {
            char1 = file1.get();
            if (char1 == '/') {
                while (char1 != '\n') {
                    char1 = file1.get();
                }
                char1 = file1.get();
            }
            else {
                file1.seekg(-1,std::ios::cur);
            }
        }

        if (num1 > 9999999) {
            printf("read_boundary: wrong format or too many points\n");
            // printf("%s\n",file_str.c_str());
            exit(0);
        }

        file_str = file_str + char1;
        num1++;
    }
    file1.close();
    // std::cout<<file_str<<std::endl;
    // exit(0);  

    row = std::count(file_str.begin(),file_str.end(),'\n');

    double **array1 = new double *[row];
    for (int i = 0; i < row; i++) {
        array1[i] = new double [4];
        for (int j = 0; j < 4; j++) {
            array1[i][j] = 0;
        }
    }

    string str1, str2; 
    for (int i = 0, j = 0, k = 0; ; ) {
        str2 = file_str[k];
        while (str2.find_first_of("0123456789Ee+-.") == string::npos) {
            if (str2 == "\n") {
                i++;
            }
            k++;
            str2 = file_str[k];
        }

        while (str2.find_first_of("0123456789Ee+-.") != string::npos) {
            str1 += str2;
            k++;
            str2 = file_str[k];
        }

        if (str1.find_first_of("0123456789Ee+-.") == string::npos) {
            break;
        }
        else if (k >= file_str.size()) {
            break;
        }
        array1[i][j] = stod(str1);
        str1.clear(); 
        
        j++;
        if (j != 0 && j%4 == 0) {
            j = 0;
        }
    }

    // for (int i = 0; i < row; i++) {
    //     printf("%0.9f %0.9f %0.9f %0.9f\n", array1[i][0], array1[i][1], array1[i][2], array1[i][3]);
    // }
    // exit(0);    

    double array2[3] = {array1[0][0], array1[0][1], array1[0][2]};
    int num2 = 0, *array3 = new int [1];
    array3[0] = 0;
    for (int i = 1; i < row; i++) {
        if (array2[0] == array1[i][0] && array2[1] == array1[i][1] && array2[2] == array1[i][2]) {
            array3[num2] = i;
            num2++;
            array3 = enlarge_1d(array3, num2, num2 + 1);
            if (i < row - 1) {
                array2[0] = array1[i+1][0];
                array2[1] = array1[i+1][1];
                array2[2] = array1[i+1][2];
                i++;
            }
        }
    }
    // printf("num2 = %d\n", num2);
    // for (int i = 0; i < num2; i++) {
    //     printf("array3[%d] = %d\n", i, array3[i]);
    // }
    // exit(0);

    int num3 = array3[0];
    for (int i = 1; i < num2; i++) {
        if (array3[i] - array3[i-1] > num3) {
            num3 = array3[i] - array3[i-1];
        }
    }
    // rearrangement
    int num5 = 0;
    double ***array4 = new double **[num2];
    for (int i = 0; i < num2; i++) {
        array4[i] = new double *[num3];
        for (int j = 0; j < num3; j++) {
            array4[i][j] = new double [4] (); 
        }

        if (i == 0) {
            num5 = array3[0] + 1;
        }
        else {
            num5 = array3[i] - array3[i-1];
        }
        for (int j = 0; j < num5; j++) {
            array4[i][j][0] = array1[array3[i] - num5 + 1 + j][0];
            array4[i][j][1] = array1[array3[i] - num5 + 1 + j][1];
            array4[i][j][2] = array1[array3[i] - num5 + 1 + j][2];
            array4[i][j][3] = array1[array3[i] - num5 + 1 + j][3];
            if (i == 0) {
                if (j == array3[0]) {
                    break;
                }
            }
            else {
                if (j == array3[i] - 1) {
                    break;
                }
            }
        }
    }

    // for (int k = 0; k < num2; k++) {
    //     for (int i = 0; i < num3; i++) {
    //         printf("%0.9f %0.9f %0.9f %0.9f\n", array4[k][i][0],array4[k][i][1],array4[k][i][2],array4[k][i][3]);
    //     }
    //     printf("\n\n\n");
    // }
    // exit(0);

    // catmull-rom spline
    double ***array5 = new double **[num2];
    int num4 = (inflation+1) * (num3-1) + 1;
    double c0[3], c1[3], c2[3], c3[3];
    double p0[3], p1[3], p2[3], p3[3];
    double step = 1.0 / (inflation+1);
    for (int i = 0; i < num2; i++) {
        array5[i] = new double *[num4];
        for (int j = 0; j < num4; j++) {
            array5[i][j] = new double [4] ();
        }

        if (i == 0) {
            num5 = array3[0] + 1;
        }
        else {
            num5 = array3[i] - array3[i-1];
        }
        for (int j = 0; j < num5 - 1; j++) {
            if (j == 0) {
                for (int k = 0; k < 3; k++) {
                    p0[k] = array1[array3[i]-1][k];
                    p1[k] = array4[i][j][k];
                    p2[k] = array4[i][1][k];
                    p3[k] = array4[i][2][k];
                }
            }
            else if (j == num5 - 1) {
                for (int k = 0; k < 3; k++) {
                    p0[k] = array4[i][j-1][k];
                    p1[k] = array4[i][j][k];
                    p2[k] = array4[i][1][k];
                    p3[k] = array4[i][2][k];
                }
            }
            else if (j == num5 - 2) {
                for (int k = 0; k < 3; k++) {
                    p0[k] = array4[i][j-1][k];
                    p1[k] = array4[i][j][k];
                    p2[k] = array4[i][0][k];
                    p3[k] = array4[i][1][k];
                }
            }
            else if (j == num5 - 3) {
                for (int k = 0; k < 3; k++) {
                    p0[k] = array4[i][j-1][k];
                    p1[k] = array4[i][j][k];
                    p2[k] = array4[i][j+1][k];
                    p3[k] = array4[i][0][k];
                }
            }
            else {
                for (int k = 0; k < 3; k++) {
                    p0[k] = array4[i][j-1][k];
                    p1[k] = array4[i][j][k];
                    p2[k] = array4[i][j+1][k];
                    p3[k] = array4[i][j+2][k];
                }
            }

            for (int k = 0; k < 3; k++) {
                c0[k] = p1[k];
                c1[k] = (-tao) * p0[k] + (tao) * p2[k];
                c2[k] = (2.0*tao) * p0[k] + (tao-3.0) * p1[k] + (3.0-2.0*tao) * p2[k] + (-tao) * p3[k];
                c3[k] = (-tao) * p0[k] + (2.0-tao) * p1[k] + (tao-2.0) * p2[k] + (tao) * p3[k];
                for (int l = 0; l < inflation + 1; l++) {
                    array5[i][(inflation+1)*j+l][k] =  c0[k] + c1[k]*power_nonnegative(step*l,1) + c2[k]*power_nonnegative(step*l,2) + c3[k]*power_nonnegative(step*l,3);
                }    
            }
            for (int l = 0; l < inflation + 1; l++) {
                array5[i][(inflation+1)*j+l][3] = array4[i][j][3];
            }
            
        }
        for (int k = 0; k < 4; k++) {
            array5[i][(inflation+1)*(num5-1)][k] = array4[i][num5-1][k];
        }
    }

    array5[0][0] = enlarge_1d(array5[0][0], 4, 6);
    array5[0][0][4] = num2;
    array5[0][0][5] = num4;

    // printf("array5[0][0][4] = %d, array5[0][0][4] = %d\n", (int) array5[0][0][4], (int) array5[0][0][5]);
    // for (int k = 0; k < num2; k++) {
    //     for (int i = 0; i < num4; i++) {
    //         printf("%0.9f %0.9f %0.9f %0.9f\n", array5[k][i][0],array5[k][i][1],array5[k][i][2],array5[k][i][3]);
    //     }
    //     // printf("\n\n\n");
    // }
    // exit(0);

    for (int i = 0; i < row; i++) {
        delete [] array1[i];
    }
    delete [] array1;

    delete [] array3;

    for (int i = 0; i < num2; i++) {
        for (int j = 0; j < num3; j++) {
            delete [] array4[i][j];
        }
        delete [] array4[i];
    }    
    delete [] array4;

    return array5;
}

