//
//  Matrix.hpp
//  ACG-Assignment2
//
//  Created by Bernhard Fritz on 07/12/15.
//  Copyright © 2015 Bernhard Fritz. All rights reserved.
//

#ifndef Matrix_h
#define Matrix_h

class Matrix {
public:
    double m[16];
    
    Matrix()
    {
        double tmp[16] =
        {
            1.0, 0.0, 0.0, 0.0,
            0.0, 1.0, 0.0, 0.0,
            0.0, 0.0, 1.0, 0.0,
            0.0, 0.0, 0.0, 1.0
        };
        
        memcpy(m, tmp, 16 * sizeof(double));
    }
    
    Matrix(double* _m)
    {
        memcpy(m, _m, 16 * sizeof(double));
    }
    
    double &operator[](int i)
    {
        return m[i];
    }
    
    Vector operator*(Vector v) const
    {
        double tmp[4] =
        {
            v.x, v.y, v.z, 1.0
        };
        
        return Vector(m[0] * tmp[0] + m[1] * tmp[1] + m[2] * tmp[2] + m[3] * tmp[3],
                      m[4] * tmp[0] + m[5] * tmp[1] + m[6] * tmp[2] + m[7] * tmp[3],
                      m[8] * tmp[0] + m[9] * tmp[1] + m[9] * tmp[2] + m[10] * tmp[3]);
    }
    
    Matrix operator*(Matrix _m) const
    {
        double tmp[16];
        
        tmp[0] = m[0] * _m[0] + m[1] * _m[4] + m[2] * _m[8] + m[3] * _m[12];
        tmp[1] = m[0] * _m[1] + m[1] * _m[5] + m[2] * _m[9] + m[3] * _m[13];
        tmp[2] = m[0] * _m[2] + m[1] * _m[6] + m[2] * _m[10] + m[3] * _m[14];
        tmp[3] = m[0] * _m[3] + m[1] * _m[7] + m[2] * _m[11] + m[3] * _m[15];
        
        tmp[4] = m[4] * _m[0] + m[5] * _m[4] + m[6] * _m[8] + m[7] * _m[12];
        tmp[5] = m[4] * _m[1] + m[5] * _m[5] + m[6] * _m[9] + m[7] * _m[13];
        tmp[6] = m[4] * _m[2] + m[5] * _m[6] + m[6] * _m[10] + m[7] * _m[14];
        tmp[7] = m[4] * _m[3] + m[5] * _m[7] + m[6] * _m[11] + m[7] * _m[15];
        
        tmp[8] = m[8] * _m[0] + m[9] * _m[4] + m[10] * _m[8] + m[11] * _m[12];
        tmp[9] = m[8] * _m[1] + m[9] * _m[5] + m[10] * _m[9] + m[11] * _m[13];
        tmp[10] = m[8] * _m[2] + m[9] * _m[6] + m[10] * _m[10] + m[11] * _m[14];
        tmp[11] = m[8] * _m[3] + m[9] * _m[7] + m[10] * _m[11] + m[11] * _m[15];
        
        tmp[12] = m[12] * _m[0] + m[13] * _m[4] + m[14] * _m[8] + m[15] * _m[12];
        tmp[13] = m[12] * _m[1] + m[13] * _m[5] + m[14] * _m[9] + m[15] * _m[13];
        tmp[14] = m[12] * _m[2] + m[13] * _m[6] + m[14] * _m[10] + m[15] * _m[14];
        tmp[15] = m[12] * _m[3] + m[13] * _m[7] + m[14] * _m[11] + m[15] * _m[15];
        
        return Matrix(tmp);
    }
    
};

#endif /* Matrix_h */
