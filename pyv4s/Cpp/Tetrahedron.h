#include <vector>

struct Tetrahedron {
    Vector H1;
    Vector H2;
    Vector L1;
    Vector L2;
};

float determinant_3x3(float matrix[3][3]) {
    return matrix[0][0]*matrix[1][1]*matrix[2][2] + matrix[0][1]*matrix[1][2]*matrix[2][0] + matrix[0][2]*matrix[1][0]*matrix[2][1]
         - matrix[0][2]*matrix[1][1]*matrix[2][0] - matrix[0][0]*matrix[1][2]*matrix[2][1] - matrix[0][1]*matrix[1][0]*matrix[2][2];
}

Vector CramersRule(float matrix[3][3], float m_independents[3]) {
        float A= determinant_3x3(matrix);
        float Ax_matrix[3][3]= {
            { m_independents[0] , matrix[0][1] , matrix[0][2] },
            { m_independents[1] , matrix[1][1] , matrix[1][2] },
            { m_independents[2] , matrix[2][1] , matrix[2][2] }
        };
        float Ax= determinant_3x3(Ax_matrix);
        float Ay_matrix[3][3]= {
            { matrix[0][0] , m_independents[0] , matrix[0][2] },
            { matrix[1][0] , m_independents[1] , matrix[1][2] },
            { matrix[2][0] , m_independents[2] , matrix[2][2] }
        };
        float Ay= determinant_3x3(Ay_matrix);
        float Az_matrix[3][3]= {
            { matrix[0][0] , matrix[0][1] , m_independents[0] },
            { matrix[1][0] , matrix[1][1] , m_independents[1] },
            { matrix[2][0] , matrix[2][1] , m_independents[2] }
        };
        float Az= determinant_3x3(Az_matrix);
        Vector output;
        output.x= Ax/A;
        output.y= Ay/A;
        output.z= Az/A;
        return output;
}

Vector* checkPBC(Vector O, Vector H1, Vector H2, Vector bounds) {
    Vector* output= new Vector[3];
    output[0].x= O.x;
    output[0].y= O.y;
    output[0].z= O.z;
    Vector H_i[2]= {H1, H2};
    for(int i= 0; i < 2; i++) {
        float dx= H_i[i].x - O.x;
        float dy= H_i[i].y - O.y;
        float dz= H_i[i].z - O.z;
        dx-= round(dx / bounds.x) * bounds.x;
        dy-= round(dy / bounds.y) * bounds.y;
        dz-= round(dz / bounds.z) * bounds.z;
        output[i+1].x= O.x + dx;
        output[i+1].y= O.y + dy;
        output[i+1].z= O.z + dz;
    }
    return output;
}

std::vector<Vector> getPerfectTetrahedron(Vector O_real, Vector H1_real, Vector H2_real, Vector bounds) {
    const float R= 1; //Anstrongs between O and perfect vertices
    const float theta= acos(-1./3.)/2; //Perfect angle for tetrahedron
    const float phy= getAngle(H1_real,O_real,H2_real,bounds) / 2.0; //I need the angle between OH and b

    Vector* pbc_vectors= checkPBC(O_real, H1_real, H2_real, bounds);
    Vector O= pbc_vectors[0];
    Vector H1= pbc_vectors[1];
    Vector H2= pbc_vectors[2];
    delete(pbc_vectors);

    Vector OH1= substracVectors(H1,O);
    Vector OH2= substracVectors(H2,O);
    Vector b= addVectors(productByScalar(OH1,1./magnitudeOfVector(OH1)),productByScalar(OH2,1./magnitudeOfVector(OH2)));
    Vector nu_OH_not_norm= crossProduct(OH1,OH2);
    Vector nu_OH= productByScalar(nu_OH_not_norm,1./magnitudeOfVector(nu_OH_not_norm));
    Vector hydrogens[2]= {OH1,OH2};
    Vector h[2];
    for(int i= 0; i < 2; i++) {
        float A_matriz[3][3]= {
            {      b.x       ,      b.y       ,      b.z       },
            { hydrogens[i].x , hydrogens[i].y , hydrogens[i].z },
            {    nu_OH.x     ,    nu_OH.y     ,    nu_OH.z     }
        };
        float m_indep[3]= {R*magnitudeOfVector(b)*float(cos(theta)), R*magnitudeOfVector(hydrogens[i])*float(cos(theta-phy)), 0.};
        h[i]= addVectors(CramersRule(A_matriz,m_indep),O);
    }
    Vector m_H= productByScalar(addVectors(h[0],h[1]),0.5);
    float delta= magnitudeOfVector(substracVectors(h[0],m_H));
    Vector m_L= substracVectors(productByScalar(O,2.),m_H);
    Vector L1= addVectors(productByScalar(nu_OH,delta),m_L);
    Vector L2= addVectors(productByScalar(nu_OH,-delta),m_L);

    std::vector<Vector> tetrahedron(4);
    tetrahedron[0]= h[0];
    tetrahedron[1]= h[1];
    tetrahedron[2]= L1;
    tetrahedron[3]= L2;
    return tetrahedron;
}
